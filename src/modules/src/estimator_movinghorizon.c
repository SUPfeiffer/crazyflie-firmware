/**
 * Model Predictive Localization with UWB based on S.Li(2019)
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

#include "stabilizer.h"
#include "estimator_movinghorizon.h"
#include "sensfusion6.h"
#include "sensors.h"
#include "position_estimator.h"
#include "log.h"

#include <stdlib.h>
#include <math.h>
#include "math3d.h"
#include "cf_math.h"

#include "stm32f4xx.h"
#include "FreeRTOS.h"
#include "queue.h"
#include "task.h"


// Measurements of a UWB Tx/Rx
static xQueueHandle tdoaDataQueue;
#define UWB_QUEUE_LENGTH (10)

static inline bool stateEstimatorHasTDOAPacket(tdoaMeasurement_t *uwb) {
  return (pdTRUE == xQueueReceive(tdoaDataQueue, uwb, 0));
}


#define ATTITUDE_UPDATE_RATE RATE_250_HZ
#define ATTITUDE_UPDATE_DT 1.0f/ATTITUDE_UPDATE_RATE

#define POS_UPDATE_RATE RATE_100_HZ
#define POS_UPDATE_DT 1.0f/POS_UPDATE_RATE

#define MOVING_HORIZON_MIN_WINDOW_SIZE 4
#define MOVING_HORIZON_MAX_WINDOW_SIZE 20
#define NEWTON_RAPHSON_THRESHOLD 0.001f
#define RANSAC_ITERATIONS 5
#define RANSAC_ERROR_THRESHOLD 0.02f    // Maximum accounted error per timestep when comparing RANSAC iterations

#define CONST_G 9.81f
#define CONST_K_AERO 0.35f

// shift variables in the window by one step
void movingHorizonShiftWindow(float *window, int windowSize);

// for position correction with any number of tdoa measurements (less accurate)
void projectOnTDOA(tdoaMeasurement_t *tdoa, point_t *original, point_t *projection);

// position calculation with > 3 tdoa measurements (returns false otherwise) using Newton-Raphson
// for multilateration
bool positionFromTDOA(point_t prediction, point_t *measurement);

// calculate a new error model [delta_x,delta_vx] based on the positioning errors in errorWindow
// at the times in timeWindow
void updateErrorModel(float *errorWindow, float *timeWindow, int8_t windowSize, float *errorModel);

// Solves the least squares problem y=mx+b
// N: Length of Vectors x and y
// params[2] = [b,m]
void linearLeastSquares(float *x, float *y, float N, float *params);

static bool isInit = false;
static point_t loc_prediction;
static point_t loc_measurement;
static velocity_t vel_prediction;

// variables in moving windows
static float time_w[MOVING_HORIZON_MAX_WINDOW_SIZE];
static float dx_w[MOVING_HORIZON_MAX_WINDOW_SIZE];
static float dy_w[MOVING_HORIZON_MAX_WINDOW_SIZE];
 
static uint8_t windowSize;
static float samplingRatio = 0.5;

// Error estimations
static float errorModel_x[2];
static float errorModel_y[2];

void estimatorMovingHorizonInit(void)
{
    if(!isInit){
        tdoaDataQueue = xQueueCreate(UWB_QUEUE_LENGTH, sizeof(tdoaMeasurement_t));
    }
    else{
        xQueueReset(tdoaDataQueue);
    }

    sensfusion6Init();

    // reset variables
    loc_prediction.x = 0.0f; loc_prediction.y = 0.0f; loc_prediction.z = 0.0f;
    loc_measurement.x = 0.0f; loc_measurement.y = 0.0f; loc_measurement.z = 0.0f;
    vel_prediction.x = 0.0f; vel_prediction.y = 0.0f; vel_prediction.z = 0.0f;
    int8_t i;
    for (i=0;i<MOVING_HORIZON_MAX_WINDOW_SIZE;i++){
        time_w[i] = 0.0f;
        dx_w[i] = 0.0f;
        dy_w[i] = 0.0f;
    }
    windowSize = 0;
    errorModel_x[0] = 0.0f;
    errorModel_x[1] = 0.0f;
    errorModel_y[0] = 0.0f;
    errorModel_y[1] = 0.0f;

    isInit = true;
}

bool estimatorMovingHorizonTest(void)
{
    bool pass = true;
    pass &= sensfusion6Test();
    return pass;
}

void estimatorMovingHorizon(state_t *state, sensorData_t *sensorData, control_t *control, const uint32_t tick)
{
    // generate attitude estimate with complementary filter
    sensorsAcquire(sensorData, tick); // Read sensors at full rate (1000Hz)
    if (RATE_DO_EXECUTE(ATTITUDE_UPDATE_RATE, tick)) {
        sensfusion6UpdateQ(sensorData->gyro.x, sensorData->gyro.y, sensorData->gyro.z,
                       sensorData->acc.x, sensorData->acc.y, sensorData->acc.z,
                       ATTITUDE_UPDATE_DT);

        // Save attitude, adjusted for the legacy CF2 body coordinate system
        sensfusion6GetEulerRPY(&state->attitude.roll, &state->attitude.pitch, &state->attitude.yaw);

        // Save quaternion, hopefully one day this could be used in a better controller.
        // Note that this is not adjusted for the legacy coordinate system
        sensfusion6GetQuaternion(
        &state->attitudeQuaternion.x,
        &state->attitudeQuaternion.y,
        &state->attitudeQuaternion.z,
        &state->attitudeQuaternion.w);

        state->acc.z = sensfusion6GetAccZWithoutGravity(sensorData->acc.x,
                                                        sensorData->acc.y,
                                                        sensorData->acc.z);
        // keep position/velocity updates in z for now (velocity update not necessary here?)
        // positionUpdateVelocity(state->acc.z, ATTITUDE_UPDATE_DT);
    }

    if (RATE_DO_EXECUTE(POS_UPDATE_RATE, tick)) {
        // update z position
        positionEstimate(state, sensorData, POS_UPDATE_DT, tick);
        
        // predict current state
        loc_prediction.x += vel_prediction.x * POS_UPDATE_DT;
        loc_prediction.y += vel_prediction.y * POS_UPDATE_DT;
        vel_prediction.x += (-CONST_G*tanf(state->attitude.pitch) - CONST_K_AERO*vel_prediction.x) * POS_UPDATE_DT;
        vel_prediction.x += (CONST_G*tanf(state->attitude.roll) - CONST_K_AERO*vel_prediction.y) * POS_UPDATE_DT;

        time_w[windowSize] = xTaskGetTickCount(); // Do we need to account for wrapping?
        
        // Only execute filter if we can calculate a new position from tdoa measurements      
        if (positionFromTDOA(loc_prediction,&loc_measurement)){

            if (windowSize == MOVING_HORIZON_MAX_WINDOW_SIZE){
                movingHorizonShiftWindow(dx_w, windowSize);
                movingHorizonShiftWindow(dy_w, windowSize);
                movingHorizonShiftWindow(time_w, windowSize);
            }
            else{ 
                windowSize += 1;
            }
            
            dx_w[windowSize] = loc_measurement.x - loc_prediction.x;
            dy_w[windowSize] = loc_measurement.y - loc_prediction.y;
            
            if (windowSize >= MOVING_HORIZON_MIN_WINDOW_SIZE){
                updateErrorModel(dx_w, time_w, windowSize, errorModel_x);
                updateErrorModel(dy_w, time_w, windowSize, errorModel_y);
            }
        }
        // correction of prediction
        state->position.x = loc_prediction.x + errorModel_x[0] + errorModel_x[1] * (time_w[windowSize]-time_w[0]);
        state->position.y = loc_prediction.y + errorModel_y[0] + errorModel_y[1] * (time_w[windowSize]-time_w[0]);
        state->velocity.x = vel_prediction.x + errorModel_x[1];
        state->velocity.y = vel_prediction.y + errorModel_y[1];
    }
}


void movingHorizonShiftWindow(float *window, int windowSize){
    uint8_t i;
    for (i=0; i<windowSize; i++){
        window[i] = window[i+1];
    }
}

// same function already exist in estimator_kalman.c
// might make sense to transfer it to estimator.c to make it easily available to other estimators (avoids dependencies on kalman filter)
static bool stateEstimatorEnqueueExternalMeasurement(xQueueHandle queue, void *measurement)
{
  portBASE_TYPE result;
  bool isInInterrupt = (SCB->ICSR & SCB_ICSR_VECTACTIVE_Msk) != 0;

  if (isInInterrupt) {
    portBASE_TYPE xHigherPriorityTaskWoken = pdFALSE;
    result = xQueueSendFromISR(queue, measurement, &xHigherPriorityTaskWoken);
    if(xHigherPriorityTaskWoken == pdTRUE)
    {
      portYIELD();
    }
  } else {
    result = xQueueSend(queue, measurement, 0);
  }
  return (result==pdTRUE);
}

bool estimatorMovingHorizonEnqueueTDOA(const tdoaMeasurement_t *uwb){
      return stateEstimatorEnqueueExternalMeasurement(tdoaDataQueue, (void *)uwb);
}

void projectOnTDOA(tdoaMeasurement_t *tdoa, point_t *origin, point_t *projection){
    float tdoa_measured = tdoa->distanceDiff;

    float x1 = tdoa->anchorPosition[1].x, y1 = tdoa->anchorPosition[1].y, z1 = tdoa->anchorPosition[1].z;
    float x0 = tdoa->anchorPosition[0].x, y0 = tdoa->anchorPosition[0].y, z0 = tdoa->anchorPosition[0].z;
    float x = origin->x, y = origin->y, z = origin->z;
        
    float dx1 = x - x1;
    float dy1 = y - y1;
    float dz1 = z - z1;

    float dy0 = y - y0;
    float dx0 = x - x0;
    float dz0 = z - z0;

    float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
    float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2));
    vector_t grad = {
        .x = (dx1/d1 - dx0/d0),
        .y = (dy1/d1 - dy0/d0), 
        .z = (dz1/d1 - dz0/d0),};
        
    // Newton-Raphson with tdoa gradient constant (value @originalPoint)
    // Initial Guess: original point (t=0)
    float error = d1 - d0 - tdoa_measured;
    float error_der;
    float t = 0;
    while (fabsf(error) >= 0.001f)
    {
        error_der = (grad.x * (dx1 + t*grad.x) + grad.y * (dy1 + t*grad.y) + grad.z * (dz1 + t*grad.z))/d1
                    - (grad.x * (dx0 + t*grad.x) + grad.y * (dy0 + t*grad.y) + grad.z * (dz0 + t*grad.z))/d0;

        t -= error/error_der;

        d1 = sqrtf(powf(dx1+t*grad.x,2)+powf(dy1+t*grad.y,2)+powf(dz1+t*grad.z,2));
        d0 = sqrtf(powf(dx0+t*grad.x,2)+powf(dy0+t*grad.y,2)+powf(dz0+t*grad.z,2));
         
        error = d1 - d0 - tdoa_measured;
    }
        
    projection->x = x + t*grad.x;
    projection->y = y + t*grad.y;
    projection->z = z + t*grad.z;
}

bool positionFromTDOA(point_t prediction, point_t *measurement){
    
    if (uxQueueMessagesWaiting(tdoaDataQueue) <3){ 
        return false;
    }

    // get all measurements from queue
    tdoaMeasurement_t tdoa[UWB_QUEUE_LENGTH];
    int8_t tdoa_count = 0;
    while (stateEstimatorHasTDOAPacket(&tdoa[tdoa_count]) && tdoa_count<UWB_QUEUE_LENGTH){
        tdoa_count++;
    }

    
    float uwb_position[3] = {prediction.x, prediction.y, prediction.z};
    float threshold = NEWTON_RAPHSON_THRESHOLD;
    
    float delta[3];
    arm_matrix_instance_f32 delta_m = {3, 1, (float *)delta};
    
    float S_d[tdoa_count];
    arm_matrix_instance_f32 S_m = {tdoa_count, 1, (float *)S_d};

    float J_d[tdoa_count*3];
    arm_matrix_instance_f32 J_m = {tdoa_count, 3, (float *)J_d};
    
    float JT_d[tdoa_count*3];
    arm_matrix_instance_f32 JT_m = {3, tdoa_count, (float *)JT_d};
    
    float JTJ_d[3*3];
    arm_matrix_instance_f32 JTJ_m = {3, 3, (float *)JTJ_d};

    float JTJinv_d[3*3];
    arm_matrix_instance_f32 JTJinv_m ={3, 3, (float *)JTJinv_d};

    float JTJpinv_d[tdoa_count*3];
    arm_matrix_instance_f32 JTJpinv_m = {3, tdoa_count, (float *)JTJpinv_d};

    int8_t i;
    do{
        for (i=0;i<tdoa_count;i++){
            float x0 = tdoa[i].anchorPosition[0].x, y0 = tdoa[i].anchorPosition[0].y, z0 = tdoa[i].anchorPosition[0].z;
            float x1 = tdoa[i].anchorPosition[1].x, y1 = tdoa[i].anchorPosition[1].y, z1 = tdoa[i].anchorPosition[1].z;
            
            float dx1 = uwb_position[0] - x1;
            float dy1 = uwb_position[1] - y1;
            float dz1 = uwb_position[2] - z1;

            float dx0 = uwb_position[0] - x0;
            float dy0 = uwb_position[1] - y0;
            float dz0 = uwb_position[2] - z0;  
            
            float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
            float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2));

            S_d[i] = d1 - d0 - tdoa[i].distanceDiff;
            J_d[i+0] = (dx1/d1 - dx0/d0);
            J_d[i+1] = (dy1/d1 - dy0/d0);
            J_d[i+2] = (dz1/d1 - dz0/d0);
        }

        arm_mat_trans_f32(&J_m, &JT_m);                 // J'
        arm_mat_mult_f32(&JT_m,&J_m,&JTJ_m);            // (J'J)
        arm_mat_inverse_f32(&JTJ_m,&JTJinv_m);          // (J'J)-1
        arm_mat_mult_f32(&JTJinv_m,&JT_m,&JTJpinv_m);   // (J'J)-1 *J'
        arm_mat_mult_f32(&JTJpinv_m,&S_m,&delta_m);     // delta = [(J'J)-1 J]*S 
        
        uwb_position[0] -= delta[0];
        uwb_position[1] -= delta[1];
        uwb_position[2] -= delta[2];
        
    } while(delta[0] > threshold || delta[1] > threshold || delta[2] > threshold);

    measurement->x = uwb_position[0];
    measurement->y = uwb_position[1];
    measurement->z = uwb_position[2];
    
    return true;
}

void updateErrorModel(float *errorWindow, float *timeWindow, int8_t windowSize, float *errorModel){
    uint8_t i;
    float dt_w[windowSize];
    for (i=0;i<windowSize;i++){
        dt_w[i] = timeWindow[i]-timeWindow[0]; 
    }

    // RANSAC
    uint8_t N_samples = (uint8_t)floor(windowSize*samplingRatio);
    if (N_samples < 2) { N_samples = 2; }
    
    uint8_t ransac_idx[N_samples];
    float totalError, stepError;
    float bestIteration = windowSize * RANSAC_ERROR_THRESHOLD;
    float sampleErrorModel[2];
    for (i=0;i<RANSAC_ITERATIONS;i++){
        // select N_sample unique samples
        float dt_s[N_samples], dx_s[N_samples];
        uint8_t j = 0;
        while(j<N_samples){
            ransac_idx[j] = (uint8_t)floor(windowSize * rand()/RAND_MAX );
            bool isUnique = true;
            for(int8_t k=0;k<j;k++){
                isUnique &= (ransac_idx[j] != ransac_idx[k]);
            }
            if (isUnique){
                dt_s[j] = dt_w[ransac_idx[j]];
                dx_s[j] = errorWindow[ransac_idx[j]];
                j++;
            }
        }
        
        linearLeastSquares(dt_s, dx_s, N_samples, sampleErrorModel);

        totalError = 0;
        for (j=0;j<windowSize;j++){
            stepError = abs(sampleErrorModel[1] * dt_w[j] + sampleErrorModel[0] - dx_w[j]);
            if (stepError > RANSAC_ERROR_THRESHOLD){ stepError = RANSAC_ERROR_THRESHOLD; }
            
            totalError += stepError;
        }

        if (totalError < bestIteration){
            bestIteration = totalError;
            errorModel[0] = sampleErrorModel[0];
            errorModel[1] = sampleErrorModel[1];
        }
    }
}

void linearLeastSquares(float *x, float *y, float N, float *params){
    // Linear Least Squares (no RANSAC)
    float sumX = 0.0f;
    float sumX2 = 0.0f;
    float sumY = 0.0f;
    float sumXY = 0.0f;
        
    uint8_t i;
    for (i=0;i<N;i++){
        sumX += x[i];
        sumX2 += powf(x[i],2);
        sumY += y[i];
        sumXY += x[i] * y[i];
    }

    float denom = N * sumX2 - powf(sumX,2);
    
    params[0] = (sumX2 * sumY - sumXY * sumX) / denom;
    params[1] = (N * sumXY - sumX * sumY) / denom;
}


// Log group
LOG_GROUP_START(mhe)
    LOG_ADD(LOG_FLOAT, uwbX, &loc_measurement.x)
    LOG_ADD(LOG_FLOAT, uwbY, &loc_measurement.y)
    LOG_ADD(LOG_FLOAT, uwbZ, &loc_measurement.z)
    LOG_ADD(LOG_FLOAT, dx, &errorModel_x[0])
    LOG_ADD(LOG_FLOAT, dvx, &errorModel_x[1])
    LOG_ADD(LOG_FLOAT, dy, &errorModel_y[0])
    LOG_ADD(LOG_FLOAT, dvy, &errorModel_y[1])
LOG_GROUP_STOP(mhe)