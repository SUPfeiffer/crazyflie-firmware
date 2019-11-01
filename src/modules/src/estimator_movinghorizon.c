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
#include <math.h>

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
#define ATTITUDE_UPDATE_DT 1.0/ATTITUDE_UPDATE_RATE

#define POS_UPDATE_RATE RATE_100_HZ
#define POS_UPDATE_DT 1.0/POS_UPDATE_RATE

#define MOVING_HORIZON_MAX_WINDOW_SIZE 10

#define CONST_G 9.81f
#define CONST_K_AERO 0.35f

struct stateWindow_s {
    float x[MOVING_HORIZON_MAX_WINDOW_SIZE];
    float y[MOVING_HORIZON_MAX_WINDOW_SIZE];
    float vx[MOVING_HORIZON_MAX_WINDOW_SIZE];
    float vy[MOVING_HORIZON_MAX_WINDOW_SIZE];
};

void movingHorizonShiftWindow(float *window[], int windowSize);
void projectOnTDOA(tdoaMeasurement_t *tdoa, point_t *original, point_t *projection);

point_t loc_prediction ={
    .x = 0.0f,
    .y = 0.0f,
    .z = 0.0f,
};
point_t loc_measurement ={
    .x = 0.0f,
    .y = 0.0f,
    .z = 0.0f,
};
velocity_t vel_prediction ={
    .x = 0.0f,
    .y = 0.0f,
    .z = 0.0f,
};

struct stateWindow_s delta ={
    .x = {0.0f},
    .y = {0.0f},
    .vx = {0.0f},
    .vy = {0.0f},
};

uint8_t currentWindowSize = 0;

float ls_sumX = 0;
float ls_sumX2 = 0;
float lsx_sumY = 0;
float lsx_sumXY = 0;
float lsy_sumY = 0;
float lsy_sumXY = 0;
float ls_denom = 0;

float errorEst_x = 0;
float errorEst_vx = 0;
float errorEst_y = 0;
float errorEst_vy = 0;

void estimatorMovingHorizonInit(void)
{
    sensfusion6Init();
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
    
    // Moving Horizon predictor for x,y
    if (currentWindowSize == MOVING_HORIZON_MAX_WINDOW_SIZE){
        // remove discarded timestep from Least Squares Sums
        lsx_sumY -= delta.x[0];
        lsy_sumY -= delta.y[0];
        // Shift window
        movingHorizonShiftWindow(delta.x,currentWindowSize);
        movingHorizonShiftWindow(delta.y,currentWindowSize);
    }
    else{
        currentWindowSize += 1;
        ls_sumX += POS_UPDATE_DT;
        ls_sumX2 += POS_UPDATE_DT*POS_UPDATE_DT;
        ls_denom = currentWindowSize * LS_sumX2 - LS_sumX * LS_sumX;
    }

    // predict next state
    uint8_t k = currentWindowSize;
    loc_prediction.x += vel_prediction.x * POS_UPDATE_DT;
    loc_prediction.y += vel_prediction.y * POS_UPDATE_DT;
    vel_prediction.x += (-CONST_G*tan(state->attitude.pitch) - CONST_K_AERO*vel_prediction.x) * POS_UPDATE_DT;
    vel_prediction.x += (CONST_G*tan(state->attitude.roll) - CONST_K_AERO*vel_prediction.y) * POS_UPDATE_DT;
    
    // Get UWB position estimates
    tdoaMeasurement_t tdoa;
    point_t projection[UWB_QUEUE_LENGTH]
    vector_t sum_for_avg = {.x=0,.y=0,.z=0};
    uint8_t tdoa_count = 0;
    while (stateEstimatorHasTDOAPacket(&tdoa) && tdoa_count < UWB_QUEUE_LENGTH){
        projectOnTDOA(&tdoa, &loc_prediction, &projection[tdoa_count]);
        tdoa_count++;
        sum_for_avg.x += projection[tdoa_count].x;
        sum_for_avg.y += projection[tdoa_count].y;
        sum_for_avg.z += projection[tdoa_count].z;
    }
    // No outlier rejection yet
    loc_measurement.x = sum_for_avg.x/tdoa_count;
    loc_measurement.y = sum_for_avg.x/tdoa_count;
    loc_measurement.z = sum_for_avg.x/tdoa_count;

    // estimate error
    delta.x[k] = loc_measurement.x - loc_prediction.x;
    delta.y[k] = loc_measurement.y - loc_prediction.y;
    
    // update least squares sums
    lsx_sumXY += (- POS_UPDATE_DT*LSx_sumY) + currentWindowSize*POS_UPDATE_DT*delta.x[k];
    lsx_sumY += delta.x[k];
    lsy_sumXY += (- POS_UPDATE_DT*LSy_sumY) + currentWindowSize*POS_UPDATE_DT*delta.y[k];
    lsy_sumY += delta.y[k];

    errorEst_x = (ls_sumX2 * lsx_sumY - lsx_sumXY * ls_sumX) / ls_denom;
    errorEst_vx = (currentWindowSize * lsx_sumXY - ls_sumX * lsx_sumY) / ls_denom;
    errorEst_y = (ls_sumX2 * lsy_sumY - lsy_sumXY * ls_sumX) / ls_denom;
    errorEst_vy = (currentWindowSize * lsy_sumXY - ls_sumX * lsy_sumY) / ls_denom;

    // correction of prediction
    state->position.x = loc_prediction.x + errorEst_x + errorEst_vx * k * POS_UPDATE_DT;
    state->position.y = loc_prediction.y + errorEst_y + errorEst_vy * k * POS_UPDATE_DT;
    state->velocity.x = vel_prediction.x + errorEst_vx;
    state->velocity.y = vel_prediction.y + errorEst_vy;
 
  }
}


void movingHorizonShiftWindow(float *window[], int windowSize){
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
    float x = origin.x, y = origin.y, z = origin.z;
        
    float dx1 = x - x1;
    float dy1 = y - y1;
    float dz1 = z - z1;

    float dy0 = y - y0;
    float dx0 = x - x0;
    float dz0 = z - z0;

    float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
    float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2))
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

        d1 = sqrtf(powf(dx1+t*grad.x)+powf(dy1+t*grad.y)+powf(dz1+t*grad.z));
        d0 = sqrtf(powf(dx0+t*grad.x)+powf(dy0+t*grad.y)+powf(dz0+t*grad.z));
         
        error = d1 - d0 - tdoa_measured;
    }
        
    projection.x = x + t*grad.x;
    projection.y = y + t*grad.y;
    projection.z = z + t*grad.z;
}