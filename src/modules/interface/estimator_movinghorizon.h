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
#ifndef __ESTIMATOR_MOVINGHORIZON_H__
#define __ESTIMATOR_MOVINGHORIZON_H__

#include "stabilizer_types.h"

void estimatorMovingHorizonInit(void);
bool estimatorMovingHorizonTest(void);
void estimatorMovingHorizon(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick);

bool estimatorMovingHorizonEnqueueTDOA(const tdoaMeasurement_t *uwb)


#endif