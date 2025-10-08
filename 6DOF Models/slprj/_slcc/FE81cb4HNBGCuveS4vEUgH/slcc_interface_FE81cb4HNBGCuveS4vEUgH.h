#include "customcode_FE81cb4HNBGCuveS4vEUgH.h"
#ifdef __cplusplus
extern "C" {
#endif


/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
DLL_EXPORT_CC extern const char_T *get_dll_checksum_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern real32_T pressureToAltitude_FE81cb4HNBGCuveS4vEUgH(real32_T in0);
DLL_EXPORT_CC extern Quaternion QuatStep_FE81cb4HNBGCuveS4vEUgH(Quaternion q, Vector w, real32_T dt);
DLL_EXPORT_CC extern real32_T vdot_FE81cb4HNBGCuveS4vEUgH(Vector v1, Vector v2);
DLL_EXPORT_CC extern Vector vcross_FE81cb4HNBGCuveS4vEUgH(Vector v1, Vector v2);
DLL_EXPORT_CC extern Vector vscale_FE81cb4HNBGCuveS4vEUgH(Vector v, real32_T k);
DLL_EXPORT_CC extern Vector vadd_FE81cb4HNBGCuveS4vEUgH(Vector v1, Vector v2);
DLL_EXPORT_CC extern Vector QuatRot_FE81cb4HNBGCuveS4vEUgH(Vector v, Quaternion q);
DLL_EXPORT_CC extern Vector TrapInt_FE81cb4HNBGCuveS4vEUgH(Vector vint, Vector v, Vector vprev, real32_T dt);
DLL_EXPORT_CC extern void UpdatePose_FE81cb4HNBGCuveS4vEUgH(Vector v_a, Vector *v_ac, Vector *v_v, Vector *v_d, Vector v_w, Quaternion *q, real32_T dt);
DLL_EXPORT_CC extern Vector zeroVec_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern Vector newVec_FE81cb4HNBGCuveS4vEUgH(real32_T x, real32_T y, real32_T z);
DLL_EXPORT_CC extern Vector createVectorFromStruct_FE81cb4HNBGCuveS4vEUgH(void *ptr);
DLL_EXPORT_CC extern real32_T vnorm_FE81cb4HNBGCuveS4vEUgH(Vector v);
DLL_EXPORT_CC extern SensorData sensorFrame2SensorData_FE81cb4HNBGCuveS4vEUgH(SensorFrame frame);
DLL_EXPORT_CC extern void fp_init_FE81cb4HNBGCuveS4vEUgH(FlightPhase *s_flight_phase, StateEst *current_state, Vector *imu_up, Vector *high_g_up, AverageBuffer *acc_buffer, AverageBuffer *baro_buffer);
DLL_EXPORT_CC extern void fp_update_FE81cb4HNBGCuveS4vEUgH(SensorFrame *data, GPS_Fix_TypeDef *gps, FlightPhase *s_flight_phase, StateEst *current_state, Vector *imu_up, Vector *high_g_up, AverageBuffer *acc_buffer, AverageBuffer *baro_buffer);
DLL_EXPORT_CC extern StateEst zeroState_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern void push_to_average_FE81cb4HNBGCuveS4vEUgH(real32_T new_value, AverageBuffer *buffer);
DLL_EXPORT_CC extern FlightPhase simulink_test_input_FE81cb4HNBGCuveS4vEUgH(int32_T initialize, real32_T a_x, real32_T a_y, real32_T a_z, real32_T pressure, real32_T time, real32_T gps_pos, real32_T gps_vel);
DLL_EXPORT_CC extern int32_T get_deploy_main_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern int32_T get_deploy_drogue_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern real32_T get_velNED_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern real32_T get_accNED_FE81cb4HNBGCuveS4vEUgH(void);
DLL_EXPORT_CC extern real32_T get_posNED_FE81cb4HNBGCuveS4vEUgH(void);

/* Function Definitions */
DLL_EXPORT_CC const uint8_T *get_checksum_source_info(int32_T *size);
#ifdef __cplusplus
}
#endif

