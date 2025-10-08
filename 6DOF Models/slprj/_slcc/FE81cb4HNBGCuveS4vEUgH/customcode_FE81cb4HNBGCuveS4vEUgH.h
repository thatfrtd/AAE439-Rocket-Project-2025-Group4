#ifndef __customcode_FE81cb4HNBGCuveS4vEUgH_h__
#define __customcode_FE81cb4HNBGCuveS4vEUgH_h__

/* Include files */
#include "mex.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tmwtypes.h"


/* Helper definitions for DLL support */
#if defined _WIN32 
  #define DLL_EXPORT_CC __declspec(dllexport)
#else
  #if __GNUC__ >= 4
    #define DLL_EXPORT_CC __attribute__ ((visibility ("default")))
  #else
    #define DLL_EXPORT_CC
  #endif
#endif
/* Custom Code from Simulation Target dialog */
#include "pressure_altitude.h"
#include "flight_estimation.h"
#include "data.h"
#include "max_m10s.h"
#include "i2c/i2c.h"

/* Function Declarations */
#ifdef __cplusplus
extern "C" {
#endif
#define customcode_FE81cb4HNBGCuveS4vEUgH_initializer()

#define customcode_FE81cb4HNBGCuveS4vEUgH_terminator()
#ifdef __cplusplus
}
#endif

#endif

