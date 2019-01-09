//
// Created by skkim9226 on 2017-12-23.
//

#ifndef PTGPU_NATIVE_LIB_H
#define PTGPU_NATIVE_LIB_H

#ifdef __ANDROID__
#include <android/asset_manager.h>
#endif
#include "../../assets/sdcard/include/camera.h"

#if defined(__ANDROID__)
#include <android/log.h>

#define LOG_TAG __FILE__
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)
#else
#define LOG_TAG __FILE__
#define  LOGI(...)  printf(__VA_ARGS__)
#define  LOGE(...)  printf(__VA_ARGS__)
#endif

#define MAX_FN 255
#define DIFF_LEFTRIGHTEYE 10.0f
//#define CURRENT_SAMPLE

extern int useGPU;
extern int forceWorkSize;
extern char kernelFileName[];
extern char bvhFileName[];
extern int pixelCount;
extern Camera cameraLeft, cameraRight;
#ifdef __ANDROID__
extern AAssetManager *mgr;
extern char *strResPath;
#endif

extern void AddWallLight();
extern void SetUpOpenCL();
extern void BuildBVH();
extern void BuildKDtree();
extern unsigned int *DrawFrame();
extern unsigned int *DrawFrameVR(short bleft);

#endif //PTGPU_NATIVE_LIB_H
