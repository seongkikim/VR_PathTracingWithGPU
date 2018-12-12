#include <jni.h>
#include <string>
#include <string.h>
#include <android/asset_manager.h>
#include <android/asset_manager_jni.h>

#include "include/displayfunc.h"
#include "include/native-lib.h"
#include "../assets/sdcard/include/geom.h"

#ifdef __cplusplus
//extern "C" {
#endif

AAssetManager   *mgr;
char *strResPath;

extern "C" JNIEXPORT jfloatArray
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_initSmallPtGPU(JNIEnv* env, jobject /* this */, jint u, jint f, jstring k, jint w, jint h, jstring s, jstring r, jobject assetManager, jboolean bvr)
{
    bool walllight = false;
    srand(time(NULL));

    useGPU = u;
    forceWorkSize = f;
#ifndef EXP_KERNEL
    strcpy(kernelFileName, env->GetStringUTFChars(k, NULL));
#endif

    width = w;
    height = h;

    mgr = AAssetManager_fromJava(env, assetManager);

    char *strFilePath = (char *)env->GetStringUTFChars(s, NULL);
    strResPath = (char *)env->GetStringUTFChars(r, NULL);
    char *strFullPath = (char *)malloc(strlen(strFilePath) + strlen(strResPath) + 1);

    strcpy(strFullPath, strResPath);
    strcat(strFullPath, "/");
    strcat(strFullPath, strFilePath);

    if (Read(strFullPath, &walllight, bvr)) {
        if (walllight) AddWallLight();

#if (ACCELSTR == 0)
        SetUpOpenCL();
#elif (ACCELSTR == 1)
        SetUpOpenCL();
		BuildBVH();
#elif (ACCELSTR == 2)
		BuildKDtree();
		SetUpOpenCL();
#endif
		if (bvr) ReInitVR(0);
        else ReInit(0);
    }

    jfloat cameraOrigTarg[6];
    if (bvr) {
        cameraOrigTarg[0] = (cameraLeft.orig.s[0] + cameraRight.orig.s[0]) / 2.0f, cameraOrigTarg[1] = (cameraLeft.orig.s[1] + cameraRight.orig.s[1]) / 2.0f, cameraOrigTarg[2] = (cameraLeft.orig.s[2] + cameraRight.orig.s[2]) / 2.0f;
        cameraOrigTarg[3] = (cameraLeft.target.s[0] + cameraRight.target.s[0]) / 2.0f, cameraOrigTarg[4] = (cameraLeft.target.s[1] + cameraRight.target.s[1]) / 2.0f, cameraOrigTarg[5] = (cameraLeft.target.s[2] + cameraRight.target.s[2]) / 2.0f;
    }
    else {
        cameraOrigTarg[0] = cameraLeft.orig.s[0], cameraOrigTarg[1] = cameraLeft.orig.s[1], cameraOrigTarg[2] = cameraLeft.orig.s[2];
        cameraOrigTarg[3] = cameraLeft.target.s[0], cameraOrigTarg[4] = cameraLeft.target.s[1], cameraOrigTarg[5] = cameraLeft.target.s[2];
    }

    jfloatArray jf_array = env->NewFloatArray(6);
    env->SetFloatArrayRegion(jf_array, 0, 6, (const jfloat *)cameraOrigTarg);

    return jf_array;
}

extern "C" JNIEXPORT void
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_reinitCamera(JNIEnv* env, jobject /* this */, jfloat origx, jfloat origy, jfloat origz, jfloat targx, jfloat targy, jfloat targz)
{
    cameraLeft.orig.s[0] = origx,cameraLeft.orig.s[1] = origy, cameraLeft.orig.s[2] = origz;
    cameraLeft.target.s[0] = targx, cameraLeft.target.s[1] = targy, cameraLeft.target.s[2] = targz;

    ReInit(0);
}

extern "C" JNIEXPORT void
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_reinitCameraVR(JNIEnv* env, jobject /* this */, jfloat origx, jfloat origy, jfloat origz, jfloat targx, jfloat targy, jfloat targz)
{
    cameraLeft.orig.s[0] = origx - DIFF_LEFTRIGHTEYE,cameraLeft.orig.s[1] = origy, cameraLeft.orig.s[2] = origz;
    cameraLeft.target.s[0] = targx, cameraLeft.target.s[1] = targy, cameraLeft.target.s[2] = targz;

    cameraRight.orig.s[0] = origx + DIFF_LEFTRIGHTEYE,cameraRight.orig.s[1] = origy, cameraRight.orig.s[2] = origz;
    cameraRight.target.s[0] = targx, cameraRight.target.s[1] = targy, cameraRight.target.s[2] = targz;

    ReInitVR(0);
}

extern "C" JNIEXPORT void
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_touchFunc(JNIEnv* env, jobject /* this */, int deltax, int deltay, jboolean bvr)
{
    touchFunc(deltax, deltay, bvr);
}

extern "C" JNIEXPORT jintArray
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_updateRendering(JNIEnv* env, jobject /* this */) {
    unsigned int *pixels = DrawFrame();

    // Get JNI Env for all function calls
    jintArray ji_array = env->NewIntArray(pixelCount);
    env->SetIntArrayRegion(ji_array, 0, pixelCount, (const jint *)pixels);

    return ji_array;
}

extern "C" JNIEXPORT jintArray
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_updateRenderingVR(JNIEnv* env, jobject /* this */, jboolean bleft) {
    unsigned int *pixels = DrawFrameVR(bleft);

    // Get JNI Env for all function calls
    jintArray ji_array = env->NewIntArray(pixelCount);
    env->SetIntArrayRegion(ji_array, 0, pixelCount, (const jint *)pixels);

    return ji_array;
}

extern "C" JNIEXPORT void
JNICALL Java_gamemobile_kmu_ac_kr_vrapp3_VRApp3Renderer_finishRendering(JNIEnv* env, jobject /* this */) {
    // Get JNI Env for all function calls
}

#ifdef __cplusplus
//}
#endif