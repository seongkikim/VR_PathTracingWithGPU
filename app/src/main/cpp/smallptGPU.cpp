/*
Copyright (c) 2009 David Bucciarelli (davibu@interfree.it)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
 * Based on smallpt, a Path Tracer by Kevin Beason, 2008
 * Modified by David Bucciarelli to show the output via OpenGL/GLUT, ported
 * to C, work with float, fixed RR, ported to OpenCL, etc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

// Jens's patch for MacOS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include "include/CL/cl.h"
#endif

#include "include/CLBVH.h"
#include "include/KDTree.h"
#include "../assets/sdcard/include/camera.h"
#include "include/scene.h"
#include "include/displayfunc.h"
#include "../assets/sdcard/include/geom.h"
#include "include/geomfunc.h"
#include "include/native-lib.h"

#ifndef __ANDROID__
#include <GL/glut.h>
#define M_PI       3.14159265358979323846 

bool Read(char *fileName, bool *walllight);
#else 
#include <android/asset_manager.h>
#include <pthread.h>
#endif

#ifdef EXP_KERNEL
Ray *ray;
Vec *throughput;
Result *result;
char *specularBounce, *terminated;

static cl_mem rayBuffer, throughputBuffer, specularBounceBuffer, terminatedBuffer, resultBuffer;
static cl_kernel kernelGen, kernelRadiance, kernelFill;
char kernelFileName[MAX_FN] = "kernels/rendering_kernel_exp.cl";
#else
static cl_kernel kernel;
char kernelFileName[MAX_FN] = "kernels/rendering_kernel.cl";
#endif

#ifdef CPU_PARTRENDERING
static cl_kernel kernelBox;
#ifdef EXP_KERNEL
static cl_kernel kernelGenBox, kernelRadianceBox, kernelFillBox;

char **specularBounce_boxes, **terminated_boxes;
unsigned int **pixels_boxes, **seeds_boxes;
Vec **colors_boxes, **throughputs_boxes;
Ray **rays_boxes;
Result **results_boxes;

cl_mem *colorboxBuffer, *pixelboxBuffer, *rayboxBuffer, *seedsboxBuffer, *throughputboxBuffer,  *specularBounceboxBuffer, *terminatedboxBuffer, *resultboxBuffer;
char *cpuProcessed;
#endif
#endif

#if (ACCELSTR == 1)
static cl_mem btnBuffer;
static cl_mem btlBuffer;
BVHNodeGPU *btn, *btl;
static cl_kernel kernelRad, kernelBvh, kernelOpt;
char bvhFileName[MAX_FN] = "kernels/BVH.cl";
#elif (ACCELSTR == 2)
static cl_mem kngBuffer;
static cl_mem knBuffer;
KDNodeGPU *pkngbuf;
int *pknbuf; 
short kngCnt;
short knCnt;
#endif

#define MAX_STYPE 255
#define MAX_INCLUDE 255
#define MAX_ERROR 255
#define MAX_LOG 255
//#define PTX_ERROR

/* OpenCL variables */
static cl_context context;
static cl_mem colorBufferLeft, colorBufferRight, pixelBufferLeft, pixelBufferRight, seedBuffer, shapeBuffer, fhiBuffer, cameraBufferLeft, cameraBufferRight;
static cl_command_queue commandQueue;
static cl_program program;

static Vec *colorsLeft, *colorsRight;
static unsigned int *seeds;
static short currentSampleLeft, currentSampleRight = 0;

/* Options */
int useGPU = 1;
int forceWorkSize = 0;

unsigned int workGroupSize = 1;
short shapeCnt = 0, lightCnt = 0;
int pixelCount;
Camera cameraLeft, cameraRight;
Shape *shapes;
FirstHitInfo *fhi;

#if (ACCELSTR == 1)
void BuildBVH();
#elif (ACCELSTR == 2)
void BuildKDtree();
#endif

#define clErrchk(ans) { clAssert((ans), __FILE__, __LINE__); }

inline void clAssert(cl_int code, const char *file, int line)
{
	if (code != CL_SUCCESS)
	{
		LOGI("Error: %d in %s (%d)\n", code, file, line);
	}
}

void FreeBuffers() {
#ifdef EXP_KERNEL
	if (resultBuffer)
		clErrchk(clReleaseMemObject(resultBuffer));

	if (terminatedBuffer)
		clErrchk(clReleaseMemObject(terminatedBuffer));

	if (specularBounceBuffer)
		clErrchk(clReleaseMemObject(specularBounceBuffer));

	if (throughputBuffer)
		clErrchk(clReleaseMemObject(throughputBuffer));

	if (rayBuffer)
		clErrchk(clReleaseMemObject(rayBuffer));

	if (fhiBuffer)
		clErrchk(clReleaseMemObject(fhiBuffer));
#ifdef CPU_PARTRENDERING
    int cntBuffers = (width * height) / (BWIDTH * BHEIGHT);

    for(int i = 0; i < cntBuffers; i++)
    {
        if (resultboxBuffer[i])
            clErrchk(clReleaseMemObject(resultboxBuffer[i]));

        if (terminatedboxBuffer[i])
            clErrchk(clReleaseMemObject(terminatedboxBuffer[i]));

        if (specularBounceboxBuffer[i])
            clErrchk(clReleaseMemObject(specularBounceboxBuffer[i]));

        if (throughputboxBuffer[i])
            clErrchk(clReleaseMemObject(throughputboxBuffer[i]));

        if (seedsboxBuffer[i])
            clErrchk(clReleaseMemObject(seedsboxBuffer[i]));

        if (rayboxBuffer[i])
            clErrchk(clReleaseMemObject(rayboxBuffer[i]));

        if (pixelboxBuffer[i])
            clErrchk(clReleaseMemObject(pixelboxBuffer[i]));

        if (colorboxBuffer[i])
            clErrchk(clReleaseMemObject(colorboxBuffer[i]));
    }

    if (resultboxBuffer) free(resultboxBuffer);
    if (terminatedboxBuffer) free(terminatedboxBuffer);
    if (specularBounceboxBuffer) free(specularBounceboxBuffer);
    if (throughputboxBuffer) free(throughputboxBuffer);
    if (seedsboxBuffer) free(seedsboxBuffer);
    if (rayboxBuffer) free(rayboxBuffer);
    if (pixelboxBuffer) free(pixelboxBuffer);
    if (colorboxBuffer) free(colorboxBuffer);

    free(cpuProcessed);
#endif
#endif
	if (seedBuffer)
	{
		clErrchk(clReleaseMemObject(seedBuffer));
	}	
	if (pixelBufferRight)
	{
		clErrchk(clReleaseMemObject(pixelBufferRight));
	}
    if (pixelBufferLeft)
    {
        clErrchk(clReleaseMemObject(pixelBufferLeft));
    }
	if (colorBufferRight)
	{
		clErrchk(clReleaseMemObject(colorBufferRight));
	}
    if (colorBufferLeft)
    {
        clErrchk(clReleaseMemObject(colorBufferLeft));
    }
#ifdef EXP_KERNEL
	if (result) free(result);
	if (terminated) free(terminated);
	if (specularBounce) free(specularBounce);
	if (throughput) free(throughput);
	if (ray) free(ray);
	if (fhi) free(fhi);

#ifdef CPU_PARTRENDERING
    for(int i = 0; i < cntBuffers; i++) {
        if (results_boxes[i]) free(results_boxes[i]);
        if (terminated_boxes[i]) free(terminated_boxes[i]);
        if (specularBounce_boxes[i]) free(specularBounce_boxes[i]);
        if (throughputs_boxes[i]) free(throughputs_boxes[i]);
        if (seeds_boxes[i]) free(seeds_boxes[i]);
        if (rays_boxes[i]) free(rays_boxes[i]);
        if (pixels_boxes[i]) free(pixels_boxes[i]);
        if (colors_boxes[i]) free(colors_boxes[i]);
    }

    if (results_boxes) free(results_boxes);
    if (terminated_boxes) free(terminated_boxes);
    if (specularBounce_boxes) free(specularBounce_boxes);
    if (throughputs_boxes) free(throughputs_boxes);
    if (seeds_boxes) free(seeds_boxes);
    if (rays_boxes) free(rays_boxes);
    if (pixels_boxes) free(pixels_boxes);
    if (colors_boxes) free(colors_boxes);
#endif
#endif
	if (seeds) free(seeds);
	if (colorsRight) free(colorsRight);
    if (colorsLeft) free(colorsLeft);
    if (pixelsRight) free(pixelsRight);
	if (pixelsLeft) free(pixelsLeft);
}

void AllocateBuffers() {
	pixelCount = width * height;
	int i;
    
	colorsLeft = (Vec *)malloc(sizeof(Vec) * pixelCount);
    colorsRight = (Vec *)malloc(sizeof(Vec) * pixelCount);

    seeds = (unsigned int *)malloc(sizeof(unsigned int) * pixelCount * 2);

	for (i = 0; i < pixelCount * 2; i++) {
		seeds[i] = rand();
		if (seeds[i] < 2) seeds[i] = 2;
	}

	pixelsLeft = (unsigned int *)malloc(sizeof(unsigned int) * pixelCount);
	pixelsRight = (unsigned int *)malloc(sizeof(unsigned int) * pixelCount);

	cl_int status;
	cl_uint sizeBytes = sizeof(Vec) * width * height;

    colorBufferLeft = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeBytes, NULL, &status);
	clErrchk(status);

    colorBufferRight = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeBytes, NULL, &status);
    clErrchk(status);

	sizeBytes = sizeof(unsigned char[4]) * width * height;
    pixelBufferLeft = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeBytes, NULL, &status);
	clErrchk(status);

    pixelBufferRight = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeBytes, NULL, &status);
    clErrchk(status);

    sizeBytes = sizeof(unsigned int) * width * height * 2;
	seedBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeBytes, NULL, &status);
	clErrchk(status);

	status = clEnqueueWriteBuffer(commandQueue, colorBufferLeft, CL_TRUE, 0, sizeBytes, colorsLeft, 0, NULL, NULL);
	clErrchk(status);

    status = clEnqueueWriteBuffer(commandQueue, colorBufferRight, CL_TRUE, 0, sizeBytes, colorsRight, 0, NULL, NULL);
    clErrchk(status);

	status = clEnqueueWriteBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeBytes, seeds, 0, NULL, NULL);
	clErrchk(status);

#ifdef EXP_KERNEL
	ray = (Ray *)malloc(sizeof(Ray) * width * height);

	rayBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Ray) *  width * height, NULL, &status);
	clErrchk(status);

	fhi = (FirstHitInfo *)malloc(sizeof(FirstHitInfo) * width * height);

	fhiBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(FirstHitInfo) * width * height, NULL, &status);
	clErrchk(status);

	throughput = (Vec *)malloc(sizeof(Vec) * width * height);

	throughputBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Vec) *  width * height, NULL, &status);
	clErrchk(status);

	specularBounce = (char *)malloc(sizeof(char) * width * height);

	specularBounceBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char) *  width * height, NULL, &status);
	clErrchk(status);

	terminated = (char *)malloc(sizeof(char) * width * height);

	terminatedBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char) *  width * height, NULL, &status);
	clErrchk(status);

	result = (Result *)malloc(sizeof(Result) * width * height);

	resultBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Result) *  width * height, NULL, &status);
	clErrchk(status);

#ifdef CPU_PARTRENDERING
    int cntBuffers = (width * height) / (BWIDTH * BHEIGHT);
	
	cpuProcessed = (char *)malloc(sizeof(char) * cntBuffers);
    colors_boxes = (Vec **) malloc(sizeof(Vec *) * cntBuffers);
    pixels_boxes = (unsigned int **) malloc(sizeof(unsigned int *) * cntBuffers);
    rays_boxes = (Ray **) malloc(sizeof(Ray *) * cntBuffers);
    seeds_boxes =  (unsigned int **) malloc(sizeof(unsigned int *) * cntBuffers);
    throughputs_boxes =  (Vec **) malloc(sizeof(Vec *) * cntBuffers);
    specularBounce_boxes =  (char **) malloc(sizeof(char *) * cntBuffers);
    terminated_boxes =  (char **) malloc(sizeof(char *) * cntBuffers);
    results_boxes =  (Result **) malloc(sizeof(Result *) * cntBuffers);

    colorboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    pixelboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    rayboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    seedsboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    throughputboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    specularBounceboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    terminatedboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);
    resultboxBuffer = (cl_mem *) malloc(sizeof(cl_mem) * cntBuffers);	
	    
	for(int i = 0; i < cntBuffers; i++) {
        cpuProcessed[i] = -1;

        int wratio = width / BWIDTH, x = (i % wratio) * BWIDTH, y = (i / wratio) * BHEIGHT;

        pixels_boxes[i] = (unsigned int *)malloc(sizeof(unsigned int) * BWIDTH * BHEIGHT);

        pixelboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);

        rays_boxes[i] = (Ray *)malloc(sizeof(Ray) * BWIDTH * BHEIGHT);

        rayboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Ray) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);

        seeds_boxes[i] = (unsigned int *)malloc(sizeof(unsigned int) * BWIDTH * BHEIGHT * 2);
        colors_boxes[i] = (Vec *)malloc(sizeof(Vec) * BWIDTH * BHEIGHT);

        for(int j = 0; j < BHEIGHT; j++) {
            memcpy((char *) seeds_boxes[i] + j * sizeof(unsigned int) * BWIDTH * 2,
                   (char *) seeds + x * sizeof(unsigned int) * 2 + (y + j) * sizeof(unsigned int) * width * 2,
                   sizeof(unsigned int) * BWIDTH * 2);

            memcpy((char *) colors_boxes[i] + j * sizeof(Vec) * BWIDTH,
                   (char *) colors + x * sizeof(Vec) + (y + j) * width * sizeof(Vec),
                   sizeof(Vec) * BWIDTH);
        }

        seedsboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int) * BWIDTH * BHEIGHT * 2, NULL, &status);
        clErrchk(status);

        clErrchk(clEnqueueWriteBuffer(commandQueue, seedsboxBuffer[i], CL_TRUE, 0, sizeof(unsigned int) * BWIDTH * BHEIGHT * 2, seeds_boxes[i], 0, NULL, NULL));

        colorboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Vec) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);

        clErrchk(clEnqueueWriteBuffer(commandQueue, colorboxBuffer[i], CL_TRUE, 0, sizeof(Vec) * BWIDTH * BHEIGHT, colors_boxes[i], 0, NULL, NULL));

        throughputs_boxes[i] = (Vec *) malloc(sizeof(Vec) * BWIDTH * BHEIGHT);

        throughputboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Vec) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);

        specularBounce_boxes[i] = (char *) malloc(sizeof(char) * BWIDTH * BHEIGHT);

        specularBounceboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);

        terminated_boxes[i] = (char *) malloc(sizeof(char) * BWIDTH * BHEIGHT);

        terminatedboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);

        results_boxes[i] = (Result *) malloc(sizeof(Result) * BWIDTH * BHEIGHT);

        resultboxBuffer[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Result) * BWIDTH * BHEIGHT, NULL, &status);
        clErrchk(status);
    }
#endif
#endif
}

char *ReadSources(const char *fileName) {
	AAsset *as = AAssetManager_open(mgr, fileName, AASSET_MODE_UNKNOWN);
	if (NULL == as) {
		LOGE("Failed to open asset '%s'\n", fileName);
		return NULL;
	}

	off_t size = AAsset_getLength(as);
	char *src = (char *)malloc(sizeof(char) * size + 1);
	if (!src) {
		LOGE("Failed to allocate memory for file '%s'\n", fileName);
		return NULL;
	}

	LOGI("Reading file '%s' (size %ld bytes)\n", fileName, size);
	int res = AAsset_read (as, src, size);
	if (res != sizeof(char) * size) {
		LOGE("Failed to read file '%s' (read %lu)\n Content: %s\n", fileName, res, src);
		return NULL;
	}

	AAsset_close(as);

	src[res] = '\0'; /* NULL terminated */

	return src;
}

void SetUpOpenCL() {
	cl_device_type dType;

	if (useGPU) dType = CL_DEVICE_TYPE_GPU;
	else dType = CL_DEVICE_TYPE_CPU;

	// Select the platform
    cl_uint numPlatforms;
	cl_platform_id platform = NULL;

	clErrchk(clGetPlatformIDs(0, NULL, &numPlatforms));

	if (numPlatforms > 0) {
		cl_platform_id *platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id) * numPlatforms);

		clErrchk(clGetPlatformIDs(numPlatforms, platforms, NULL));

		unsigned int i;
		for (i = 0; i < numPlatforms; ++i) {
			char pbuf[100];

			clErrchk(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(pbuf), pbuf, NULL));
			clErrchk(clGetPlatformIDs(numPlatforms, platforms, NULL));

			LOGI("OpenCL Platform %d: %s\n", i, pbuf);
		}

		platform = platforms[0];
		free(platforms);
	}

	// Select the device
	cl_device_id devices[32];
	cl_uint deviceCount;

	clErrchk(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 32, devices, &deviceCount));

	int deviceFound = 0;
	cl_device_id selectedDevice;
	unsigned int i;

	for (i = 0; i < deviceCount; ++i) {
		cl_device_type type = 0;

		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL));

		char stype[MAX_STYPE];
		switch (type) {
			case CL_DEVICE_TYPE_ALL:
				strcpy(stype, "TYPE_ALL");
				break;
			case CL_DEVICE_TYPE_DEFAULT:
				strcpy(stype, "TYPE_DEFAULT");
				break;
			case CL_DEVICE_TYPE_CPU:
				strcpy(stype, "TYPE_CPU");
				if (!useGPU && !deviceFound) {
					selectedDevice = devices[i];
					deviceFound = 1;
				}
				break;
			case CL_DEVICE_TYPE_GPU:
				strcpy(stype, "TYPE_GPU");
				if (useGPU && !deviceFound) {
					selectedDevice = devices[i];
					deviceFound = 1;
				}
				break;
			default:
				strcpy(stype, "TYPE_UNKNOWN");
				break;
		}

		LOGI("OpenCL Device %d: Type = %s\n", i, stype);

		char buf[256];
		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(char[256]), &buf, NULL));

		LOGI("OpenCL Device %d: Name = %s\n", i, buf);

		cl_uint units = 0;
		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &units, NULL));

		LOGI("OpenCL Device %d: Compute units = %u\n", i, units);

		size_t gsize = 0;
		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &gsize, NULL));

		LOGI("OpenCL Device %d: Max. work group size = %d\n", i, (unsigned int)gsize);
	}

	if (!deviceFound) {
		LOGE("Unable to select an appropriate device\n");
		exit(0);
		return ;
	}

	// Create the context
	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platform, 0 };
	cl_context_properties *cprops = (NULL == platform) ? NULL : cps;
	cl_int status;

	context = clCreateContext(cprops, 1, &selectedDevice, NULL, NULL, &status);
	clErrchk(status);

    /* Get the device list data */
	size_t deviceListSize;
	clErrchk(clGetContextInfo(context, CL_CONTEXT_DEVICES, 32, devices, &deviceListSize));

	/* Print devices list */
	for (i = 0; i < deviceListSize / sizeof(cl_device_id); ++i) {
		cl_device_type type = 0;
		char stype[MAX_STYPE];

		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL));

		switch (type) {
			case CL_DEVICE_TYPE_ALL:
				strcpy(stype, "TYPE_ALL");
				break;
			case CL_DEVICE_TYPE_DEFAULT:
				strcpy(stype, "TYPE_DEFAULT");
				break;
			case CL_DEVICE_TYPE_CPU:
				strcpy(stype, "TYPE_CPU");
				break;
			case CL_DEVICE_TYPE_GPU:
				strcpy(stype, "TYPE_GPU");
				break;
			default:
				strcpy(stype, "TYPE_UNKNOWN");
				break;
		}

		LOGI("[SELECTED] OpenCL Device %d: Type = %s\n", i, stype);

		char buf[256];
		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(char[256]), &buf, NULL));

		LOGI("[SELECTED] OpenCL Device %d: Name = %s\n", i, buf);

		cl_uint units = 0;
		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &units, NULL));

		LOGI("[SELECTED] OpenCL Device %d: Compute units = %u\n", i, units);

		size_t gsize = 0;
		clErrchk(clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &gsize, NULL));

		LOGI("[SELECTED] OpenCL Device %d: Max. work group size = %d\n", i, (unsigned int)gsize);
	}

	cl_command_queue_properties prop = 0;

	commandQueue = clCreateCommandQueue(context, devices[0], prop, &status);
	clErrchk(status);

	/*------------------------------------------------------------------------*/
#ifdef __ANDROID__
	shapeBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Shape) * shapeCnt, NULL, &status);
#else
    shapeBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Shape) * shapeCnt, NULL, &status);
#endif
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, shapeBuffer, CL_TRUE, 0, sizeof(Shape) * shapeCnt, shapes, 0, NULL, NULL));

#if (ACCELSTR == 1)
	/*------------------------------------------------------------------------*/
	btnBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(BVHNodeGPU) * (shapeCnt-1), NULL, &status);
	clErrchk(status);

	btlBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(BVHNodeGPU) * shapeCnt, NULL, &status);
	clErrchk(status);
#elif (ACCELSTR == 2)
	/*------------------------------------------------------------------------*/
	kngBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(KDNodeGPU) * (kngCnt), NULL, &status);
	clErrchk(status);

	knBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int) * knCnt, NULL, &status);
	clErrchk(status);
#endif
	cameraBufferLeft = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Camera), NULL, &status);
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, cameraBufferLeft, CL_TRUE, 0, sizeof(Camera), &cameraLeft, 0, NULL, NULL));

	cameraBufferRight = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Camera), NULL, &status);
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, cameraBufferRight, CL_TRUE, 0, sizeof(Camera), &cameraRight, 0, NULL, NULL));

	AllocateBuffers();

	/*------------------------------------------------------------------------*/
	/* Create the kernel program */
	const char *sources = ReadSources(kernelFileName);
	program = clCreateProgramWithSource(context, 1, &sources, NULL, &status);
	clErrchk(status);

	char strInclude[MAX_INCLUDE];
	strcpy(strInclude, "-DGPU_KERNEL -I. -I");
	strcat(strInclude, strResPath);
	strcat(strInclude, "/include");

	status = clBuildProgram(program, 1, devices, strInclude, NULL, NULL);
	clErrchk(status);

	if (status != CL_SUCCESS) {
		LOGE("Failed to build OpenCL kernel: %d\n", status);

        size_t retValSize;
		clErrchk(clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &retValSize));

        char *buildLog = (char *)malloc(retValSize + 1);
		clErrchk(clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, retValSize, buildLog, NULL));

        buildLog[retValSize] = '\0';
#ifdef PTX_ERROR
		// Query binary (PTX file) size
		size_t bin_sz;
		status = clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL);

		// Read binary (PTX file) to memory bufferf
		unsigned char *bin = (unsigned char *)malloc(bin_sz);
		status = clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL);

		char ptxFileName[255];
		strcpy(ptxFileName, kernelFileName);
		strcat(ptxFileName, ".ptx");

		// Save PTX to add_vectors_ocl.ptx
		FILE *fp = fopen(ptxFileName, "wb");
		fwrite(bin, sizeof(char), bin_sz, fp);
		fclose(fp);

		free(bin);
#endif
		LOGE("OpenCL Programm Build Log: %s\n", buildLog);

		char strError[MAX_ERROR];
		strcpy(strError, strResPath);
		strcat(strError, "/error_renderingkernel.txt");

		FILE *fp = fopen(strError, "wt");
		fwrite(buildLog, sizeof(char), retValSize + 1, fp);
		fclose(fp);

		free(buildLog);
		return ;
    }
#ifdef EXP_KERNEL
	kernelGen = clCreateKernel(program, "GenerateCameraRay_exp", &status);
	clErrchk(status);

	kernelRadiance = clCreateKernel(program, "RadiancePathTracing_exp", &status);
	clErrchk(status);

	kernelFill = clCreateKernel(program, "FillPixel_exp", &status);
	clErrchk(status);
#ifdef CPU_PARTRENDERING
    kernelGenBox = clCreateKernel(program, "GenerateCameraRay_expbox", &status);
    clErrchk(status);

    kernelRadianceBox = clCreateKernel(program, "RadiancePathTracing_expbox", &status);
    clErrchk(status);

    kernelFillBox = clCreateKernel(program, "FillPixel_expbox", &status);
    clErrchk(status);
#endif
#else
	kernel = clCreateKernel(program, "RadianceGPU", &status);
	clErrchk(status);
#endif

#if (ACCELSTR == 1)
	/* Create the kernel program */
	const char *sourcesBvh = ReadSources(bvhFileName);
	program = clCreateProgramWithSource(context, 1, &sourcesBvh, NULL, &status);
	clErrchk(status);

    //char strInclude[MAX_INCLUDE];
    strcpy(strInclude, "-DGPU_KERNEL -I. -I");
    strcat(strInclude, strResPath);
    strcat(strInclude, "/include");

	status = clBuildProgram(program, 1, devices, strInclude, NULL, NULL);
	clErrchk(status);

	if (status != CL_SUCCESS) {
		LOGE("Failed to build OpenCL kernel (BVH): %d\n", status);

		size_t retValSize;
		clErrchk(clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &retValSize));

		char *buildLog = (char *)malloc(retValSize + 1);
		clErrchk(clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, retValSize, buildLog, NULL));

		buildLog[retValSize] = '\0';
#if 0
		// Query binary (PTX file) size
		size_t bin_sz;
		clErrchk(clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL));

		// Read binary (PTX file) to memory bufferf
		unsigned char *bin = (unsigned char *)malloc(bin_sz);
		clErrchk(clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL));

		char ptxFileName[255];
		strcpy(ptxFileName, bvhFileName);
		strcat(ptxFileName, ".ptx");

		// Save PTX to add_vectors_ocl.ptx
		FILE *fp = fopen(ptxFileName, "wb");
		fwrite(bin, sizeof(char), bin_sz, fp);
		fclose(fp);

		free(bin);
#endif
		LOGE("OpenCL Programm Build Log: %s\n", buildLog);

		char strError[MAX_ERROR];
		strcpy(strError, strResPath);
		strcat(strError, "/error_BVH.txt");

		FILE *fp = fopen(strError, "wt");
		fwrite(buildLog, sizeof(char), retValSize + 1, fp);
		fclose(fp);

		free(buildLog);
		return;
	}

	kernelRad = clCreateKernel(program, "kernelConstructRadixTree", &status);
	clErrchk(status);

	kernelBvh = clCreateKernel(program, "kernelConstructBVHTree", &status);
	clErrchk(status);

	kernelOpt = clCreateKernel(program, "kernelOptimize", &status);
	clErrchk(status);
#endif

	// LordCRC's patch for better workGroupSize
	size_t gsize = 0;
#ifdef EXP_KERNEL
	clErrchk(clGetKernelWorkGroupInfo(kernelGen, devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &gsize, NULL));
#else
	clErrchk(clGetKernelWorkGroupInfo(kernel, devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &gsize, NULL));
#endif

	workGroupSize = (unsigned int) gsize;
	LOGI("OpenCL Device 0: kernel work group size = %d\n", workGroupSize);

	if (forceWorkSize > 0) {
		LOGI("OpenCL Device 0: forced kernel work group size = %d\n", forceWorkSize);
		workGroupSize = forceWorkSize;
	}
}

#ifdef EXP_KERNEL
void ExecuteKernel(cl_kernel p_kernel, int cnt_kernel) {
    /* Enqueue a kernel run call */
    size_t globalThreads[1];

    globalThreads[0] = cnt_kernel;

    if (globalThreads[0] % workGroupSize != 0) globalThreads[0] = (globalThreads[0] / workGroupSize + 1) * workGroupSize;

    size_t localThreads[1];

    localThreads[0] = workGroupSize;

    clErrchk(clEnqueueNDRangeKernel(commandQueue, p_kernel, 1, NULL, globalThreads, localThreads, 0, NULL, NULL));
}
#else
void ExecuteKernel() {
	/* Enqueue a kernel run call */
	size_t globalThreads[1];

	globalThreads[0] = width * height;

	if (globalThreads[0] % workGroupSize != 0) globalThreads[0] = (globalThreads[0] / workGroupSize + 1) * workGroupSize;

	size_t localThreads[1];

	localThreads[0] = workGroupSize;

	clErrchk(clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, globalThreads, localThreads, 0, NULL, NULL));
}
#endif

#ifdef CPU_PARTRENDERING
#ifdef EXP_KERNEL
void DrawBoxCPU(short xstart, short ystart, short bwidth, short bheight, short twidth, short theight, double &cpuTotalTime, double &rwTotalTime, short kindex) {
    const float invWidth = 1.f / twidth;
    const float invHeight = 1.f / theight;

    double cpuStartTime = WallClockTime();

#pragma omp parallel for
    for (int y = ystart; y < ystart + bheight; y++) { /* Loop over image rows */
        for (int x = xstart; x < xstart + bwidth; x++) { /* Loop cols */
            const int i = (y - ystart + 1) > 0 ? (y - ystart) * bwidth + (x - xstart) : (x - xstart);//(bheight - (y - ystart) - 1) * bwidth + (x - xstart);
            const int i2 = i << 1;

            const float r1 = GetRandom(&seeds_boxes[kindex][i2], &seeds_boxes[kindex][i2 + 1]) - .5f;
            const float r2 = GetRandom(&seeds_boxes[kindex][i2], &seeds_boxes[kindex][i2 + 1]) - .5f;
            const float kcx = (x + r1) * invWidth - .5f;
            const float kcy = (y + r2) * invHeight - .5f;

            Vec rdir;
            vinit(rdir,
                  camera.x.s[0] * kcx + camera.y.s[0] * kcy + camera.dir.s[0],
                  camera.x.s[1] * kcx + camera.y.s[1] * kcy + camera.dir.s[1],
                  camera.x.s[2] * kcx + camera.y.s[2] * kcy + camera.dir.s[2]);

            Vec rorig;
            vsmul(rorig, 0.1f, rdir);
            vadd(rorig, rorig, camera.orig)

            vnorm(rdir);
            const Ray ray = { rorig, rdir };

            Vec r;
            vinit(r, 1.0f, 1.0f, 1.0f);

            RadiancePathTracing(shapes, shapeCnt, lightCnt,
#if (ACCELSTR == 1)
                    btn, btl,
#elif (ACCELSTR == 2)
                    pkngbuf, kngCnt, pknbuf, knCnt,
#endif
                    &ray, &seeds_boxes[kindex][i2], &seeds_boxes[kindex][i2 + 1], &r);

            if (currentSample == 0)
                colors_boxes[kindex][i] = r;
            else {
                const float k1 = currentSample;
                const float k2 = 1.f / (k1 + 1.f);

                vsmul(colors_boxes[kindex][i], k1, colors_boxes[kindex][i]);
                vadd(colors_boxes[kindex][i], colors_boxes[kindex][i], r);
                vsmul(colors_boxes[kindex][i], k2, colors_boxes[kindex][i]);
            }

            pixels_boxes[kindex][i] = (toInt((colors_boxes[kindex][i]).s[0]) << 16) |
                                     (toInt((colors_boxes[kindex][i]).s[1]) << 8) |
                                     toInt((colors_boxes[kindex][i]).s[2]) | 0xff000000;
        }
    }
    cpuTotalTime += (WallClockTime() - cpuStartTime);
}

void DrawBoxExpKernel(short xstart, short ystart, short bwidth, short bheight, short twidth, short theight, double &setTotalTime, double &kernelTotalTime, double &rwTotalTime, short kindex) {
    int pindex = 0;
    double setStartTime, kernelStartTime;
    int rayCnt = bwidth * bheight;

#if 0
    /* Set kernel arguments */
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&cameraBuffer));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&seedsboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *)&xstart));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *)&ystart));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *)&bwidth));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *)&bheight));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *)&twidth));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *)&theight));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&rayboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&throughputboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&specularBounceboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&terminatedboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *)&resultboxBuffer[kindex]));

    ExecuteKernel(kernelGenBox, rayCnt);
    //clFinish(commandQueue);
#endif
    //clErrchk(clEnqueueWriteBuffer(commandQueue, cameraBuffer, CL_TRUE, 0, sizeof(Camera), &camera, 0, NULL, NULL));
    //clErrchk(clEnqueueWriteBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) * width * height * 2, seeds, 0, NULL, NULL));
    clErrchk(clEnqueueWriteBuffer(commandQueue, rayboxBuffer[kindex], CL_TRUE, 0, sizeof(Ray) *  bwidth * bheight, rays_boxes[kindex], 0, NULL, NULL));
    clErrchk(clEnqueueWriteBuffer(commandQueue, throughputboxBuffer[kindex], CL_TRUE, 0, sizeof(Vec) *  bwidth * bheight, throughputs_boxes[kindex], 0, NULL, NULL));
    clErrchk(clEnqueueWriteBuffer(commandQueue, specularBounceboxBuffer[kindex], CL_TRUE, 0, sizeof(char) *  bwidth * bheight, specularBounce_boxes[kindex], 0, NULL, NULL));
    clErrchk(clEnqueueWriteBuffer(commandQueue, terminatedboxBuffer[kindex], CL_TRUE, 0, sizeof(char) *  bwidth * bheight, terminated_boxes[kindex], 0, NULL, NULL));
    clErrchk(clEnqueueWriteBuffer(commandQueue, resultboxBuffer[kindex], CL_TRUE, 0, sizeof(Result) *  bwidth * bheight, results_boxes[kindex], 0, NULL, NULL));

    for (int i = 0; i < MAX_DEPTH; i++)
    {
        pindex = 0;

        setStartTime = WallClockTime();
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&shapeBuffer));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&shapeCnt));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&lightCnt));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&xstart));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&ystart));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&bwidth));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&bheight));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&twidth));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&theight));
#if (ACCELSTR == 1)
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&btnBuffer));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&btlBuffer));
#elif (ACCELSTR == 2)
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&kngBuffer));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&kngCnt));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&knBuffer));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(short), (void *)&knCnt));
#endif
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&rayboxBuffer[kindex]));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&seedsboxBuffer[kindex]));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&throughputboxBuffer[kindex]));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&specularBounceboxBuffer[kindex]));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&terminatedboxBuffer[kindex]));
        clErrchk(clSetKernelArg(kernelRadianceBox, pindex++, sizeof(cl_mem), (void *)&resultboxBuffer[kindex]));
        setTotalTime += (WallClockTime() - setStartTime);

        kernelStartTime = WallClockTime();
        ExecuteKernel(kernelRadianceBox, rayCnt);
        //clFinish(commandQueue);
        kernelTotalTime += (WallClockTime() - kernelStartTime);
    }

    pindex = 0;

    setStartTime = WallClockTime();
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(short), (void *)&xstart));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(short), (void *)&ystart));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(short), (void *)&bwidth));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(short), (void *)&bheight));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(short), (void *)&currentSample));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(cl_mem), (void *)&colorboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(cl_mem), (void *)&resultboxBuffer[kindex]));
    clErrchk(clSetKernelArg(kernelFillBox, pindex++, sizeof(cl_mem), (void *)&pixelboxBuffer[kindex]));
    setTotalTime += (WallClockTime() - setStartTime);

    kernelStartTime = WallClockTime();
    ExecuteKernel(kernelFillBox, rayCnt);
    //clFinish(commandQueue);
    kernelTotalTime += (WallClockTime() - kernelStartTime);
}

unsigned int *DrawAllBoxesExpKernel(int bwidth, int bheight, float *rCPU, bool bFirst) {
	int startSampleCount = currentSample, nGPU = 0, nCPU = 1, index, wratio = width / bwidth;
	bool cpuTurn = false;

	double startTime = WallClockTime(), setTotalTime = 0.0, kernelTotalTime = 0.0, rwTotalTime = 0.0, cpuTotalTime = 0.0;

    index = 0;

    for (int y = 0; y < height; y += bheight) {
        for (int x = 0; x < width; x += bwidth) {
            for (int i = 0; i < MAX_SPP; i++) {
                if (cpuTurn) {
                    if (cpuProcessed[index] == 0) {
                        double rwStartTime = WallClockTime();
                        clErrchk(clEnqueueReadBuffer(commandQueue, colorboxBuffer[index], CL_TRUE, 0, sizeof(Vec) * bwidth * bheight, colors_boxes[index], 0, NULL, NULL));
                        clErrchk(clEnqueueReadBuffer(commandQueue, seedsboxBuffer[index], CL_TRUE, 0, sizeof(unsigned int) * bwidth * bheight * 2, seeds_boxes[index], 0, NULL, NULL));
                        rwTotalTime += (WallClockTime() - rwStartTime);
                    }
                    DrawBoxCPU(x, y, bwidth, bheight, width, height, cpuTotalTime, rwTotalTime, index);

                    nCPU++;
                } else {
                    if (cpuProcessed[index] == 1) {
                        double rwStartTime = WallClockTime();
                        clErrchk(clEnqueueWriteBuffer(commandQueue, colorboxBuffer[index], CL_TRUE, 0, sizeof(Vec) * bwidth * bheight, colors_boxes[index], 0, NULL, NULL));
                        clErrchk(clEnqueueWriteBuffer(commandQueue, seedsboxBuffer[index], CL_TRUE, 0, sizeof(unsigned int) * bwidth * bheight * 2, seeds_boxes[index], 0, NULL, NULL));
                        rwTotalTime += (WallClockTime() - rwStartTime);
                    }
                    DrawBoxExpKernel(x, y, bwidth, bheight, width, height, setTotalTime, kernelTotalTime, rwTotalTime, index);

                    nGPU++;
                }
            }
            if (cpuTurn) cpuProcessed[index] = 1;
            else cpuProcessed[index] = 0;

            if ((float) nGPU / nCPU >= *rCPU) cpuTurn = true;
            else cpuTurn = false;

            index++;
        }
    }

    for(int i = 0; i < index; i++) {
        if (cpuProcessed[i] == 0) {
            double rwStartTime = WallClockTime();
            clErrchk(clEnqueueReadBuffer(commandQueue, pixelboxBuffer[i], CL_TRUE, 0, sizeof(unsigned int) * bwidth * bheight, pixels_boxes[i], 0, NULL, NULL));
            rwTotalTime += (WallClockTime() - rwStartTime);
        }

        const int x = (i % wratio) * bwidth, y = (i / wratio) * bheight;

        for (register int j = 0; j < bheight; j++) {
            memcpy((char *) pixels + (y + j) * width * sizeof(unsigned int) + x * sizeof(unsigned int), (char *) pixels_boxes[i] + j * sizeof(unsigned int) * bwidth, sizeof(unsigned int) * bwidth);
        }
    }

	if (bFirst) *rCPU = (cpuTotalTime / (nCPU - 1)) / ((setTotalTime + kernelTotalTime + rwTotalTime) / nGPU);

	currentSample += MAX_SPP;

	/*------------------------------------------------------------------------*/
	const double elapsedTime = WallClockTime() - startTime;
	const int samples = currentSample - startSampleCount;
	const double sampleSec = samples * height * width / elapsedTime;

	LOGI("Set time %.5f msec, Kernel time %.5f msec, CPU time %.5f msec, RW time %.5f msec, Total time %.5f msec (pass %d)  Sample/sec  %.1fK\n",
		 setTotalTime, kernelTotalTime, cpuTotalTime, rwTotalTime, elapsedTime, currentSample, sampleSec / 1000.f);

	return pixels;
}
#else
void ExecuteBoxKernel(int x, int y, int bwidth, int bheight) {
	/* Enqueue a kernel run call */
	size_t globalThreads[1];

	globalThreads[0] = bwidth * bheight;

	if (globalThreads[0] % workGroupSize != 0)
		globalThreads[0] = (globalThreads[0] / workGroupSize + 1) * workGroupSize;

	size_t localThreads[1];

	localThreads[0] = workGroupSize;

	clErrchk(clEnqueueNDRangeKernel(commandQueue, kernelBox, 1, NULL, globalThreads, localThreads, 0, NULL, NULL));
}

void DrawBox(int xstart, int ystart, int bwidth, int bheight, int twidth, int theight, double &cpuTotalTime, double &rwTotalTime) {
	const float invWidth = 1.f / twidth;
	const float invHeight = 1.f / theight;

	double rwStartTime = WallClockTime();

	clErrchk(clEnqueueReadBuffer(commandQueue, colorBuffer, CL_TRUE, 0, sizeof(Vec) * twidth * theight, colors, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, pixelBuffer, CL_TRUE, 0, sizeof(unsigned int) * twidth * theight, pixels, 0, NULL, NULL));

	rwTotalTime += (WallClockTime() - rwStartTime);

	double cpuStartTime = WallClockTime();

#pragma omp parallel for
	for (int y = ystart; y < ystart + bheight; y++) { /* Loop over image rows */
		for (int x = xstart; x < xstart + bwidth; x++) { /* Loop cols */
			const int i = (theight - y) > 0 ? (theight - y - 1) * twidth + x : x;
			const int i2 = i << 1;

			const float r1 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
			const float r2 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
			const float kcx = (x + r1) * invWidth - .5f;
			const float kcy = (y + r2) * invHeight - .5f;

			Vec rdir;
			vinit(rdir,
				camera.x.s[0] * kcx + camera.y.s[0] * kcy + camera.dir.s[0],
				camera.x.s[1] * kcx + camera.y.s[1] * kcy + camera.dir.s[1],
				camera.x.s[2] * kcx + camera.y.s[2] * kcy + camera.dir.s[2]);

			Vec rorig;
			vsmul(rorig, 0.1f, rdir);
			vadd(rorig, rorig, camera.orig)

			vnorm(rdir);
			const Ray ray = { rorig, rdir };

			Vec r;
			r.s[0] = r.s[1] = r.s[2] = 1.0f;

			RadiancePathTracing(shapes, shapeCnt, lightCnt,
#if (ACCELSTR == 1)
				btn, btl,
#elif (ACCELSTR == 2)
				pkngbuf, kngCnt, pknbuf, knCnt,
#endif
				&ray, &seeds[i2], &seeds[i2 + 1], &r);

			if (currentSample == 0)
				colors[i] = r;
			else {
				const float k1 = currentSample;
				const float k2 = 1.f / (k1 + 1.f);
				colors[i].s[0] = (colors[i].s[0] * k1 + r.s[0]) * k2;
				colors[i].s[1] = (colors[i].s[1] * k1 + r.s[1]) * k2;
				colors[i].s[2] = (colors[i].s[2] * k1 + r.s[2]) * k2;
			}

			pixels[y * twidth + x] = toInt(colors[i].s[0]) |
				(toInt(colors[i].s[1]) << 8) |
				(toInt(colors[i].s[2]) << 16);
		}
	}
	cpuTotalTime += (WallClockTime() - cpuStartTime);

	rwStartTime = WallClockTime();

	clErrchk(clEnqueueWriteBuffer(commandQueue, pixelBuffer, CL_TRUE, 0, sizeof(unsigned int) * twidth * theight, pixels, 0, NULL, NULL));
	clErrchk(clEnqueueWriteBuffer(commandQueue, colorBuffer, CL_TRUE, 0, sizeof(Vec) * twidth * theight, colors, 0, NULL, NULL));

	rwTotalTime += (WallClockTime() - rwStartTime);
}

unsigned int *DrawAllBoxes(int bwidth, int bheight, float *rCPU, bool bFirst) {
	int startSampleCount = currentSample, nGPU = 0, nCPU = 1, index = 0;
	bool cpuTurn = false;
	double startTime = WallClockTime(), setStartTime, kernelStartTime;
	double setTotalTime = 0.0, kernelTotalTime = 0.0, rwTotalTime = 0.0;
	double cpuTotalTime = 0.0;

	cl_int status;

	for (int y = 0; y < height; y += bheight) {
		for (int x = 0; x < width; x += bwidth) {
			if (cpuTurn) {
				DrawBox(x, y, bwidth, bheight, width, height, cpuTotalTime, rwTotalTime);

				nCPU++;

				if ((float)nGPU / nCPU >= *rCPU) cpuTurn = true;
				else cpuTurn = false;
			}
			else
			{
				index = 0;
				setStartTime = WallClockTime();

				/* Set kernel arguments */
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&colorBuffer));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&seedBuffer));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&cameraBuffer));
                clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&shapeBuffer));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&shapeCnt));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&lightCnt));
#if (ACCELSTR == 1)
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&btnBuffer));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&btlBuffer));
#elif (ACCELSTR == 2)
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&kngBuffer));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&kngCnt));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&knBuffer));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&knCnt));
#endif
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&x));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&y));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&bwidth));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&bheight));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&width));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&height));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(short), (void *)&currentSample));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&pixelBuffer));
#ifdef DEBUG_INTERSECTIONS
				int *debug1 = (int *)malloc(sizeof(int) * shapeCnt);
				memset(debug1, 0, sizeof(int) * shapeCnt);

				cl_mem debugBuffer1 = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int) * shapeCnt, NULL, &status);
				clErrchk(status);

				clErrchk(clEnqueueWriteBuffer(commandQueue, debugBuffer1, CL_TRUE, 0, sizeof(int) * shapeCnt, debug1, 0, NULL, NULL));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&debugBuffer1));

				float *debug2 = (float *)malloc(6 * sizeof(float) * shapeCnt);
				memset(debug2, 0, sizeof(float) * shapeCnt);

				cl_mem debugBuffer2 = clCreateBuffer(context, CL_MEM_READ_WRITE, 6 * sizeof(float) * shapeCnt, NULL, &status);
				clErrchk(status);

				clErrchk(clEnqueueWriteBuffer(commandQueue, debugBuffer2, CL_TRUE, 0, 6 * sizeof(float) * shapeCnt, debug2, 0, NULL, NULL));
				clErrchk(clSetKernelArg(kernelBox, index++, sizeof(cl_mem), (void *)&debugBuffer2));
#endif
				setTotalTime += (WallClockTime() - setStartTime);

				//--------------------------------------------------------------------------
#ifdef CURRENT_SAMPLE
				if (currentSample < 20) {
#endif
					kernelStartTime = WallClockTime();
					ExecuteBoxKernel(x, y, bwidth, bheight);
					//clFinish(commandQueue);
					kernelTotalTime += (WallClockTime() - kernelStartTime);
#ifdef CURRENT_SAMPLE
				}
				else {
					/* After first 20 samples, continue to execute kernels for more and more time */
					const float k = min(currentSample - 20, 100) / 100.f;
					const float tresholdTime = 0.5f * k * 1000.0f;
					for (;;) {
						kernelStartTime = WallClockTime();
						ExecuteBoxKernel();
						//clFinish(commandQueue);
						kernelTotalTime += (WallClockTime() - kernelStartTime);

						currentSample++;
						const float elapsedTime = WallClockTime() - startTime;
						if (elapsedTime > tresholdTime)
							break;
					}
				}
#endif
				//readStartTime = WallClockTime();

				//clFinish(commandQueue);
				//readTotalTime += (WallClockTime() - readStartTime);
#ifdef DEBUG_INTERSECTIONS
				status = clEnqueueReadBuffer(commandQueue, debugBuffer1, CL_TRUE, 0, sizeof(int) * shapeCnt, debug1, 0, NULL, NULL);
				clErrchk(status);

				clErrchk(clEnqueueReadBuffer(commandQueue, debugBuffer2, CL_TRUE, 0, 6 * sizeof(float) * shapeCnt, debug2, 0, NULL, NULL));

				FILE *f = fopen("images\\intersections.txt", "wt"); // Write image to PPM file.
				for (int i = 0; i < shapeCnt; i++) {
					fprintf(f, "%d, ", debug1[i]);
				}
				fprintf(f, "\n");
				for (int i = 0; i < 6 * shapeCnt; i++) {
					fprintf(f, "%f, ", debug2[i]);
				}
				fclose(f);

				clErrchk(clReleaseMemObject(debugBuffer2));
				free(debug2);

				clErrchk(clReleaseMemObject(debugBuffer1));
				free(debug1);
#endif
				nGPU++;

				if ((float)nGPU / nCPU >= *rCPU) cpuTurn = true;
				else cpuTurn = false;
			}
		}
	}

	double rwStartTime = WallClockTime();
	clErrchk(clEnqueueReadBuffer(commandQueue, pixelBuffer, CL_TRUE, 0, sizeof(unsigned int) * width * height, pixels, 0, NULL, NULL));
	rwTotalTime += (WallClockTime() - rwStartTime);

	if (bFirst) *rCPU = (cpuTotalTime / (nCPU - 1)) / ((setTotalTime + kernelTotalTime + rwTotalTime) / nGPU);

	currentSample++;

	/*------------------------------------------------------------------------*/
	const double elapsedTime = WallClockTime() - startTime;
	const int samples = currentSample - startSampleCount;
	const double sampleSec = samples * height * width / elapsedTime;

	LOGI("Set time %.5f msec, Kernel time %.5f msec, CPU time %.5f msec, RW time %.5f msec, Total time %.5f msec (pass %d)  Sample/sec  %.1fK\n",
		setTotalTime, kernelTotalTime, cpuTotalTime, rwTotalTime, elapsedTime, currentSample, sampleSec / 1000.f);

	return pixels;
}
#endif
unsigned int *DrawFrame()
{
    static float rCPU = 1.0f;
    static bool first = true;

#ifdef EXP_KERNEL
    unsigned int *pPixels = DrawAllBoxesExpKernel(BWIDTH, BHEIGHT, &rCPU, first);
#else
    unsigned int *pPixels = DrawAllBoxes(160, 120, &rCPU, first);
#endif
    first = false;

    return pPixels;
}
#else //NOT CPU_PARTRENDERING
#define NUM_THREADS 7
pthread_t threads[NUM_THREADS];

void *perform_work(void *arguments){
    int index = *((int *)arguments);
    double startTime = WallClockTime();

    const float invWidth = 1.f / (float)width;
    const float invHeight = 1.f / (float)height;

    for(int y = index; y < (int) height; y += NUM_THREADS) {
        for(int x = 0; x < (int) width; x++) {
            const int sgid = y > 0 ? (y - 1) * width + x : x;
            const int sgid2 = sgid << 1;

            const float r1 = GetRandom(&seeds[sgid2], &seeds[sgid2 + 1]) - .5f;
            const float r2 = GetRandom(&seeds[sgid2], &seeds[sgid2 + 1]) - .5f;

            const float kcx = (x + r1) * invWidth - .5f;
            const float kcy = (y + r2) * invHeight - .5f;

            Vec rdir;
            vinit(rdir,
                  cameraLeft.x.s[0] * kcx + cameraLeft.y.s[0] * kcy + cameraLeft.dir.s[0],
				  cameraLeft.x.s[1] * kcx + cameraLeft.y.s[1] * kcy + cameraLeft.dir.s[1],
				  cameraLeft.x.s[2] * kcx + cameraLeft.y.s[2] * kcy + cameraLeft.dir.s[2]);

            Vec rorig;
            vsmul(rorig, 0.1f, rdir);
            vadd(rorig, rorig, cameraLeft.orig);

            vnorm(rdir);
            rinit(ray[sgid], rorig, rdir);
        }
    }

    LOGI("Thread %d, RayChange time %.5f msec\n", index, WallClockTime() - startTime);

    return NULL;
}

void ChangeRays() {
    int thread_args[NUM_THREADS];

    //create all threads one by one
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_args[i] = i;
        pthread_create(&threads[i], NULL, perform_work, &thread_args[i]);
    }

    //wait for each thread to complete
    //for (int i = 0; i < NUM_THREADS; i++) {
    //    pthread_join(threads[i], NULL);
    //}
}

unsigned int *DrawFrame() {
	int len = pixelCount * sizeof(unsigned int), index = 0;
	double startTime = WallClockTime(), setStartTime, kernelStartTime, readStartTime;
	double setTotalTime = 0.0, kernelTotalTime = 0.0, readTotalTime = 0.0;
	int startSampleCount = currentSampleLeft;

#ifdef EXP_KERNEL
	for (int i = 0; i < MAX_SPP; i++)
	{
		int rayCnt = width * height;

        for (int j = 0; j < MAX_DEPTH; j++)
		{
			index = 0;

			setStartTime = WallClockTime();
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&shapeBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *)&shapeCnt));
            clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *)&lightCnt));
            clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *)&width));
            clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *)&height));
#if (ACCELSTR == 1)
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&btnBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&btlBuffer));
#elif (ACCELSTR == 2)
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&kngBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *)&kngCnt));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&knBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *)&knCnt));
#endif
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&rayBuffer));
            clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&seedBuffer));
            clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&throughputBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&specularBounceBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&terminatedBuffer));
			clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&resultBuffer));
			setTotalTime += (WallClockTime() - setStartTime);

			kernelStartTime = WallClockTime();
			ExecuteKernel(kernelRadiance, rayCnt);
			//clFinish(commandQueue);
			kernelTotalTime += (WallClockTime() - kernelStartTime);
#if 0
            //clErrchk(clEnqueueReadBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  rayCnt, ray, 0, NULL, NULL));
            //clErrchk(clEnqueueReadBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) *  rayCnt * 2, seeds, 0, NULL, NULL));
            //clErrchk(clEnqueueReadBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  rayCnt, throughput, 0, NULL, NULL));
            //clErrchk(clEnqueueReadBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, specularBounce, 0, NULL, NULL));
            clErrchk(clEnqueueReadBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, terminated, 0, NULL, NULL));
            //clErrchk(clEnqueueReadBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));

            int rc = 0, lastSwapped = rayCnt - 1;

            for(int k = 0; k < rayCnt; k++) {
                if (terminated[k] == 1)
                {
                    for(int l = lastSwapped; l > k; l--)
                    {
                        if (terminated[l] == 0) {
                            Ray tRay = ray[k]; ray[k] = ray[l]; ray[l] = tRay;

                            unsigned int tSeeds1 = seeds[2 * k]; seeds[2 * k] = seeds[2 * l]; seeds[2 * l] = tSeeds1;
                            unsigned int tSeeds2 = seeds[2 * k + 1]; seeds[2 * k + 1] = seeds[2 * l + 1]; seeds[2 * l + 1] = tSeeds2;

                            Vec tThroughput = throughput[k]; throughput[k] = throughput[l]; throughput[l] = tThroughput;
                            char tSB = specularBounce[k]; specularBounce[k] = specularBounce[l]; specularBounce[l] = tSB;
                            char tT = terminated[k]; terminated[k] = terminated[l]; terminated[l] = tT;
                            //result[k] = result[l];

                            lastSwapped = l - 1;
                            break;
                        }
                    }
                    rc++;
                }
            }

            //rayCnt = rc;
            if (rc == rayCnt)
            {
                LOGI("Loop finished\n");
                break;
            }
			//LOGI("Ray Count: %d", rc);

            //clErrchk(clEnqueueWriteBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  rayCnt, ray, 0, NULL, NULL));
            //clErrchk(clEnqueueWriteBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) *  rayCnt * 2, seeds, 0, NULL, NULL));
            //clErrchk(clEnqueueWriteBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  rayCnt, throughput, 0, NULL, NULL));
            //clErrchk(clEnqueueWriteBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, specularBounce, 0, NULL, NULL));
            //clErrchk(clEnqueueWriteBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, terminated, 0, NULL, NULL));
            //clErrchk(clEnqueueWriteBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));
#endif
        }

		index = 0;

		setStartTime = WallClockTime();
		clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *)&width));
		clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *)&height));
		clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *)&currentSampleLeft));
        clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *)&colorBufferLeft));
        clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *)&resultBuffer));
		clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *)&pixelBufferLeft));
		setTotalTime += (WallClockTime() - setStartTime);

		kernelStartTime = WallClockTime();
		ExecuteKernel(kernelFill, width * height);
		//clFinish(commandQueue);
		kernelTotalTime += (WallClockTime() - kernelStartTime);

#if 1
		index = 0;

		/* Set kernel arguments */
		setStartTime = WallClockTime();
		clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&cameraBufferLeft));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&seedBuffer));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(short), (void *)&width));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(short), (void *)&height));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&rayBuffer));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&throughputBuffer));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&specularBounceBuffer));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&terminatedBuffer));
        clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&resultBuffer));
		setTotalTime += (WallClockTime() - setStartTime);

		kernelStartTime = WallClockTime();
		ExecuteKernel(kernelGen, rayCnt);
		//clFinish(commandQueue);
		kernelTotalTime += (WallClockTime() - kernelStartTime);
#else
		ChangeRays();

		//clErrchk(clEnqueueWriteBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  width * height, ray, 0, NULL, NULL));
		clErrchk(clEnqueueWriteBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  width * height, throughput, 0, NULL, NULL));
		clErrchk(clEnqueueWriteBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  width * height, specularBounce, 0, NULL, NULL));
		clErrchk(clEnqueueWriteBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  width * height, terminated, 0, NULL, NULL));
		clErrchk(clEnqueueWriteBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));
        clErrchk(clEnqueueWriteBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  width * height, ray, 0, NULL, NULL));
#endif
	}

	//--------------------------------------------------------------------------
	/* Enqueue readBuffer */
    readStartTime = WallClockTime();
    clErrchk(clEnqueueReadBuffer(commandQueue, pixelBufferLeft, CL_TRUE, 0, len, pixelsLeft, 0, NULL, NULL));
	//clFinish(commandQueue);
	readTotalTime += (WallClockTime() - readStartTime);

    currentSampleLeft += MAX_SPP;
#else
	setStartTime = WallClockTime();

	/* Set kernel arguments */
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&colorBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&seedBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &cameraBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&shapeBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *) &shapeCnt));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&lightCnt));
#if (ACCELSTR == 1)
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&btnBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&btlBuffer));
#elif (ACCELSTR == 2)
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&kngBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&kngCnt));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&knBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&knCnt));
#endif
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *) &width));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *) &height));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&currentSample));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &pixelBuffer));
#ifdef DEBUG_INTERSECTIONS
    cl_int status;
    int *debug1 = (int *)malloc(sizeof(int) * shapeCnt);
    memset(debug1, 0, sizeof(int) * shapeCnt);

    cl_mem debugBuffer1 = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int) * shapeCnt, NULL, &status);
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, debugBuffer1, CL_TRUE, 0, sizeof(int) * shapeCnt, debug1, 0, NULL, NULL));

	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &debugBuffer1));

	float *debug2 = (float *)malloc(6 * sizeof(float) * shapeCnt);
	memset(debug2, 0, sizeof(float) * shapeCnt);

	cl_mem debugBuffer2 = clCreateBuffer(context, CL_MEM_READ_WRITE, 6 * sizeof(float) * shapeCnt, NULL, &status);
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, debugBuffer2, CL_TRUE, 0, 6 * sizeof(float) * shapeCnt, debug2, 0, NULL, NULL));

	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &debugBuffer2));
#endif
	setTotalTime += (WallClockTime() - setStartTime);

	//--------------------------------------------------------------------------
#ifdef CURRENT_SAMPLE
	if (currentSample < 20) {
#endif
	kernelStartTime = WallClockTime();
	ExecuteKernel();
	//clFinish(commandQueue);
	kernelTotalTime += (WallClockTime() - kernelStartTime);

#ifdef CURRENT_SAMPLE
    }
    else {
        /* After first 20 samples, continue to execute kernels for more and more time */
        const float k = min(currentSample - 20, 100) / 100.f;
        const float tresholdTime = 0.5f * k * 1000.0f;
        for (;;) {
            kernelStartTime = WallClockTime();
            ExecuteKernel();
            //clFinish(commandQueue);
            kernelTotalTime += (WallClockTime() - kernelStartTime);

            currentSample++;
            const float elapsedTime = WallClockTime() - startTime;
            if (elapsedTime > tresholdTime)
                break;
        }
    }
#endif
    readStartTime = WallClockTime();
    //--------------------------------------------------------------------------
    /* Enqueue readBuffer */
	clErrchk(clEnqueueReadBuffer(commandQueue, pixelBuffer, CL_TRUE, 0, len, pixels, 0, NULL, NULL));
	//clFinish(commandQueue);
	readTotalTime += (WallClockTime() - readStartTime);
#ifdef DEBUG_INTERSECTIONS
	clErrchk(clEnqueueReadBuffer(commandQueue, debugBuffer1, CL_TRUE, 0, sizeof(int) * shapeCnt, debug1, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, debugBuffer2, CL_TRUE, 0, 6 * sizeof(float) * shapeCnt, debug2, 0, NULL, NULL));

	char strLog[MAX_LOG];
	strcpy(strLog, strResPath);
	strcat(strLog, "/intersections.txt");

    FILE *f = fopen(strLog, "wt"); // Write image to PPM file.
	for(int i = 0; i < shapeCnt; i++) {
		fprintf(f, "%d, ", debug1[i]);
	}
	fprintf(f, "\n");
	for(int i = 0; i < 6 * shapeCnt; i++) {
		fprintf(f, "%f, ", debug2[i]);
	}
	fclose(f);

	clErrchk(clReleaseMemObject(debugBuffer2));
    free(debug2);

	clErrchk(clReleaseMemObject(debugBuffer1));
	free(debug1);
#endif
    currentSample++;
#endif

	/*------------------------------------------------------------------------*/
	const double elapsedTime = WallClockTime() - startTime;
	const int samples = currentSampleLeft - startSampleCount;
	const double sampleSec = samples * height * width / elapsedTime;
	LOGI("Set time %.5f msec, Kernel time %.5f msec, Read time %.5f msec, Total time %.5f msec (pass %d)  Sample/sec  %.1fK\n",
		setTotalTime, kernelTotalTime, readTotalTime, elapsedTime, currentSampleLeft, sampleSec / 1000.f);

	return pixelsLeft;
}

bool RectangleIntersect(const Vec n, const Vec p0, const Vec l0, const Vec l, float *t)
{
    // assuming vectors are all normalized
    float denom = vdot(n, l);
    if (denom > 1e-6) {
        Vec p0l0;
        vsub(p0l0, p0, l0); //= p0 - l0;
        *t = vdot(p0l0, n) / denom;
        return (*t >= 0);
    }

    return false;
}

unsigned int *DrawFrameVR(short bleft) {
	int len = pixelCount * sizeof(unsigned int), index = 0;
	double startTime = WallClockTime(), setStartTime, kernelStartTime, readStartTime, setTotalTime = 0.0, kernelTotalTime = 0.0, readTotalTime = 0.0;

	int currentSample;

	if (bleft) currentSample = currentSampleLeft;
	else currentSample = currentSampleRight;

	const int startSampleCount = currentSample;

#ifdef EXP_KERNEL
	if (bleft) {
        for (int i = 0; i < MAX_SPP; i++) {
            int rayCnt = width * height;
#if 1
            index = 0;

            /* Set kernel arguments */
            setStartTime = WallClockTime();
            if (bleft) { clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &cameraBufferLeft)); }
            else clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &cameraBufferRight));

            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &seedBuffer));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(short), (void *) &width));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(short), (void *) &height));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &rayBuffer));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &throughputBuffer));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &specularBounceBuffer));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &terminatedBuffer));
            clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *) &resultBuffer));
            setTotalTime += (WallClockTime() - setStartTime);

            kernelStartTime = WallClockTime();
            ExecuteKernel(kernelGen, rayCnt);
            //clFinish(commandQueue);
            kernelTotalTime += (WallClockTime() - kernelStartTime);
#else
            ChangeRays();

            //clErrchk(clEnqueueWriteBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  width * height, ray, 0, NULL, NULL));
            clErrchk(clEnqueueWriteBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  width * height, throughput, 0, NULL, NULL));
            clErrchk(clEnqueueWriteBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  width * height, specularBounce, 0, NULL, NULL));
            clErrchk(clEnqueueWriteBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  width * height, terminated, 0, NULL, NULL));
            clErrchk(clEnqueueWriteBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));
            clErrchk(clEnqueueWriteBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  width * height, ray, 0, NULL, NULL));
#endif
            for (short j = 0; j < MAX_DEPTH; j++) {
                index = 0;

                setStartTime = WallClockTime();
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &shapeBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &shapeCnt));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &lightCnt));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &width));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &height));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &j));
#if (ACCELSTR == 1)
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&btnBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *)&btlBuffer));
#elif (ACCELSTR == 2)
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &kngBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &kngCnt));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &knBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(short), (void *) &knCnt));
#endif
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &rayBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &seedBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &throughputBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &specularBounceBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &terminatedBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &resultBuffer));
                clErrchk(clSetKernelArg(kernelRadiance, index++, sizeof(cl_mem), (void *) &fhiBuffer));
                setTotalTime += (WallClockTime() - setStartTime);

                kernelStartTime = WallClockTime();
                ExecuteKernel(kernelRadiance, rayCnt);
                //clFinish(commandQueue);
                kernelTotalTime += (WallClockTime() - kernelStartTime);
#if 0
                //clErrchk(clEnqueueReadBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  rayCnt, ray, 0, NULL, NULL));
                //clErrchk(clEnqueueReadBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) *  rayCnt * 2, seeds, 0, NULL, NULL));
                //clErrchk(clEnqueueReadBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  rayCnt, throughput, 0, NULL, NULL));
                //clErrchk(clEnqueueReadBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, specularBounce, 0, NULL, NULL));
                clErrchk(clEnqueueReadBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, terminated, 0, NULL, NULL));
                //clErrchk(clEnqueueReadBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));

                int rc = 0, lastSwapped = rayCnt - 1;

                for(int k = 0; k < rayCnt; k++) {
                    if (terminated[k] == 1)
                    {
                        for(int l = lastSwapped; l > k; l--)
                        {
                            if (terminated[l] == 0) {
                                Ray tRay = ray[k]; ray[k] = ray[l]; ray[l] = tRay;

                                unsigned int tSeeds1 = seeds[2 * k]; seeds[2 * k] = seeds[2 * l]; seeds[2 * l] = tSeeds1;
                                unsigned int tSeeds2 = seeds[2 * k + 1]; seeds[2 * k + 1] = seeds[2 * l + 1]; seeds[2 * l + 1] = tSeeds2;

                                Vec tThroughput = throughput[k]; throughput[k] = throughput[l]; throughput[l] = tThroughput;
                                char tSB = specularBounce[k]; specularBounce[k] = specularBounce[l]; specularBounce[l] = tSB;
                                char tT = terminated[k]; terminated[k] = terminated[l]; terminated[l] = tT;
                                //result[k] = result[l];

                                lastSwapped = l - 1;
                                break;
                            }
                        }
                        rc++;
                    }
                }

                //rayCnt = rc;
                if (rc == rayCnt)
                {
                    LOGI("Loop finished\n");
                    break;
                }
                //LOGI("Ray Count: %d", rc);

                //clErrchk(clEnqueueWriteBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  rayCnt, ray, 0, NULL, NULL));
                //clErrchk(clEnqueueWriteBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) *  rayCnt * 2, seeds, 0, NULL, NULL));
                //clErrchk(clEnqueueWriteBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  rayCnt, throughput, 0, NULL, NULL));
                //clErrchk(clEnqueueWriteBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, specularBounce, 0, NULL, NULL));
                //clErrchk(clEnqueueWriteBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  rayCnt, terminated, 0, NULL, NULL));
                //clErrchk(clEnqueueWriteBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));
#endif
            }

            index = 0;

            setStartTime = WallClockTime();
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *) &width));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *) &height));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *) &currentSample));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(short), (void *) &bleft));

            if (bleft) { clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &colorBufferLeft)); }
            else clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &colorBufferRight));

            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &resultBuffer));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &pixelBufferLeft));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &pixelBufferRight));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &cameraBufferRight));
            clErrchk(clSetKernelArg(kernelFill, index++, sizeof(cl_mem), (void *) &fhiBuffer));
            setTotalTime += (WallClockTime() - setStartTime);

            kernelStartTime = WallClockTime();
            ExecuteKernel(kernelFill, width * height);
            //clFinish(commandQueue);
            kernelTotalTime += (WallClockTime() - kernelStartTime);
#if 1
            clEnqueueReadBuffer(commandQueue, fhiBuffer, CL_TRUE, 0, sizeof(FirstHitInfo) * width * height, fhi, 0, NULL, NULL);
            clEnqueueReadBuffer(commandQueue, colorBufferLeft, CL_TRUE, 0, sizeof(Vec) * width * height, colorsLeft, 0, NULL, NULL);

            for(int yLeft = 0; yLeft < height; yLeft++) {
                for (int xLeft = 0; xLeft < width; xLeft++) {
                    const int locPixelLeft = yLeft * width + xLeft;
                    const Vec l0 = fhi[locPixelLeft].firstHitPoint;

                    Vec p0;
                    vadd(p0, cameraRight.start, cameraRight.end);
                    vsmul(p0, 0.5f, p0);

                    Vec n;
                    vsmul(n, -1, cameraRight.dir)

                    Vec vl;
                    vsub(vl, cameraRight.orig, fhi[locPixelLeft].firstHitPoint);
					vnorm(vl);

                    float t;
                    bool bint = RectangleIntersect(n, p0, l0, vl, &t);

                    if (bint) {
                        Vec p;
                        vsmul(p, t, vl);
                        vadd(p, p, l0);
                        //vsmad(p, t, vl, l0);

						if (p.s[0] >= cameraRight.start.s[0] && p.s[0] < cameraRight.end.s[0] && p.s[1] >= cameraRight.start.s[1] && p.s[1] < cameraRight.end.s[1])
						{
							int xRight = ((p.s[0] - cameraRight.start.s[0]) / (cameraRight.end.s[0] - cameraRight.start.s[0])) * (float)width + .5f;
							int yRight = ((p.s[1] - cameraRight.start.s[1]) / (cameraRight.end.s[1] - cameraRight.start.s[1])) * (float)height + .5f;

							const int locPixelRight = (height - yRight - 1) * width + xRight;

							colorsRight[locPixelRight] = colorsLeft[locPixelLeft];

							pixelsRight[locPixelRight] = (toInt(colorsLeft[locPixelLeft].s[0]) << 16) |
														 (toInt(colorsLeft[locPixelLeft].s[1]) << 8) |
														 (toInt(colorsLeft[locPixelLeft].s[2])) | 0xff000000;
						}
                    }
                }
            }

            const float invWidth = 1.f / (float)width;
            const float invHeight = 1.f / (float)height;

            for(int yRight = 0; yRight < height; yRight++) {
                for (int xRight = 0; xRight < width; xRight++) {
                    const int i = yRight * width + xRight;
                    if (!viszero(colorsRight[i])) continue;

                    const int locPixelRight = (height - yRight - 1) * width + xRight;
                    const int i2 = i << 1;

                    const float r1 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
                    const float r2 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
                    const float kcx = (xRight + r1) * invWidth - .5f;
                    const float kcy = (yRight + r2) * invHeight - .5f;

                    Vec rdir;
                    vinit(rdir,
                          cameraRight.x.s[0] * kcx + cameraRight.y.s[0] * kcy + cameraRight.dir.s[0],
                          cameraRight.x.s[1] * kcx + cameraRight.y.s[1] * kcy + cameraRight.dir.s[1],
                          cameraRight.x.s[2] * kcx + cameraRight.y.s[2] * kcy + cameraRight.dir.s[2]);

                    Vec rorig;
                    vsmul(rorig, 0.1f, rdir);
                    vadd(rorig, rorig, cameraRight.orig)

                    vnorm(rdir);
                    const Ray ray = { rorig, rdir };

                    Vec r;
                    vinit(r, 1.0f, 1.0f, 1.0f);

                    RadiancePathTracing(shapes, shapeCnt, lightCnt,
#if (ACCELSTR == 1)
                            btn, btl,
#elif (ACCELSTR == 2)
                            pkngbuf, kngCnt, pknbuf, knCnt,
#endif
                            &ray, &seeds[i2], &seeds[i2 + 1], &r);

                    if (currentSample == 0)
                        colorsRight[i] = r;
                    else {
                        const float k1 = currentSample;
                        const float k2 = 1.f / (k1 + 1.f);

                        vsmul(colorsRight[i], k1, colorsRight[i]);
                        vadd(colorsRight[i], colorsRight[i], r);
                        vsmul(colorsRight[i], k2, colorsRight[i]);
                    }

                    pixelsRight[locPixelRight] = (toInt((colorsRight[i]).s[0]) << 16) |
                                              (toInt((colorsRight[i]).s[1]) << 8) |
                                              toInt((colorsRight[i]).s[2]) | 0xff000000;
                }
            }
#else
            //clEnqueueReadBuffer(commandQueue, fhiBuffer, CL_TRUE, 0, sizeof(FirstHitInfo) * width * height, fhi, 0, NULL, NULL);
            //clEnqueueReadBuffer(commandQueue, colorBufferLeft, CL_TRUE, 0, sizeof(Vec) * width * height, colorsLeft, 0, NULL, NULL);
            clEnqueueReadBuffer(commandQueue, pixelBufferLeft, CL_TRUE, 0, sizeof(unsigned char[4]) * width * height, pixelsLeft, 0, NULL, NULL);

            //memset(pixelsRight, 0, sizeof(unsigned char[4]) * width * height);

            for(int yLeft = 1; yLeft <= height; yLeft++) {
                for (int xLeft = DIFF_LEFTRIGHTEYE; xLeft < width; xLeft++) {
                	//if (xLeft > 2 * DIFF_LEFTRIGHTEYE) {
						const int locPixelLeft = (yLeft - 1) * width + xLeft;
						const int locPixelRight = (yLeft - 1) * width + (xLeft - DIFF_LEFTRIGHTEYE);

						pixelsRight[locPixelRight] = pixelsLeft[locPixelLeft];
                	//}
                }
            }
#endif
        }
    }
	//--------------------------------------------------------------------------
	/* Enqueue readBuffer */
	readStartTime = WallClockTime();

	if (bleft) { clErrchk(clEnqueueReadBuffer(commandQueue, pixelBufferLeft, CL_TRUE, 0, len, pixelsLeft, 0, NULL, NULL)); }
	//else clErrchk(clEnqueueReadBuffer(commandQueue, pixelBufferRight, CL_TRUE, 0, len, pixelsRight, 0, NULL, NULL));

	//clFinish(commandQueue);
	readTotalTime += (WallClockTime() - readStartTime);

    currentSample += MAX_SPP;

	if (bleft) currentSampleLeft = currentSample;
	else currentSampleRight = currentSample;
#else
	setStartTime = WallClockTime();

	/* Set kernel arguments */
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&colorBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&seedBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &cameraBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&shapeBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *) &shapeCnt));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&lightCnt));
#if (ACCELSTR == 1)
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&btnBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&btlBuffer));
#elif (ACCELSTR == 2)
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&kngBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&kngCnt));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *)&knBuffer));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&knCnt));
#endif
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *) &width));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *) &height));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(short), (void *)&currentSample));
	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &pixelBuffer));
#ifdef DEBUG_INTERSECTIONS
    cl_int status;
    int *debug1 = (int *)malloc(sizeof(int) * shapeCnt);
    memset(debug1, 0, sizeof(int) * shapeCnt);

    cl_mem debugBuffer1 = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int) * shapeCnt, NULL, &status);
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, debugBuffer1, CL_TRUE, 0, sizeof(int) * shapeCnt, debug1, 0, NULL, NULL));

	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &debugBuffer1));

	float *debug2 = (float *)malloc(6 * sizeof(float) * shapeCnt);
	memset(debug2, 0, sizeof(float) * shapeCnt);

	cl_mem debugBuffer2 = clCreateBuffer(context, CL_MEM_READ_WRITE, 6 * sizeof(float) * shapeCnt, NULL, &status);
	clErrchk(status);

	clErrchk(clEnqueueWriteBuffer(commandQueue, debugBuffer2, CL_TRUE, 0, 6 * sizeof(float) * shapeCnt, debug2, 0, NULL, NULL));

	clErrchk(clSetKernelArg(kernel, index++, sizeof(cl_mem), (void *) &debugBuffer2));
#endif
	setTotalTime += (WallClockTime() - setStartTime);

	//--------------------------------------------------------------------------
#ifdef CURRENT_SAMPLE
	if (currentSample < 20) {
#endif
	kernelStartTime = WallClockTime();
	ExecuteKernel();
	//clFinish(commandQueue);
	kernelTotalTime += (WallClockTime() - kernelStartTime);

#ifdef CURRENT_SAMPLE
    }
    else {
        /* After first 20 samples, continue to execute kernels for more and more time */
        const float k = min(currentSample - 20, 100) / 100.f;
        const float tresholdTime = 0.5f * k * 1000.0f;
        for (;;) {
            kernelStartTime = WallClockTime();
            ExecuteKernel();
            //clFinish(commandQueue);
            kernelTotalTime += (WallClockTime() - kernelStartTime);

            currentSample++;
            const float elapsedTime = WallClockTime() - startTime;
            if (elapsedTime > tresholdTime)
                break;
        }
    }
#endif
    readStartTime = WallClockTime();
    //--------------------------------------------------------------------------
    /* Enqueue readBuffer */
	clErrchk(clEnqueueReadBuffer(commandQueue, pixelBuffer, CL_TRUE, 0, len, pixels, 0, NULL, NULL));
	//clFinish(commandQueue);
	readTotalTime += (WallClockTime() - readStartTime);
#ifdef DEBUG_INTERSECTIONS
	clErrchk(clEnqueueReadBuffer(commandQueue, debugBuffer1, CL_TRUE, 0, sizeof(int) * shapeCnt, debug1, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, debugBuffer2, CL_TRUE, 0, 6 * sizeof(float) * shapeCnt, debug2, 0, NULL, NULL));

	char strLog[MAX_LOG];
	strcpy(strLog, strResPath);
	strcat(strLog, "/intersections.txt");

    FILE *f = fopen(strLog, "wt"); // Write image to PPM file.
	for(int i = 0; i < shapeCnt; i++) {
		fprintf(f, "%d, ", debug1[i]);
	}
	fprintf(f, "\n");
	for(int i = 0; i < 6 * shapeCnt; i++) {
		fprintf(f, "%f, ", debug2[i]);
	}
	fclose(f);

	clErrchk(clReleaseMemObject(debugBuffer2));
    free(debug2);

	clErrchk(clReleaseMemObject(debugBuffer1));
	free(debug1);
#endif
    currentSample++;
#endif

	/*------------------------------------------------------------------------*/
	const double elapsedTime = WallClockTime() - startTime;
	const int samples = currentSample - startSampleCount;

	const double sampleSec = samples * height * width / elapsedTime;
	LOGI("Set time %.5f msec, Kernel time %.5f msec, Read time %.5f msec, Total time %.5f msec (pass %d)  Sample/sec  %.1fK, Left %d\n",
		 setTotalTime, kernelTotalTime, readTotalTime, elapsedTime, currentSample, sampleSec / 1000.f, bleft);

	if (bleft) return pixelsLeft;
	else return pixelsRight;
}
#endif

void ReInitScene() {
    currentSampleLeft = 0;
	currentSampleRight = 0;

    // Redownload the scene
	clErrchk(clEnqueueWriteBuffer(commandQueue, shapeBuffer, CL_TRUE, 0, sizeof(Shape) * shapeCnt, shapes, 0, NULL, NULL));
}

void ReInit(const int reallocBuffers) {
    // Check if I have to reallocate buffers
    if (reallocBuffers) {
        FreeBuffers();
        UpdateCamera(false);
        AllocateBuffers();
    } else {
        UpdateCamera(false);
	}

	clErrchk(clEnqueueWriteBuffer(commandQueue, cameraBufferLeft, CL_TRUE, 0, sizeof(Camera), &cameraLeft, 0, NULL, NULL));
#if (ACCELSTR == 1)
	clErrchk(clEnqueueWriteBuffer(commandQueue, btnBuffer, CL_TRUE, 0, sizeof(BVHNodeGPU) * (shapeCnt - 1), btn, 0, NULL, NULL));
	clErrchk(clEnqueueWriteBuffer(commandQueue, btlBuffer, CL_TRUE, 0, sizeof(BVHNodeGPU) * (shapeCnt), btl, 0, NULL, NULL));
#elif (ACCELSTR == 2)
	clErrchk(clEnqueueWriteBuffer(commandQueue, kngBuffer, CL_TRUE, 0, sizeof(KDNodeGPU) * (kngCnt), pkngbuf, 0, NULL, NULL));
	clErrchk(clEnqueueWriteBuffer(commandQueue, knBuffer, CL_TRUE, 0, sizeof(int) * (knCnt), pknbuf, 0, NULL, NULL));
#endif

    currentSampleLeft = 0;
	currentSampleRight = 0;

#ifdef EXP_KERNEL
#ifdef CPU_PARTRENDERING
	int pindex = 0, kindex = 0;
	short bwidth = BWIDTH, bheight = BHEIGHT, twidth = width, theight = height;

    for(int ystart = 0; ystart < height; ystart += bheight) {
        for(int xstart = 0; xstart < width; xstart += bwidth) {
            for(int i = 0; i< bheight; i++) {
                memcpy((char *) colors_boxes[kindex] + i * sizeof(Vec) * bwidth,
                       (char *) colors + xstart * sizeof(Vec) + (ystart + i) * width * sizeof(Vec),
                       sizeof(Vec) * bwidth);
            }

            clErrchk(clEnqueueWriteBuffer(commandQueue, colorboxBuffer[kindex], CL_TRUE, 0, sizeof(Vec) * bwidth * bheight, colors_boxes[kindex], 0, NULL, NULL));

			pindex = 0;

            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &cameraBuffer));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &seedsboxBuffer[kindex]));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *) &xstart));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *) &ystart));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *) &bwidth));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *) &bheight));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *) &twidth));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(short), (void *) &theight));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &rayboxBuffer[kindex]));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &throughputboxBuffer[kindex]));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &specularBounceboxBuffer[kindex]));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &terminatedboxBuffer[kindex]));
            clErrchk(clSetKernelArg(kernelGenBox, pindex++, sizeof(cl_mem), (void *) &resultboxBuffer[kindex]));

            ExecuteKernel(kernelGenBox, bwidth * bheight);
            //clFinish(commandQueue);

            //clErrchk(clEnqueueReadBuffer(commandQueue, cameraBuffer, CL_TRUE, 0, sizeof(Camera), &camera, 0, NULL, NULL));
            //clErrchk(clEnqueueReadBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) * width * height * 2, seeds, 0, NULL, NULL));
            clErrchk(clEnqueueReadBuffer(commandQueue, rayboxBuffer[kindex], CL_TRUE, 0, sizeof(Ray) *  bwidth * bheight, rays_boxes[kindex], 0, NULL, NULL));
            clErrchk(clEnqueueReadBuffer(commandQueue, throughputboxBuffer[kindex], CL_TRUE, 0, sizeof(Vec) *  bwidth * bheight, throughputs_boxes[kindex], 0, NULL, NULL));
            clErrchk(clEnqueueReadBuffer(commandQueue, specularBounceboxBuffer[kindex], CL_TRUE, 0, sizeof(char) *  bwidth * bheight, specularBounce_boxes[kindex], 0, NULL, NULL));
            clErrchk(clEnqueueReadBuffer(commandQueue, terminatedboxBuffer[kindex], CL_TRUE, 0, sizeof(char) *  bwidth * bheight, terminated_boxes[kindex], 0, NULL, NULL));
            clErrchk(clEnqueueReadBuffer(commandQueue, resultboxBuffer[kindex], CL_TRUE, 0, sizeof(Result) *  bwidth * bheight, results_boxes[kindex], 0, NULL, NULL));

            kindex++;
        }
    }
#else
	int index = 0;

	/* Set kernel arguments */
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&cameraBufferLeft));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&seedBuffer));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(short), (void *)&width));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(short), (void *)&height));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&rayBuffer));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&throughputBuffer));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&specularBounceBuffer));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&terminatedBuffer));
	clErrchk(clSetKernelArg(kernelGen, index++, sizeof(cl_mem), (void *)&resultBuffer));

	ExecuteKernel(kernelGen, width * height);
	clFinish(commandQueue);

	//clErrchk(clEnqueueReadBuffer(commandQueue, cameraBuffer, CL_TRUE, 0, sizeof(Camera), &camera, 0, NULL, NULL));
	//clErrchk(clEnqueueReadBuffer(commandQueue, seedBuffer, CL_TRUE, 0, sizeof(unsigned int) * width * height * 2, seeds, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, rayBuffer, CL_TRUE, 0, sizeof(Ray) *  width * height, ray, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, throughputBuffer, CL_TRUE, 0, sizeof(Vec) *  width * height, throughput, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, specularBounceBuffer, CL_TRUE, 0, sizeof(char) *  width * height, specularBounce, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, terminatedBuffer, CL_TRUE, 0, sizeof(char) *  width * height, terminated, 0, NULL, NULL));
	clErrchk(clEnqueueReadBuffer(commandQueue, resultBuffer, CL_TRUE, 0, sizeof(Result) *  width * height, result, 0, NULL, NULL));
#endif
#endif
}

void ReInitVR(const int reallocBuffers) {
	// Check if I have to reallocate buffers
	if (reallocBuffers) {
		FreeBuffers();
		UpdateCamera(true);
		AllocateBuffers();
	} else {
		UpdateCamera(true);
	}

	clErrchk(clEnqueueWriteBuffer(commandQueue, cameraBufferLeft, CL_TRUE, 0, sizeof(Camera), &cameraLeft, 0, NULL, NULL));
	clErrchk(clEnqueueWriteBuffer(commandQueue, cameraBufferRight, CL_TRUE, 0, sizeof(Camera), &cameraRight, 0, NULL, NULL));
#if (ACCELSTR == 1)
	clErrchk(clEnqueueWriteBuffer(commandQueue, btnBuffer, CL_TRUE, 0, sizeof(BVHNodeGPU) * (shapeCnt - 1), btn, 0, NULL, NULL));
	clErrchk(clEnqueueWriteBuffer(commandQueue, btlBuffer, CL_TRUE, 0, sizeof(BVHNodeGPU) * (shapeCnt), btl, 0, NULL, NULL));
#elif (ACCELSTR == 2)
	clErrchk(clEnqueueWriteBuffer(commandQueue, kngBuffer, CL_TRUE, 0, sizeof(KDNodeGPU) * (kngCnt), pkngbuf, 0, NULL, NULL));
	clErrchk(clEnqueueWriteBuffer(commandQueue, knBuffer, CL_TRUE, 0, sizeof(int) * (knCnt), pknbuf, 0, NULL, NULL));
#endif
	currentSampleLeft = 0;
	currentSampleRight = 0;
}

#ifdef WIN32
static int mouseX = 0, mouseY = 0;
static int mouseButton = 0;

#define TWO_PI 6.28318530717958647693f
#define PI_OVER_TWO 1.57079632679489661923f

#define MOVE_STEP 10.0f
#define ROTATE_STEP (1.f * M_PI / 180.f)

void idleFunc(void) {
	glutPostRedisplay();
}

void displayFunc(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glRasterPos2i(0, 0);
#ifdef CPU_PARTRENDERING
	static float rCPU = 1.0f;
	static bool first = true;

	DrawAllBoxes(160, 120, &rCPU, first);
	first = false;
#else
	DrawFrame();
#endif
	glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
	glutSwapBuffers();
}

void reshapeFunc(int newWidth, int newHeight) {
	width = newWidth;
	height = newHeight;

	glViewport(0, 0, width, height);
	glLoadIdentity();
	glOrtho(0.f, width - 1.f, 0.f, height - 1.f, -1.f, 1.f);

	ReInit(1);

	glutPostRedisplay();
}

/// gets the current mouse position and compares to the last one to check if
/// we need to change the pitch/yaw of the camera
/// it only changes the camera if the mouse button is pressed
void motionFunc(int x, int y) {
	int deltaX = mouseX - x;
	int deltaY = mouseY - y;

	if (deltaX != 0 || deltaY != 0) {
		// rotate the camera using pitch (nodding movement) and yaw (nonono movement)
		if (mouseButton == GLUT_LEFT_BUTTON) {
			camera.yaw += deltaX * 0.01;
			camera.yaw = camera.yaw - TWO_PI * floor(camera.yaw / TWO_PI);
			camera.pitch += -deltaY * 0.01;
			camera.pitch = clamp(camera.pitch, -PI_OVER_TWO, PI_OVER_TWO);
		}

		glutSetCursor(GLUT_CURSOR_CROSSHAIR);

		mouseX = x;
		mouseY = y;
		ReInit(0);
	}
	else {
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}

/// simply records the mouse state
void mouseFunc(int button, int state, int x, int y) {
	mouseButton = button;
	mouseX = x;
	mouseY = y;

	motionFunc(x, y);
}

void keyFunc(unsigned char key, int x, int y) {
	switch (key) {
	case 'p': {
		FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.

		if (!f) {
			fprintf(stderr, "Failed to open image file: image.ppm\n");
		}
		else {
			fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);

			int x, y;

			for (y = height - 1; y >= 0; --y) {
				unsigned char *p = (unsigned char *)(&pixels[y * width]);
				for (x = 0; x < width; ++x, p += 4)
					fprintf(f, "%d %d %d ", p[0], p[1], p[2]);
			}

			fclose(f);
		}
		break;
	}
	case 27: /* Escape key */
		fprintf(stderr, "Done.\n");
		exit(0);
		break;
	case ' ': /* Refresh display */
		ReInit(1);
		break;
	case 'a': {
		Vec dir = camera.x;
		vnorm(dir);
		vsmul(dir, -MOVE_STEP, dir);
		vadd(camera.orig, camera.orig, dir);
		vadd(camera.target, camera.target, dir);
		ReInit(0);
		break;
	}
	case 'd': {
		Vec dir = camera.x;
		vnorm(dir);
		vsmul(dir, MOVE_STEP, dir);
		vadd(camera.orig, camera.orig, dir);
		vadd(camera.target, camera.target, dir);
		ReInit(0);
		break;
	}
	case 'w': {
		Vec dir = camera.dir;
		vsmul(dir, MOVE_STEP, dir);
		vadd(camera.orig, camera.orig, dir);
		vadd(camera.target, camera.target, dir);
		ReInit(0);
		break;
	}
	case 's': {
		Vec dir = camera.dir;
		vsmul(dir, -MOVE_STEP, dir);
		vadd(camera.orig, camera.orig, dir);
		vadd(camera.target, camera.target, dir);
		ReInit(0);
		break;
	}
	case 'r':
		camera.orig.y += MOVE_STEP;
		camera.target.y += MOVE_STEP;
		ReInit(0);
		break;
	case 'f':
		camera.orig.y -= MOVE_STEP;
		camera.target.y -= MOVE_STEP;
		ReInit(0);
		break;
	default:
		break;
	}
}

void specialFunc(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_UP: {
		camera.pitch += -0.01;
		camera.pitch = clamp(camera.pitch, -PI_OVER_TWO, PI_OVER_TWO);
		ReInit(0);
		break;
	}
	case GLUT_KEY_DOWN: {
		camera.pitch += 0.01;
		camera.pitch = clamp(camera.pitch, -PI_OVER_TWO, PI_OVER_TWO);
		ReInit(0);
		break;
	}
	case GLUT_KEY_LEFT: {
		camera.yaw += 0.01;
		camera.yaw = camera.yaw - TWO_PI * floor(camera.yaw / TWO_PI);
		ReInit(0);
		break;
	}
	case GLUT_KEY_RIGHT: {
		camera.yaw += -0.01;
		camera.yaw = camera.yaw - TWO_PI * floor(camera.yaw / TWO_PI);
		ReInit(0);
		break;
	}
	case GLUT_KEY_PAGE_UP:
		camera.target.y += MOVE_STEP;
		ReInit(0);
		break;
	case GLUT_KEY_PAGE_DOWN:
		camera.target.y -= MOVE_STEP;
		ReInit(0);
		break;
	default:
		break;
	}
}

void InitGlut(int argc, char *argv[], char *windowTittle) {
	glutInitWindowSize(width, height);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInit(&argc, argv);

	glutCreateWindow(windowTittle);

	glutReshapeFunc(reshapeFunc);
	glutMouseFunc(mouseFunc);
	glutMotionFunc(motionFunc);
	glutKeyboardFunc(keyFunc);
	glutSpecialFunc(specialFunc);
	glutDisplayFunc(displayFunc);
	glutIdleFunc(idleFunc);

	glViewport(0, 0, width, height);
	glLoadIdentity();
	glOrtho(0.f, width - 1.f, 0.f, height - 1.f, -1.f, 1.f);
}

int main(int argc, char *argv[]) {
    amiSmallptCPU = 0;

    if (argc == 7) {
		bool walllight = true;
		srand(time(NULL));

        useGPU = atoi(argv[1]);
        forceWorkSize = atoi(argv[2]);
		//strcpy(forceWorkSize, argv[2]);
#ifndef EXP_KERNEL
		strcpy(kernelFileName, argv[3]);
#endif
        width = atoi(argv[4]);
        height = atoi(argv[5]);

		Reaad(argv[6], &walllight);
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
		UpdateCamera();
    } else if (argc == 1) {
		srand(time(NULL));

		shapeCnt = sizeof(CornellSpheres) / sizeof(Sphere);
		shapes = (Shape *)malloc(sizeof(Shape) * shapeCnt);

		for(int i = 0; i < shapeCnt; i++)
		{
			shapes[i].type = SPHERE;
			shapes[i].s = CornellSpheres[i];
			shapes[i].e = t[i].e;
			shapes[i].c = t[i].c;
			shapes[i].refl = t[i].refl;
		}

        vinit(camera.orig, 50.f, 45.f, 205.6f);
        vinit(camera.target, 50.f, 45 - 0.042612f, 204.6);
	}
	else
	{
		LOGE("Usage: %s\n", argv[0]);
		LOGE("Usage: %s <use CPU/GPU device (0=CPU or 1=GPU)> <workgroup size (0=default value or anything > 0 and power of 2)> <kernel file name> <window width> <window height> <scene file>\n", argv[0]);

		exit(-1);
	}

    /*------------------------------------------------------------------------*/
	InitGlut(argc, argv, (char *)"SmallPTGPU (Added the BVH and KDTree for the intersection tests)");

	LOGI("Acceleration for intersection test: ");
#if (ACCELSTR == 0)
	LOGI("No Accel\n");
#elif (ACCELSTR == 1)
	LOGI("BVH\n");
#elif (ACCELSTR == 2)
	LOGI("KDTree\n");
#endif

	glutMainLoop();

    return 0;
}
#endif

void AddWallLight()
{
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { WALL_RAD,{ WALL_RAD + 25.0f, 0.0f, 0.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .75f, .25f, .25f }; shapes[shapeCnt++].refl = DIFF; /* Left */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { WALL_RAD,{ -WALL_RAD - 25.0f, 0.0f, 0.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .25f, .25f, .75f }; shapes[shapeCnt++].refl = DIFF; /* Rght */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { WALL_RAD,{ 0.0f, 0.0f, WALL_RAD - 25.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .75f, .75f, .75f }; shapes[shapeCnt++].refl = DIFF; /* Back */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { WALL_RAD,{ 0.0f, 0.0f, -WALL_RAD + 100.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { 0.f, 0.f, 0.f }; shapes[shapeCnt++].refl = DIFF; /* Frnt */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { WALL_RAD,{ 0.0f, WALL_RAD + 25.0f, 0.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .75f, .75f, .75f }; shapes[shapeCnt++].refl = DIFF; /* Botm */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { WALL_RAD,{ 0.0f, -WALL_RAD - 25.0f, 0.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .75f, .75f, .75f }; shapes[shapeCnt++].refl = DIFF; /* Top */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { 5.0f,{ 10.0f, -17.0f, 0.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .9f, .9f, .9f }; shapes[shapeCnt++].refl = SPEC; /* Mirr */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { 5.0f,{ -10.0f, -17.0f, 0.0f } };  shapes[shapeCnt].e = { 0.f, 0.f, 0.f }; shapes[shapeCnt].c = { .9f, .9f, .9f }; shapes[shapeCnt++].refl = REFR; /* Glas */
    shapes[shapeCnt].type = SPHERE; shapes[shapeCnt].s = { 2.f,{ 10.0f, 15.0f, 0.0f } };  shapes[shapeCnt].e = { 12.f, 12.f, 12.f }; shapes[shapeCnt].c = { 0.f, 0.f, 0.f }; shapes[shapeCnt++].refl = DIFF; /* Lite */
}

#if (ACCELSTR == 1)
void BuildBVH()
{
	CLBVH *pCB = new CLBVH(shapes, shapeCnt, commandQueue, context, kernelRad, kernelBvh, kernelOpt);

	pCB->buildRadixTree();
	pCB->buildBVHTree();
	//pCB->optimize();

	pCB->getTrees(&btn, &btl);
	/*
    pCB->makeNaiveBVHTree();
    pCB->getTree(&nbtn, &nbtnCnt);
    */
}
#elif (ACCELSTR == 2)
void BuildKDtree()
{
    std::vector<Shape *> vs;

    for (int i = 0; i < shapeCnt; i++) {
        shapes[i].index = i;
        //shapes[i].morton_code = 0;

        vs.push_back(&shapes[i]);
    }

    KDTree *kdTree = new KDTree();
    KDTreeNode *rootNode = kdTree->build(vs, 0);
    //kdTree->printNode(rootNode, 0);

    kdTree->getTrees(rootNode, &pkngbuf, &kngCnt, &pknbuf, &knCnt);
}
#endif