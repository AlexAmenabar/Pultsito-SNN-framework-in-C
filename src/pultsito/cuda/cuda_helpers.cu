#include <stdio.h>
#include <cuda_runtime.h>

#include "snn_library.h"
#include "cuda/cuda_helpers.cuh"


// function to get the memory available in the GPU
extern "C" int get_memory_info(cuda_info_t *cuda_info, int dev) {
    
    size_t free_mem, total_mem;

    // Get memory info from the currently active device
    cudaError_t err = cudaMemGetInfo(&free_mem, &total_mem);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
        return 1;
    }

    printf("  > GPU memory usage:\n");
    printf("    > Total memory: %.2f MB (%.2f GB)\n", total_mem / 1024.0 / 1024.0, total_mem / 1024.0 / 1024.0 / 1024.0);
    printf("    > Free memory : %.2f MB\n", free_mem / 1024.0 / 1024.0);
    printf("    > Used memory : %.2f MB\n", (total_mem - free_mem) / 1024.0 / 1024.0);


    // store free memory in device
    cuda_info->gpu_total_mem[dev] = (double)total_mem;
    cuda_info->gpu_free_mem[dev] = (double)free_mem;  

    return 0;
}


extern "C" cuda_info_t* getGPUProperties(){

    // define varibales
    cuda_info_t *cuda_info;
    int nDevices;

    // allocate memory to store cuda data
    cuda_info = (cuda_info_t*)malloc(sizeof(cuda_info_t));
    cudaGetDeviceCount(&nDevices); // get number of available devices
  
    // store the number of devices
    cuda_info->nDevices = (size_t)nDevices;
    cuda_info->gpu_total_mem = (double*)malloc(nDevices * sizeof(double));
    cuda_info->gpu_free_mem = (double*)malloc(nDevices * sizeof(double));
    cuda_info->shared_memory_mem = (double*)malloc(nDevices * sizeof(double));


    printf("Number of devices: %d\n", nDevices);
    
    // loop over devices and print information
    for (size_t i = 0; i < (size_t)nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf(" > Device Number: %zu\n", i);
        printf("  > Device name: %s\n", prop.name);
        printf("  > Total global memory (Gbytes) %.1f\n",(float)(prop.totalGlobalMem)/1024.0/1024.0/1024.0);
        printf("  > Shared memory per block (Kbytes) %.1f\n",(float)(prop.sharedMemPerBlock)/1024.0);
        printf("  > minor-major: %d-%d\n", prop.minor, prop.major);
        printf("  > Warp-size: %d\n", prop.warpSize);
        printf("  > Concurrent kernels: %s\n", prop.concurrentKernels ? "yes" : "no");
        printf("  > Max Threads per Multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);
        printf("  > Multiprocessor count: %d\n", prop.multiProcessorCount);
        printf("  > Max blocks dim x: %d\n", prop.maxGridSize[0]);
        printf("  > Max blocks dim y: %d\n", prop.maxGridSize[1]);
        printf("  > Max blocks dim z: %d\n", prop.maxGridSize[2]);
        printf("  > Max thredas dim x: %d\n", prop.maxThreadsDim[0]);
        printf("  > Max thredas dim y: %d\n", prop.maxThreadsDim[1]);
        printf("  > Max thredas dim z: %d\n", prop.maxThreadsDim[2]);
        printf("  > Max threads per block: %d\n", prop.maxThreadsPerBlock);

        // set i device as active and get memory info
        cudaSetDevice(i);
        get_memory_info(cuda_info, i);

        // store some data // TODO: in the future more data should be stored
        cuda_info->shared_memory_mem[i] = prop.sharedMemPerBlock; // Bytes
    }

    cudaSetDevice(0);

    return cuda_info;
}