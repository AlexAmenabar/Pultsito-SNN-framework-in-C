#include <stdio.h>
#include <cuda_runtime.h>

#include "snn_library.h"
#include "cuda/cuda_helpers.cuh"



extern "C" int get_memory_info(cuda_info_t *cuda_info, int dev) {
    
    size_t free_mem, total_mem;

    // Get memory info from the currently active device
    cudaError_t err = cudaMemGetInfo(&free_mem, &total_mem);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
        return 1;
    }

    printf("GPU memory usage:\n");
    printf("  Total memory: %.2f MB (%.2f GB)\n", total_mem / 1024.0 / 1024.0, total_mem / 1024.0 / 1024.0 / 1024.0);
    printf("  Free memory : %.2f MB\n", free_mem / 1024.0 / 1024.0);
    printf("  Used memory : %.2f MB\n", (total_mem - free_mem) / 1024.0 / 1024.0);


    // store free memory in device
    cuda_info->gpu_mem[dev] = (double)free_mem;
    cuda_info->gpu_usable_mem = (double)free_mem;  

    return 0;
}


extern "C" cuda_info_t* getProperties(){

    cuda_info_t *cuda_info;
    int nDevices;


    // allocate memory to store cuda info
    cuda_info = (cuda_info_t*)malloc(sizeof(cuda_info_t));
    cudaGetDeviceCount(&nDevices);
  
    // store the number of devices
    cuda_info->nDevices = nDevices;
    cuda_info->gpu_mem = (double*)malloc(nDevices * sizeof(double));


    printf("Number of devices: %d\n", nDevices);
    
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Memory Clock Rate (MHz): %d\n",
            prop.memoryClockRate/1024);
        printf("  Memory Bus Width (bits): %d\n",
            prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): %.1f\n",
            2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
        printf("  Total global memory (Gbytes) %.1f\n",(float)(prop.totalGlobalMem)/1024.0/1024.0/1024.0);
        printf("  Shared memory per block (Kbytes) %.1f\n",(float)(prop.sharedMemPerBlock)/1024.0);
        printf("  minor-major: %d-%d\n", prop.minor, prop.major);
        printf("  Warp-size: %d\n", prop.warpSize);
        printf("  Concurrent kernels: %s\n", prop.concurrentKernels ? "yes" : "no");
        printf("  Concurrent computation/communication: %s\n\n",prop.deviceOverlap ? "yes" : "no");

        // set i device as active and get memory info
        cudaSetDevice(i);
        get_memory_info(cuda_info, i);
    }

    cudaSetDevice(0);

    return cuda_info;
}