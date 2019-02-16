//============================================================================
// Name        : inverse.cpp
// Author      : Erman Gurses
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cuda_runtime.h>

#define lookup2D(ptr,x,y,matY) (ptr)[ (x) * (matY) + (y)]

#define lookup2DK(ptr,threadId,x,y,matY) (ptr)[ (threadId)*(matY)*(matY)+ (x) * (matY) + (y)]

#define lookup3D(ptr,i,j,k,Nx,Ny,Nz) &(ptr)[( (Nx) * (Ny) * (k)   +           \
                                              (Nx) * (j ) + (i) ) * matY*matY ]

#define lookup5D(ptr,i,j,k,Nx,Ny,Nz,x,y,matY) (ptr)[((Nx) * (Ny)    * (k)   +            \
                                                     (Nx) *  (j)    + (i))  * matSize +  \
                                                      (x) *  (matY) + (y)]

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
   if (code != cudaSuccess){
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ int getGlobalIdx_3D_3D(){ 

int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z   ) + 
                         (threadIdx.z * (blockDim.x * blockDim.y)) + 
                         (threadIdx.y * blockDim.x) + threadIdx.x;

  return threadId; 
}

__global__ void inverse(double *A, double *B, double *C, int matY){

int blockId = blockIdx.x + blockIdx.y * gridDim.x + 
               gridDim.x * gridDim.y * blockIdx.z;

int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z) + 
                         (threadIdx.z * (blockDim.x * blockDim.y)) + 
                         (threadIdx.y * blockDim.x) + threadIdx.x;  


//printf("threadId %d\n",threadId);
  
 __shared__ double temp[3];
      
  for (int j = 0; j < matY; j++){
    for (int i = 0; i < matY; i++){
      lookup2DK(B,threadId,j,i,matY) = 0,0;
      lookup2DK(C,threadId,j,i,matY) =lookup2DK(A,threadId,j,i,matY);
    }
    lookup2DK(B,threadId,j,j,matY) = 1.0;
  }

  for (int k = 0; k < matY-1; k++){
    for (int j = 0; j < matY; j++){
      lookup2DK(B,threadId,k,j,matY) = lookup2DK(B,threadId,k,j,matY) /
			               lookup2DK(C,threadId,k,k,matY);
    }
    for (int j = 0; j < matY; j++){
      temp[j] = lookup2DK(C,threadId,k,k,matY);
    }
    for (int j = 0; j < matY; j++){
      lookup2DK(C,threadId,k,j,matY) = lookup2DK(C,threadId,k,j,matY) / temp[j];
    }
    for (int i = k + 1; i < matY; i++){
      for (int l = 0; l < matY; l++){
        lookup2DK(B,threadId,i,l,matY) = 
                  lookup2DK(B,threadId,i,l,matY) -
	          lookup2DK(C,threadId,i,k,matY) * 
                  lookup2DK(B,threadId,k,l,matY);

         temp[l] = lookup2DK(C,threadId,i,l,matY) -
                   lookup2DK(C,threadId,i,k,matY) * 
                   lookup2DK(C,threadId,k,l,matY);
      }//l
      for (int l = 0; l < matY; l++){
        lookup2DK(C,threadId,i,l,matY)   = temp[l];
      }//l
    } //i
  }//k
  for (int i = 0; i < matY; i++){
    lookup2DK(B,threadId,matY-1,i,matY) = lookup2DK(B,threadId,matY-1,i,matY) /
			                  lookup2DK(B,threadId,matY-1,matY-1,matY);

    lookup2DK(C,threadId,matY-1,i,matY) = lookup2DK(C,threadId,matY-1,i,matY) /
			                  lookup2DK(C,threadId,matY-1,matY-1,matY);
  }

  for (int k = matY-1; k >= 1; k--){
    for (int i = 0; i <= k-1; i++){
      for (int l = 0; l < matY; l++){
        lookup2DK(B,threadId,i,l,matY) = 
                  lookup2DK(B,threadId,i,l,matY) - 
                  lookup2DK(C,threadId,i,k,matY) * 
                  lookup2DK(B,threadId,k,l,matY);

        temp[l] = lookup2DK(C,threadId,i,l,matY) -
                  lookup2DK(C,threadId,i,k,matY) *
	          lookup2DK(C,threadId,k,l,matY);
      }//l
      for (int l = 0; l < matY; l++){
        lookup2DK(C,threadId,i,l,matY) = temp[l];
      }//l
    }//i
  }//k 
}

int main() {
  int k, j, i;
  int Nx = 4;
  int Ny = 4;
  int Nz = 4;
  int matX= 3;
  int matY= 3;
  int matSize = matX * matY;
  int requiredSpaceBlock = sizeof(double) * matSize;
  int requiredSpace = requiredSpaceBlock * Nx * Ny * Nz;
  double *B, *M, *C, *B_Device, *M_Device, *C_Device;

  M = (double*) malloc(requiredSpace);
  B = (double*) malloc(requiredSpace);  
  C = (double*) malloc(requiredSpace);

  for(k = 0; k < Nz; k++){
    for(j = 0; j < Ny; j++){
      for(i = 0; i < Nx; i++){
        lookup5D(M,i,j,k,Nx,Ny,Nz,0,0,matY) =  1;
        lookup5D(M,i,j,k,Nx,Ny,Nz,0,1,matY) = -1;
        lookup5D(M,i,j,k,Nx,Ny,Nz,0,2,matY) = -1;
        lookup5D(M,i,j,k,Nx,Ny,Nz,1,0,matY) = -1;
        lookup5D(M,i,j,k,Nx,Ny,Nz,1,1,matY) =  2;
        lookup5D(M,i,j,k,Nx,Ny,Nz,1,2,matY) =  3;
        lookup5D(M,i,j,k,Nx,Ny,Nz,2,0,matY) =  1;
        lookup5D(M,i,j,k,Nx,Ny,Nz,2,1,matY) =  1;   
        lookup5D(M,i,j,k,Nx,Ny,Nz,2,2,matY) =  4;
      }
    }
  } 
   
  gpuErrchk(cudaMalloc((void**)&B_Device,requiredSpace));
  gpuErrchk(cudaMalloc((void**)&M_Device,requiredSpace));
  gpuErrchk(cudaMalloc((void**)&C_Device,requiredSpace));
  
  gpuErrchk(cudaMemcpy(B_Device,B,requiredSpace,cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(M_Device,M,requiredSpace,cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(C_Device,C,requiredSpace,cudaMemcpyHostToDevice));
              
  int blockX = 2; 
  int blockY = 2; 
  int blockZ = 2; 

  int gridX = ceil(Nx/blockX); 
  int gridY = ceil(Ny/blockY);
  int gridZ = ceil(Nz/blockZ);
   
  dim3 blockSize(blockX,blockY,blockZ);
  dim3 gridSize(gridX,gridY,gridZ);
 
  // Kernel Call       
  inverse<<<gridSize,blockSize>>>(M_Device,B_Device,C_Device,matY);
 
  gpuErrchk( cudaPeekAtLastError() ); 
  gpuErrchk(cudaDeviceSynchronize());
      
  gpuErrchk(cudaMemcpy(B,B_Device,requiredSpace,cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(M,M_Device,requiredSpace,cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(C,C_Device,requiredSpace,cudaMemcpyDeviceToHost));
  

  for(int k = 0; k < Nz; k++){
    for(int j = 0; j < Ny; j++){
      for(int i = 0; i < Nx; i++){
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,0,0,matY));
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,0,1,matY));
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,0,2,matY));
        printf("\n");

        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,1,0,matY));
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,1,1,matY));
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,1,2,matY));
        printf("\n");

        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,2,0,matY));
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,2,1,matY));
        printf("%9f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,2,2,matY));
        printf("\n\n");
      }
    }
  }

 return 0;
}

