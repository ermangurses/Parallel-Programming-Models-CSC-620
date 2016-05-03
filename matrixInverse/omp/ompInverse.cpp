//============================================================================
// Name        : inverse.cpp
// Author      : Erman Gurses
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>

#define lookup2D(ptr,x,y,blockY) (ptr)[ (x) * (blockY) + (y)]

#define lookup3D(ptr,i,j,k,Nx,Ny,Nz) &(ptr)[( (Nx) * (Ny) * (k)   +            \
                                              (Nx) *  (j) + (i) ) * blockSize ]

#define lookup5D(ptr,i,j,k,Nx,Ny,Nz,x,y,blockY) (ptr)[((Nx) *  (Ny) * (k)   +            \
                                                       (Nx) *  (j) + (i)) * blockSize +  \
                                                       (x) *  (blockY) + (y)]

using namespace std;

void inverse(double *A, double *B, int neqs);
void prt (double *A, int neqs);

int main() {
  int k, j, i;
  int Nx = 50;
  int Ny = 50;
  int Nz = 50;
  int blockX= 3;
  int blockY= 3;
  int blockSize = blockX * blockY;

  int requiredSpaceBlock = sizeof(double) * blockSize;

  int requiredSpace = requiredSpaceBlock * Nx * Ny * Nz;


  double *B,*M;

  M = (double*) malloc(requiredSpace);
  B = (double*) malloc(requiredSpace);
 
  omp_set_num_threads(8); 
 
  #pragma omp parallel for private(k,j,i) schedule(static) collapse(3)
  for( k = 0; k < Nz; k++){
    for( j = 0; j < Ny; j++){
      for( i = 0; i < Nx; i++){
              lookup5D(M,i,j,k,Nx,Ny,Nz,0,0,blockY) =  1;
              lookup5D(M,i,j,k,Nx,Ny,Nz,0,1,blockY) = -1;
              lookup5D(M,i,j,k,Nx,Ny,Nz,0,2,blockY) = -1;
              lookup5D(M,i,j,k,Nx,Ny,Nz,1,0,blockY) = -1;
              lookup5D(M,i,j,k,Nx,Ny,Nz,1,1,blockY) =  2;
              lookup5D(M,i,j,k,Nx,Ny,Nz,1,2,blockY) =  3;
              lookup5D(M,i,j,k,Nx,Ny,Nz,2,0,blockY) =  1;
              lookup5D(M,i,j,k,Nx,Ny,Nz,2,1,blockY) =  1;
              lookup5D(M,i,j,k,Nx,Ny,Nz,2,2,blockY) =  4;
      }
    }
  } 
 double start_time = omp_get_wtime();
   
  #pragma omp parallel for shared(start_time, Nx, Ny, Nz, M, B) private(k,j,i) schedule(static) collapse(3)
  for( k = 0; k < Nz; k++){
    for( j = 0; j < Ny; j++){
      for( i = 0; i < Nx; i++){
        inverse(lookup3D(M,i,j,k,Nx,Ny,Nz),lookup3D(B,i,j,k,Nx,Ny,Nz),blockY); 
      }
    }
  }

  double end_time = omp_get_wtime();
  double time =  (end_time - start_time);
  printf("\nElapsedTime:%f\n",time);

/*   for(int k = 0; k < Nz; k++){
     for(int j = 0; j < Ny; j++){
       for(int i = 0; i < Nx; i++){

              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,0,0,blockY));
              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,0,1,blockY));
              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,0,2,blockY));

              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,1,0,blockY));
              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,1,1,blockY));
              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,1,2,blockY));

              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,2,0,blockY));
              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,2,1,blockY));
              printf("%f  ",lookup5D(B,i,j,k,Nx,Ny,Nz,2,2,blockY));
              printf("\n");
        }
      }
    }
*/
 return 0;
}

void inverse(double *A, double *B, int blockY){

  double * temp, *C;
  temp = (double*) malloc(blockY);
  C    = (double*) malloc(blockY*blockY);

  for (int j = 0;  j < blockY; j++){
    for ( int i = 0; i < blockY; i++){
    	lookup2D(B,j,i,blockY) = 0,0;
        lookup2D(C,j,i,blockY) =lookup2D(A,j,i,blockY);
    }
    lookup2D(B,j,j,blockY) = 1.0;
  }

  for (int k = 0;  k < blockY-1; k++){
    for (int j = 0;  j < blockY; j++){
       lookup2D(B,k,j,blockY) = lookup2D(B,k,j,blockY) /
				lookup2D(C,k,k,blockY);
	}
	for (int j = 0;  j < blockY; j++){
		temp[j] = lookup2D(C,k,k,blockY);
	}
	for (int j = 0;  j < blockY; j++){
		lookup2D(C,k,j,blockY) = lookup2D(C,k,j,blockY) / temp[j];
	}
	for (int i = k+1;  i < blockY; i++){
	  for (int l = 0;  l < blockY; l++){
	    lookup2D(B,i,l,blockY) = lookup2D(B,i,l,blockY) -
				     lookup2D(C,i,k,blockY) *
				     lookup2D(B,k,l,blockY);

	    temp[l] = lookup2D(C,i,l,blockY) -
        	      lookup2D(C,i,k,blockY) *
        	      lookup2D(C,k,l,blockY);
	  }//l
          for (int l = 0;  l < blockY; l++){
    	    lookup2D(C,i,l,blockY)   = temp[l];
	  }//l
	} //i
   }//k
  for (int i = 0;  i < blockY; i++){
	  lookup2D(B,blockY-1,i,blockY)  =  lookup2D(B,blockY-1,i,blockY) /
			                    lookup2D(B,blockY-1,blockY-1,blockY);

	  lookup2D(C,blockY-1,i,blockY)  =  lookup2D(C,blockY-1,i,blockY) /
			                    lookup2D(C,blockY-1,blockY-1,blockY);
  }

  for (int k = blockY-1;  k >= 1; k--){
    for (int i = 0;  i <=k-1; i++){
      for (int l = 0;  l < blockY; l++){
        lookup2D(B,i,l,blockY) = lookup2D(B,i,l,blockY) -
				 lookup2D(C,i,k,blockY) *
			         lookup2D(B,k,l,blockY);
        temp[l] = lookup2D(C,i,l,blockY) -
                  lookup2D(C,i,k,blockY) *
	          lookup2D(C,k,l,blockY);
      }//l
      for (int l = 0;  l < blockY; l++){
    	  lookup2D(C,i,l,blockY) = temp[l];
      }//l
    }//i
  }//k
  //free(temp);
}
void prt (double *A, int blockY){

  printf("\n");
  for (int i = 0;  i < blockY; i++){
    for ( int j = 0; j < blockY; j++){
      printf("%5.2f   ",lookup2D(A,i,j,blockY));
    }
    printf("\n");
  }
  printf("\n");

}
