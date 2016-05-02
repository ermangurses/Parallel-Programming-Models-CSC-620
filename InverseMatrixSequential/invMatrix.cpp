//============================================================================
// Name        : inverse.cpp
// Author      : Erman Gurses
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#define lookup2D(ptr,x,y,blockY) (ptr)[ (x) * (blockY) + (y)]

#define lookup3D(ptr,i,j,k,Nx,Ny,Nz,x,y) (ptr)[((Nx) * (Ny) * (k)   +            \
                                                (Nx) *  (j) + (i)) * blockSize ]

#define lookup5D(ptr,i,j,k,Nx,Ny,Nz,x,y,blockY) (ptr)[((Nx) *  (Ny) * (k)   +            \
                                                       (Nx) *  (j) + (i)) * blockSize +  \
                                                        (x) *  (blockY) + (y)]

using namespace std;

void inverse(double *A, double *B, double *C, int neqs);
void prt (double *A, int neqs);

int main() {

  int Nx = 2;
  int Ny = 3;
  int Nz = 4;
  int blockX= 3;
  int blockY= 3;
  int blockSize = blockX * blockY;

  int requiredSpaceBlock = sizeof(double) * blockSize;

  int requiredSpace = requiredSpaceBlock * Nx * Ny * Nz;


  double *A,*B,*C,*M;

  M = (double*) malloc(requiredSpace);
  A = (double*) malloc(requiredSpaceBlock);
  B = (double*) malloc(requiredSpaceBlock);
  C = (double*) malloc(requiredSpaceBlock);


  /* for(int k = 0; k < Nz; k++){
	   for(int j = 0; j < Ny; j++){
		   for(int i = 0; i < Nx; i++){
			   for(int y = 0; y < blockY; y++){
				   for(int x = 0; x < blockX; x++){
					   lookup5D(M,i,j,k,Nx,Ny,Nz,x,y,blockY) = k;
				   }
			   }
		   }
	   }
   }*/

   for(int k = 0; k < Nz; k++){
	   for(int j = 0; j < Ny; j++){
		   for(int i = 0; i < Nx; i++){
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
   for(int k = 0; k < Nx*Ny*Nz*blockX*blockY; k++){
		   printf("%3.f   ",M[k]);
		   if((k+1)%9==0){
			   printf("\n");
		   }
   }

  lookup2D(A,0,0,blockY) =   1;
  lookup2D(A,0,1,blockY) =  -1;
  lookup2D(A,0,2,blockY) =  -1;
  lookup2D(A,1,0,blockY) =  -1;
  lookup2D(A,1,1,blockY) =   2;
  lookup2D(A,1,2,blockY) =   3;
  lookup2D(A,2,0,blockY) =   1;
  lookup2D(A,2,1,blockY) =   1;
  lookup2D(A,2,2,blockY) =   4;

  inverse(A,B,C,blockY);

   free (A);
   free (B);
   free (C);
   free (M);

 return 0;
}

void inverse(double *A, double *B, double *C, int blockY){

  double * temp;
  temp = (double*) malloc(blockY);

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
  prt(B,blockY);
  prt(C,blockY);
}
void prt (double *A, int blockY){

  printf("\n");
  for (int i = 0;  i < 3; i++){
    for ( int j = 0; j < 3; j++){
      printf("%5.2f   ",lookup2D(A,i,j,blockY));
    }
    printf("\n");
  }
  printf("\n");

}

