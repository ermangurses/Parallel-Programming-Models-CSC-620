NVCC=nvcc

#CUDA code generation flags
GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30
GENCODE_SM35    := -gencode arch=compute_35,code=sm_35
GENCODE_FLAGS   := $(GENCODE_SM20) $(GENCODE_SM30) $(GENCODE_SM35)

OPTIONS_CUDA   := -O4 -ftz=false  $(GENCODE_FLAGS)

all: cudaInverse.cu 
	$(NVCC) $^ -o  $@ $(OPTIONS_CUDA)
