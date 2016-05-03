CXX=g++
COMMON_DIR=../common

OPTIONS= -O3 -fopenmp 


VARIANT = inverse




$(VARIANT): % : %.cpp 
	$(CXX)  $^ -o $@ $(OPTIONS)
run:
	./inverse
clean:
	rm inverse
