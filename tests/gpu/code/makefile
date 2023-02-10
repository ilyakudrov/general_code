INCLUDE_C=../../../lib/cpu/include
SRC_C=../../../lib/cpu/src
INCLUDE_CUDA=../../../lib/gpu/include
SRC_CUDA=../../../lib/gpu/src
OBJ=../../../lib/gpu/obj
CXX=g++
CUDA=nvcc
CXXFLAGS=-Ofast,-fopenmp,-march=native,-mtune=native
FLAGS_CUDA=-arch=sm_61 -use_fast_math -O3 --ptxas-options=-v -Xptxas -dlcm=cg
CLIME_PATH=../../../lib/cpu/src/c-lime
ARCH=test
OBJFILES_C=$(OBJ)/link_$(ARCH).o $(OBJ)/data_$(ARCH).o $(OBJ)/decomposition_$(ARCH).o $(OBJ)/monopoles_$(ARCH).o $(OBJ)/loop_$(ARCH).o $(OBJ)/matrix_$(ARCH).o
OBJFILES_CLIME=$(OBJ)/lime_fseeko_$(ARCH).o  $(OBJ)/lime_header_$(ARCH).o  $(OBJ)/lime_reader_$(ARCH).o  $(OBJ)/lime_utils_$(ARCH).o  $(OBJ)/lime_writer_$(ARCH).o
# OBJFILES_CUDA=$(OBJ)/observables/plaket.o $(OBJ)/general/matrix.o

TARGETS=plaket

$(TARGETS): %: %.cu $(OBJFILES) $(OBJFILES_C) $(OBJFILES_CUDA) $(OBJFILES_CLIME)
	$(CUDA) -ccbin $(CXX) -Xcompiler $(CXXFLAGS) $(FLAGS_CUDA) $@.cu $(OBJFILES_C) $(OBJFILES_CUDA) $(OBJFILES_CLIME) -o $@_$(ARCH)

$(OBJFILES_C): $(OBJ)/%_$(ARCH).o: $(SRC_C)/%.cpp
	$(CUDA) -ccbin $(CXX) -Xcompiler $(CXXFLAGS) $(FLAGS_CUDA) -x cu -dc -c $< -o $@

$(OBJFILES_CLIME): $(OBJ)/%_$(ARCH).o: $(CLIME_PATH)/%.c
	$(CUDA) -ccbin $(CXX) -Xcompiler $(CXXFLAGS) $(FLAGS_CUDA) -x cu -dc -c $< -o $@

$(OBJFILES_CUDA): $(OBJ)/%.o : $(SRC_CUDA)/%.cu
	@mkdir -p $(@D)
	$(CUDA) -ccbin $(CXX) -Xcompiler $(CXXFLAGS) $(FLAGS_CUDA) -x cu -dc -c $< -o $@

clean:
	rm -rf $(OBJ)/*.o
