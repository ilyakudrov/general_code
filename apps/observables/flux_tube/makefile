INCLUDE_PATH=../../../lib/cpu/include
OBJ_PATH=../../../lib/cpu/obj
GCC=g++
MATRIX_PLAKET=su2
MATRIX_WILSON=su2
MATRIX=su2
FLAGS=-Ofast -fopenmp -march=native -mtune=native -std=c++17
C_PATH=../../../lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJ_FILES=$(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/polyakov_loops_$(ARCH).o $(OBJ_PATH)/smearing_$(ARCH).o \
	$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o $(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/flux_tube_$(ARCH).o $(OBJ_PATH)/mag_$(ARCH).o
HEADERFILES=$(INCLUDE_PATH)/matrix.h $(INCLUDE_PATH)/indexing.h
EIGH_PATH="/home/ilya/soft/source/eigen-master"

flux_tube_wilson: flux_tube_wilson.cpp $(OBJ_FILES) $(HEADERFILES)
	$(GCC) flux_tube_wilson.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -DMATRIX_PLAKET=$(MATRIX_PLAKET) -DMATRIX_WILSON=$(MATRIX_WILSON) -I $(EIGH_PATH) -o flux_tube_wilson_$(MATRIX_PLAKET)_$(MATRIX_WILSON)_$(ARCH)
flux_tube_schwinger: flux_tube_schwinger.cpp $(OBJ_FILES) $(HEADERFILES)
	$(GCC) flux_tube_schwinger.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -DMATRIX=$(MATRIX) -I $(EIGH_PATH) -o flux_tube_schwinger_$(MATRIX)_$(ARCH)