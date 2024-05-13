OBJ_PATH=../../lib/cpu/obj
GCC=g++
MATRIX=su2
FLAGS=-Ofast -fopenmp -march=native -mtune=native -DMATRIX=$(MATRIX)
C_PATH=../../lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c $(C_PATH)/lime_header.c $(C_PATH)/lime_reader.c $(C_PATH)/lime_utils.c $(C_PATH)/lime_writer.c
OBJ_FILES=$(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/matrix_$(ARCH).o $(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o \
	$(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o

transform_to_qc2dstag: transform_to_qc2dstag.cpp $(OBJ_FILES)
	$(GCC) transform_to_qc2dstag.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -o transform_to_qc2dstag_$(MATRIX)_$(ARCH)