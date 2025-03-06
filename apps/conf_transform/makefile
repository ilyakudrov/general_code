INCLUDE_PATH=../../lib/cpu/include
OBJ_PATH=../../lib/cpu/obj
GCC=g++
MATRIX=su2
FLAGS=-Ofast -fopenmp -march=native -mtune=native -std=c++17 -DMATRIX=$(MATRIX)
C_PATH=../../lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c $(C_PATH)/lime_header.c $(C_PATH)/lime_reader.c $(C_PATH)/lime_utils.c $(C_PATH)/lime_writer.c
OBJ_FILES=$(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o \
	$(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o
HEADERFILES=$(INCLUDE_PATH)/matrix.h $(INCLUDE_PATH)/indexing.h
EIGH_PATH="/home/ilya/soft/source/eigen-master"

transform_to_qc2dstag: transform_to_qc2dstag.cpp $(OBJ_FILES) $(HEADERFILES)
	$(GCC) transform_to_qc2dstag.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -I $(EIGH_PATH) -o transform_to_qc2dstag_$(MATRIX)_$(ARCH)