INCLUDE_PATH=../../lib/cpu/include
OBJ_PATH=../../lib/cpu/obj
GCC=g++
C_PATH=../../lib/cpu/src/c-lime
ARCH=test
MATRIX=su2
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o \
	$(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o $(OBJ_PATH)/smearing_$(ARCH).o
HEADERFILES=$(INCLUDE_PATH)/matrix.h $(INCLUDE_PATH)/indexing.h
EIGH_PATH="/home/ilya/soft/source/eigen-master"

polyakov_loop: polyakov_loop.cpp $(OBJFILES) $(HEADERFILES)
	$(GCC) $@.cpp $(OBJFILES) $(C_LIME) -DMATRIX=$(MATRIX) -I $(EIGH_PATH) -fopenmp -o polyakov_loop_$(MATRIX)_$(ARCH)
