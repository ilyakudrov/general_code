OBJ_PATH=../../lib/cpu/obj
GCC=g++
C_PATH=../../lib/cpu/src/c-lime
ARCH=test
MATRIX=su2
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/matrix_$(ARCH).o \
	$(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o

polyakov_loop_correlator: polyakov_loop_correlator.cpp $(OBJFILES)
	$(GCC) $@.cpp $(OBJFILES) $(C_LIME) -DMATRIX=$(MATRIX) -fopenmp -o polyakov_loop_correlator_$(MATRIX)_$(ARCH)
