OBJ_PATH=../lib/cpu/obj
GCC=g++
C_PATH=../lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/flux_tube_$(ARCH).o $(OBJ_PATH)/matrix_$(ARCH).o \
	$(OBJ_PATH)/smearing_$(ARCH).o $(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/eigen_$(ARCH).o $(OBJ_PATH)/abelian_projection_su3_$(ARCH).o \
	$(OBJ_PATH)/mag_$(ARCH).o

TARGETS=mag_fixation

$(TARGETS): %: %.cpp $(OBJFILES)
	$(GCC) $@.cpp $(OBJFILES) $(C_LIME) -fopenmp -o $@_$(ARCH)
