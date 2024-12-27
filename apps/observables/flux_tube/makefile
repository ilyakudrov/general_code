OBJ_PATH=../../../lib/cpu/obj
GCC=g++
MATRIX_PLAKET=su2
MATRIX_WILSON=su2
FLAGS=-Ofast -fopenmp -march=native -mtune=native -DMATRIX_PLAKET=$(MATRIX_PLAKET) -DMATRIX_WILSON=$(MATRIX_WILSON)
C_PATH=../../../lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJ_FILES=$(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/smearing_$(ARCH).o $(OBJ_PATH)/matrix_$(ARCH).o \
	$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o $(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/flux_tube_$(ARCH).o

flux_tube_wilson: flux_tube_wilson.cpp $(OBJ_FILES)
	$(GCC) flux_tube_wilson.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -o flux_tube_wilson_$(MATRIX_PLAKET)_$(MATRIX_WILSON)_$(ARCH)
flux_tube_schwinger: flux_tube_schwinger.cpp $(OBJ_FILES)
	$(GCC) flux_tube_schwinger.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -o flux_tube_schwinger_$(MATRIX_PLAKET)_$(MATRIX_WILSON)_$(ARCH)