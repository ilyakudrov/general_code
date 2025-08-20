INCLUDE_PATH=../../../lib/cpu/include
OBJ_PATH=../../../lib/cpu/obj
GCC=g++
MATRIX=su2
FLAGS=-Ofast -fopenmp -march=native -mtune=native -std=c++17 -DMATRIX=$(MATRIX)
C_PATH=../../../lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJ_FILES=$(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/smearing_$(ARCH).o \
	$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o $(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/gluon_propagator_$(ARCH).o \
	$(OBJ_PATH)/Landau_U1_$(ARCH).o $(OBJ_PATH)/Landau_su2_$(ARCH).o
HEADERFILES=$(INCLUDE_PATH)/matrix.h $(INCLUDE_PATH)/indexing.h
EIGH_PATH="/home/ilya/soft/source/eigen-master"

gluon_propagator: gluon_propagator.cpp $(OBJ_FILES) $(HEADERFILES)
	$(GCC) gluon_propagator.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -I $(EIGH_PATH) -o gluon_propagator_$(MATRIX)_$(ARCH)

gluon_propagator_abelian: gluon_propagator_abelian.cpp $(OBJ_FILES) $(HEADERFILES)
	$(GCC) gluon_propagator_abelian.cpp $(OBJ_FILES) $(C_LIME) $(FLAGS) -I $(EIGH_PATH) -o gluon_propagator_abelian_$(ARCH)