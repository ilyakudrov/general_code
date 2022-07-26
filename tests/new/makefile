INCLUDE_PATH=../../lib/cpu/include
SRC_PATH=../../lib/cpu/src
OBJ_PATH=../../lib/cpu/obj
GCC=g++
FLAGS=-Ofast -fopenmp -march=native -mtune=native
C_PATH=../../lib/cpu/src/c-lime
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
ARCH=test
OBJFILES=$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/flux_tube_$(ARCH).o $(OBJ_PATH)/matrix_$(ARCH).o \
	$(OBJ_PATH)/smearing_$(ARCH).o $(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/eigen_$(ARCH).o $(OBJ_PATH)/abelian_projection_su3_$(ARCH).o \
	$(OBJ_PATH)/mag_$(ARCH).o $(OBJ_PATH)/Landau_U1_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o

SRC_TARGETS=link data basic_observables flux_tube matrix smearing monopoles loop eigen abelian_projection_su3

indexing_test: indexing_test.cpp src/indexing.cpp $(OBJFILES) $(C_LIME)
	$(GCC) $(FLAGS) $^ -o indexing_test_$(ARCH)
wilson_loop_indexing_test: wilson_loop_indexing_test.cpp src/wilson_loop_indexing.cpp $(OBJFILES) $(C_LIME)
	$(GCC) $(FLAGS) $^ -o wilson_loop_indexing_test_$(ARCH)
smearing_indexing_test: smearing_indexing_test.cpp src/smearing_indexing.cpp src/wilson_loop_indexing.cpp  $(OBJFILES) $(C_LIME)
	$(GCC) $(FLAGS) $^ -o smearing_indexing_test_$(ARCH)
$(OBJFILES): $(OBJ_PATH)/%_$(ARCH).o: $(SRC_PATH)/%.cpp
	$(GCC) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_PATH)/*.o
