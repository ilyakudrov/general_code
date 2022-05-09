INCLUDE_PATH=../../lib/cpu/include
SRC_PATH=../../lib/cpu/src
OBJ_PATH=../../lib/cpu/obj
GCC=g++
FLAGS=-O3 -ffast-math -fopenmp
C_PATH=../../lib/cpu/src/c-lime
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/basic_observables.o $(OBJ_PATH)/flux_tube.o \
$(OBJ_PATH)/matrix.o $(OBJ_PATH)/smearing.o $(OBJ_PATH)/monopoles.o $(OBJ_PATH)/loop.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/abelian_projection_su3.o
SRC_TARGETS=link data basic_observables flux_tube matrix smearing monopoles loop eigen abelian_projection_su3

indexing_test: indexing_test.cpp src/indexing.cpp $(OBJFILES) $(C_LIME)
	$(GCC) $(FLAGS) $^ -o indexing_test
wilson_loop_indexing_test: wilson_loop_indexing_test.cpp src/wilson_loop_indexing.cpp $(OBJFILES) $(C_LIME)
	$(GCC) $(FLAGS) $^ -o wilson_loop_indexing_test
smearing_indexing_test: smearing_indexing_test.cpp src/smearing_indexing.cpp src/wilson_loop_indexing.cpp  $(OBJFILES) $(C_LIME)
	$(GCC) $(FLAGS) $^ -o smearing_indexing_test
$(OBJFILES): $(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp
	$(GCC) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_PATH)/*.o
