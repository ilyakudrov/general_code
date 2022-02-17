INCLUDE_PATH=../../lib/cpu/include
SRC_PATH=../../lib/cpu/src
OBJ_PATH=../../lib/cpu/obj
GCC=g++
FLAGS=-O3 -march=native -ffast-math
C_PATH=../../lib/cpu/src/c-lime
C_FILES=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/basic_observables.o $(OBJ_PATH)/flux_tube.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/smearing.o $(OBJ_PATH)/monopoles.o $(OBJ_PATH)/loop.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/abelian_projection_su3.o
indexing_test: indexing_test.cpp src/indexing.cpp $(OBJFILES)
	$(GCC) $(FLAGS) indexing_test.cpp src/indexing.cpp $(OBJFILES) $(C_FILES) -o indexing_test
$(OBJ_PATH)/link.o: $(SRC_PATH)/link.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/link.cpp -o $(OBJ_PATH)/link.o
$(OBJ_PATH)/data.o: $(SRC_PATH)/data.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/data.cpp -o $(OBJ_PATH)/data.o
$(OBJ_PATH)/basic_observables.o: $(SRC_PATH)/basic_observables.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/basic_observables.cpp -o $(OBJ_PATH)/basic_observables.o
$(OBJ_PATH)/flux_tube.o: $(SRC_PATH)/flux_tube.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/flux_tube.cpp -o $(OBJ_PATH)/flux_tube.o
$(OBJ_PATH)/matrix.o: $(SRC_PATH)/matrix.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/matrix.cpp -o $(OBJ_PATH)/matrix.o
$(OBJ_PATH)/smearing.o: $(SRC_PATH)/smearing.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/smearing.cpp -o $(OBJ_PATH)/smearing.o
$(OBJ_PATH)/eigen.o: $(SRC_PATH)/eigen.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/eigen.cpp -o $(OBJ_PATH)/eigen.o
$(OBJ_PATH)/monopoles.o: $(SRC_PATH)/monopoles.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/monopoles.cpp -o $(OBJ_PATH)/monopoles.o
$(OBJ_PATH)/loop.o: $(SRC_PATH)/loop.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/loop.cpp -o $(OBJ_PATH)/loop.o
$(OBJ_PATH)/abelian_projection_su3.o: $(SRC_PATH)/abelian_projection_su3.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/abelian_projection_su3.cpp -o $(OBJ_PATH)/abelian_projection_su3.o
clean:
	rm -rf $(OBJ_PATH)/*.o
