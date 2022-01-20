INCLUDE_PATH=../../../lib/cpu/include
SRC_PATH=../../../lib/cpu/src
OBJ_PATH=../../../lib/cpu/obj
GCC=g++
FLAGS=-O3 -march=native -mtune=native -ffast-math
C_PATH=../../../lib/cpu/src/c-lime
C_FILES=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/basic_observables.o $(OBJ_PATH)/flux_tube.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/smearing.o $(OBJ_PATH)/monopoles.o $(OBJ_PATH)/loop.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/abelian_projection_su3.o
test_general: $(OBJ_PATH)/test_general.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_general.o $(OBJFILES) $(C_FILES) -o test_general
test_smearing: $(OBJ_PATH)/test_smearing.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_smearing.o $(OBJFILES) $(C_FILES) -o test_smearing
test_eigenvectors: $(OBJ_PATH)/test_eigenvectors.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_eigenvectors.o $(OBJFILES) $(C_FILES) -o test_eigenvectors
test_flux_tube: $(OBJ_PATH)/test_flux_tube.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_flux_tube.o $(OBJFILES) $(C_FILES) -o test_flux_tube
test_monopoles: $(OBJ_PATH)/test_monopoles.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_monopoles.o $(OBJFILES) $(C_FILES) -o test_monopoles
monopole_decomposition: $(OBJ_PATH)/monopole_decomposition.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/monopole_decomposition.o $(OBJFILES) $(C_FILES) -o monopole_decomposition
$(OBJ_PATH)/test_general.o: test_general.cpp
	$(GCC) $(FLAGS) -c test_general.cpp -o $(OBJ_PATH)/test_general.o
$(OBJ_PATH)/test_smearing.o: test_smearing.cpp
	$(GCC) $(FLAGS) -c test_smearing.cpp -o $(OBJ_PATH)/test_smearing.o
$(OBJ_PATH)/test_eigenvectors.o: test_eigenvectors.cpp
	$(GCC) $(FLAGS) -c test_eigenvectors.cpp -o $(OBJ_PATH)/test_eigenvectors.o
$(OBJ_PATH)/test_flux_tube.o: test_flux_tube.cpp
	$(GCC) $(FLAGS) -c test_flux_tube.cpp -o $(OBJ_PATH)/test_flux_tube.o
$(OBJ_PATH)/test_monopoles.o: test_monopoles.cpp
	$(GCC) $(FLAGS) -c test_monopoles.cpp -o $(OBJ_PATH)/test_monopoles.o
$(OBJ_PATH)/monopole_decomposition.o: monopole_decomposition.cpp
	$(GCC) $(FLAGS) -c monopole_decomposition.cpp -o $(OBJ_PATH)/monopole_decomposition.o
$(OBJ_PATH)/link.o: $(SRC_PATH)/link.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/link.cpp -o $(OBJ_PATH)/link.o
$(OBJ_PATH)/data.o: $(SRC_PATH)/data.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/data.cpp -o $(OBJ_PATH)/data.o
$(OBJ_PATH)/basic_observables.o: $(SRC_PATH)/basic_observables.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/basic_observables.cpp -o $(OBJ_PATH)/basic_observables.o
$(OBJ_PATH)/flux_tube.o: $(SRC_PATH)/flux_tube.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/flux_tube.cpp -o $(OBJ_PATH)/flux_tube.o
$(OBJ_PATH)/result.o: $(SRC_PATH)/result.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/result.cpp -o $(OBJ_PATH)/result.o
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
