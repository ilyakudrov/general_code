INCLUDE_PATH=../../../lib/cpu/include
SRC_PATH=../../../lib/cpu/src
OBJ_PATH=../../../lib/cpu/obj
GCC=g++
FLAGS=-O3 -D DOUBLE -g
OBJFILES=$(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/smearing.o
test_general: $(OBJ_PATH)/test_general.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_general.o $(OBJFILES) -lprofiler -o test_general
test_smearing: $(OBJ_PATH)/test_smearing.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_smearing.o $(OBJFILES) -o test_smearing
test_eigenvectors: $(OBJ_PATH)/test_eigenvectors.o $(OBJFILES)
	$(GCC) $(OBJ_PATH)/test_eigenvectors.o $(OBJFILES) -o test_eigenvectors
$(OBJ_PATH)/test_general.o: test_general.cpp
	$(GCC) $(FLAGS) -c test_general.cpp -o $(OBJ_PATH)/test_general.o
$(OBJ_PATH)/test_smearing.o: test_smearing.cpp
	$(GCC) $(FLAGS) -c test_smearing.cpp -o $(OBJ_PATH)/test_smearing.o
$(OBJ_PATH)/test_eigenvectors.o: test_eigenvectors.cpp
	$(GCC) $(FLAGS) -c test_eigenvectors.cpp -o $(OBJ_PATH)/test_eigenvectors.o
$(OBJ_PATH)/link.o: $(SRC_PATH)/link.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/link.cpp -o $(OBJ_PATH)/link.o
$(OBJ_PATH)/data.o: $(SRC_PATH)/data.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/data.cpp -o $(OBJ_PATH)/data.o
$(OBJ_PATH)/observables.o: $(SRC_PATH)/observables.cpp
	$(GCC) $(FLAGS) -c $(SRC_PATH)/observables.cpp -o $(OBJ_PATH)/observables.o
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
clean:
	rm -rf $(OBJ_PATH)/*.o
