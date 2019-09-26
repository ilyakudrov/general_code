INCLUDE_PATH=../../../lib/cpu/include
SRC_PATH=../../../lib/cpu/src
OBJ_PATH=../../../lib/cpu/obj
GCC=g++
test_general: $(OBJ_PATH)/test_general.o $(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/monopoles.o
	$(GCC) -I$(INCLUDE_PATH) -O3 $(OBJ_PATH)/test_general.o $(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/monopoles.o -o test_general
test_smearing: $(OBJ_PATH)/test_smearing.o $(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/monopoles.o
	$(GCC) -I$(INCLUDE_PATH) -O3 $(OBJ_PATH)/test_smearing.o $(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/monopoles.o -o test_smearing
test_eigenvectors: $(OBJ_PATH)/test_eigenvectors.o $(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/monopoles.o
	$(GCC) -I$(INCLUDE_PATH) -O3 $(OBJ_PATH)/test_eigenvectors.o $(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/observables.o $(OBJ_PATH)/result.o $(OBJ_PATH)/matrix.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/monopoles.o -o test_eigenvectors
$(OBJ_PATH)/test_general.o: test_general.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c test_general.cpp -o $(OBJ_PATH)/test_general.o
$(OBJ_PATH)/test_smearing.o: test_smearing.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c test_smearing.cpp -o $(OBJ_PATH)/test_smearing.o
$(OBJ_PATH)/test_eigenvectors.o: test_eigenvectors.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c test_eigenvectors.cpp -o $(OBJ_PATH)/test_eigenvectors.o
$(OBJ_PATH)/link.o: $(SRC_PATH)/link.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/link.cpp -o $(OBJ_PATH)/link.o
$(OBJ_PATH)/data.o: $(SRC_PATH)/data.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/data.cpp -o $(OBJ_PATH)/data.o
$(OBJ_PATH)/observables.o: $(SRC_PATH)/observables.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/observables.cpp -o $(OBJ_PATH)/observables.o
$(OBJ_PATH)/result.o: $(SRC_PATH)/result.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/result.cpp -o $(OBJ_PATH)/result.o
$(OBJ_PATH)/matrix.o: $(SRC_PATH)/matrix.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/matrix.cpp -o $(OBJ_PATH)/matrix.o
$(OBJ_PATH)/eigen.o: $(SRC_PATH)/eigen.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/eigen.cpp -o $(OBJ_PATH)/eigen.o
$(OBJ_PATH)/monopoles.o: $(SRC_PATH)/monopoles.cpp
	$(GCC) -I$(INCLUDE_PATH) -O3 -c $(SRC_PATH)/monopoles.cpp -o $(OBJ_PATH)/monopoles.o
clean:
	rm -rf $(OBJ_PATH)/*.o
