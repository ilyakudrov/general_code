INCLUDE_PATH=../../../lib/cpu/include
SRC_PATH=../../../lib/cpu/src
OBJ_PATH=../../../lib/cpu/obj
GCC=g++
FLAGS=-O3 #-ffast-math
C_PATH=../../../lib/cpu/src/c-lime
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/basic_observables.o $(OBJ_PATH)/flux_tube.o $(OBJ_PATH)/matrix.o \
	$(OBJ_PATH)/smearing.o $(OBJ_PATH)/monopoles.o $(OBJ_PATH)/loop.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/abelian_projection_su3.o \
	$(OBJ_PATH)/Landau_gauge.o

TARGETS=test_general test_smearing test_eigenvectors test_flux_tube test_monopoles monopole_decomposition test_Landau_gauge

$(TARGETS): %: %.cpp $(OBJFILES)
	$(GCC) $@.cpp $(OBJFILES) $(C_LIME) -o $@

$(OBJFILES): $(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp
	$(GCC) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_PATH)/*.o
