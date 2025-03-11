INCLUDE_PATH=./lib/cpu/include
SRC_PATH=./lib/cpu/src
OBJ_PATH=./lib/cpu/obj
GCC=g++
FLAGS=-Ofast -fopenmp -march=native -mtune=native -std=c++17
C_PATH=./lib/cpu/src/c-lime
ARCH=test
C_LIME=$(C_PATH)/lime_fseeko.c  $(C_PATH)/lime_header.c  $(C_PATH)/lime_reader.c  $(C_PATH)/lime_utils.c  $(C_PATH)/lime_writer.c
OBJFILES=$(OBJ_PATH)/link_$(ARCH).o $(OBJ_PATH)/data_$(ARCH).o $(OBJ_PATH)/basic_observables_$(ARCH).o $(OBJ_PATH)/polyakov_loops_$(ARCH).o $(OBJ_PATH)/flux_tube_$(ARCH).o \
	$(OBJ_PATH)/smearing_$(ARCH).o $(OBJ_PATH)/monopoles_$(ARCH).o $(OBJ_PATH)/loop_$(ARCH).o $(OBJ_PATH)/eigen_$(ARCH).o $(OBJ_PATH)/abelian_projection_su3_$(ARCH).o \
	$(OBJ_PATH)/mag_$(ARCH).o $(OBJ_PATH)/Landau_U1_$(ARCH).o $(OBJ_PATH)/decomposition_$(ARCH).o $(OBJ_PATH)/gluon_propagator_$(ARCH).o
HEADERFILES=$(INCLUDE_PATH)/matrix.h $(INCLUDE_PATH)/indexing.h
EIGH_PATH="/home/ilya/soft/source/eigen-master"

all: $(OBJFILES)

$(OBJFILES): $(OBJ_PATH)/%_$(ARCH).o: $(SRC_PATH)/%.cpp $(HEADERFILES)
	$(GCC) $(FLAGS) -I $(EIGH_PATH) -c $< -o $@

clear:
	rm -rf $(OBJ_PATH)/*.o