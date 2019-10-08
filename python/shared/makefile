PATH_SRC=../../lib/cpu/src
PATH_INCLUDE=../../lib/cpu/include

result: $(PATH_SRC)/result.cpp
	g++ -O4 -shared -o result.so -fPIC -I$(PATH_INCLUDE) $(PATH_SRC)/result.cpp
