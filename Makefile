CC=g++ -std=c++11
INCL=-I"/Users/philip/Documents/coding/terminal graphics/include" -I/opt/homebrew/opt/openblas/include -I"/opt/homebrew/Cellar/eigen/3.4.0_1/include"
LIB=-L/opt/homebrew/opt/openblas/lib -lopenblas

main: src/main.cpp src/utils.cpp 
	$(CC) $(INCL) $(LIB) src/main.cpp src/utils.cpp -o bin/main

test: src/test.cpp
	$(CC) $(INCL) $(LIB) src/test.cpp -o bin/test



