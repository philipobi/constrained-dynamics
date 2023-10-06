CC=g++

main: src/main.cpp src/utils.cpp 
	$(CC) -L/opt/homebrew/opt/openblas/lib -I"/Users/philip/Documents/coding/terminal graphics/include" -I/opt/homebrew/opt/openblas/include -lopenblas -std=c++11 src/main.cpp src/utils.cpp -o bin/main

test: src/test.cpp

	$(CC) -L/opt/homebrew/opt/openblas/lib -I"/Users/philip/Documents/coding/terminal graphics/include" -I/opt/homebrew/opt/openblas/include -lopenblas -std=c++11 src/test.cpp -o bin/test



