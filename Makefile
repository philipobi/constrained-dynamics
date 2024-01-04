CC=g++ -std=c++17
INCL=-I"/Users/philip/Documents/coding/terminal graphics/include" -I/opt/homebrew/Cellar/eigen/3.4.0_1/include

main: src/main.cpp src/utils.cpp 
	$(CC) $(INCL) src/main.cpp src/utils.cpp -o bin/main

donut: src/donut.cpp src/utils.cpp 
	$(CC) $(INCL) src/donut.cpp src/utils.cpp -o bin/donut

debug_main: src/main.cpp src/utils.cpp 
	$(CC) $(INCL) src/main.cpp src/utils.cpp -g -o bin/main

test: src/test.cpp

	$(CC) $(INCL) src/test.cpp src/utils.cpp -o bin/test


