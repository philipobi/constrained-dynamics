INCL=-Iinclude

sparse: src/sparse_linalg.c src/linalg.c
	gcc $(INCL) -c -o bin/sparse_linalg.o src/sparse_linalg.c
	gcc $(INCL) -c -o bin/linalg.o src/linalg.c
	gcc -o bin/sparse bin/sparse_linalg.o bin/linalg.o 
	

main: src/sparse_linalg.c src/linalg.c src/solver.c src/main.c
	gcc $(INCL) -c -o bin/sparse_linalg.o src/sparse_linalg.c
	gcc $(INCL) -c -o bin/linalg.o src/linalg.c
	gcc $(INCL) -c -o bin/solver.o src/solver.c
	gcc $(INCL) -c -o bin/main.o src/main.c
	gcc -o bin/main bin/sparse_linalg.o bin/linalg.o bin/solver.o bin/main.o