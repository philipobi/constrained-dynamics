INCL=-Iinclude

sparse: src/sparse_linalg.c src/linalg.c
	gcc $(INCL) -c -o bin/sparse_linalg.o src/sparse_linalg.c
	gcc $(INCL) -c -o bin/linalg.o src/linalg.c
	gcc -o bin/sparse bin/sparse_linalg.o bin/linalg.o 
	