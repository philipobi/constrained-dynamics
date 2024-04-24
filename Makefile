COMPFLAGS := -g -Iinclude
LINKFLAGS := -g

components := main simulation constraint linalg sparse_linalg solver
files := $(foreach component,$(components),bin/$(component).o)

build: clean $(files)
	gcc $(LINKFLAGS) -o bin/main $(files)

bin/%.o:
	gcc $(COMPFLAGS) -c -o $@ src/$(*F).c

clean:
	rm -f bin/*.o

