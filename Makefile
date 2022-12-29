CFLAGS  := -Wall -Wextra -Wpedantic -ggdb
INC_DIR := include/

all: solve inverse

inverse: inverse.o matrix.o
	$(CC) -o $@ $^

inverse.o: src/inverse.c include/matrix.h
	$(CC) -c -I $(INC_DIR) $(CFLAGS) -o $@ $<

solve: solve.o matrix.o
	$(CC) -o $@ $^

solve.o: src/solve.c include/matrix.h
	$(CC) -c -I $(INC_DIR) $(CFLAGS) -o $@ $<

matrix.o: src/matrix.c include/matrix.h
	$(CC) -c -I $(INC_DIR) $(CFLAGS) -o $@ $<

clean:
	$(RM) *.o solve
