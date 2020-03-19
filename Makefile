all:
	gcc -Wall -Wextra -L/usr/lib fem.c -o fem -lgsl -lgslcblas -lm
run:
	./fem
