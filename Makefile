all:
	gcc -Wall -Wextra -L/usr/lib fem.c -o fem -lgsl -lgslcblas -lm
run:
	./fem
allmesh:
	gcc -Wall -Wextra mesh.c -o mesh 
runmesh:
	./mesh

