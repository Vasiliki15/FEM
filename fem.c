#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#define NT 24
#define N 21
// #define READ_FROM_FILE_INPUT

typedef enum boundary_cond {
	INNER_NODE, DIRICHLET
} boundary_cond;

double find_det(const double array_b[][2]);
void construct_b_T(const double x[3], const double y[3], double array_b_T[][2]);

void local_stiffness(const double array_b_T[][2],
		     const double det,
		     double array_local[][3],
		     double vector_local[3]);
void find_cord_of_nodes(const int elem,
			double x[3],
			double y[3],
			int global_node[3]);
int gauss_elimination_gsl(const double *A, const double *b, double *u, int dim);
void construct_global_matrix(double *A_global, double *b, int dim);
gsl_vector *construct_gsl_vector(const double *b, int dim);
gsl_matrix *construct_gsl_matrix(const double *A, int dim);
void mark_boundary_nodes(boundary_cond *cond, int dim);
void impose_boundary_conditions(boundary_cond *cond,
				double *A,
				double *b,
				int dim);
void read_from_file(const char *file_elements,
		    const char *file_nodes,
		    int connectivity[][3],
		    double nodes[][2]);


#ifdef READ_FROM_FILE_INPUT
static int connectivity[NT][3];
static double nodes[N][2];
#else
double nodes[21][2] = {
	{-1.0, -1.0}, {-0.5, -1.0}, {0, -1.0}, {0, -0.5}, {0, 0}, {0.5, 0},
	{1.0, 0},
	{1.0, 0.5}, {1.0, 1.0}, {0.5, 1.0}, {0, 1.0}, {-0.5, 1.0}, {-1.0, 1.0},
	{-1.0, 0.5}, {-1.0, 0}, {-1.0, -0.5}, {-0.5, -0.5}, {-0.5, 0},
	{-0.5, 0.5}, {0, 0.5}, {0.5, 0.5}
};
int connectivity[24][3] = {
	{2, 16, 1}, {16, 2, 17}, {3, 17, 2}, {17, 3, 4}, {17, 15, 16},
	{15, 17, 18},
	{4, 18, 17}, {18, 4, 5}, {19, 15, 18}, {15, 19, 14}, {20, 18, 5},
	{18, 20, 19}, {12, 14, 19}, {14, 12, 13}, {11, 19, 20}, {19, 11, 12},
	{6, 20, 5}, {20, 6, 21}, {7, 21, 6}, {21, 7, 8}, {21, 11, 20},
	{11, 21, 10}, {8, 10, 21}, {10, 8, 9}
};
 #endif

int main()
{
	boundary_cond *cond;
	int j;
	double *A_global;
	double *b;
	double u[N];
	int res;
	FILE *fp;


	A_global = calloc(N * N, sizeof(double));
	if (A_global == NULL) {
		perror("Allocate A_global\n");
		return 1;
	}
	cond = calloc(N, sizeof(double));
	if (cond == NULL) {
		perror("Allocate cond\n");
		return 1;
	}
	b = calloc(N, sizeof(double));
	if (cond == NULL) {
		perror("Allocate b\n");
		return 1;
	}

	#ifdef READ_FROM_FILE_INPUT
	read_from_file("connectivety_elems.dat",
		       "nodes_cord.dat",
		       connectivety,
		       nodes);
	#endif


	mark_boundary_nodes(cond, N);
	construct_global_matrix(A_global, b, N);
	impose_boundary_conditions(cond, A_global, b, N);
	res = gauss_elimination_gsl(A_global, b, u, N);
	if (res == 0) {
		printf("System solved succesfully\n");
	} else {
		printf("there was an error in solving the system\n");
		return -1;
	}
	fp = fopen("output.dat", "w+");
	for (j = 0; j < N; j++) {
		fprintf(fp, "%lf %lf %lf\n", nodes[j][0], nodes[j][1], u[j]);
	}
	fclose(fp);

	return 0;
}

/*The function reads data from files and stores them to corresponding connectivety and nodes arrays*/
void read_from_file(const char *file_elements,
		    const char *file_nodes,
		    int connectivity[][3],
		    double nodes[][2])
{
	int i, j;
	FILE *fp;

	fp = fopen(file_elements, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
	}
	for (i = 0; i < NT; i++) {
		for (j = 0; j < 3; j++) {
			fscanf(fp, "%d", &connectivity[i][j]);
		}
	}
	fclose(fp);
	fp = fopen(file_nodes, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < 2; j++) {
			fscanf(fp, "%lf", &nodes[i][j]);
		}
	}
	fclose(fp);
}

/*impose boundary condition*/
void impose_boundary_conditions(boundary_cond *cond,
				double *A,
				double *b,
				int dim)
{
	int i, j;

	for (i = 0; i < dim; i++) {
		if (cond[i] == DIRICHLET) {
			b[i] = 0;
			for (j = 0; j < dim; j++) {
				A[i * dim + j] = 0;
			}
			A[i * dim + i] = 1;
		}
	}
}

/* mark nodes on the boundary,to impose Dirichlet boundary condition */
void mark_boundary_nodes(boundary_cond *cond, int dim)
{
	int i;

	for (i = 0; i < dim; i++) {
		if ( (nodes[i][1] == -1) || (nodes[i][0] == -1) ||
		     (nodes[i][1] == 1) || (nodes[i][0] == 1) ||
		     (nodes[i][1] == 0 && nodes[i][0] >= 0) ||
		     (nodes[i][0] == 0 && nodes[i][1] < 0) ) {
			cond[i] = DIRICHLET;
		} else {
			cond[i] = INNER_NODE;
		}
	}
}

/* Find the det of B matrix*/
double find_det(const double array_b_T[][2])
{
	return (fabs(array_b_T[0][0] * array_b_T[1][1] - array_b_T[1][0] *
		     array_b_T[0][1]));
}

/*Construct transpose of tranformation matrix B */
void construct_b_T(const double x[3], const double y[3], double array_b_T[][2])
{
	array_b_T[0][0] = x[0] - x[2];
	array_b_T[0][1] = y[0] - y[2];
	array_b_T[1][0] = x[1] - x[2];
	array_b_T[1][1] = y[1] - y[2];
}

/**
 * @brief Find the x,y coordinates of every node
 *
 *
 * @param[in] element of which the nodes are concerned
 * @param[out] a vector of x coordinates of the three nodes of every element
 * @param[out] a vector of y coordinates of the three nodes of every element
 * @param[out] a vector of nodes of th element
 *
 * @return none
 */
void find_cord_of_nodes(const int elem,
			double x[3],
			double y[3],
			int global_node[3])
{
	int j, sel;
	for (j = 0; j < 3; j++) {
		sel = connectivity[elem][j];
		global_node[j] = connectivity[elem][j];
		x[j] = nodes[sel - 1][0];
		y[j] = nodes[sel - 1][1];
	}
}

/**
 * @brief Assemply of local stiffness matrix
 *
 * Parameter array should be the transpose of transformation mastrix B
 *
 *
 * @Note This function constructs a 3x3 local array
 *
 * @param[in] transpose B
 * @param[in] det of B
 * @param[out] x The resulting local array and corresponding local right hand side vector
 *
 * @return none
 */

void local_stiffness(const double array_b_T[][2],
		     const double det,
		     double array_local[][3],
		     double vector_local[3])
{
	int i, j, res;
	double area = 0.5;
	double x1[2], x2[2], b_gsl[2];
	double inner_prod;
	double d_lambda[2][3] = {
		{1, 0, -1}, {0, 1, -1}
	};


	for (i = 0; i < 3; i++) {
		b_gsl[0] = d_lambda[0][i];
		b_gsl[1] = d_lambda[1][i];
		res = gauss_elimination_gsl((const double *)array_b_T,
					    b_gsl,
					    x1,
					    2);
		if (res != 0) {
			printf("there was an error in solving the system\n");
			return;
		}

		for (j = 0; j < 3; j++) {
			b_gsl[0] = d_lambda[0][j];
			b_gsl[1] = d_lambda[1][j];
			res = gauss_elimination_gsl((const double *)array_b_T,
						    b_gsl,
						    x2,
						    2);
			if (res != 0) {
				printf(
					"there was an error in solving the system\n");
				return;
			}
			inner_prod = x1[0] * x2[0] + x1[1] * x2[1];
			array_local[i][j] = area * det * inner_prod;
		}
		vector_local[i] = 0.5;
	}
}

/**
 * @brief Assemply of global stiffness matrix
 *
 * @param[in] dimension of the matrix
 * @param[out] global stiffness matrix
 * @param[out] global stiffness matrix right hand side vector
 *
 * @return none
 */

void construct_global_matrix(double *A_global, double *b, int dim)
{
	int i, j, k, global_node[3];
	double x[3], y[3], det;
	double array_b_T[2][2], array_local[3][3],
	       vector_local[3];


	for (k = 0; k < NT; k++) {
		find_cord_of_nodes(k, x, y, global_node);
		construct_b_T(x, y, array_b_T);
		det = find_det(array_b_T);
		local_stiffness(array_b_T, det, array_local, vector_local);
		for (i = 0; i < 3; i++) {
			b[global_node[i] - 1] += vector_local[i];
			for (j = 0; j < 3; j++) {
				A_global[(global_node[i] - 1) * dim +
					 (global_node[j] -
					  1)] +=  array_local[i][j];
			}
		}
	}
}

/*Construct matrix in gsl format*/
gsl_matrix *construct_gsl_matrix(const double *A, int dim)
{
	int i, j;

	gsl_matrix *gsl_m = gsl_matrix_alloc(dim, dim);
	if (gsl_m == NULL) {
		return NULL;
	}

	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			gsl_matrix_set(gsl_m, i, j, A[i * dim + j]);
		}
	}

	return gsl_m;
}

/*Construct vector in gsl format*/
gsl_vector *construct_gsl_vector(const double *b, int dim)
{
	int i;

	gsl_vector *gsl_v = gsl_vector_alloc(dim);
	if (gsl_v == NULL) {
		return NULL;
	}

	for (i = 0; i < dim; i++) {
		gsl_vector_set(gsl_v, i, b[i]);
	}

	return gsl_v;
}

/*Solves linear system*/
int gauss_elimination_gsl(const double *A, const double *b, double *u, int dim)
{
	int s, i;

	gsl_matrix *A_m = construct_gsl_matrix(A, dim);
	gsl_vector *b_v = construct_gsl_vector(b, dim);
	gsl_vector *u_v = gsl_vector_alloc(dim);
	gsl_permutation *p = gsl_permutation_alloc(dim);

	if (A_m == NULL || b_v == NULL || u_v == NULL || p == NULL) {
		printf("Cannot construct matrix or/and vectors\n");
		return -1;
	}

	/* Solve the system*/
	gsl_linalg_LU_decomp(A_m, p, &s);
	gsl_linalg_LU_solve(A_m, p, b_v, u_v);
	printf("\nu = \n");
	gsl_vector_fprintf(stdout, u_v, "%lf");

	for (i = 0; i < dim; i++) {
		u[i] = gsl_vector_get(u_v, i);
	}
	return 0;
}
