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
#define d 3


typedef enum boundary_cond {
    inner_node,DIRICHLET
}boundary_cond;


/*static int connectivity_up[96][3];
static double nodes_up[65][2];*/
	

double find_det(const double array_b[][2]);
void gauss_elimination(double A[][3], double vector_x[2]);
void construct_b(const double x[3],const double y[3], double array_b_T[][2]);
void construct_b(const double x[3],const double y[3], double array_b[][2]);
void construct_global_matrix(double A_global[][N], double b[N]);
void local_stiffness(const double array_b_T[][2],const double array_b[][2],const double det, double array_local[][3], double vector_local[3]);
void find_cord_of_nodes(const int elem, double x[3], double y[3],int global_node[3]);
void phi_value(const double x[3], double phi[3]);
void construct_global_vector(double b[N], boundary_cond cond[N]);
int gauss_elimination_gsl(const double A[][N], const double b[N], double u[N]);
gsl_matrix *construct_gsl_matrix(const double A[][N]);
gsl_vector *construct_gsl_vector(const double b[N]);
gsl_matrix *construct_gsl_matrix_2(const double A[][2]);
gsl_vector *construct_gsl_vector_2(const double b[2]);
int gauss_elimination_gsl_2(const double A[][2], const double b[2], double u[2]);


double nodes_up[21][2]={{-1.0,-1.0},{-0.5,-1.0},{0,-1.0},{0,-0.5},{0,0},{0.5,0},{1.0,0},{1.0,0.5},{1.0,1.0},{0.5,1.0},{0,1.0},{-0.5,1.0},{-1.0,1.0},{-1.0,0.5},{-1.0,0},{-1.0,-0.5},{-0.5,-0.5},{-0.5,0},{-0.5,0.5},{0,0.5},{0.5,0.5}};
int connectivity_up[24][3] = {{2,16,1},{16,2,17},{3,17,2},{17,3,4},{17,15,16},{15,17,18},{4,18,17},{18,4,5},{19,15,18},{15,19,14},{20,18,5},{18,20,19},{12,14,19},{14,12,13},{11,19,20},{19,11,12},{6,20,5},{20,6,21},{7,21,6},{21,7,8},{21,11,20},{11,21,10},{8,10,21},{10,8,9}};

int main(){
	enum boundary_cond cond_d[N];
	int i,j;
    static double A_global[N][N]={0};
    static double b[N]={0};
	double u[N];
	int res;
	FILE *fp;
	

	/*fp=fopen("con.dat", "r");
	if (fp == NULL){
	    printf("cannot open file\n");
	 }
	 
	 for(i=0; i<96; i++){
	    for(j=0; j<3; j++){
	          fscanf(fp,"%d" , &connectivity_up[i][j]);
	     }
	  }
	 fclose(fp);
	 fp=fopen("nodes.dat", "r");
	 if (fp == NULL){
	    printf("cannot open file\n");
	 }
	 for (i=0; i<65; i++){
	    for(j=0; j<2; j++) {
	        fscanf(fp, "%lf", &nodes_up[i][j]);
	    }
	 } 
	 fclose(fp);*/

/* mark nodes on the boundary,to impose Dirichlet boundary condition */
	for (i=0; i<N; i++){
		if ( (nodes_up[i][1]==-1) || (nodes_up[i][0]==-1) || (nodes_up[i][1]==1) || (nodes_up[i][0]==1) || (nodes_up[i][1]==0 && nodes_up[i][0]>=0) || (nodes_up[i][0]==0 && nodes_up[i][1]<0) ){
		    cond_d[i]=DIRICHLET;
		}
		else {
			cond_d[i]=inner_node;
		}
	} 
	construct_global_matrix(A_global,b);
    /*impose boundary condition*/
   for (i=0; i<N; i++){
        if (cond_d[i]==DIRICHLET) {
            b[i]=0;
            for (j=0; j<N; j++){
                A_global[i][j]=0;
            }
            A_global[i][i]=1;
        }
    }
   
    res=gauss_elimination_gsl(A_global,b,u);
    if (res==0){
        printf("System solved succesfully\n");
    }
    else{
        printf("there was an error in solving the system\n");
        return -1;
    }
    fp=fopen("output21.dat", "w+");
    for (j=0; j<N; j++){
        fprintf(fp,"%lf %lf %lf\n", nodes_up[j][0], nodes_up[j][1], u[j]);
    }
    fclose(fp);
 
	return 0;
}



double find_det(const double array_b[][2]){
	return (fabs(array_b[0][0]*array_b[1][1]-array_b[1][0]*array_b[0][1]));
}



void construct_b_T(const double x[3],const double y[3], double array_b_T[][2]){
	array_b_T[0][0]=x[0]-x[2];
	array_b_T[0][1]=y[0]-y[2];
	array_b_T[1][0]=x[1]-x[2];
	array_b_T[1][1]=y[1]-y[2];
}

void construct_b(const double x[3],const double y[3], double array_b[][2]){
	array_b[0][0]=x[0]-x[2];
	array_b[0][1]=x[1]-x[2];
	array_b[1][0]=y[0]-y[2];
	array_b[1][1]=y[1]-y[2];
}


void find_cord_of_nodes(const int elem, double x[3], double y[3], int global_node[3]){
	int j,sel;
		for (j=0; j<3; j++){
			sel=connectivity_up[elem][j];
			global_node[j]=connectivity_up[elem][j];
			x[j]=nodes_up[sel-1][0];
			y[j]=nodes_up[sel-1][1];		
		}
}


/**
 * @brief Perform Gauss Elimination Method on Ax=b
 *
 * Parameter array should be A augmented with the values of vector b. The result
 * lies in vector x
 *
 * @Note This function works on arrays of order 2.
 *
 * @param[in] array Augmented array (A + b)
 * @param[out] x The resulting vector
 *
 * @return none
 */

void local_stiffness(const double array_b_T[][2],const double array_b[][2],const double det, double array_local[][3], double vector_local[3]){
	int i,j,res;
	double area=0.5;
	double x1[2], x2[2], x3[2], phi[3], b_gsl[2];
    double inner_prod;
	double d_lambda[2][3]={ {1,0,-1} , {0,1,-1} };
	double lambda[2][3]={ {1,0,0} , {0,1,0} };

	for (i=0; i<3; i++){
	          b_gsl[0]=d_lambda[0][i];
	          b_gsl[1]=d_lambda[1][i];
	          res=gauss_elimination_gsl_2(array_b_T, b_gsl, x1);
	          if (res!=0){
              printf("there was an error in solving the system\n");
              return ;
              }
	          
	          for (j=0; j<3; j++){ 
	                b_gsl[0]=d_lambda[0][j];
	                b_gsl[1]=d_lambda[1][j];
	                res=gauss_elimination_gsl_2(array_b_T, b_gsl,x2);
	                if (res!=0){
                    printf("there was an error in solving the system\n");
                    return ;
                    }
	                inner_prod=x1[0]*x2[0]+x1[1]*x2[1];
	                array_local[i][j]=area*det*inner_prod;
	          }
	          b_gsl[0]=lambda[0][i];
	          b_gsl[1]=lambda[1][i];
	          res= gauss_elimination_gsl_2(array_b, b_gsl, x3);
	          if (res!=0){
              printf("there was an error in solving the system\n");
              return ;
              }
	          phi_value(x3,phi);
	          vector_local[i]=area*det*phi[i];       
	 }
}

void phi_value(double phi[3]){
/* The right-hand side is computed using Midpoint Quadrature on Triangles*/
    double coef=1/6;
    double f=1;
    
    phi[0]= f*coef*(1/2+1/2+0); /*Mindpoints of sides of a unit triangle for x*/
    phi[1]= f*coef*(0+1/2+1/2); /*Mindpoints of sides of a unit triangle for y*/
    phi[2]= f*coef*(1/2+0+1/2); /*Mindpoints of sides of a unit triangle for 1-x-y*/
}

void construct_global_matrix(double A_global[][N], double b[N]){
    int i,j,k,global_node[3];
    double x[3],y[3],det;
    double array_b[2][2],array_b_T[2][2], array_local[3][3], vector_local[3];
    
 
   for (k=0; k<NT; k++){
        find_cord_of_nodes(k,x,y,global_node);
        construct_b(x,y,array_b);
        construct_b_T(x,y,array_b_T);
        det=find_det(array_b);
        local_stiffness(array_b_T,array_b,det,array_local,vector_local);    
        for (i=0; i<3; i++){
            b[global_node[i]-1]=b[global_node[i]-1]+vector_local[i];
            for (j=0; j<3; j++){
                     A_global[global_node[i]-1][global_node[j]-1]= A_global[global_node[i]-1][global_node[j]-1] + array_local[i][j];
             }
        }
   }
}

gsl_matrix *construct_gsl_matrix(const double A[][N]){
    int i,j;

    gsl_matrix *gsl_m = gsl_matrix_alloc(N,N);
    if (gsl_m == NULL) {
        return NULL;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            gsl_matrix_set(gsl_m, i, j, A[i][j]);
        }
    }

    return gsl_m;
}

gsl_vector *construct_gsl_vector(const double b[N]){
    int i;

    gsl_vector *gsl_v = gsl_vector_alloc(N);
    if (gsl_v == NULL) {
        return NULL;
    }

    for (i = 0; i < N; i++) {
        gsl_vector_set(gsl_v, i, b[i]);
    }

    return gsl_v;
}
gsl_matrix *construct_gsl_matrix_2(const double A[][2]){
    int i,j;

    gsl_matrix *gsl_m = gsl_matrix_alloc(2,2);
    if (gsl_m == NULL) {
        return NULL;
    }

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            gsl_matrix_set(gsl_m, i, j, A[i][j]);
        }
    }

    return gsl_m;
}
gsl_vector *construct_gsl_vector_2(const double b[2]){
    int i;

    gsl_vector *gsl_v = gsl_vector_alloc(2);
    if (gsl_v == NULL) {
        return NULL;
    }

    for (i = 0; i < 2; i++) {
        gsl_vector_set(gsl_v, i, b[i]);
    }

    return gsl_v;
}


int gauss_elimination_gsl(const double A[][N], const double b[N], double u[N]){
    int s,i;

     gsl_matrix *A_m = construct_gsl_matrix(A);
     gsl_vector *b_v = construct_gsl_vector(b);
     gsl_vector *u_v = gsl_vector_alloc(N);
     gsl_permutation *p = gsl_permutation_alloc(N);

    if (A_m == NULL || b_v == NULL || u_v == NULL || p == NULL) {
        printf("Cannot construct matrix or/and vectors\n");
        return -1;
    }
    
    /* Solve the system*/
    gsl_linalg_LU_decomp(A_m, p, &s);
    gsl_linalg_LU_solve(A_m, p, b_v, u_v);
    printf("\nu = \n");
    gsl_vector_fprintf (stdout, u_v, "%lf");
 
    for (i = 0; i < N; i++) {
        u[i]=gsl_vector_get(u_v,i);
    }
    return 0;
}
int gauss_elimination_gsl_2(const double A[][2], const double b[2], double u[2]){
    int s,i;

     gsl_matrix *A_m = construct_gsl_matrix_2(A);
     gsl_vector *b_v = construct_gsl_vector_2(b);
     gsl_vector *u_v = gsl_vector_alloc(2);
     gsl_permutation *p = gsl_permutation_alloc(2);

    if (A_m == NULL || b_v == NULL || u_v == NULL || p == NULL) {
        printf("Cannot construct matrix or/and vectors\n");
        return -1;
    }
    
    /* Solve the system*/
    gsl_linalg_LU_decomp(A_m, p, &s);
    gsl_linalg_LU_solve(A_m, p, b_v, u_v);
   /* printf("\nu = \n");
    gsl_vector_fprintf (stdout, u_v, "%lf");*/
 
    for (i = 0; i < 2; i++) {
        u[i]=gsl_vector_get(u_v,i);
    }
    return 0;
}



