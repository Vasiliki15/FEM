#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#define NT 24
#define N 21
#define d 3


enum boundary_cond {inner_node,Dirichlet};


double find_det(const double array_b[][2]);
void gauss_elimination(double A[][3], double vector_x[2]);
void construct_b(double x[3], double y[3], double array_b[][2]);
void construct_global_matrix(double A_global[][N]);
void local_stiffness(const double array_b[][2], double det, double array_local[][3]); 
void find_cord_of_nodes(int elem, double x[3], double y[3],int global_node[3]);

int main(){
	enum boundary_cond cond_d[21];
	int i,j;
	double A_global[N][N]={0};
	double b[N]={0};
	
	
/* mark nodes on the boundary,to impose Dirichlet boundary condition */
	for (i=0; i<N; i++){
		if (i<16) {
			cond_d[i]=Dirichlet;
		}
		else {
			cond_d[i]=inner_node;
		}
	}

    for (i=0; i<N; i++){
        if (cond_d[i]==Dirichlet) {
            
            for (j=0; j<N; j++){
                A_global[i][j]=0;
            }
            A_global[i][i]=1;
        }
       }
                
   
      
	construct_global_matrix(A_global);
	/*for (i=0; i<N; i++){
		for(j=0; j<N; j++){
		    
			   if ( i==j){
			      if (A_global[i][j]==0){
			        printf("fuck fuck fuckkkkkkkk\n");
			       }
			      }
			
		}
	}*/
	return 0;
}

double nodes[21][2]={{-1.0,-1.0},{-0.5,-1.0},{0,-1.0},{0,-0.5},{0,0},{0.5,0},{1.0,0},{1.0,0.5},{1.0,1.0},{0.5,1.0},{0,1.0},{-0.5,1.0},{-1.0,1.0},{-1.0,0.5},{-1.0,0},{-1.0,-0.5},{-0.5,-0.5},{-0.5,0},{-0.5,0.5},{0,0.5},{0.5,0.5}};
int connectivity[24][3] = {{2,16,1},{16,2,17},{3,17,2},{17,3,4},{17,15,16},{15,17,18},{4,18,17},{18,4,5},{19,15,18},{15,19,14},{20,18,5},{18,20,19},{12,14,19},{14,12,13},{11,19,20},{19,11,12},{6,20,5},{20,6,21},{7,21,6},{21,7,8},{21,11,20},{11,21,10},{8,10,21},{10,8,9}};

double find_det(const double array_b[][2]){
	return (fabs(array_b[0][0]*array_b[1][1]-array_b[1][0]*array_b[0][1]));
}



void construct_b(double x[3], double y[3], double array_b[][2]){

	array_b[0][0]=x[0]-x[2];
	array_b[0][1]=y[0]-y[2];
	array_b[1][0]=x[1]-x[2];
	array_b[1][1]=y[1]-y[2];
}


void find_cord_of_nodes(int elem, double x[3], double y[3], int global_node[3]){
	int i,j,sel;
	

		for (j=0; j<3; j++){
			sel=connectivity[elem][j];
			global_node[j]=connectivity[elem][j]-1;
			x[j]=nodes[sel-1][0];
			y[j]=nodes[sel-1][1];
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
void gauss_elimination(double array[][3], double x[2])
{
	int i,j,k,n;
	int sum=0;
	double c;
	double p[2][2]={ {0 ,1}, {1, 0} };
	double array_temp[2][2];
	
	
    /* Find the order */
	n = sizeof(array[0]) / sizeof(array[0][0]) - 1;
	
	for (i=0; i<2; i++){
	    for (j=0; j<3; j++){
	        printf("%lf\n", array[i][j]);
	      }
	      }
   
    if (array[0][0] == 0){
            for (i=0; i<2; i++){
                  for (j=0; j<3; j++){
                        array_temp[i][j]=p[i][0]*array[0][j]+p[i][1]*array[1][j];
                                
                   }
             }
              for (i=0; i<2; i++){
                   for(j=0; j<3; j++){
                        array[i][j]=array_temp[i][j];
                   }
              }
            x[1]=array[1][2]/array[1][1];
            x[0]=(array[0][2]-array[0][1]*x[1])/array[0][0];
           printf("x0=%lf, x1=%lf\n", x[0], x[1]);
               return;                
   }
    
    /* loop for the generation of upper triangular matrix*/
	for(j=0; j<n; j++) {
        for(i=0; i<n; i++){
		    if(i>j){
                c=array[i][j]/array[j][j];
                
                for(k=0; k<n+1; k++){
                    array[i][k]=array[i][k]-c*array[j][k];
                }
            }
        }
    }

    x[n-1]=array[n-1][n]/array[n-1][n-1];

    /* this loop is for backward substitution*/
    for(i=n-2; i>=0; i--){
        sum=0;
        for(j=i; j<n; j++){
            sum=sum+array[i][j]*x[j];
        }
        x[i]=(array[i][n]-sum)/array[i][i];
    }
    
    for(i=0; i<n; i++){
        printf("\nx%d=%lf\t\n",i,x[i]);
    }
    
}

void local_stiffness(const double array_b[][2], double det, double array_local[][3]){
	int i,j;
	double area=0.5;
	double x1[2], x2[2];
    double inner_prod;
	double lambda[2][3]={ {1,0,-1} , {0,1,-1} };
	double array_b_aug[2][3]= { {array_b[0][0], array_b[0][1], 0}, {array_b[1][0], array_b[1][1], 0} };

	
	
	for (i=0; i<3; i++){
	          array_b_aug[0][2]=lambda[0][i];
	          array_b_aug[1][2]=lambda[1][i];
	          gauss_elimination(array_b_aug,x1);
	          for (j=0; j<3; j++){ 
	                array_b_aug[0][2]=lambda[0][j];
	                array_b_aug[1][2]=lambda[1][j];
	                gauss_elimination(array_b_aug,x2);
	                inner_prod=x1[0]*x2[0]+x1[1]*x2[1];
	                array_local[i][j]=area*det*inner_prod;
	            
	          }
	 }

}

void construct_global_matrix(double A_global[][N]){
    int i,j,k,global_node[3];
    double x[3],y[3],det;
    double array_b[2][2], array_local[3][3];
    
   /*for (k=0; k<NT; k++){*/
        find_cord_of_nodes(20,x,y,global_node);
        construct_b(x,y,array_b);
        det=find_det(array_b);
        local_stiffness(array_b,det,array_local);      
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                     A_global[global_node[i]][global_node[j]]= A_global[global_node[i]][global_node[j]] + array_local[i][j];
             }
        }
   
}





