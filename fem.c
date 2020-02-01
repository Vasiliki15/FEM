#include <stdio.h>
#define NT 24
#define N 21
#define d 3


enum boundary_cond {inner_node,Dirichlet};


int find_det(const int array_b[][2]);
void local_stiffness(const int array_b[][2], int inv_array_b[][2], float A_local[][3]);
void find_inv_matrix(const int array_b[][2], int inv_array_b[][2]);
/*void construct_b(int x1, int x2, int x3);*/


int main(){
	enum boundary_cond cond_d[21];
	int i,j,det;
	float A_global[N][N]={0};
	float A_local[3][3];
	float *local_array;
	float b[N]={0};
	int array_b[2][2]={{2,1},{3,2}};
	int inv_array_b[2][2];
	

/* mark nodes on the boundary,to impose Dirichlet boundary condition */
	for (i=0; i<N; i++){
		if (i<16) {
			cond_d[i]=Dirichlet;
		}
		else {
			cond_d[i]=inner_node;
		}
	}
	/*local=&A_local[0][1];*/
	find_inv_matrix(array_b, inv_array_b);
	local_stiffness(array_b,inv_array_b,A_local);
	for (i=0; i<3; i++){
		for(j=0; j<3; j++){
			printf("%f\n", A_local[i][j]);
		}
	}
	return 0;
}

float nodes[21][2]={{-1.0,-1.0},{-0.5,-1.0},{0,-1.0},{0,-0.5},{0,0},{0.5,0},{1.0,0},{1.0,0.5},{1.0,1.0},{0.5,1.0},{0,1.0},{-0.5,1.0},{-1.0,1.0},{-1.0,0.5},{-1.0,0},{-1.0,-0.5},{-0.5,-0.5},{-0.5,0},{-0.5,0.5},{0,0.5},{0.5,0.5}};
int elem[24][3] = {{2,16,1},{16,2,17},{3,17,2},{17,3,4},{17,15,16},{15,17,18},{4,18,17},{18,4,5},{19,15,18},{15,19,14},{20,18,5},{18,20,19},{12,14,19},{14,12,13},{11,19,20},{19,11,12},{6,20,5},{20,6,21},{7,21,6},{21,7,8},{21,11,20},{11,21,10},{8,10,21},{10,8,9}};

int find_det(const int array_b[][2]){
	int det;
	det= array_b[0][0]*array_b[1][1]-array_b[1][0]*array_b[0][1];
	return det;
}

void find_inv_matrix(const int array_b[][2], int inv_array_b[][2]){
	int det;
	
	det=find_det(array_b);

	inv_array_b[0][0] = array_b[1][1]/det;
	inv_array_b[0][1] = -array_b[1][0]/det;
	inv_array_b[1][0] = -array_b[0][1]/det;
	inv_array_b[1][1] = array_b[0][0]/det;
}

void local_stiffness(const int array_b[][2], int inv_array_b[][2],float A_local[][3]) {
	int det,i,j;
	float area=0.5;
	int phib[2][3]={ {inv_array_b[0][0], inv_array_b[0][1], -inv_array_b[0][0]-inv_array_b[0][1]}, {inv_array_b[1][0], inv_array_b[1][1], -inv_array_b[1][0]-inv_array_b[1][1]} };
	
	det=find_det(array_b);
	for(i=0; i<3; i++){
		for(j=0;j<3;j++){
			A_local[i][j]= area*det*(phib[0][i]*phib[0][j] + phib[1][i]*phib[1][j]);
		}
	}
}

/*void construct_b(int x1, int x2, int x3){*/


