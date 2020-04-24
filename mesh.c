#include <stdio.h>

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


int main()
{
	int i, j, k;
	int side1[24][2], side2[24][2], side3[24][2];
	int side4[48][2], side5[48][2], side6[48][2];
	double nodes_up[95][2];
	int res1, res2;
	int nodes_ptr[21][21] = {0};
	int flag;
	int connect_up[96][3] = {0};
	FILE *fp;


	for (i = 0; i < 95; i++) {
		for (j = 0; j < 2; j++) {
			nodes_up[i][j] = 100;
		}
	}

	for (i = 0; i < 24; i++) {
		for (j = 0; j < 2; j++) {
			side1[i][j] = connectivity[i][j];
		}
		side2[i][0] = connectivity[i][1];
		side2[i][1] = connectivity[i][2];
		side3[i][0] = connectivity[i][2];
		side3[i][1] = connectivity[i][0];
	}

	for (i = 0; i < 21; i++) {
		for (j = 0; j < 2; j++) {
			nodes_up[i][j] = nodes[i][j];
		}
	}


	j = 21;
	for (i = 0; i < 24; i++) {
		/*First side*/
		res1 = side1[i][0];
		res2 = side1[i][1];
		flag = 0;
		if (i == 0) {
			nodes_up[j][0] =
				(nodes[res1 - 1][0] + nodes[res2 - 1][0]) / 2;
			nodes_up[j][1] =
				(nodes[res1 - 1][1] + nodes[res2 - 1][1]) / 2;
			j++;
			nodes_ptr[res1 - 1][res2 - 1] = j;
			side4[i][0] = res1;
			side4[i][1] = j;
			side4[47 - i][0] = j;
			side4[47 - i][1] = res2;
		} else {
			for (k = 0; k < i; k++) {
				if ( (side1[k][0] == res2 &&
				      side1[k][1] == res1)  ) {
					flag = 1;
				} else
				if ( (side2[k][0] == res2 &&
				      side2[k][1] == res1)  ) {
					printf("found on side2 line 62\n");
				} else
				if ( (side3[k][0] == res2 &&
				      side3[k][1] == res1)  ) {
					printf("found on side3 line 65\n");
				}
			}
			if (flag == 1) {
				side4[i][0] = res1;
				side4[i][1] = nodes_ptr[res2 - 1][res1 - 1];
				side4[47 -
				      i][0] = nodes_ptr[res2 - 1][res1 - 1];
				side4[47 - i][1] = res2;
			} else {
				nodes_up[j][0] =
					(nodes[res1 - 1][0] +
					 nodes[res2 - 1][0]) / 2;
				nodes_up[j][1] =
					(nodes[res1 - 1][1] +
					 nodes[res2 - 1][1]) / 2;
				j++;
				nodes_ptr[res1 - 1][res2 - 1] = j;
				side4[i][0] = res1;
				side4[i][1] = j;
				side4[47 - i][0] = j;
				side4[47 - i][1] = res2;
			}
		}

		/*Side two */
		res1 = side2[i][0];
		res2 = side2[i][1];
		flag = 0;
		if (i == 0) {
			nodes_up[j][0] =
				(nodes[res1 - 1][0] + nodes[res2 - 1][0]) / 2;
			nodes_up[j][1] =
				(nodes[res1 - 1][1] + nodes[res2 - 1][1]) / 2;
			j++;
			nodes_ptr[res1 - 1][res2 - 1] = j;
			side5[i][0] = res1;
			side5[i][1] = j;
			side5[47 - i][0] = j;
			side5[47 - i][1] = res2;
		} else {
			for (k = 0; k < i; k++) {
				if ( (side2[k][0] == res2 &&
				      side2[k][1] == res1)  ) {
					flag = 1;
				} else
				if ( (side3[k][0] == res2 &&
				      side3[k][1] == res1)  ) {
					flag = 1;
				}
			}
			if (flag == 1) {
				side5[i][0] = res1;
				side5[i][1] = nodes_ptr[res2 - 1][res1 - 1];
				side5[47 -
				      i][0] = nodes_ptr[res2 - 1][res1 - 1];
				side5[47 - i][1] = res2;
			} else {
				nodes_up[j][0] =
					(nodes[res1 - 1][0] +
					 nodes[res2 - 1][0]) / 2;
				nodes_up[j][1] =
					(nodes[res1 - 1][1] +
					 nodes[res2 - 1][1]) / 2;
				j++;
				nodes_ptr[res1 - 1][res2 - 1] = j;
				side5[i][0] = res1;
				side5[i][1] = j;
				side5[47 - i][0] = j;
				side5[47 - i][1] = res2;
			}
		}

		/*Third side */
		res1 = side3[i][0];
		res2 = side3[i][1];
		flag = 0;
		if (i == 0) {
			nodes_up[j][0] =
				(nodes[res1 - 1][0] + nodes[res2 - 1][0]) / 2;
			nodes_up[j][1] =
				(nodes[res1 - 1][1] + nodes[res2 - 1][1]) / 2;
			j++;
			nodes_ptr[res1 - 1][res2 - 1] = j;
			side6[i][0] = res1;
			side6[i][1] = j;
			side6[47 - i][0] = j;
			side6[47 - i][1] = res2;
		} else {
			for (k = 0; k < i; k++) {
				if ( (side3[k][0] == res2 &&
				      side3[k][1] == res1)  ) {
					flag = 1;
				}
			}
			if (flag == 1) {
				side6[i][0] = res1;
				side6[i][1] = nodes_ptr[res2 - 1][res1 - 1];
				side6[47 -
				      i][0] = nodes_ptr[res2 - 1][res1 - 1];
				side6[47 - i][1] = res2;
			} else {
				nodes_up[j][0] =
					(nodes[res1 - 1][0] +
					 nodes[res2 - 1][0]) / 2;
				nodes_up[j][1] =
					(nodes[res1 - 1][1] +
					 nodes[res2 - 1][1]) / 2;
				j++;
				nodes_ptr[res1 - 1][res2 - 1] = j;
				side6[i][0] = res1;
				side6[i][1] = j;
				side6[47 - i][0] = j;
				side6[47 - i][1] = res2;
			}
		}
	}
	j = 0;
	for (i = 0; i < 24; i++) {
		connect_up[j][0] = side4[i][0];
		connect_up[j][1] = side4[i][1];
		connect_up[j][2] = side6[i][1];
		j++;
		connect_up[j][0] = side4[47 - i][0];
		connect_up[j][1] = side4[47 - i][1];
		connect_up[j][2] = side5[47 - i][0];
		j++;
		connect_up[j][0] = side5[47 - i][0];
		connect_up[j][1] = side5[47 - i][1];
		connect_up[j][2] = side6[47 - i][0];
		j++;
		connect_up[j][0] = side4[i][1];
		connect_up[j][1] = side5[i][1];
		connect_up[j][2] = side6[i][1];
		j++;
	}

	fp = fopen("connectivety_elems.dat", "w+");
	for (j = 0; j < 96; j++) {
		fprintf(fp,
			"%d %d %d\n",
			connect_up[j][0],
			connect_up[j][1],
			connect_up[j][2]);
	}
	fclose(fp);
	fp = fopen("nodes_cord.dat", "w+");
	for (j = 0; j < 65; j++) {
		fprintf(fp, "%lf %lf\n", nodes_up[j][0], nodes_up[j][1]);
	}
	fclose(fp);

	return 0;
}
