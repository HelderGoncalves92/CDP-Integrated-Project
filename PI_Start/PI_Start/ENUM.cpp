//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "ENUM.h"
#include "math.h"
int* ENUM(double* c, double** mu, int ini, int fim){
	double cL, *cT, *y, *v;
	int *delta, s = ini, t = ini, *d, i, h;
	int *u, *uT;

	//Alocate memory for every vector
	cT = (double*)_mm_malloc(((dimVector + 1)) * sizeof(double), 64);
	y = (double*)_mm_malloc((dimVector + 1) * sizeof(double), 64);
	v = (double*)_mm_malloc((dimVector + 1) * sizeof(double), 64);
	delta = (int*)_mm_malloc((dimVector + 1) * sizeof(int), 64);
	d = (int*)_mm_malloc((dimVector + 1) * sizeof(int), 64);
	u = (int*)_mm_malloc((dimVector + 1) * sizeof(int), 64);
	uT = (int*)_mm_malloc((dimVector + 1) * sizeof(int), 64);

	cL = c[ini];
	u[ini] = uT[ini] = 1;
	y[ini] = delta[ini] = 0;
	d[ini] = 1;

	for(i = ini + 1; i <= fim; i++){
		cT[i] = 0;
		u[i] = 0;
		uT[i] = 0;
		y[i] = 0;
		delta[i] = 0;
		d[i] = 1;
	}

	while(t <= fim){
		cT[t] = cT[t + 1] + pow(y[t] + uT[t],2) * c[t];
		if (cT[t] < cL){
			if (t > ini){
				t--;
				y[t] = 0;
				for (i = t + 1; i <= s; i++){
					y[t] += uT[i], mu[i][t];
				}
				uT[t] = v[t] = round(-y[t]);
				delta[t] = 0;
				if (uT[t] > -y[t]){
					d[t] = -1;
				}
				else{
					d[t] = 1;
				}
			}
			else{
				cL = cT[ini];
				for (i = ini; i < fim; i++){
					u[i] = uT[i];
				}
			}
		}
		else{
			t++;
			s = maxn(s,t);
			if(t < s){
				delta[t] = -delta[t];
			}
			if(delta[t]*d[t] >= 0){
				delta[t] += d[t];
			}
			uT[t] = v[t] + delta[t];
		}
	}
	
	return uT;
}

