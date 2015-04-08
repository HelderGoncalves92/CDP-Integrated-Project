//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "math.h"
int* ENUM(double** b, double* c, double** mu, int ini, int fim){
	int k = fim - ini;
	double cL, cT[k+1], y[k+1], v[k+1];
	int delta[k+1], s = ini, t = ini, d[k+1], i, h;
	int u[k+1], uT[k+1];
	cL = c[ini];
	u[0] = uT[0] = 1;
	y[0] = delta[0] = 0;
	d[0] = 1;

	for(i = 1; i <= k; i++){
		cT[i] = 0;
		u[i] = 0;
		uT[i] = 0;
		y[i] = 0;
		delta[i] = 0;
		d[i] = 1;
	}

	while(t <= fim){
		cT[t - ini] = cT[t - ini + 1] + pow(y[t - ini] + uT[t - ini],2) * c[t];
		if (cT[t - ini] < cL){
			if (t > ini){
				t--;
				y[t - 1] = 0;
				for (i = t + 1; i <= s; i++){
					y[t - ini] += uT[i - ini], mu[i][t];
				}
				uT[t - ini] = v[t - ini] = round(-y[t - ini]);
				delta[t - ini] = 0;
				if (uT[t - ini] > -y[t - ini]){
					d[t - ini] = -1;
				}
				else{
					d[t - ini] = 1;
				}
			}
			else{
				cL = cT[0];
				for (i = 0; i < k; i++){
					u[i] = uT[i];
				}
			}
		}
		else{
			t++;
			s = fmax(s,t);
			if(t < s){
				delta[t - ini] = -delta[t - ini];
			}
			if(delta[t - ini]*d[t - ini] >= 0){
				delta[t - ini] += d[t - ini];
			}
			uT[t - ini] = v[t - ini] + delta[t - ini];
		}
	}
	
	double bnew = 0, aux = 0, aux2, retval;
	for(i = 0; i < k; i++){
		bnew += u[i]*b[i];
	}
	for(i = 0; i < k; i++){
		aux2 = 0;
		for(h = i; h < k; h++){
			aux2 += u[h] * mu[h][i];
		}
		aux2 = pow(aux2,2);
		aux += aux2;
	}
	//Supostamente retorna este valor mais um tal "minimal place" que nao tou bem a ver o que e
	retval = fmin(bnew, aux);

	return retval;
}

