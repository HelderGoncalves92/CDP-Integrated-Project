//NOTA: Esta a receber k como se fosse k-j e assume-se que j=0
//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "math.h"

double* ENUM(double* b[], double c[], double* mu[], int k){
	double cL, cT[k+1], y[k+1], v[k+1];
	int delta[k+1], s = 0, t = 0, d[k+1], i, h;
	int u[k+1], uT[k+1];
	cL = c[0];
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

	while(t <= k){
		cT[t] = cT[t+1] + (y[t] + uT[t])^2 * c[t];
		if(cT[t] < cL){
			if(t > 0){
				t--;
				y[t] = 0;
				for(i = t+1; i <= s; i++){
					y[t] += uT[i], mu[i][t];
				}
				uT = v[t] = round(-y[t]);
				delta[t] = 0;
				if(uT[t] > -y[t]){
					d[t] = -1;
				}else{
					d[t] = 1;
				}
			}else{
				t++;
				s = fmax(s,t);
				if(t < s){
					delta[t] = -delta[t];
				}
				if(delta[t]*d[t] >= 0){
					delta[t] += d[t];
				}
				uT[t] = v[t] + delta[t];
			}
		}
	}
	double bnew = 0, aux = 0, aux2, retval;
	for(i = 0; i < k; i++){
		bnew += u[i]*b[i];
	}
	for(i = 0; i < k; i++){
		aux2 = 0;
		for(h = i; h < k; h++){
			aux2 += u[h], mu[h][i];
		}
		aux2 = pow(aux2,2);
		aux += aux2;
	}
	//Supostamente retorna este valor mais um tal "minimal place" que nao tou bem a ver o que e
	retval = fmin(bnew, aux);

	return retval;
}

