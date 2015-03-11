double* EnumWPrun (int n, double* mu[], double b[]){
	double C = b[0];
	//O i começa a 0 pois e o primeiro indice
	int i = 0, j, k;
	double dist[n+1], c[n];
	double e[n+1][n];
	int delta[n], d[n], u[n], uL[n];
	int last_nonzero = 1;

	u[0] = uL[0] = 1;

	for(j = 1; j < n; j++){
		u[j] = 0;
		uL[j] = 0;
	}
	for(j = 0; j < n; j++){
		c[j] = 0;
		delta[j] = 0;
		d[j] = j + 1;
	}
	for(j = 0; j < n + 1; j++){
		for(k = 0; k < n; k++){
			e[j][k] = 0;
		}
	}

	while(true){
		dist[i] = dist[i+1] + (u[i] - c[i])^2 * b[i];
		if(dist[i] < C){
			if(i != 0){
				//move down
				i--;
				d[i - 1] = max(d[i - 1], d[i]);
				for(j = d[i]; j <= i+1; j--){
					e[j][i] = e[j+1][i] + u[j] * mu[j][i];
				}
				c[i] = -e[i+1][i];
				u[i] = round(c[i]);
				delta[i] = 1;
			}else{
				//update best vector
				C = dist[i];
				for(j = 0; j < n; j++){
					uL[j] = u[j];
				}
			}
		}else{
			if(i = n){
				return uL;
			}
			//move up
			i++;
			d[i - 1] = i;
			if(i >= last_nonzero){
				last_nonzero = i;
				u[i]++;
			}else{
				if(u[i] > c[i]){
					u[i] = u[i] - delta[i];
				}else{
					u[i] = u[i] + delta[i];
				}
				delta[i]++;
			}
		}
	}

}

