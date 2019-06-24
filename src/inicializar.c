#include "general.h"
#include "inicializar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Definicion de funciones
int set_box(float *x,int N,double L){
	int n=cbrt(N),i,j,k;
	double dl=L/n;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			for(k=0;k<n;k++){
				*(x+3*i+(3*n)*j+(3*n*n)*k)=dl*(i+0.5);
				*(x+3*i+(3*n)*j+(3*n*n)*k+1)=dl*(j+0.5);
				*(x+3*i+(3*n)*j+(3*n*n)*k+2)=dl*(k+0.5);
			}
		}
}

return 0;
}

int set_v(float *v,int N,double T){
	int i,k;
	double sigma=sqrt(T);
	for(i=0;i<3*N;i++){
		*(v+i)=gaussiana(0.0,sigma);
	}
	double Vcm[3]={0,0,0};
	for(i=0;i<N;i++){
		for(k=0;k<3;k++){
			Vcm[k]+=*(v+3*i+k)/N;
		}
	}
	for(i=0;i<N;i++){
		for(k=0;k<3;k++){
			*(v+3*i+k)-=Vcm[k];
		}
	}
return 0;
}

