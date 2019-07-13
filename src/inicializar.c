#include "general.h"
#include "inicializar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int set_box(double *x,int N,double L);
int set_v(double *velocidad,int N,double T);

// Definicion de funciones
int set_box(double *x,int N,double L){
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

int set_v(double *velocidad,int N,double T){
	int i,k;
	double sigma=sqrt(T);
	for(i=0;i<3*N;i++){
		*(velocidad+i)=gaussiana(0.0,sigma);
	}
	double Vcm[3]={0.0,0.0,0.0};
	for(i=0;i<N;i++){
		for(k=0;k<3;k++){
			Vcm[k]+=*(velocidad+3*i+k)/N;
		}
	}
	for(i=0;i<N;i++){
		for(k=0;k<3;k++){
			*(velocidad+3*i+k)-=Vcm[k];
		}
	}

return 0;
}


int save_checkpoint(char *nombre, double *x, double *velocidad,int N){
	int i;

	FILE *ftermalizado= fopen(nombre, "w");

	for(i=0;i<3*N;i++){
		fprintf(ftermalizado, "%lf %lf\n" ,*(x+i),*(velocidad+i));
	}

fclose(ftermalizado);
return 0;
}


int load_checkpoint(char *nombre, double *x, double *velocidad,int N){
//	int i;

	FILE *ftermalizado= fopen(nombre, "r");

	//for(i=0;i<3*N;i++){
		if(!fscanf(ftermalizado, "%lf%lf", x,velocidad)){printf("Error en load checkpoin");}
	//	if(!fscanf(ftermalizado, "%lf", (+2*i+1))){break;}
	//}

fclose(ftermalizado);
return 0;
}

