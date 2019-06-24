#include "general.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Definicion de funciones
double aleatorio(void){
	return ((double)rand())/((double)RAND_MAX);
}


double gaussiana(double mu, double sigma){
	int i,n=10;
	double z=0.0;
	for(i=0;i<n;i++){
		z+=aleatorio();
	}
	z=sqrt(12+n)*(z/((double)n)-0.5);
return z*sigma+mu;
}

double norma(double *v, int longitud){
	int i;
	double acum=0.0;
	for(i=0;i<longitud;i++){
		acum+=*(v+i)*(*(v+i));
	}
return acum;
}

int pruebageneral(int N, double mu, double sigma,double radio){
	printf("Hasta aca");	
	int i;
	double *v,fase;
	
	v=(double*)malloc(2*sizeof(double));
	FILE *fgen= fopen("pruebagen", "w");
	for(i=0;i<N;i++){
		fase=0.005*((double)i);
		*v=radio*cos(fase);
		*(v+1)=radio*sin(fase);
		fprintf(fgen, "%lf %lf %lf\n",aleatorio(),gaussiana(mu,sigma),norma(v,2));
	}
free(v);
fclose(fgen);
return 0;
}
