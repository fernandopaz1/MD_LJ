#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
#include "interaccion.h"
#include "avanzar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Definicion de funciones


int velocity_verlet(double *x, double *velocidad, double *fuerza, double *fuerza_vieja, double h, int N){
	int i;

	for(i=0; i<3*N; i++){
		*(x+i)+=*(velocidad+i)*h+0.5*h*h*(*(fuerza+i));
		*(velocidad+i)+=0.5*h*((*(fuerza+i))+(*(fuerza_vieja+i)));
	}
return 0;
}


int paso(double *x, double *velocidad,double *fuerza, double *fuerza_vieja, double *potencial,int N, double rc, double V0, double L, double h){
	int i, j;

	

	for(i=0;i<N;i++){	
		*(fuerza_vieja+3*i)=*(fuerza+3*i);
		*(fuerza+3*i)=0.0;
		*(potencial+i)=0.0;
	}																											

		
	
	for(i=0;i<(N-1);i++){
		for(j=i+1;j<N;j++){
			Lenard_Jones(fuerza, potencial, 1.0, 1.0, rc, x, i, j, V0,L);
		}	
	}


	velocity_verlet(x,velocidad,fuerza,fuerza_vieja,h,N);


return 0;
}

int simulacion(int N, int N_frames,double rc, double L, double h, double T){
	int l,i,M;
	double *x, *velocidad, *fuerza, *fuerza_vieja, *potencial,*E_cinetica,V0,*histo,*min,*max,a,b;
	

	M=3000; //cantidad de columnas del histograma 

	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	E_cinetica=(double*)malloc(3*sizeof(double));
	histo=(double*)malloc(2*M*sizeof(double));
	min=(double*)malloc(sizeof(double));
	max=(double*)malloc(sizeof(double));

	
	char filename[255];
  	sprintf(filename, "prueba.lammpstrj");


	FILE *fhisto= fopen("Histograma", "w");

	set_box(x, N, L);
	set_v(velocidad,N,T);
	
	normalizacion_velocidad(velocidad, E_cinetica,T, N);

	V0=V0_LJ(rc,1.0, 1.0);

	for(l=0;l<N_frames;l++){
		

		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);
		
		minmax(velocidad,3*N, min, max);		
		
		a=*min;
		b=*max;

		histograma(histo,velocidad,3*N,a,b,M);
		for(i=0;i<3*N;i++){ 
			sprintf(fhisto, "%lf" ,*(histo+i));}
		sprintf(fhisto, "\n" );
		save_lammpstrj(filename, x, velocidad, N, L, l);	

	}
	

free(x);
free(velocidad);
free(fuerza);
free(fuerza_vieja);
free(potencial);
free(E_cinetica);
free(fhisto);
return 0;
}


