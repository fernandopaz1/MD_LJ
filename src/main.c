#include "general.h"
#include "interaccion.h"
#include "inicializar.h"
#include "avanzar.h"
#include "visualizacion.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

int problema3(int N, int N_frames,int N_correlacion,double rc, double L, double h, double T);

int main(int argc, char *argv[]){

	double total_time;
	clock_t start, end;
	start = clock();


	
	srand(time(NULL));
	

	int N = 125;
	double rho=0.1;
	double L=cbrt(N/rho);
	double T=1.5;
	int N_frames = 60001;
	//int N_frames = 100;
	double h=0.00005;
	double rc=2.56;
//	int N_prom=10;
	int N_correlacion=500;

//	problema1(N, N_frames,rc, L, h, T);
//	problema2(N, N_frames,rc, L, h, T);
//	problema2(N,N_frames,N_prom,N_correlacion,rc,L,h,T);
	problema3(N, N_frames,N_correlacion,rc, L, h, T);
	

/*
	char filename[255];
  	sprintf(filename, "prueba");
	
	double *x = (double*)malloc(3*N*sizeof(double));
	double *velocidad = (double*)malloc(3*N*sizeof(double));
	
	set_box(x, N, L);
	set_v(velocidad,N,T);
	save_checkpoint(filename,x,velocidad,N);

*/
	end = clock();								//time count stops 


	total_time = ((double) (((double)(end - start)) / (double)CLOCKS_PER_SEC)/60.0);	
	printf("\nEl tiempo (minutos) requerido es:  %lf \n", total_time);		//calulate total time
  return 0;
}


int problema3(int N, int N_frames,int N_correlacion,double rc, double L, double h, double T){
	int l;
	double *x, *velocidad, *fuerza, *fuerza_vieja, *potencial,V0,densidad,*p;


	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	double *E_cinetica=(double*)malloc(3*sizeof(double));
	p = (double *)malloc(sizeof(double));

	char filename[255];
  	sprintf(filename, "problema3 densidad=.lammpstrj");

	V0=V0_LJ(rc,1.0, 1.0);

	densidad=N/(L*L*L);
	

	printf("%lf\n\n",densidad);
	
	L=cbrt(N/densidad);
	

	set_box(x, N, L);
	set_v(velocidad,N,T);

	

//	load_checkpoint(totalLine2, x, velocidad,N);

	
	normalizacion_velocidad(velocidad, E_cinetica,T, N);

	//printf("\nHasta aca\n");	

	for(l=0;l<6000;l++){
		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h,p);		
	}
	
	int j=0;
	for(l=0;l<N_frames;l++){
		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h,p);
		j++;
		if(j==N_correlacion){
			j=0;
			save_lammpstrj(filename, x, velocidad, N, L, l);
		}
	}
		
free(x);
free(velocidad);
free(fuerza);
free(fuerza_vieja);
free(potencial);
free(E_cinetica);
return 0;
}

