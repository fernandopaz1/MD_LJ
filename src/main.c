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

int problema3(int N, int N_frames,double rc, double L, double h, double T);

int main(int argc, char *argv[]){

	double total_time;
	clock_t start, end;
	start = clock();


	
	srand(time(NULL));
	

	int N = 125;
	double rho=0.6;
	double L=cbrt(N/rho);
	double T=1.5;
	int N_frames = 100000;
	//int N_frames = 100;
	double h=0.0001;
	double rc=2.56;


	problema1(N, N_frames,rc, L, h, T);
//	problema2(N, N_frames,rc, L, h, T);
//	problema3(N, N_frames,rc, L, h, T);
	

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


int problema3(int N, int N_frames,double rc, double L, double h, double T){
	int l,Q=400,i,N_sum=3;
	double *x, *velocidad, *fuerza, *fuerza_vieja, *potencial,V0,densidad,*E_cinetica,dr;


	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	E_cinetica=(double*)malloc(3*sizeof(double));
	double *distrad = (double *)malloc(Q*sizeof(double));
	double *distsum = (double *)malloc(Q*sizeof(double));
	
	char filename[255];
  	sprintf(filename, "prueba.lammpstrj");

	V0=V0_LJ(rc,1.0, 1.0);

	densidad=N/(L*L*L);
	

	for(densidad=0.6000;densidad<1.1;densidad+=0.2){
	printf("%lf\n\n",densidad);
	
	L=cbrt(N/densidad);
	
//////////////////////////////////////////////////////////////////////    Esta parte es para nombrar archivos
	char *dens=" densidad= ";	
	char *buffer2=malloc(256);
	snprintf(buffer2, 256, "%lf", densidad);
	
	size_t len20 = strlen(buffer2);
	size_t len30 = strlen(dens);
	
	char *totalLine20 = malloc(len30 + len20 + 1);
	if (!totalLine20) abort();
	memcpy(totalLine20,        dens, len30);
	memcpy(totalLine20 + len30, buffer2, len20);
	totalLine20[len30 + len20] = '\0';
	char *nombre_energia= "Energia T= ";
	char *termalizado="termalizacion T= ";
	char *buffer=malloc(256);
	snprintf(buffer, 256, "%lf", T);
	size_t len1 = strlen(nombre_energia);
	size_t len2 = strlen(buffer);
	size_t len3 = strlen(termalizado);
	size_t len4 = strlen(totalLine20);
	char *totalLine = malloc(len1 + len2 +len4+ 1);
	char *totalLine2 = malloc(len3 + len2 +len4+ 1);
	if (!totalLine) abort();
	if (!totalLine2) abort();
	memcpy(totalLine,        nombre_energia, len1);
	memcpy(totalLine + len1, buffer, len2);
	memcpy(totalLine + len1+len2, totalLine20, len4);
	totalLine[len1 + len2+len4] = '\0';
	memcpy(totalLine2,        termalizado, len3);
	memcpy(totalLine2 + len3, buffer, len2);
	memcpy(totalLine2 + len3+len2, totalLine20, len4);
	totalLine2[len3 + len2+len4] = '\0';
/////////////////////////////////////////////////////////////////////	

	set_box(x, N, L);
	set_v(velocidad,N,T);

	

//	load_checkpoint(totalLine2, x, velocidad,N);

	
	normalizacion_velocidad(velocidad, E_cinetica,T, N);

	printf("\nHasta aca\n");	

	for(l=0;l<100000;l++){
		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);		
	}
	
	

	
//	save_checkpoint(totalLine2,x,velocidad,N);
	
	dr=L/((double)(Q + 2));
	FILE *fg= fopen(totalLine20, "w");
	i=0;
	while(i<N_sum){
		for(l=0;l<15000;l++){
			paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);		
		}
		distrib_radial(distrad, x,N,L,densidad,Q);
		for(i=0;i<Q;i++){
			*(distsum+i)+=*(distrad+i);
		}
		i++;
	}
	for(i=0;i<Q;i++){
		fprintf(fg, "%lf %lf\n" ,dr*(1.0*i+1.0),*(distsum+i)/((double)N_sum));
	}
	
	fclose(fg);
	free(totalLine);
	free(totalLine2);
	free(totalLine20);	
	}
		
free(x);
free(velocidad);
free(fuerza);
free(fuerza_vieja);
free(potencial);
free(E_cinetica);
free(distrad);
return 0;
}

