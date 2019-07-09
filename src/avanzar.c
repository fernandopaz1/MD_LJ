#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
#include "interaccion.h"
#include "avanzar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

// Definicion de funciones


int velocity_verlet(double *x, double *velocidad, double *fuerza, double *fuerza_vieja, double h, int N,double L){
	int i;

	for(i=0; i<3*N; i++){
		*(x+i)+=*(velocidad+i)*h+0.5*h*h*(*(fuerza+i));
		if(*(x+i)<0.0){
			*(x+i)+=L;
		}
		if(*(x+i)>L){
			*(x+i)-=L;
		}
		*(velocidad+i)+=0.5*h*((*(fuerza+i))+(*(fuerza_vieja+i)));
	}
return 0;
}


int paso(double *x, double *velocidad,double *fuerza, double *fuerza_vieja, double *potencial,int N, double rc, double V0, double L, double h){
	int i, j;

	

	for(i=0;i<N;i++){
		for(j=0;j<3;j++){	
			*(fuerza_vieja+3*i+j)=*(fuerza+3*i+j);
			*(fuerza+3*i+j)=0.0;}
		*(potencial+i)=0.0;
	}																											

		
	
	for(i=0;i<(N-1);i++){
		for(j=i+1;j<N;j++){
			Lenard_Jones(fuerza, potencial, 1.0, 1.0, rc, x, i, j, V0,L);
		}	
	}


	velocity_verlet(x,velocidad,fuerza,fuerza_vieja,h,N,L);


return 0;
}


int descorrelacionar(double *x, double *velocidad,double *fuerza, double *fuerza_vieja, double *potencial,int N, double rc, double V0, double L, double h){
	double a;

	a=1.0/sqrt(N);
	while(coeficiente_verlet(L,N,x)>a){
			paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);
	}

return 0;
}



double coeficiente_verlet(double L, int N, double *x){
    	int i, m;
    	double a, lambda_x, lambda_y, lambda_z, lambda, c;
    	double PI;
    	PI=3.14;
	m=cbrt(N);
    	a=L/m; 
    	c=(2*PI)/a; 
    	lambda_x=0.0;
    	lambda_y=0.0;
    	lambda_z=0.0;
       
    	for (i=0; i<N;i++){
    	      lambda_x+=cos(c*(*(x+3*i)-a/2));
    	      lambda_y+=cos(c*(*(x+3*i+1)-a/2));
    	      lambda_z+=cos(c*(*(x+3*i+2)-a/2));
    	}

    	lambda=(lambda_x+lambda_y+lambda_z)/(3*N); 

return lambda;

}


int simulacion(int N, int N_frames,double rc, double L, double h, double T){
	int l;
	double *x, *velocidad, *fuerza, *fuerza_vieja, *potencial,V0, Ec, Ep, *E_cinetica,*E_potencial,coef_verlet;


	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	E_cinetica=(double*)malloc(3*sizeof(double));
	E_potencial=(double*)malloc(sizeof(double));
	
	char filename[255];
  	sprintf(filename, "prueba.lammpstrj");

	V0=V0_LJ(rc,1.0, 1.0);


	for(T=0.728;T>0.2;T-=0.05){
	printf("%lf\n\n",T);
	set_box(x, N, L);
	set_v(velocidad,N,T);
	
	normalizacion_velocidad(velocidad, E_cinetica,T, N);

//	descorrelacionar(x,velocidad,fuerza,fuerza_vieja,potencial,N,rc,V0,L,h);
//////////////////////////////////////////////////////////////////////    Esta parte es para nombrar archivos
	char *nombre_energia= "Energia T= ";
	char *termalizado="termalizacion T= ";

	char *buffer=malloc(256);
	snprintf(buffer, 256, "%lf", T);

	size_t len1 = strlen(nombre_energia);
	size_t len2 = strlen(buffer);
	size_t len3 = strlen(termalizado);
	

	char *totalLine = malloc(len1 + len2 + 1);
	char *totalLine2 = malloc(len3 + len2 + 1);
	if (!totalLine) abort();
	if (!totalLine2) abort();

	memcpy(totalLine,        nombre_energia, len1);
	memcpy(totalLine + len1, buffer, len2);
	totalLine[len1 + len2] = '\0';

	memcpy(totalLine2,        termalizado, len3);
	memcpy(totalLine2 + len3, buffer, len2);
	totalLine[len3 + len2] = '\0';

/////////////////////////////////////////////////////////////////////	


	FILE *fe= fopen(totalLine, "w");


//	load_checkpoint(totalLine2, x, velocidad,N);
	
	

	for(l=0;l<N_frames;l++){
		

		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);
		

		Ec=energia_cinetica(velocidad,E_cinetica, N);

		Ep=energia_potencial(potencial,E_potencial, N, 1.0);
				
		coef_verlet=coeficiente_verlet(L,N,x);
		
		fprintf(fe, "%lf %lf %lf %lf\n" ,Ec,Ep,Ec+Ep,coef_verlet);

		save_lammpstrj(filename, x, velocidad, N, L, l);	

	}
	save_checkpoint(totalLine2,x,velocidad,N);
//	free(buffer);	
	free(totalLine);
	free(totalLine2);
	fclose(fe);
	}
		
free(x);
free(velocidad);
free(fuerza);
free(fuerza_vieja);
free(potencial);
free(E_cinetica);
return 0;
}


