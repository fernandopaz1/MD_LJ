#include "general.h"
#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Definicion de funciones

double delta(float *x, int i, int j, double L,double *delta_x){
	int k;
	double distancia_cuad;

	for(k=0;k<3;k++){
		*(delta_x+k)=(double)(*(x+3*i+k)-*(x+3*j+k));
		*(delta_x+k)-=L*floor((*(delta_x+k)/L)*2);   // si la distancia es mayor a la de la caja interactua con la imagen
	}
	distancia_cuad=norma(delta_x,3);
return distancia_cuad;
}


int Lenard_Jones(double *fuerzas, double *potencial, double e, double sigma, double rc, float *x, int i, int j, double V0,double L){		
	double distancia_cuad, r_sext, r_doc,*delta_x,a;
	int k;

	delta_x=(double*)malloc(3*sizeof(double));

	distancia_cuad=delta(x,i,j,L,delta_x);

	if(distancia_cuad<rc){
		r_sext=pow(sigma/distancia_cuad,3);
		r_doc=pow(r_sext,2);
		*potencial=4.0*e*(r_doc-r_sext);
		for(k=0;k<3;k++){
			a=6.0*e*(*(delta_x+k)/distancia_cuad)*((*potencial));
			*(fuerzas+3*i+k)+=a;
			*(fuerzas+3*j+k)-=a;	
		}
		*potencial-=V0;
	}
free(delta_x);
return 0;
}

double V0_LJ(double rc,double sigma, double e){
	double r_sext, r_doc;
	
	r_sext=pow(sigma/(rc*rc),3);
	r_doc=pow(r_sext,2);

return 4.0*e*(r_doc-r_sext);
}

int tabla(double e, double sigma, double paso){
	double r,r_cuad, r_sext, r_doc,V;
	
	FILE *finterp= fopen("tabla", "w");
	
	r=0.001;
	for(r=0.001;r<2.56;r+=paso){
		r_cuad=r*r;
		r_sext=pow(sigma/r_cuad,3);
		r_doc=pow(r_sext,2);
		V=4.0*e*(r_doc-r_sext);
		fprintf(finterp, "%lf %lf %lf\n",r,r_cuad,V);
	}

fclose(finterp);
return 0;
}




int Lenard_Jones_interp(double *fuerzas, double *potencial, double e, double sigma, double rc, float *x, int i, int j, double V0, double *finterp, double paso, double L){		
	double distancia_cuad,*delta_x,a;
	int indice,k;	

	delta_x=(double*)malloc(3*sizeof(double));

	distancia_cuad=delta(x,i,j,L,delta_x);

	if(distancia_cuad<rc){	
		
		indice=(int)((distancia_cuad-(*finterp+1))/paso);
	
		*potencial=*(finterp+3*indice+2);
		for(k=0;k<3;k++){
			a=6.0*e*(*(delta_x+k)/distancia_cuad)*(*potencial);
			*(fuerzas+3*i+k)+=a;
			*(fuerzas+3*j+k)-=a;	
		}
		*potencial-=V0;
	}
free(delta_x);
return 0;
}


