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
	z=sqrt(12*n)*(z/((double)n)-0.5);
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

void histograma(double *y,double *x,int n,double a,double b,int m)
{

  /* 
     Histograma del vector de datos x[0]...x[n-1] 
     
     M  = cantidad de columnas del histograma
     a,b= limites del del histograma 
     
     y[0]...y[m-1]  valores de las columnas (normalizados)
     y[m]...y[2m-1] valores de las marcas de clase

  */

  int    i,j;
  double h,hh,s;

  s=1.0/(double)n;
  h=(b-a)/(double)m;
  hh=h/2.0;

  for(i=0;i<m;i++)
    {
      *(y+i)=0.0;
      *(y+m+i)=(double)i*h+a+hh;
    }

  for(i=0;i<n;i++)
    {
      j=(int)floor((*(x+i)-a)/h);

      if (j<0) j=0;
      if (j>m) j=m-1;

      y[j]=y[j]+s;
    }
}

int minmax(double *vector,int N, double *min, double *max){
	int i;
	
	*min=*vector;
	*max=*vector;
	
	for(i=0;i<N;i++){
		if(*(vector+i)<(*min)){
			*min=*(vector+i);
		}
		if(*(vector+i)>(*max)){
			*max=*(vector+i);
		}
	}	
	
return 0;
}



double energia_cinetica(double *velocidad,double *E_cinetica, int N){
	int i,k;

	for(k=0;k<3;k++){
			*(E_cinetica+k)=0.0;
		}

	for(i=0;i<N;i+=3){
		for(k=0;k<3;k++){
			*(E_cinetica+k)+=0.5*pow(*(velocidad+3*i+k),2);
		}
	}	
	
	double a=(*(E_cinetica)+(*(E_cinetica+1))+(*(E_cinetica+2)))/((double)N);
return a;
}

double energia_potencial(double *potencial,double *E_potencial, int N, double masa){
	int i;

	*(E_potencial)=0.0;

	for(i=0;i<N;i++){
		*(E_potencial)+=0.5*(*(potencial+i));
	}	

return (*E_potencial)/((double)N);
}



double energia(double *velocidad, double *potencial,double *E_cinetica,double *E_potencial, int N, double masa){
	double Ec=energia_cinetica(velocidad,E_cinetica, N);
	double Ep=energia_potencial(potencial,E_potencial, N, masa);

return (Ec+Ep);
}

int normalizacion_velocidad(double *velocidad, double *E_cinetica,double T, int N){
	int i;
	double Tmedida=energia_cinetica(velocidad,E_cinetica, N);
	double factor=sqrt(T/Tmedida);
	for(i=0;i<3*N;i++){
		*(velocidad+i)=factor*(*(velocidad+i));
		
	}
return 0;
}


int pruebageneral(int N, double mu, double sigma,double radio){
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
