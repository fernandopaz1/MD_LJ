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
    	c=(4*PI)/a; 
    	lambda_x=0.0;
    	lambda_y=0.0;
    	lambda_z=0.0;
       
    	for (i=0; i<N;i++){
    	      lambda_x+=cos(c*(*(x+3*i)));
    	      lambda_y+=cos(c*(*(x+3*i+1)));
    	      lambda_z+=cos(c*(*(x+3*i+2)));
    	}

    	lambda=(lambda_x+lambda_y+lambda_z)/(3*N); 

return lambda;

}


int problema1(int N, int N_frames,double rc, double L, double h, double T){
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
	T=0.4;
	//int i=0;
	for(T=0.4000;T<2.0;T+=0.1){
//	if(i%5!=0 || i>0.728){
	printf("%lf\n\n",T);
	set_box(x, N, L);
	set_v(velocidad,N,T);
	
//	normalizacion_velocidad(velocidad, E_cinetica,T, N);

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
	free(totalLine);	free(totalLine2);
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




int problema2(int N, int N_frames,double rc, double L, double h, double T){
	int i,l,Q=100;
	double *x, *velocidad, *fuerza, *fuerza_vieja, *potencial,V0, Ec, Ep, *E_cinetica,*E_potencial,coef_verlet,densidad,P_exceso,*drij;


	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	E_cinetica=(double*)malloc(3*sizeof(double));
	E_potencial=(double*)malloc(sizeof(double));
	double *distrad = (double *)malloc(Q*sizeof(double));
	drij=(double*)malloc(N*N*sizeof(double));
	
	char filename[255];
  	sprintf(filename, "prueba.lammpstrj");

	V0=V0_LJ(rc,1.0, 1.0);

	densidad=N/(L*L*L);

	
////////////////////////////////////////////////////////////////////////// Nombrando la densidad
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

	set_box(x, N, L);
	set_v(velocidad,N,T);

/////////////////////////////////////////////////////////////////////////
	//int i=0;
	for(T=2.0;T>0.4;T-=0.05){
	printf("%lf\n\n",T);

	
	normalizacion_velocidad(velocidad, E_cinetica,T, N);
	
	for(i=0;i<15000;i++){
		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);	
	}
//////////////////////////////////////////////////////////////////////    Esta parte es para nombrar archivos
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
	totalLine[len3 + len2+len4] = '\0';

/////////////////////////////////////////////////////////////////////	


	FILE *fe= fopen(totalLine, "w");


	//load_checkpoint(totalLine2, x, velocidad,N);
	
	

	for(l=0;l<N_frames;l++){
		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);
		

		Ec=energia_cinetica(velocidad,E_cinetica, N);

		Ep=energia_potencial(potencial,E_potencial, N, 1.0);
				
		coef_verlet=coeficiente_verlet(L,N,x);
	
		P_exceso=presion(Ec,L,N,fuerza,x,rc);
		
		fprintf(fe, "%lf %lf %lf %lf %lf %lf\n" ,Ec,Ep,Ec+Ep,coef_verlet,P_exceso,lyderman(L,N,x,drij));
		if(l%10==0){
			save_lammpstrj(filename, x, velocidad, N, L, l);	
		}
	}
	save_checkpoint(totalLine2,x,velocidad,N);
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
free(distrad);
free(drij);
return 0;
}


double presion(double Ec,double L, int N,double *fuerza,double *x,double rc){
	int i,j,k;
	double *delta_x,p_exceso,distancia_cuad;

	delta_x=(double*)malloc(3*sizeof(double));
	double rc_cuad=rc*rc;
	p_exceso=0.0;

	for(i=0;i<(N-1);i++){
		for(j=i+1;j<N;j++){
			distancia_cuad=delta(x,i,j,L,delta_x);
			if(distancia_cuad<rc_cuad){
				for(k=0;k<3;k++){
					p_exceso+=*(delta_x+k)*(*(fuerza+3*i+k));
				}	
			}

		}
	
	}
free(delta_x);
return p_exceso/((double)N);
}

double distrib_radial(double *distrad, double *x, double N, double L, double densidad, double Q) {
    int bin;
    double rij, rij2, dR1;
    double dr[3];

    dR1 = L/((double)(Q + 2)); // longitud de un bin

    for(int i=0; i<N-1; i++) {

        for(int j=i+1; j<N; j++) {

            rij2 = 0;

            for(int k=0; k<3; k++) {
                // calcula los dk con k = (x, y, z)
                dr[k] = x[i*3+k] - x[j*3+k];

                // condiciones de contorno para dk
                if(dr[k] > 0.5*L){
                    dr[k] -= L;
                }
                else if(dr[k] < -(0.5*L)){
                    dr[k] += L;
                }

                // suma las diferencias cuadradas
                rij2 += dr[k] * dr[k];
            }

            // calcula el módulo de la distancia
            rij = sqrt(rij2);

            bin = floor(rij/dR1);
            distrad[bin] += 1.0 / (4 * M_PI * rij * rij * dR1 * densidad);

        }
    }

    return 0;
}



double lyderman(double L, int N,double *x, double *drij){
	int i,j;
	double *delta_x,distancia_cuad,r,ly,media,media_cuad;

	delta_x=(double*)malloc(3*sizeof(double));
	
	media=0.0;
	media_cuad=0.0;

	for(i=0;i<(N-1);i++){
		for(j=i+1;j<N;j++){
			distancia_cuad=delta(x,i,j,L,delta_x);
			r=sqrt(distancia_cuad);
			media+=r;
			media_cuad+=distancia_cuad;	
		}
	
	}
	ly=sqrt(media_cuad-(media)*(media))/((media)*((double)(N*(N-1))));
free(delta_x);
return ly;
}

/*

int correlacion(int N,int N_frames,double rc,double L,double h,double T){
	int i, j, k, N_prom,N_div,N_total;
	double *E_cuad,*delta_E,*E_acum,porcion,*rho;
	double 	*E,*x,*velocidad,*fuerza,*fuerza_vieja,*potencial,*E_cinetica,*E_potencial;
	
	FILE *fterm= fopen("correlacion1a", "w");
	

	N_prom=dim*dim;
		
	N_div=100000;	
	N_total=10*N_div;
	porcion=N_div/20;

	E=(double*)malloc(sizeof(double));
	rho=(double*)malloc(porcion*sizeof(double));
	E_acum=(double*)malloc(N_total*sizeof(double));
	E_cuad=(double*)malloc(sizeof(double));
	E=(double*)malloc(sizeof(double));
		
	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	E_cinetica=(double*)malloc(3*sizeof(double));
	E_potencial=(double*)malloc(sizeof(double));

	V0=V0_LJ(rc,1.0, 1.0);

	//int i=0;
	for(T=0.4000;T<2.0;T+=0.01){
	printf("%lf\n\n",T);
	set_box(x, N, L);
	set_v(velocidad,N,T);
	
	normalizacion_velocidad(velocidad, E_cinetica,T, N);

	descorrelacionar(x,velocidad,fuerza,fuerza_vieja,potencial,N,rc,V0,L,h);


	for(l=0;l<N_frames;l++){
		

		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);
		

		Ec=energia_cinetica(velocidad,E_cinetica, N);

		Ep=energia_potencial(potencial,E_potencial, N, 1.0);
		
		*(E+i)+=Ec+Ep;
	}


	T=0.4;
	
	

	for(J=0.1;J<0.6;J=J+0.1){

		fprintf(fterm, "%lf",J);	
		p_int(J,p_i);
		
		for(i=0;i<10*N_div;i++){
			flip_spin(red,dim,spin,delta_E,J, 0.0, p_l,p_i,m,H);    //termzlizacion
				}
		for(i=0;i<10;i++){
			media_cuad=mean_cuad(m_acum, i*N_div,N_div)/a_cuad;
			desvio=std_dev(m_acum,i*N_div, N_div,media_cuad)/a_cuad;
			for(j=0;j<porcion;j++){
				*(rho_m+j)=*(rho_m+j)+correl(m_acum,i*N_div, j, N_div,porcion, media_cuad, desvio)/10.0;
			}			
		}		
		i=0;
		for(j=0;j<porcion;j++){
				fprintf(fterm, " %lf",*(rho_m+j));
			}
			
		fprintf(fterm,"\n");

	}

fclose(fterm);
return 0;
}

*/
