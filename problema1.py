#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 27 02:04:18 2019

@author: paz
"""


import numpy as np
import math as mt
import matplotlib.pyplot as plt

plt.ion()

N=125
h=0.0001
a=1.0/np.sqrt(N)


T_ini=0.4
T_fin=2.0
delta_T=0.1

len_T=int((T_fin-T_ini)/delta_T)

T_array=np.linspace(T_ini,T_fin,len_T)
tau=np.zeros(len(T_array))
term=np.zeros(len(T_array))

energia_media=np.zeros(len_T)
energia_cuad_med=np.zeros(len_T)
fluctuaciones=np.zeros(len_T)

def corrnp(X,i):
    """
    llama a la funcion corrcoef de numpy, pero
    la hice para que le des la notacion simplificada
    de Correlacion (la casera). Tiene la ventana
    de que tarda 10 veces menos que la casera."""
    if i == 0:
        return np.corrcoef(X,X)[1,0]
    else:
        return np.corrcoef(X[:-i],X[i:])[1,0]

def Correlador(A):
    promedios = 20#int(input("\nCuantos promedios queres hacer?\n"))
    len_data = len(A[:,0])  # longitud de la tira de datos, original
    # redefino A sin esos valores
    A = A[int(len_data/10):]
    len_data = len(A[:,0]) # redefino len_data sin el 10% inicial
    # En estas condiciones, a los datos restantes quiero calcularles la corre.
    # para eso voy a querer armar una matriz en la que cada fila sea un cacho
    # del vector original, para luego promediar las correlaciones.
    # 
    # quiero promediar 20 veces, asi que divido el vector restante en 20
    # voy a armar una matriz de NxM (N columnas, M filas)
    
    M = promedios #esto podría ser numero de promedios
    N = int(len_data/M)
    matriz = np.zeros([M,N])
    for i in range(M):
        for j in range(N):
            matriz[i,j] = A[i*M + j, 1]
    c_k_matrix = np.zeros([M,int(N/2)])
    for k in range(M):
        c_k = []
        for i in range(int(N/2)):
            c_k.append(corrnp(matriz[k, :],i))
        c_k_matrix[k] = c_k
    
    return np.mean(c_k_matrix, axis = 0)


def Tau(A):
    valor=0
    for i in range(len(A)-1):
        if abs(A[i]) > 0.1 and abs(A[i+1]) < 0.1:
            valor = i+1
        
    return valor

T=T_ini
j=0
while(T<T_fin):

	#energias=np.loadtxt('/home/paz/MD_LJ/archivos/problema1/Energia T= {}'.format('%.6f' %T))
	energias=np.loadtxt('Energia T= {}'.format('%.6f' %T))

	Ec=energias[:,0]
	Ep=energias[:,1]
	E=energias[:,2]

	verlet=energias[:,3]

	A=np.zeros((len(E),2))

	A[:,0]=E[:]
	A[:,1]=E*E

	
	i=0
	while(verlet[i]>a):
		indice_term=i
		i=i+1
	term[j]=indice_term


	###################################################
	#
	#	Fluctuaciones temporales de la energía
	#
	###################################################


	energia_media[j]=np.mean(E[indice_term:])
	energia_cuad_med[j]=0

	for i in range(indice_term,len(E)):
		energia_cuad_med[j]=energia_cuad_med[j]+E[i]*E[i]
	energia_cuad_med[j]=energia_cuad_med[j]/(len(E)-indice_term)

	fluctuaciones[j]=np.sqrt(energia_cuad_med[j]-energia_media[j]*energia_media[j])

	print("Las fluctuaciones de la energía son:  ",fluctuaciones[j])


	###################################################
	#
	#	Energías
	#
	###################################################

	c=Correlador(A)
	tau[j]=Tau(c)

	t=np.linspace(0,h*len(Ec),len(Ec))

	if(j%2==0):
		plt.figure(2)
		line1,=plt.plot(verlet,label='T={}, Pasos={}'.format('%.2f' %T, '%d' %indice_term))
		#line2,=plt.plot([t[indice_term],t[indice_term]],[min(verlet),verlet[indice_term]])
		line3,=plt.plot([0,len(E)],[verlet[indice_term],verlet[indice_term]],linestyle='--')
		line2,=plt.plot([0,len(E)],[-verlet[indice_term],-verlet[indice_term]],linestyle='--')
		plt.xlabel('Tiempo')
		plt.ylabel('Coeficiente de Verlet')
		plt.legend()
		
		plt.figure(1)
		line1,=plt.plot(t,E,label='Energía $\\sigma_E={}$'.format('%.4f' %fluctuaciones[j]))
		line2,=plt.plot(t,Ec,label='Energía Cinética')
		line3,=plt.plot(t,Ep,label='Energía Potencial')
		line4,=plt.plot([t[indice_term],t[len(E)-1]],[energia_media[j],energia_media[j]],label='Energía Media = {}'.format('%.4f' %energia_media[j]))
		plt.xlabel('Tiempo')
		plt.ylabel('Energia')
		plt.legend()
		plt.show()
		

		plt.figure(7)
		line1=plt.plot(c,label='Correlacion T$^*={}$  $\\tau={}$'.format('%.2f' %T,'%d' %tau[j]))
		plt.xlabel('k')
		plt.ylabel('$\\rho(k)$')
		#plt.xlim([-100,3000])
		#plt.ylim([-2,2])
		plt.legend()
		plt.show()


	plt.show()
	T=T+delta_T
	j=j+1

tau_max=np.max(tau)

plt.figure(5)
plt.scatter(T_array,tau,label='$max(\\tau)$={}'.format('%d' %tau_max))
plt.xlabel('Temperatura')
plt.ylabel('Pasos')
plt.legend()
plt.show

term_max=np.max(term)

plt.figure(5)
plt.scatter(T_array,term,label='Pasos termalizacion max={}'.format('%d' %term_max))
plt.xlabel('Temperatura')
plt.ylabel('Pasos')
plt.legend()
plt.show

p=np.polyfit(T_array,energia_media,1)
plt.figure(3)
line1=plt.scatter(T_array,energia_media,label='Energía por partícula')
line2=plt.plot(T_array,np.polyval(p,T_array),label='$c_v={}$'.format('%.4f' %p[0]))
plt.xlabel('Temperatura')
plt.ylabel('Energía')
plt.legend()
plt.show()

plt.figure(4)
line1=plt.scatter(T_array,fluctuaciones,label='Energía por partícula')
plt.xlabel('Temperatura')
plt.ylabel('fluctuaciones')
plt.legend()
plt.show()



