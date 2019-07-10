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
h=0.00005
a=1.0/np.sqrt(N)


T_ini=0.4
T_fin=2.0
delta_T=0.05

densidad=0.6

len_T=int((T_fin-T_ini)/delta_T)

T_array=np.linspace(T_ini,T_fin,len_T)

energia_media=np.zeros(len_T)
energia_cuad_med=np.zeros(len_T)
fluctuaciones=np.zeros(len_T)


T=T_ini
j=0
while(T<T_fin):

#	energias=np.loadtxt('/home/paz/MD_LJ/archivos/problema1/Energia T= {}'.format('%.6f' %T))
	energias=np.loadtxt('Energia T= {} densidad= {}'.format('%.6f' %T,'%.6f' %densidad))

	Ec=energias[:,0]
	Ep=energias[:,1]
	E=energias[:,2]

	verlet=energias[:,3]


	
	i=0
	indice_term=0
	while(verlet[i]>a):
		indice_term=i
		i=i+1



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

	#print("Las fluctuaciones de la energía son:  ",fluctuaciones[j])


	###################################################
	#
	#	Energías
	#
	###################################################

	

	t=np.linspace(0,h*len(Ec),len(Ec))

	if(j%10==0):
		plt.figure(2)
		line1,=plt.plot(t,verlet,label='T={}, Tiempo={}'.format('%.2f' %T, '%.2f' %t[indice_term]))
		line2,=plt.plot([t[indice_term],t[indice_term]],[min(verlet),verlet[indice_term]])
		line3,=plt.plot([0,h*len(E)],[verlet[indice_term],verlet[indice_term]])
		plt.xlabel('Tiempo')
		plt.ylabel('Coeficiente de Verlet')
		plt.legend()
		plt.show()
	
		plt.figure(1)
		line1,=plt.plot(t,E,label='Energía $\\sigma_E={}$'.format('%.4f' %fluctuaciones[j]))
		line2,=plt.plot(t,Ec,label='Energía Cinética')
		line3,=plt.plot(t,Ep,label='Energía Potencial')
		line4,=plt.plot([t[indice_term],t[len(E)-1]],[energia_media[j],energia_media[j]],label='Energía Media = {}'.format('%.4f' %energia_media[j]))
		plt.xlabel('Tiempo')
		plt.ylabel('Energia')
		plt.legend()
		plt.show()

	
	T=T+delta_T
	j=j+1





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
