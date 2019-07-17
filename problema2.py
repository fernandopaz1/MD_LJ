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


T_ini=0.4+0.05
T_fin=2.0
delta_T=0.05

densidad=0.6

for i in range(0,3):
	densidad=0.6+i*0.2
 

	V=N/densidad

	len_T=int(np.floor(T_fin-T_ini)/delta_T)

	T_array=np.zeros(len_T-1)

	cv=np.zeros(len_T-1)

	energia_media=np.zeros(len_T-1)
	energia_cuad_med=np.zeros(len_T-1)
	fluctuaciones=np.zeros(len_T-1)
	p_exceso_medio=np.zeros(len_T-1)

	presion=np.zeros(len_T-1)

	indice_term=0
	T=T_ini
	j=0
	while(j<len_T-1):

		energias=np.loadtxt('Energia T= {} densidad= {}'.format('%.6f' %T,'%.6f' %densidad))
	#	energias=np.loadtxt('/home/paz/MD_LJ/archivos/problema2/Energia T= {} densidad= {}'.format('%.6f' %T,'%.6f' %densidad))
		T_array[j]=T
		Ec=energias[:,0]
		Ep=energias[:,1]
		E=energias[:,2]
		P_exceso=energias[:,3]
		p_exceso_medio[j]=np.mean(P_exceso)








		###################################################
		#
		#	Fluctuaciones temporales de la energía
		#
		###################################################


		energia_media[j]=np.mean(E[indice_term:])
		energia_cuad_med[j]=np.mean(E[indice_term:]*E[indice_term:])
		auxiliar=np.mean(1/Ec[indice_term:])
		aux2=np.mean(Ec[indice_term:])
		T_array[j]=(3/2)*aux2

		cv[j]=1/(N-T_array[j]*((1.5*N-1))*auxiliar)
		presion[j]=(2/3)*densidad*np.mean(Ec[indice_term:])+(1/(3*V))*p_exceso_medio[j]

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
	line1=plt.scatter(T_array,energia_media,label='$\\rho={}$'.format('%.2f' %densidad))
	line2=plt.plot(T_array,np.polyval(p,T_array),label='$c_v={}$'.format('%.4f' %p[0]))
	plt.xlabel('Temperatura')
	plt.ylabel('Energía')
	plt.legend()
	plt.show()

	plt.figure(4)
	line1=plt.scatter(T_array,fluctuaciones,label='$\\rho={}$'.format('%.2f' %densidad))
	plt.xlabel('Temperatura')
	plt.ylabel('fluctuaciones')
	plt.legend()
	plt.show()


	plt.figure(5)
	line1=plt.scatter(T_array,presion/densidad,label='$\\rho={}$'.format('%.2f' %densidad))
	plt.xlabel('Temperatura')
	plt.ylabel('presion')
	plt.legend()
	plt.show()


	plt.figure(6)
	line1=plt.scatter(T_array,p_exceso_medio/(3*V),label='$\\rho={}$'.format('%.2f' %densidad))
	plt.xlabel('Temperatura')
	plt.ylabel('Presion de exceso')
	plt.legend()
	plt.show()



	plt.figure(7)
	line1=plt.scatter(T_array,cv,label='cv $\\rho={}$'.format('%.2f' %densidad))
	plt.xlabel('Temperatura')
	plt.ylabel('cv')
	plt.legend()
	plt.show()
