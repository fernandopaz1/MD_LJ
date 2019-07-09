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

N=512
h=0.00005
#h=1

energias=np.loadtxt('Energia T= 0.728000')

Ec=energias[:,0]
Ep=energias[:,1]
E=energias[:,2]

verlet=energias[:,3]


a=1.0/np.sqrt(N)
i=0
while(verlet[i]>a):
	indice_term=i
	i=i+1



###################################################
#
#	Fluctuaciones temporales de la energía
#
###################################################


energia_media=np.mean(E[indice_term:])
energia_cuad_med=0

for i in range(indice_term,len(E)):
	energia_cuad_med=energia_cuad_med+E[i]*E[i]
energia_cuad_med=energia_cuad_med/(len(E)-indice_term)

fluctuaciones=np.sqrt(energia_cuad_med-energia_media*energia_media)

print("Las fluctuaciones de la energía son:  ",fluctuaciones)


###################################################
#
#	Energías
#
###################################################



t=np.linspace(0,h*len(Ec),len(Ec))

plt.figure(1)
line1,=plt.plot(t,E,label='Energía $\\sigma_E={}$'.format('%.4f' %fluctuaciones))
line2,=plt.plot(t,Ec,label='Energía Cinética')
line3,=plt.plot(t,Ep,label='Energía Potencial')
line4,=plt.plot([t[indice_term],t[len(E)-1]],[energia_media,energia_media],label='Energía Media = {}'.format('%.4f' %energia_media))
plt.xlabel('Tiempo')
plt.ylabel('Energia')
plt.legend()
plt.show()


plt.figure(2)
line1,=plt.plot(t,verlet,label='T=0.728  , h={}'.format('%.5f' %h))
line2,=plt.plot([t[indice_term],t[indice_term]],[min(verlet),verlet[indice_term]],label='Tiempo de Termalización = {}'.format('%.4f' %t[indice_term]))
line3,=plt.plot([0,h*len(E)],[verlet[indice_term],verlet[indice_term]],label='Coef. Verlet =$\\frac{1}{\\sqrt{512}}$')
plt.xlabel('Tiempo')
plt.ylabel('Coeficiente de Verlet')
plt.legend()
plt.show()


