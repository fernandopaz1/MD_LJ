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

datos=np.loadtxt('pruebagen')
mu=5.0
sigma=1.0



aleatorio=datos[:,0]
gaussiana=datos[:,1]
norma=datos[:,2]

plt.figure(1)
plt.hist(aleatorio,bins=500,label='probabilidad',normed=True)
plt.xlabel('probabilidad')
plt.ylabel('frecuencias')
plt.legend()
plt.show()

x=np.linspace(mu-sigma,mu+sigma,1000)

plt.figure(2)
plt.hist(gaussiana,bins=500,label='Gaussiana',normed=True)
plt.plot(x,(np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma)))/np.sqrt(2*np.pi*sigma*sigma),label='Gaussiana')
plt.xlabel('$x$')
plt.ylabel('Frecuencia')
plt.legend()
plt.show()


plt.figure(3)
plt.plot(norma,label='Norma')
plt.xlabel('paso')
plt.ylabel('norma')
plt.legend()
plt.show()

