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

T=1.5


densidad=0.6

V=N/densidad
i=1
while(i<4):
	g=np.loadtxt(' densidad= {}'.format('%.6f' %densidad))


	r=g[:,0]
	g_r=g[:,1]


	plt.figure(i)
	plt.plot(r,g_r,label='Distribución radial $\\rho$= {}'.format('%.2f' %densidad))
	plt.legend()
	plt.xlabel('r')
	plt.ylabel('g(r)')
	plt.show()	
	
	densidad=densidad+0.2
	i=i+1
