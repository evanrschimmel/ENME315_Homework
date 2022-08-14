# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 22:14:16 2022

@author: schimmer
"""

import matplotlib.pyplot as plt
import numpy as np

Tb = 350  # base temperature, degC
Tinf = 25 # ambient temperature, degC
D = 0.01 # fin diameter, m
L = 0.05 # fin length, m
h = 250 # convection heat transfer coefficient, W/m^2*K
k = 240 # thermal conductivity, W/m*K
Ac = np.pi * (D * D) / 4 # cross-sectional area, m^2
P = np.pi * D # perimeter, m
m = np.sqrt((h * P) / (k * Ac)) # m-value, m^-1
M = np.sqrt(h * P * k * Ac) * (Tb - Tinf) # M-value, W

x = np.linspace(0,L,1000)

# Infinity Boundary Condition
Q_inffin = M
T_inffin = (Tb - Tinf) * np.exp(-m*x) + Tinf

# Adiabatic Boundary Condition
Q_adia = M * np.tanh(m * L)
T_adia = (np.cosh(m * L - m * x) / np.cosh(m * L)) * (Tb - Tinf) + Tinf

# Defined Tip Temperature Boundary Condition
Ttip = 250 # tip temperature, degC
Q_def = M * (np.cosh(m * L) - ((Ttip - Tinf) / (Tb - Tinf))) / np.sinh(m * L)
T_def = ((((Ttip - Tinf) / (Tb - Tinf)) * np.sinh(m * x) + np.sinh(m * L - m * x)) / (np.sinh(m * L))) * (Tb - Tinf) + Tinf

# Convection Boundary Condition
Q_conv = M * ((np.sinh(m * L) + (h / (m * k)) * np.cosh(m * L)) / (np.cosh(m * L) + (h / (m * k)) * np.sinh(m * L)))
T_conv = ((np.cosh(m * L - m * x) + (h / (m * k)) * np.sinh(m * L - m * x)) / (np.cosh(m * L) + (h / (m * k)) * np.sinh(m * L))) * (Tb - Tinf) + Tinf


plt.plot(x*1000,T_inffin,label='Infinity BC')
plt.plot(x*1000,T_adia,label='Adiabatic BC')
plt.plot(x*1000,T_def,label="Defined Tip BC")
plt.plot(x*1000,T_conv,label="Convection BC")
plt.title('Temperature Along Cylindrical Fin')
plt.xlabel('Length [mm]')
plt.ylabel('Temperature [degC]')
plt.legend()
plt.show()