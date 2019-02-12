# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:48:40 2018

Test script for 2D model

@author: olem
"""

import model_2D as mod
import scipy as sci
import numpy as np
import matplotlib.pyplot as plt

aao = np.linspace(-np.pi,np.pi,200)
z = np.linspace(60e3,120e3,60)
VAT = np.zeros((aao.size,z.size))

atm = mod.atmosphere(aao,z,VAT)

scene = mod.scene(atm)

aao_sat = np.deg2rad(0)
z_sat = 600e3

z_tan = 60e3
sat = mod.satellite(aao_sat,z_sat,z_tan)
scene.add_sensor(sat)

z_tan = 75e3
sat = mod.satellite(aao_sat,z_sat,z_tan)
scene.add_sensor(sat)

z_tan = 80e3
sat = mod.satellite(aao_sat,z_sat,z_tan)
scene.add_sensor(sat)

[x,y] = scene.calc_path(scene.satellites[0],100)

L = scene.calc_radiance(scene.satellites[0],100)
                   
K = scene.calc_jacobian(10)