#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 14:49:27 2017

@author: anqil
"""

import math
import numpy as np

def domain_size(r_bot, r_top, r_tan, r_sat, h_field):
    d = math.sqrt(r_sat**2 - r_tan**2)
    theta_1 = math.atan(d/(r_tan-r_bot))
    theta_2 = math.asin(r_bot / r_top * math.sin(theta_1))
    theta_3 = math.pi - theta_1
    theta = math.acos(r_tan/r_sat)

    beta = math.pi - theta_2 - theta_3
    alpha = theta - beta
    gamma = math.pi - theta_1- theta_2

    y_rad_start = alpha
    y_rad_end = theta + gamma

    a = h_field**2 + d**2 + (r_bot - r_tan)**2
    b = 2*r_tan*(r_bot - r_tan) - 2*d**2
    c = d**2 + r_tan**2 - r_top**2
    k = -b/2/a + math.sqrt(b**2 - 4*a*c)/2/a
    x_rad = math.atan(k*h_field/(r_tan + k*(r_bot-r_tan)))
    
    return y_rad_start, y_rad_end, x_rad


def pathleng(heights, Xi):
    deltaz = heights[1:] - heights[:-1]
    deltaz = np.append(deltaz,deltaz[-1])
    heights = np.append(heights, heights[-1]+deltaz[-1])

    nheights = len(heights)
    Re = 6370 # km Earth radius

    if Xi==90:
        	Zt = heights
    else:
        	Zt = (Re + heights) * np.sin(Xi*np.pi/180) - Re

    pathl = np.zeros((nheights, nheights))
    for j in range(nheights):
        h = heights[j:] 
        Ztj = Zt[j]
        pathl[j,j:] = np.sqrt(h**2 + 2*Re*(h-Ztj) - Ztj**2)
        
    pathl[:,:-1] = pathl[:, 1:] - pathl[:, 0:-1]
    pathl = pathl[:-1, :-1]
    pathl = np.triu(pathl)
    heights = heights[:-1]
    nheights=nheights-1

#    if Xi>90:
#        for j in range(nheights):
#            if Zt[j] > 0:
#                I = find(heights<heights[j] & heights>Zt[j])
#                
#                if ((isempty(I))):
#                    I = max(1,j-1)
#                else:
#                    I=np.append(max(I[1]-1,1),I)
#                    
#                h=heights(I)+deltaz(I)
#                Ztj=Zt(j);
#                pathl(j,I) = sqrt(h.*h+2*Re*(h-Ztj)-Ztj*Ztj)
#                                
#                if (isempty(find(I==1, 1))):
#                    pathl(j,I)=(pathl(j,I)-pathl(j,max(I-1,1)))
#                else:
#                    J=I(I~=1)
#                    pathl(j,J)=(pathl(j,J)-pathl(j,max(J-1,1)))
#                        
#                
#            else Zt(j)<=0:
#                pathl(j,:) = zeros(size(pathl(j,:)))
#                
#
#	pathl1=fliplr(pathl)
#	nanregion=find(isnan(pathl)==1)
#	pathl2=(triu((pathl'),1))';
#	pathl2(nanregion')=zeros(size(nanregion'))
#    
#	pathlength=pathl+pathl2
#	columnpathl=[pathl1 pathl2]
#    
    return pathl


def satmove(V, x, y, z, ti, v_field, h_field, v_pix_num, h_pix_num, 
            r_bot, r_top, r_tan, r_sat):
    
    from scipy.interpolate import interpn
    G = 6.673e-20 # km3/(s2*kg2)
    Mearth = 5.98e24 # kg
    R = 6370 #km earth radius
    angular_speed = math.sqrt(G*Mearth/r_sat**3) #radians/s
    theta = math.acos(r_tan/r_sat)
    y_rad_start, y_rad_end, x_rad = domain_size(r_bot, r_top, r_tan, r_sat, h_field)
    Sat = np.array([0, 0, r_sat])  #initial satellite location (u,v,w) 

    Pu = np.linspace(-h_field, h_field, h_pix_num)
    Pv = np.linspace((r_tan-v_field)*math.sin(theta), 
                     (r_tan+v_field)*math.sin(theta), v_pix_num)
    Pw = np.linspace((r_tan-v_field)*math.cos(theta), 
                     (r_tan+v_field)*math.cos(theta), v_pix_num)

    k_start = (r_sat * math.sin(theta) - r_top * math.sin(theta 
               - y_rad_start)) / (r_sat * math.sin(theta)) 
    k_end = (r_sat * math.sin(theta) + r_top * math.sin(y_rad_end
             - theta)) / (r_sat * math.sin(theta))
    
    point_num_along_k = 3000
    
    k = np.linspace(k_start, k_end, point_num_along_k) #interval points along the ray
#    lu = np.empty((point_num_along_k,h_pix_num,v_pix_num)) # u coordinate of rays
#    lv = np.empty((point_num_along_k,h_pix_num,v_pix_num)) # v coordinate of rays
#    lw = np.empty((point_num_along_k,h_pix_num,v_pix_num)) # w coordinate of rays
#    for inter in range(point_num_along_k):
#        for i in range(h_pix_num):
#            for j in range(v_pix_num):
#                lu[inter,i,j] = Sat[0] + k[inter] * (Pu[i]-Sat[0])
#                lv[inter,i,j] = Sat[1] + k[inter] * (Pv[j]-Sat[1])
#                lw[inter,i,j] = Sat[2] + k[inter] * (Pw[j]-Sat[2])
    lu = Sat[0] + np.matmul(k[:,None], Pu[None,:] - Sat[0])
    lu = np.repeat(lu[:,:,None], v_pix_num, axis=2)
    lv = Sat[1] + np.matmul(k[:,None], Pv[None,:] - Sat[1])
    lv = np.repeat(lv[:,None,:], h_pix_num, axis=1)
    lw = Sat[2] + np.matmul(k[:,None], Pw[None,:] - Sat[2])
    lw = np.repeat(lw[:,None,:], h_pix_num, axis=1)
    
    dlu = np.diff(lu, axis=0)
    dlv = np.diff(lv, axis=0)
    dlw = np.diff(lw, axis=0)
    dl = np.sqrt(dlu**2 + dlv**2 + dlw**2) 
    dl = np.append(dl, dl[-1,None,:], axis=0)*1e5 #pathlength in cm
    

#=== change coordinate system of rays(u,v,w)->(alpha,beta,rho) ===
    alpha = np.arctan(lu / lw)
    beta = np.arctan(lv / lw)
    rho = np.sqrt(lu**2 + lv**2 + lw**2)
    
#=== Create coordinates of the query points(xq,yq,zq) === 
    xq = alpha * r_tan
    yq = beta * r_tan + ti * angular_speed * r_tan
    zq = rho - R
    
    q_points = np.array([xq.flatten(), yq.flatten(), zq.flatten()]).T
    Vq = interpn((x, y, z), V, q_points, 
                 method='nearest', bounds_error=False, fill_value=0)
    
    limb_radiance = np.sum(Vq.reshape(point_num_along_k, h_pix_num, v_pix_num) * dl,
                   axis=0)
    return limb_radiance, xq,yq,zq,dl
