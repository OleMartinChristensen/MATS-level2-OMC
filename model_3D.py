# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:46:51 2018

2D jacobian calulation on spherical geometry

@author: olem
"""

import numpy as np
import scipy as sci
from itertools import compress
import math as m

def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az

def sph2cart(r,elev,az):
    x = r*np.cos(elev)*np.cos(az)
    y = r*np.cos(elev)*np.sin(az)
    z = r*np.sin(elev)
    
    return x, y, z

class position(object):
    
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    @classmethod
    def from_spherical(cls,r,elev,az):
        
        x,y,z = sph2cart(r,elev,az)
        return cls(x,y,z)
    
    def get_cart(self):
        
        return self.x,self.y,self.z
    
    def get_spherical(self):
        
        r,elev,az = cart2sph(self.x,self.y,self.z)
        
        return self.r,self.elev,self.az
     
    def __add__(self,other):
        xx = self.x + other.x
        yy = self.y + other.y
        zz = self.z + other.z
        
        return position(xx,yy,zz)
    
    def __sub__(self,other):
        xx = self.x - other.x
        yy = self.y - other.y
        zz = self.z - other.z
        
        return position(xx,yy,zz)    
    
class look_vector(object):
    
    def __init__(self,dx,dy,dz):
        l = np.sqrt(dx**2 + dy**2 + dz**2)
        self.dx = dx/l
        self.dy = dy/l
        self.dz = dz/l
    
    @classmethod
    def from_tanpos(self,satpos,tanpos):
        lookpos = tanpos-satpos
        return look_vector(lookpos.x,lookpos.y,lookpos.z)       
    
    def from_az_el(self,satpos,azimuth,elevation):
        
        
        return look_vector(lookpos.x,lookpos.y,lookpos.z)      

class satellite(object):
    
    """Satellite object containing position and look vector of satellite"""

    def __init__(self, pos,look):
        self.aao = aao
        self.z = z_sat
        self.z_tan = z_tan

        self.spheroid = spheroid()
        self.r_tan = z_tan + self.spheroid.r

        self.r = z_sat+self.spheroid.r

        self.x = np.cos(self.aao)*self.r
        self.y = np.sin(self.aao)*self.r        
        
        self.r_tan = self.spheroid.r+z_tan
        self.aao_tan = self.aao - np.arccos(self.r_tan/self.r)

        self.x_tan = np.cos(self.aao_tan)*self.r_tan
        self.y_tan = np.sin(self.aao_tan)*self.r_tan 


        x_look_full = self.x_tan-self.x
        y_look_full = self.y_tan-self.y
        
        self.x_look = x_look_full/np.sqrt(x_look_full**2+y_look_full**2)
        self.y_look = y_look_full/np.sqrt(x_look_full**2+y_look_full**2)
        
class spheroid(object):
    
    r = 6371e3 #radius of earth in m    
    
class atmosphere(object):
    
    """A 2D grid consisting for a circular geometry defined through the angle from equator
    and altitudes above the spheriod

    Attributes:
        aao[np.array]: angle from equator along orbit axis (rad)
        z[np.array]: altitude above spheriod (m)
        
        b[np.array]: distance from equator along circumference at reference radius
        r[np.array]: distance from center of earth
        
        aao_grid: a 2D array of the grid positions
        z_grid: a 2D array of the grid positions
    """

    def __init__(self, aao, z,VAT):
        
        #Make some checks on the input
        
        self.aao = aao
        self.z = z
               
        self.spheroid = spheroid()
        self.r = z+self.spheroid.r
        self.r_ref = self.r[0]    #reference radius where distance b is calculated               
        self.b = aao*self.r_ref
        
        self.aao_grid,self.r_grid = np.meshgrid(self.aao, self.r)
        self.x_grid = np.cos(self.aao_grid)*self.r_grid
        self.y_grid = np.sin(self.aao_grid)*self.r_grid

        self.VAT = VAT+1

class scene(object):
    
    def __init__(self, atmosphere: atmosphere):
        self.atmosphere = atmosphere
        self.satellites = []
        
    def add_sensor(self,satellite: satellite):
        
        self.satellites.append(satellite)
    
    def calc_path(self,satellite: satellite,dl):
        
        z_max = np.max(self.atmosphere.z)
        r_top = z_max + self.atmosphere.spheroid.r
        #Note satellite must be above the atmosphere. 
        init_aao_pos = satellite.aao_tan + np.arccos(satellite.r_tan/r_top)
        end_aao_pos = satellite.aao_tan - np.arccos(satellite.r_tan/r_top)
        
        x_0,y_0 = pol2cart(init_aao_pos,r_top)
        x_end,y_end = pol2cart(end_aao_pos,r_top)
        
        L_tot = np.sqrt((x_end-x_0)**2 + (y_end-y_0)**2)
        L = np.arange(0,L_tot,dl)
        
        path_x = x_0 + L*satellite.x_look
        path_y = y_0 + L*satellite.y_look
        
        return path_x,path_y
        
    def calc_radiance(self,satellite: satellite,dl):   
        
        x,y = self.calc_path(satellite,dl)
        aao,r = cart2pol(x,y)
        q_points = np.array([aao,r]).T
        Vq = sci.interpolate.interpn((self.atmosphere.aao, self.atmosphere.r), self.atmosphere.VAT, q_points, 
                 method='nearest', bounds_error=True, fill_value=0)
        L_tot = np.sum(Vq)*dl
        return L_tot
        
    def calc_jacobian_row(self,satellite: satellite,dl):   
        x,y = self.calc_path(satellite,dl)
        aao,r = cart2pol(x,y)
        aaoindex = np.zeros(aao.size,dtype=int)*-1
        rindex =np.zeros(aao.size,dtype=int)*-1
        linindex = np.zeros(aao.size,dtype=int)*-1
                           
        for i in np.arange(aao.size):
            aaoindex[i] = list(compress(range(len(self.atmosphere.aao)), self.atmosphere.aao>aao[i]))[0]
            rindex[i] = list(compress(range(len(self.atmosphere.r)), self.atmosphere.r<r[i]))[-1]
            linindex[i] = np.ravel_multi_index((aaoindex[i],rindex[i]),(self.atmosphere.aao.shape[0],self.atmosphere.r.shape[0]),order='F')
        
        index,value = np.unique(linindex,return_counts=True)
        return index,value
        
    def calc_jacobian(self,dl):
        
        K = np.zeros([len(self.satellites),self.atmosphere.aao.size*self.atmosphere.r.size])
        
        for i in np.arange(len(self.satellites)):
            index,value = self.calc_jacobian_row(self.satellites[i],dl)
            K[i,index] = value
        
        return K
