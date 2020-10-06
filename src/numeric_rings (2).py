# -*- coding: utf-8 -*-
"""Module with functions for Numerical calculation of Rayleigh Integral using the Rings
"""
import math
import cmath
import numpy
import time
import logging
from multiprocessing import Pool, Manager
from datetime import timedelta
import ribs
import pickle
import scipy


def calc_numeric_ring():
#     raschet integrala Rayleigh s ploskoy i so sfericheskoy poverhnosti, garmonicheskaya volna
    #   INTEGER jtheta,jz,numtheta,numphi,numz,jr,numr, jphi
    #   REAL  a0,a0F,f0,c0,F,kF,theta0,z,hz,theta,dtheta,pi,rf,hr,r,
    #  &      phi, hphi         
    # REAL, ALLOCATABLE :: amp(:,:), intr(:)     
    #   COMPLEX :: i, pressure, e1
    
    # i=(0,1)
      # pi=4.0*atan(1.0)
# *****************************************************************
# c   pole na osi sfericheskoy chashki VG, vsyo v systeme CI
# c     piston source, velocity is given at spherical cap, Rayleigh
    F = 0.13
    a0 = 0.170
    f0 = 1000000.0
    c0 = 1500.0
    kF = 2.0*math.pi*f0*F/c0
    theta0 = math.asin(a0/F)
    numtheta = 250
    numz = 400
    numr = 100
    dtheta = theta0/numtheta
    hz = (0.2e-4)/F
    hr = (0.2e-4)/F
    # ALLOCATE(amp(0:numz,0:numr))
    # ALLOCATE(intr(0:numr))
    amp = []
    intr = []

    # OPEN(unit=50,file="intens.txt",STATUS= 'UNKNOWN')  
    for jz in range(numz):
        z = 0.05/0.07 + hz*jz
        for jr in range(-numr, numr):
            print(jz, jr)
            r = jr*hr
            pressure = 0.5*math.exp(1j*kF*math.sqrt(z*z+r*r))/math.sqrt(z*z+r*r)*dtheta
            for jtheta in range(1,numtheta-1):
                theta = jtheta*dtheta
                hphi = 2.0*math.pi/(2.0*jtheta+1)
                numphi = 2*jtheta+1
                e1 = 0.0*1j
                for jphi in range(0, numphi-1):
                    phi = hphi*jphi
                    rf = math.sqrt(math.sin(theta)**2.0+r*r-2.0*math.sin(theta)*r*math.cos(phi)+(z-1.0+math.cos(theta))**2.0)
                    e1 = e1+math.exp(1j*kF*rf)/rf
                pressure += e1*math.sin(theta)/(2.0*jtheta+1)
            amp[jz,jr] = cmath.abs(pressure)*kF*dtheta
            intr[jr] = amp[jz,jr]


if __name__ == '__main__':
    pass
    # profile_main_cycle()
    # test_slicing()
    # test case
    # trans = transducer.Transducer("Probe1", 1, 4.0e-03, 5.0e-03, 6.0e06)
    # tris = sliceTransOnTrisProper(trans, 100)

    # trans = transducer.Transducer(
    #   name = "Probe1",
    #   num_of_elements = 1,
    #   element_width = 4.0e-03,
    #   element_length = 5.0e-03,
    # frequency = 6.0e06
    # frequency = 12.0e06
    # frequency = 4.0e06
    #   frequency = 11.6e06
    #   )

    # water = medium.Medium(
    #   name = "Water",
    #   density = 1000.0,
    #   speed_of_sound = 1488.0
    #   )
