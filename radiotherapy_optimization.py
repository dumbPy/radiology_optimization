#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 14:47:13 2018
@author: dumbPy
git    : https://github.com/dumbPy
"""
#import pulp as pp
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance as dist
from math import exp
#from amplpy import AMPL
import pulp


radi     =  3  #Radius of the circular organ region   
grid_sz  =  6  # Length of the grid square area (grid size)
grid_res = 20  # 1 unit length to be devided in grid_res Number of pizels per length
beamlets = 20  # Number of beamlets per beam. All the beamlets in a beam have same angle
pos      =  5  # Total angular positions of the setup

#All length parameters will be multiplied with grid_res from hence forth
grid_sz *= grid_res
radi    *= grid_res
center   = (radi, radi)
# Get all angular positions    
positions = np.linspace(0, 360, num=pos, endpoint=False)


def IN_TUMOR(coordinates):
    """
    Uses 6x6 grid cordinates
    """
    i = coordinates[0]; j = coordinates[1]
    res = grid_res
    if i>=2.5*res and i<=3.5*res and j>=1*res and j<=5*res:
        return 1
    elif (((i>=1.5*res and i<=2.5*res) or(i>=3.5*res and i<=4.5*res)) and 
          ((j>=1*res and j<=2*res) or (j>=4*res and j<=5*res))):
        return 1
    else:
        return 0


def GET_Dijp(angle, beamlet_no):
    """
    Beamlet Number starts at 0 to n-1 for n beamlets
    Returns a D_ij_p matrix
    """
    def IN_COL(coordinates):
        i = coordinates[0]; j = coordinates[1]
        res = grid_res
#    if i>=2.975*res and i<=3.025*res:
        if i>=(2.5+beamlet_no/20)*res and i<=(2.55+beamlet_no/20)*res:
            return(exp(-1*j/(4*grid_res)))
        else:
            return 0

    col = np.asarray([[IN_COL((i, j)) 
                        for i in range(grid_sz)] for j in range(grid_sz)])

    import scipy.ndimage.interpolation as inter
    col = inter.rotate(col, angle, reshape = False)
    return col


map_risk = np.asarray([[1 if dist.euclidean((i, j), center)<= radi else 0
                        for i in range(grid_sz)] for j in range(grid_sz)])
map_tumor = np.asarray([[IN_TUMOR((i, j)) 
                        for i in range(grid_sz)] for j in range(grid_sz)])
map_risk -= map_tumor

dose_risk = np.asarray([[4 if dist.euclidean((i, j), center)<= radi else 0
                        for i in range(grid_sz)] for j in range(grid_sz)])
dose_tumor= np.asarray([[10*IN_TUMOR((i, j)) 
                        for i in range(grid_sz)] for j in range(grid_sz)])

#
dose_tumor_temp= np.asarray([[4*IN_TUMOR((i, j)) 
                        for i in range(grid_sz)] for j in range(grid_sz)])    
dose_risk -= dose_tumor_temp #Substract 4 value dosage from tumor area

#Final Dose Requirements in One Matrix
#dose in risk region = 4; dose in tumor = 10
dose = dose_risk+dose_tumor
    
"""
Optimization Code
"""

# =============================================================================
# problem.eval("""
# param theta; #A 3 dimentional theta for set T, R, N
# param SETS = 1..3;
# 
# set GRID_X;
# set GRID_Y;
# set BEAMLETS;
# set POS;  # Total number of beam angle positions to be used 
# 
# #For every positon, for every beamlet p,
# #we calculate intensity at every pixel (x, y)
# param D_ijp{POS, BEAMLETS, GRID_X, GRID_Y};
# 
# var w{POS, BEAMLETS} INTEGER >=0; #weights for each beamlet in each position
# var D_ij{POS, GRID_X, GRID_Y};
# 
# #A 3 channel T, R, N set's 1-0 mapping
# #for a channel, say R(represented by 1),
# #X, Y = 1 if (x, y) belongs to set R
# var MAP_TRN{3, GRID_X, GRID_Y}; 
# 
# #No Normal region. trn = [1, 2] for tumor and risk region
# minimize objective : sum{trn in 1..2}theta[trn]*sum{i in 1..GRID_X, j in 1..GRID_Y} MAP_TRN[trn, i, j]*sum{p in 1..POS}D_ij[p, i, j]
# 
# subject to con{i in 1..GRID_X, j in 1..GRID_Y, p in 1..POS}: D_ij[p, i, j] = sum{b in 1..BEAMLETS}D_ijp[p, b, i, j]
# 
# """)
# 
# =============================================================================

theta = [1, 10]
map_trn = [map_tumor, map_risk]
GRID_X = grid_sz*grid_res
GRID_Y = GRID_X
BEAMLETS = beamlets
POS = pos

#4D D_ijp [beam_position, beamlet, X, Y]
D_ijp = []
for beam_angle in positions:
    beam = [GET_Dijp(beam_angle, beamlet_no) for beamlet_no in range(20)]
    D_ijp.append(beam)
#dose matrix (in ampl model format)
d_ij = dose

final_map = dose_risk+dose_tumor
plt.imshow(final_map)

#prob = pp.LpProblem("Tomotherapy",pp.LpMinimize)
