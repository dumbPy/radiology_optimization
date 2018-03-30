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
import pulp as pp


radi     =  3  #Radius of the circular organ region   
grid_sz  =  6  # Length of the grid square area (grid size)
grid_res = 20  # 1 unit length to be devided in grid_res Number of pizels per length
beamlets = 20  # Number of beamlets per beam. All the beamlets in a beam have same angle
pos      =  2  # Total angular positions of the setup

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



#4D D_ijp [beam_position, beamlet, X, Y]
D_ijp = []
for beam_angle in positions:
    beam = [GET_Dijp(beam_angle, beamlet_no) for beamlet_no in range(20)]
    D_ijp.append(beam)
D_ijp = np.asarray(D_ijp)
   
"""
Optimization Code
"""

prob = pp.LpProblem("Tomotherapy", pp.LpMinimize)
# =============================================================================
# w = [[]*beamlets]*positions
# for i in range(positions):
#     for j in range(beamlets):
# =============================================================================
w = np.asarray([pp.LpVariable(f'w[{i}][{j}]', lowBound=0, upBound=1) 
                for i in range(pos) for j in range(beamlets)])

w=pp.LpVariable.dicts('weights',[(i, j) for i in range(pos)
                        for j in range(beamlets)],0, None)
    
#D_ij = np.asarray([[pp.LpVariable(f'D[{i}][{j}]', lowBound=0) for i in range(grid_sz)]for j in range(grid_sz)])


theta = [1, 10]
map_trn = [map_tumor, map_risk]
GRID_X = grid_sz
GRID_Y = GRID_X
BEAMLETS = beamlets
POS = pos

D_ij=pp.LpVariable.dicts('total intensity',[(i, j) for i in range(grid_sz)
                                      for j in range(grid_sz)],0, None, pp.LpContinuous)

y=pp.LpVariable.dicts('absolute variable',[(i, j) for i in range(grid_sz)
                                        for j in range(grid_sz)],0, None, pp.LpContinuous)
#Objective Function

prob += pp.lpSum(D_ij[(i, j)]-dose[i][j] for i in range(grid_sz) for j in range(grid_sz))
for i in range(grid_sz):
    for j in range(grid_sz):
        prob += pp.lpSum(D_ijp[beam][beamlet][i][j]*w[(beam, beamlet)] for beam in range(pos) for beamlet in range(beamlets)) == D_ij[(i, j)]


actual_dose= np.asarray([[D_ij[(i, j)].varValue 
                        for i in range(grid_sz)] for j in range(grid_sz)])    

#Objective Function

prob += np.sum(np.asarray([theta[i]*np.multiply(map_trn[i], (D_ij-dose)) for i in [0, 1]])),
prob += sum(np.asarray([D_ijp[beam][beamlet]*w[beam][beamlet] for beam in range(pos) for beamlet in range(beamlets)])) == D_ij,




#dose matrix (in ampl model format)

final_map = dose_risk+dose_tumor+actual_dose
plt.imshow(final_map)

#prob = pp.LpProblem("Tomotherapy",pp.LpMinimize)
