

param grid_sz;
param beamlets;
param theta; #A 2 dimentional theta for set T(Tumor), R(Risk)

set GRID_X := 0..grid_sz;
set GRID_Y := 0..grid_sz;
set BEAMLETS := 1..beamlets;
set POS;  # Total number of beam angle positions to be used 

#For every positon, for every beamlet p,
#we calculate intensity at every pixel (x, y)
param D_ijp{POS, BEAMLETS, GRID_X, GRID_Y};
#Expected Dose at each pixel
param d_ij{GRID_X, GRID_Y};

#weights for each beamlet in each position
var w{POS, BEAMLETS};
#Actual Dosage after optimization
var D_ij{GRID_X, GRID_Y};

#A 3 channel T, R, N set's 1-0 mapping
#for a channel, say R(represented by 1),
#X, Y = 1 if (x, y) belongs to set R
var MAP_TRN{2, GRID_X, GRID_Y} integer >= 0; 

#No Normal region. trn = [1, 2] for tumor and risk region
minimize objective : sum{trn in 1..2}theta[trn]*sum{i in GRID_X, j in GRID_Y} MAP_TRN[trn, i, j]*(D_ij[i, j]-d_ij[i, j]);

subject to con1{i in GRID_X, j in GRID_Y, p in POS}: D_ij[p, i, j] = sum{b in BEAMLETS}D_ijp[p, b, i, j]
subject to con2{p in POS, b in BEAMLETS}:              w[p, b] >= 0;
subject to con3{p in POS, b in BEAMLETS}:              w[p, b] <= 1;
