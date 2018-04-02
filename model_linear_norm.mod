

set TRN;
set GRID_X;		#set of all x cordinates
set GRID_Y;		#set of all y cordinates
set BEAMLETS;	#set of all beamlet numbers 0 to 19
set POS;  		#set of beam angle positions

#A 2 dimentional theta for set T(Tumor), R(Risk)
param theta{TRN};
#For every positon, for every beamlet p,
#we calculate intensity at every pixel (x, y)
param D_ijp{POS, BEAMLETS, GRID_X, GRID_Y};
#Desired Dose at each pixel
param d_ij{GRID_X, GRID_Y};
#A 3 channel T, R, N set's 1-0 mapping (assuming TRN is of size 3)
#(we use trn of size 2 here as we donot consider normal region)
#for a channel, say R(represented by 1),
#X, Y = 1 if (x, y) belongs to set R
param MAP_TRN{TRN, GRID_X, GRID_Y}; 
param MAP_RISK{GRID_X, GRID_Y};


#weights for each beamlet in each position
var w{POS, BEAMLETS};
#Actual Dosage after optimization
var D_ij{GRID_X, GRID_Y};
var error{GRID_X, GRID_Y};


#No Normal region. trn = [1, 2] for tumor and risk region
minimize objective : sum{trn in TRN, i in GRID_X, j in GRID_Y} theta[trn]*MAP_TRN[trn, i, j]*error[i, j];


subject to con1{i in GRID_X, j in GRID_Y}: D_ij[i, j] == sum{p in POS, b in BEAMLETS} w[p, b]*D_ijp[p, b, i, j];
subject to con2{p in POS, b in BEAMLETS}:              w[p, b] >= 0;
#subject to con3{p in POS, b in BEAMLETS}:              w[p, b] <= 1;
#Error in risk region -ve, i.e., actual dose below 4
subject to con4{i in GRID_X, j in GRID_Y}: MAP_RISK[i, j]*error[i, j] >= 0;

#error at each pixel is actual dose - max/min dose (10 or 4 for tumor and risk resp.)
#Linearised norm of error
subject to con5{i in GRID_X, j in GRID_Y}: error[i, j] >= d_ij[i, j]-D_ij[i, j];
subject to con5{i in GRID_X, j in GRID_Y}: error[i, j] >= D_ij[i, j]-d_ij[i, j];