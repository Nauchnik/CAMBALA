# z 0:0.5:3000
# zbar = 2*(z-1300)/1300;
# eps = 0.00737;
# c = 1500*(1 + eps*(zbar - 1 + exp(-zbar)));

import numpy as np 
import math 

depths = []
d = 5
while d <= 3000:
	depths.append(d)
	d += 5

ext_depths = depths
ext_depths.insert(0,0)
#print(ext_depths)
cw = []
for d in ext_depths:
	zbar = 2*(d-1300)/1300
	eps = 0.00737
	c = 1500*(1 + eps*(zbar - 1 + np.exp(-zbar)))
	cw.append(c)

with open('log_cw', 'w') as ofile:
	for c in cw:
		ofile.write(str(c) + ' ')

with open('log_c1s', 'w') as ofile:
	for c in cw[:-1]:
		ofile.write(str(c) + ' ')

with open('log_c2s', 'w') as ofile:
	for c in cw[1:]:
		ofile.write(str(c) + ' ')

rhos = [1 for d in depths]
betas = [0 for d in depths]
	
with open('log_rhos', 'w') as ofile:
	for rho in rhos:
		ofile.write(str(rho) + ' ')	
		
with open('log_betas', 'w') as ofile:
	for beta in betas:
		ofile.write(str(beta) + ' ')	

