import numpy as np
#import matplotlib.pyplot as plt
from math import exp, log, floor, fabs, pi, sqrt, pow






##################Experimental Data Generation##################
def PinkNoise(alpha, m, tn):
	pdata = np.random.rand(m, tn)	
	funky = np.vectorize(lambda x: 1/pow(x, alpha))
	return funky(pdata)


def RunPinkNoiseTest():
	data = PinkNoise(.33, 5, 200)
	data[1, :] *= .5;
	data[1, :] += data[2,:]*.5
	data[3, :] += data[4,:]*.5
	[n,tn] = np.shape(data)
	corr_list = []
	max_corrs = np.zeros(int(floor(tn/2)-1))
	max_locs = []

	for t in range(0,int(floor(tn/2))-1):
		corr_list.append(CAtLagT(data, t))
		max_locs.append(np.argmax(np.fabs(corr_list[t])))
		max_corrs[t] = np.max(np.fabs(corr_list[t]))
		
	corr_std = np.std(max_corrs)


	p_val_w_ind = [];
	#Take all the edges with statistically significant edge chance
	p_val = map( (lambda x: 1.0-ProbZ(x, tn)), max_corrs/corr_std)
	for i in range(0, len(p_val)):
		p_val_w_ind.append((i, p_val[i]))
		
		
	sorted_p_val = sorted(p_val_w_ind, key=lambda tup: tup[1])
	q = 0.15;
	k = FDRController(sorted_p_val, q)

	for i in range(0, k+1):
		print('edge: {:d}--->{:d}, p-val={:f}  t={:d}  '.format(int(floor(max_locs[sorted_p_val[i][0]]/n)), max_locs[sorted_p_val[i][0]]%n, sorted_p_val[i][1], sorted_p_val[i][0]))	


for i in range(25):
	RunPinkNoiseTest();
	print('--------');












