import numpy as np
import time
#import matplotlib.pyplot as plt
from math import exp, log, floor, fabs, pi, sqrt, pow


def LagCorrelation(v1, v2, lag):
	# cc = np.correlate(v1,np.hstack([np.zeros(lag), v2[lag:]]));
	n = v1.size
	v1 -= v1[:(n-lag)].mean()
	v2 -= v2[lag:].mean()
	v1 /= v1[:(n-lag)].std()
	v2 /= v2[lag:].std()
	return np.correlate(v1[:(n-lag)], v2[lag:])/(n-lag)

def LagCorrelationRange(v1, v2, start, stop):
	results = np.zeros(stop-start+1);
	for lag in range(start, stop):
		results[lag] = LagCorrelation(v1, v2, lag);
	print results

def FisherTransform(C):
	one_vec = np.ones(C.shape);
	return 0.5 * np.log( (one_vec + C)/(one_vec-C) )
	
def ProbZ(Zij, tn):
	an = sqrt(2*log(tn))
	bn = an - (log(log(tn)) + np.log(4*pi))/(2*an)
	return exp( -2 * exp(-an*(Zij-bn)) )

# Returns the first k sorted values that fail null hypothesis
def FDRController(p_vals, q):
	m = len(p_vals)
	k = 0;
	for i in range(m):
		# print('{:f}    {:f}'.format(p_vals[i][1], q*(1+i)/m))
		if (p_vals[i][1] > q*(i+1)/m):
			break;
		k = i;
	return k;

def InterpolateEdges(vec):
	return 0;


def Get4Neighbors(map_m, map_n, i):
	neigh = [];
	for j in [i-map_n, i+1, i+map_n, i-1]:
		if(j < 0 or j >= map_n*map_m): continue
		# left edge
		if i%map_n == 0:
			if j == i-1: continue
		# right edge 
		if i%map_n == map_n-1:
			if j == i+1: continue
		neigh.append(j);
	return neigh

def Get8Neighbors(map_m, map_n, i):
	neigh = [];
	for j in [i-map_n-1, i-map_n, i-map_n+1, i-1, i+1, i+map_n-1, i+map_n, i+map_n+1 ]:
		if(j < 0 or j >= map_n*map_m): continue
		# left edge
		if i%map_n == 0:
			if j == i-1 or j == i-map_n-1 or j == i+map_n-1: continue
		# right edge
		if i%map_n == map_n-1:
			if j == i+1 or j == i-map_n+1 or j == i+map_n+1: continue
		neigh.append(j);
	return neigh




def GetValidEdges(data):
	t_range = 1;
	[n,tn] = np.shape(data)
	corr_list = []

	map_m = 50
	map_n = 50

	Corrs = np.zeros((n,n,t_range));
	# This will change for the regions we are using 
	for t in range(0, t_range):
		for i in range(0, map_m*map_n):
			neigh = Get4Neighbors(map_m, map_n, i)
			for j in neigh:
				Corrs[i,j,t] = FisherTransform(LagCorrelation(data[i,:], data[j,:], t));

	# for t in range(0, t_range):
	# 	for i in range(0,n):
	# 		for j in range(0,n):
	# 			if i == j : continue
	# 			Corrs[i,j,t] = FisherTransform(LagCorrelation(data[i,:], data[j,:], t));

	print Corrs[:,:, 0];

	max_corrs = np.zeros(t_range)
	max_locs = []				
	# get the max locations and values
	for t in range(0, t_range):
		max_locs.append(np.argmax(np.fabs(Corrs[:,:,t])))
		max_corrs[t] = np.max(np.fabs(Corrs[:,:,t])) 
	
	corr_std = np.std(max_corrs)
	#Take all the edges with statistically significant edge chance
	p_val = map( (lambda x: 1.0-ProbZ(x, tn)), max_corrs/corr_std)
	p_val_w_ind = []
	for i in range(0, len(p_val)):
		p_val_w_ind.append((i, p_val[i]))
	sorted_p_val = sorted(p_val_w_ind, key=lambda tup: tup[1])
	q = 0.25;
	k = FDRController(sorted_p_val, q)
	edge_list = []
	for i in range(0, k+1):
		edge_list.append( (int(floor(max_locs[sorted_p_val[i][0]]/n)), max_locs[sorted_p_val[i][0]]%n, sorted_p_val[i][0]) )
		# print('edge: {:d}--->{:d}, p-val={:f}  t={:d}  '.format(int(floor(max_locs[sorted_p_val[i][0]]/n)), max_locs[sorted_p_val[i][0]]%n, sorted_p_val[i][1], sorted_p_val[i][0]))	
	
	# for edge in edge_list:
	# 	print edge
	# print max_locs
	# print '\n'
	# print max_corrs
	return 0



##################Experimental Data Generation##################
def PinkNoise(alpha, m, tn):
	pdata = np.random.rand(m, tn)	
	funky = np.vectorize(lambda x: 1/pow(x, alpha))
	return funky(pdata)



def RunPinkNoiseTest():
	data = PinkNoise(.33, 2500, 20)
	data[1, :] *= .2;
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
	# for i in range(0, k+1):
		# print('edge: {:d}--->{:d}, p-val={:f}  t={:d}  '.format(int(floor(max_locs[sorted_p_val[i][0]]/n)), max_locs[sorted_p_val[i][0]]%n, sorted_p_val[i][1], sorted_p_val[i][0]))	



# for t  in range(20, 2000, 100):
# 	data = PinkNoise(.33, t, 20)
# 	data[1, :] *= .2;
# 	data[1, :] += data[2,:]*.5
# 	data[3, :] += data[4,:]*.5
		
# 	# v1 = np.array([1., 2., 3., 4., 5., 6., 7.])
# 	# v2 = np.array([1., 1., 1., 4., 5., 6., 7.])
#bird_locs = np.genfromtxt('bird_timeseries.csv', delimiter=',');
coords = np.genfromtxt('bird_coords.csv', delimiter=',');
ts = np.genfromtxt('bird_timeseries.csv', delimiter=',');

print 
a = time.clock()
GetValidEdges(ts)
print time.clock() - a






