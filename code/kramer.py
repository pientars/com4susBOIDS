import numpy as np
#import matplotlib.pyplot as plt
from math import exp, log, floor, fabs, pi, sqrt, pow

def CrossCorrijAtLagT(data, t, i, j):
	[m,tn] = np.shape(data)
	sigi = np.std(data[i,0:(tn-t-1)])
	mui = np.mean(data[i,0:(tn-t-1)])
	sigj = np.std(data[j,t:(tn-1)])
	muj = np.mean(data[j,t:(tn-1)])


	accu = 0;
	for r in range(0, tn-t-1):
		accu += (data[i,r] - mui) * (data[j,r+t] - muj);
	return accu/(sigi*sigj*(tn-t-1));

def FisherTransform(Cijt):
	return .5 * log((1+Cijt)/(1-Cijt));


def CAtLagT(data, t):
	[n,tn] = np.shape(data)
	Corrs = np.zeros((n,n))
	# print(data);
	for i in range(0,n):
		for j in range(0,n):
			if (i == j) : continue
			# print(CrossCorrijAtLagT(data, t, i, j))
			Corrs[i,j] = FisherTransform(CrossCorrijAtLagT(data, t, i, j))
	return Corrs;	

#P{Zij <= z}
def ProbZ(Zij, tn):
	an = sqrt(2*log(tn))
	bn = an - (log(log(tn)) + np.log(4*pi))/(2*an)
	return exp( -2 * exp(-an*(Zij-bn)) )






##################Experimental Data Generation##################
def PinkNoise(alpha, m, tn):
	pdata = np.random.rand(m, tn)	
	funky = np.vectorize(lambda x: 1/pow(x, alpha))
	return funky(pdata)


#print(PinkNoise(.333, 5, 5))
#Dependent edges from 1<->2   and   3<->4
#data = np.random.rand(5,50)
data = PinkNoise(.33, 5, 200)
data[1, :] += data[2,:]*.5
data[3, :] += data[4,:]*.5
#print(data)
# data[4, 10:50] += np.ones((40))*.1

[n,tn] = np.shape(data)
corr_list = []
max_corrs = np.zeros(int(floor(tn/2)-1))
max_locs = []

for t in range(0,int(floor(tn/2))-1):
	corr_list.append(CAtLagT(data, t))
	# print(corr_list[t])
	# print("\n")
	max_locs.append(np.argmax(np.fabs(corr_list[t])))
	max_corrs[t] = np.max(np.fabs(corr_list[t]))
	
corr_std = np.std(max_corrs)

# print(max_corrs/corr_std)
# print(max_locs)


#Take all the edges with statistically significant edge chance
p_val = map( (lambda x: 1.0-ProbZ(x, tn)), max_corrs/corr_std)
for i in range(0, len(p_val)):
	if p_val[i] < 0.01:
		print('p-val={:f}  t={:d}  i={:d}  j={:d}'.format(p_val[i], i,  int(floor(max_locs[i]/n)), max_locs[i]%n))

#print(p_val)
