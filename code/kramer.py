import numpy as np
# import matplotlib.pyplot as plt
from math import exp, log, floor

def CrossCorrijAtLagT(data, t, i, j):
	[n,tn] = np.shape(data)
	sigi = np.std(data[i,:])
	mui = np.mean(data[i,:])
	sigj = np.std(data[j,:])
	muj = np.mean(data[j,:])

	accu = 0;

	for r in range(0, n-t):
		accu += (data[i,r] - mui) * (data[j,r+t] - muj);
	accu /= 1/(sigi*sigj*(n-2*t));
	return accu;

def FisherTransform(Cijt):
	return .5 * log((1+Cijt)/(1-Cijt));


def CAtLagT(data, t):
	[n,tn] = np.shape(data)
	Corrs = np.zeros((5,5))

	# print(data);
	for i in range(0,n):
		for j in range(0,n):
			if (i == j) : continue
			Corrs[i,j] = FisherTransform(CrossCorrijAtLagT(data, i, j, 2))
	return Corrs;	



data = np.random.rand(5,50)
[n,tn] = np.shape(data)
corr_list = []
tn = 4;
for t in range(0,int(floor(tn/2))):
	corr_list.append(CAtLagT(data, t))

print(corr_list);



#create the cross correlation Cij[t]

#take the Fishes transform CijF[t]

#Calculate the variance over CijF[t]

#
