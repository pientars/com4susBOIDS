import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
import numpy as np
import csv
import sys

def plot_ts(node1, node2):
	ts = np.genfromtxt('bird_timeseries.csv', delimiter=',');
	x = np.linspace(0,52,52);
	y = ts[node1,:];
	plt.plot(x, y);
	plt.hold(True);
	y = ts[node2,:];
	plt.plot(x, y);
	plt.show();
	plt.hold(False);

def my_color(val):
	return plt.cm.autumn(val)

def get_edges_from_period(period):
	data = np.genfromtxt('map_vals.csv', delimiter=',');
	edges = data[data[:, 0] == period];
	edges = edges[edges[:, 1] <= 2];
	edges = edges[:,2:];
	return edges;



def plot_birds_okay(time):
	data = np.genfromtxt('bird_coords.csv', delimiter=',');
	coords = data[:,3:5].transpose();

	# # Lambert Conformal map of USA lower 48 states
	m = Basemap(llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64,
	urcrnrlat=49, projection='lcc', lat_1=33, lat_2=45,
	lon_0=-95, resolution='c', area_thresh=10000)

	m.drawcoastlines()
	m.drawcountries(linewidth=1.5)
	m.drawstates(color='.5')

	m.drawmapboundary(fill_color='.9')
	m.fillcontinents(color='.6',lake_color='.9')
	m.drawparallels(np.arange(25,65,20),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-120,-40,20),labels=[0,0,0,1])

	#cities = []*len(coords);
	lat = coords[0,:];
	lon = coords[1,:];

	# edges will be a set of edge tuples, weighted by the third column
	# edges are denoted by node id (linear from 0..2500)
	#edges = np.array([[20, 25, 0.6], [1300,2000,0.4], [1200,2134,0.5], [300,301,1.5], [200,201,1.2], [600,700,1.8], [25,26,0.5], [27,29,0.5], [1400,1405,1.2], [407,409,1.5]]);
	print(sys.argv[1])
	time = int(sys.argv[1]);
	window = time*5;
	edges = get_edges_from_period(window);

	# get coordinates of each edge
	edge_lats_src = np.array([lat[i] for i in edges[:,0]],ndmin=2);
	edge_lons_src = np.array([lon[i] for i in edges[:,0]],ndmin=2);
	edge_lats_trg = np.array([lat[i] for i in edges[:,1]],ndmin=2);
	edge_lons_trg = np.array([lon[i] for i in edges[:,1]],ndmin=2);

	weights = edges[:,2].T;
	colors = weights;
	alphas = edges[:,3].T;
	max_weight = weights.max();
	colors = 1 - (colors/max_weight); #normalize to 0..1 and reverse

	sx, sy = m(edge_lons_src, edge_lats_src);
	tx, ty = m(edge_lons_trg, edge_lats_trg);

	x, y = m(lon, lat);
	#plt.plot(x, y, 'ro', markersize=2)

	#use line collections if this slows down
	for i in range(len(weights)-1):
		if (weights[i] > 0.1):
			plt.plot([sx[0][i], tx[0][i]], [sy[0][i], ty[0][i]], linewidth=.5, alpha=alphas[i], color=my_color(colors[i]));

	plt.savefig('bird' + str(window) + '.png')
	# #plt.show()

plot_birds_okay(0);