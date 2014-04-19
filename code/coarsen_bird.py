import numpy as np
import csv
#from numpy import genfromtxt

bird_locs = np.genfromtxt('tree_swallow_locs.csv', delimiter=',');
bird_locs = bird_locs[1::,1::];

min_lat = bird_locs[:,0].min(axis=0); max_lat = bird_locs[:,0].max(axis=0);
min_long = bird_locs[:,1].min(axis=0); max_long = bird_locs[:,1].max(axis=0);
lat_range = max_lat - min_lat; long_range = max_long - min_long;

# number of square segments along one dimension. take as param
num_div = 50; assert(num_div > 2);

# regularize the space as (ns * ns) equally sized rectangles
graph_dims = [num_div, num_div];

# a set of coordinates for each bounding box
graph = [[set() for i in xrange(numdiv)] for j in xrange(numdiv)];

#   don't think we need these ...
# graph_lats = np.linspace(min_lat, max_lat, num_div + 1);
# graph_longs = np.linspace(min_long, max_long, num_div + 1);

lat_dim = graph_lats[1] - graph_lats[0];
long_dim = graph_longs[1] - graph_longs[0];

mean_file = open('test_swallow_mean.csv', 'r');
reader = csv.reader(mean_file, delimiter = ',');

for line in reader:
	coord_id = int(line[0]) - 1;
	ts = [float(val) for val in line[1::]];
	h = get_horz_node(bird_locs[coord_id,1],min_long,long_dim);
	v = get_vert_node(bird_locs[coord_id,0],min_lat,lat_dim);
#todo: add ts to graph[v][h]

def get_horz_node(longitude, min_long, long_dim):
	return (longitude - min_long) // long_dim;

def get_vert_node(latitude, min_lat, lat_dim):
	return (latitude - min_lat) // lat_dim;