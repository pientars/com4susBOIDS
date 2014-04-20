import numpy as np
import csv

bird_locs = np.genfromtxt('tree_swallow_locs.csv', delimiter=',');
bird_locs = bird_locs[1::,1::];

min_lat = bird_locs[:,0].min(axis=0); max_lat = bird_locs[:,0].max(axis=0);
min_long = bird_locs[:,1].min(axis=0); max_long = bird_locs[:,1].max(axis=0);
lat_range = max_lat - min_lat; long_range = max_long - min_long;

# number of square segments along one dimension. take as param
num_div = 50; assert(num_div > 2);
num_times = 52;

# accumulate time-series data for each bounding box in 'graph'
graph = [[[0]*num_times for i in xrange(num_div + 1)] for j in xrange(num_div + 1)];
num_ts = [[0 for i in xrange(num_div + 1)] for j in xrange(num_div + 1)];

graph_lats = np.linspace(min_lat, max_lat, num_div + 1);
graph_longs = np.linspace(min_long, max_long, num_div + 1);
lat_dim = graph_lats[1] - graph_lats[0];
long_dim = graph_longs[1] - graph_longs[0];

# read in the file of bird means
mean_file = open('tree_swallow_mean.csv', 'r');
reader = csv.reader(mean_file, delimiter = ',');

first_line = True
for line in reader:
	if first_line:    #skip first line
		first_line = False
		continue
	coord_id = int(line[0]) - 1;
	ts = [float(val) for val in line[1::]];
	h = int(get_horz_node(bird_locs[coord_id,1],min_long,long_dim));
	v = int(get_vert_node(bird_locs[coord_id,0],min_lat,lat_dim));
	graph[v][h] = map(lambda x: x[0] + x[1], zip(graph[v][h], ts))
	#graph[v][h] = [x+y for x,y in zip(graph[v][h], ts)];
	num_ts[v][h] += 1;

for i,j in np.ndindex((num_div, num_div)):
	if (num_ts[i][j] >= 1):
		graph[i][j] = [x / num_ts[i][j] for x in graph[i][j]];

with open('bird_ts.csv', 'wb') as csvfile:
	graphwriter = csv.writer(csvfile, delimiter=',');
	for i,j in np.ndindex((num_div, num_div)):
		# place the node coordinates in the center of the bounding box
		node_lat = graph_lats[i] + (lat_dim / 2);
		node_long = graph_longs[j] + (long_dim / 2);
		data_to_write = [num_ts[i][j], i, j, node_lat, node_long, graph[i][j]];
		graphwriter.writerow(data_to_write);

def get_horz_node(longitude, min_long, long_dim):
	return int((longitude - min_long) // long_dim);

def get_vert_node(latitude, min_lat, lat_dim):
	return int((latitude - min_lat) // lat_dim);