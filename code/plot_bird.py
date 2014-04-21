import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import csv

data = np.genfromtxt('bird_coords.csv', delimiter=',');
coords = data[:,3:5].transpose();

# Lambert Conformal map of USA lower 48 states
m = Basemap(llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64,
urcrnrlat=49, projection='lcc', lat_1=33, lat_2=45,
lon_0=-95, resolution='h', area_thresh=10000)

m.drawcoastlines()
m.drawcountries(linewidth=2)
m.drawstates()

# fill the background (the oceans)
m.drawmapboundary(fill_color='aqua')
# fill the continental area
# we color the lakes like the oceans
m.fillcontinents(color='coral',lake_color='aqua')

# draw parallels and meridians
m.drawparallels(np.arange(25,65,20),labels=[1,0,0,0])
m.drawmeridians(np.arange(-120,-40,20),labels=[0,0,0,1])

cities = []*len(coords);
lat = coords[0,:];
lon = coords[1,:];

x, y = m(lon, lat)
plt.plot(x, y, 'ro', markersize=2)

# for each city,
#for city, xc, yc in zip(cities, x, y):
	# draw the city name in a yellow (shaded) box
	#plt.text(xc+250000, yc-150000, "bird", bbox=dict(facecolor='yellow', alpha=0.5))

plt.show()