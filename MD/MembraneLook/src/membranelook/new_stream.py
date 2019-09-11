# coding: utf-8
import matplotlib, matplotlib.pyplot
import numpy as np
import MDAnalysis, MDAnalysis.visualization.streamlines
import matplotlib.pyplot as plt

u1, v1, average_displacement, standard_deviation_of_displacement = MDAnalysis.visualization.streamlines.generate_streamlines('leaflet0.gro', 'l0_800_810.xtc', grid_spacing=20, MDA_selection='name PO4', start_frame=9, end_frame=10, xmin=0.00, xmax=1160.0, ymin= 0.00, ymax= 1160.0, maximum_delta_magnitude=60.0, num_cores=10)
x = np.linspace(0, 1160, 20)
y = np.linspace(0, 1160, 20)
x_c = np.linspace(0, 1160, 57)
y_c = np.linspace(0, 1160, 57)

nb_squares_x = (1160-0)/20 - 1
nb_squares_y = (1160-0)/20 - 1

xmin_in = 0 
xmax_in = 1160 

ymin_in = 0 
ymax_in = 1160 

zmin_in = 30 
zmax_in = 90 

U2 = MDAnalysis.Universe('leaflet0.gro','l0_800_810.xtc')
PO4 = U2.select_atoms('all')
PO4_coords = PO4.positions

population = np.zeros((int(nb_squares_x),int(nb_squares_y)))
z_avg_tmp = np.zeros((int(nb_squares_x),int(nb_squares_y)))

population +=  np.histogram2d(PO4_coords[:,0], PO4_coords[:,1], bins=(nb_squares_x,nb_squares_y), range=[[xmin_in, xmax_in], [ymin_in, ymax_in]])[0] 
z_avg_tmp +=  np.histogram2d(PO4_coords[:,0], PO4_coords[:,1],bins=(nb_squares_x,nb_squares_y), range=[[xmin_in, xmax_in], [ymin_in, ymax_in]], weights= PO4_coords[:,2])[0]
pop = population < 1  
z_avg_tmp = z_avg_tmp/(population + pop)
levels = matplotlib.ticker.MaxNLocator(nbins=100).tick_values(zmin_in, zmax_in)
cbar = matplotlib.pyplot.colorbar()
cbar.set_label('z_value [A]',size=15)
speed = np.sqrt(u1*u1 + v1*v1)
print(np.shape(speed))
fig = matplotlib.pyplot.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.set_xlabel('x ($\AA$)')
ax.set_ylabel('y ($\AA$)')
plt.contourf(x_c, y_c, np.transpose(z_avg_tmp), cmap=matplotlib.pyplot.cm.jet, levels=levels, vmin = zmin_in , vmax = zmax_in)	
plt.streamplot(x, y, u1, v1, density=(10,10), color=speed, linewidth=3.0*speed/speed.max(), cmap=matplotlib.pyplot.cm.jet)
plt.show()
fig.savefig('testing_streamline.png',dpi=300)
fig.colorbar(strm.lines)
