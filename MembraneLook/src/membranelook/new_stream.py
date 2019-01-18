# coding: utf-8
import matplotlib, matplotlib.pyplot
import numpy as np
import MDAnalysis, MDAnalysis.visualization.streamlines

u1, v1, average_displacement, standard_deviation_of_displacement = MDAnalysis.visualization.streamlines.generate_streamlines('min_bilayer_NP.gro', 'complete.xtc', grid_spacing=5, MDA_selection='name PO4', start_frame=100, end_frame=110, xmin=0.03000049591, xmax= 150.96008301, ymin= 0.0799999237, ymax= 150.34008789, maximum_delta_magnitude=60.0, num_cores=10)
x = np.linspace(0, 150, 30)
y = np.linspace(0, 150, 30)
speed = np.sqrt(u1*u1 + v1*v1)
fig = matplotlib.pyplot.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.set_xlabel('x ($\AA$)')
ax.set_ylabel('y ($\AA$)')
strm = ax.streamplot(x, y, u1, v1, density=(50,50), color=speed, linewidth=1*speed/speed.max())
fig.colorbar(strm.lines)
fig.savefig('testing_streamline.png',dpi=300)
