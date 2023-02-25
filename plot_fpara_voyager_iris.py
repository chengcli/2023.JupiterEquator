#! /usr/bin/env python3.6
from pylab import *
from hydrogen_cp import *

data = genfromtxt('Conrath_Gierasch_1984_fpara.csv', delimiter = ',')
lat = data[:,0]
fpara = data[:,1]
print(lat, fpara)

data_fitted = genfromtxt('Conrath_Gierasch_1984_fpara_fitted.csv', delimiter = ',')
lat2 = data_fitted[:,0]
fpara2 = data_fitted[:,1]

figure(1, figsize = (10,8))
ax = axes()
ax.scatter(lat, fpara, s = 5, c = 'k')
ax.errorbar([50], [0.325], yerr = 0.01, marker = 'o', capsize = 4)
ax.errorbar([8.2], [0.29], yerr = 0.005, marker = 'o', capsize = 4)
ax.plot(lat2, fpara2, 'C1', linewidth = 3)
ax.plot([-60., 60.], [fpara_equil(120.), fpara_equil(120.)], 'C2', linewidth = 2)
ax.set_xlabel('Latitude (degree)', fontsize = 18)
ax.set_ylabel('f$_{para}$', fontsize = 18)

savefig('figs/fpara_voyager_iris.png', bbox_inches = 'tight')
