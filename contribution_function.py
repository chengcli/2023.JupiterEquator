#! /usr/bin/env python3.6
from pylab import *
from netCDF4 import Dataset
#from snapy.harp.utils import get_ray_out

def contribution_function(case, band = '', ang = None):
  datafile = 'outputs/%s-main.nc' % case
  data = Dataset(datafile, 'r')
  pres = data['press'][0,:,0,0]/1.E5 # pa -> bar
  dlnp = log(pres[0]) - log(pres[1])
  dlnp = log(pres[1]) - log(pres[2])
  x1f = data['x1f'][:]

  inpfile = '%s.inp' % case
  #amu, aphi = get_ray_out(inpfile, band)
  amu, aphi = 0., 0.
  btau = data['%stau' % band][0,:,0,0]

  tau = cumsum(btau[::-1])[::-1]
  if ang != None:
    mu = cos(ang/180.*pi)
  else:
    mu = cos(amu/180.*pi)
  wfunc = 1./mu*exp(-tau/mu)*btau/dlnp
  wfunc /= wfunc.max()
  return pres, wfunc

def plot_contribution_function(ax, case, ang = None):
  pres, wfunc1 = contribution_function(case, 'CH1', ang)
  pres, wfunc2 = contribution_function(case, 'CH2', ang)
  pres, wfunc3 = contribution_function(case, 'CH3', ang)
  pres, wfunc4 = contribution_function(case, 'CH4', ang)
  pres, wfunc5 = contribution_function(case, 'CH5', ang)
  pres, wfunc6 = contribution_function(case, 'CH6', ang)

  ax.plot(wfunc1, pres, 'k--', color = '0.7')
  ax.plot(wfunc2, pres, 'k--', color = '0.7')
  ax.plot(wfunc3, pres, 'k--', color = '0.7')
  ax.plot(wfunc4, pres, 'k--', color = '0.7')
  ax.plot(wfunc5, pres, 'k--', color = '0.7')
  ax.plot(wfunc6, pres, 'k--', color = '0.7')
  #ax.legend(fancybox=True, framealpha=0.5)

  ax.set_ylim([20., 0.5])
  ax.set_yscale('log')
  ax.set_xlabel('contribution function', fontsize = 12)
  ax.set_ylabel('Pressure (bar)', fontsize = 12)

if __name__ == '__main__':
  case = 'saturn_vla_inversion--9.8'
  #case = '../build_disk/saturn_vla-disk-1.6'
  figure(1, figsize = (8, 8))
  ax = axes()
  plot_contribution_function(ax, case, 48.2)
  show()
  #savefig('../figs/contribution_function.png', bbox_inches = 'tight')
