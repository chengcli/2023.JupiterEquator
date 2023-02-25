#! /usr/bin/env python3.6
#
# This program compares different heat capacies of hydrogen
#
#

from pylab import *

# These data come from Table 3 of 
# Leachman, Jacobsen, Penoncello, Lemmon, 2009
# Fundamental Equations of State for Parahydrogen, Normal hydrogen and Orhtohydrogen
para_u = [4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973, 0.536078]
para_v = [499, 826.5, 970.8, 1166.2, 1341.4, 5395, 10185]

normal_u = [1.616, -0.4117, -0.792, 0.758, 1.217]
normal_v = [531, 751, 1989, 2484, 6859]

ortho_u = [2.54151, -2.3661, 1.00365, 1.22447]
ortho_v = [856, 1444, 2194, 6968]

# planck constant
h_cgs = 6.62606957e-27   

# speed of light
c_cgs = 2.99792458e+10

# boltzman constant
kBoltz_cgs = 1.3806504e-16

FACT_cgs = h_cgs*c_cgs/kBoltz_cgs

# These data come from NIST chemistry book
# 298 ~ 1000 K
shomate1 = [33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974]

# 1000 ~ 2500 K
shomate2 = [18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133]

# 2500 ~ 6000 K
shomate3 = [43.413560, -4.293079, 1.272428, -0.096876, -20.533862, -38.515158, 162.081354]

def fpara_equil(T):
  # compute partition functions
  Zpara = 0
  Zortho = 0
  JMAX = 20
  B0 = 59.322
  D0 = 4.71e-02

  DJ = zeros(JMAX)
  FJ = zeros(JMAX)

  for J in range(JMAX):
    if (J//2)*2 == J:
      DJ[J] = 2*J+1
    else:
      DJ[J] = 3*(2*J+1)
    FJ[J] = B0*J*(J + 1) - D0*(J*(J + 1))**2

  #for (size_t J = 0; J < JMAX; J += 2)
  for J in range(0, JMAX, 2):
      Zpara += DJ[J] * exp(-FJ[J]*FACT_cgs/T);

  #for (size_t J = 1; J < JMAX; J += 2)
  for J in range(1, JMAX, 2):
      Zortho += DJ[J] * exp(-FJ[J]*FACT_cgs/T);

  return Zpara/(Zpara + Zortho);

def cp_para(T):
  cp_R = 2.5
  for i in range(len(para_u)):
    cp_R += para_u[i]*(para_v[i]/T)**2*exp(para_v[i]/T)/(exp(para_v[i]/T) - 1)**2
  return cp_R

def cp_ortho(T):
  cp_R = 2.5
  for i in range(len(ortho_u)):
    cp_R += ortho_u[i]*(ortho_v[i]/T)**2*exp(ortho_v[i]/T)/(exp(ortho_v[i]/T) - 1)**2
  return cp_R

def cp_normal(T):
  cp_R = 2.5
  for i in range(len(normal_u)):
    cp_R += normal_u[i]*(normal_v[i]/T)**2*exp(normal_v[i]/T)/(exp(normal_v[i]/T) - 1)**2
  return cp_R

def cp_nist(T):
  if T < 298.:
    return NAN
  if T < 1000.:
    shomate = shomate1
  elif T < 2500.:
    shomate = shomate2
  elif T < 6000.:
    shomate = shomate3
  else:
    return NAN
  t = T/1000.
  cp = shomate[0] + shomate[1]*t + shomate[2]*t*t + shomate[3]*t*t*t + shomate[4]/(t*t)
  return cp/8.3145

def Allison_fit(T):
  return 3.298 - 0.0096*(1000/T)*(1000/T) + 0.9620*(T/1000) - 2.054*(T/1000)**2 + 1.352*(T/1000)**3

if __name__ == '__main__':
  temp_grid = linspace(100., 1000., 1001)
#for temp in temp_grid:
  figure(1, figsize = (8, 8))
  ax = axes()
  ax.plot(temp_grid, cp_para(temp_grid), 'C0', label = 'Para')
  ax.plot(temp_grid, cp_normal(temp_grid), 'C1', label = 'Normal')
  ax.plot(temp_grid, cp_ortho(temp_grid), 'C2', label = 'Ortho')
  ax.plot(temp_grid, list(map(cp_nist, temp_grid)), 'C4-', label = 'NIST', linewidth = 3)
  ax.plot(temp_grid, (1*cp_normal(temp_grid) + 0.157*2.5)/1.157, 'C3--', label = 'H$_2$:He=1:0.157')
  ax.plot(temp_grid, Allison_fit(temp_grid), 'C5--', label = 'Allison (H$_2$,He,CH4)')

  ax.legend(fontsize = 12, loc = 7)
  ax.set_xlim([100., 1000.])
  ax.set_xlabel('Temperature (K)', fontsize = 12)
  ax.set_ylabel('Cp/R', fontsize = 12)

  ax2 = ax.twinx()
  ax2.plot(temp_grid, list(map(fpara_equil, temp_grid)), 'k--')
  ax2.set_ylabel('f$_{para}$', fontsize = 12)
  print(fpara_equil(1000.))

#savefig('hydrogen_cp.png', bbox_inches = 'tight')

  show()
