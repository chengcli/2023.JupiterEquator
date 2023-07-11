#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d
from contribution_function import plot_contribution_function
import string

plt.rcParams['font.family'] = 'Serif'
plt.tick_params(axis='both', labelsize=12)

# frequencies
freq = [0.6, 1.25, 2.6, 5.2, 10., 22.]

# adiabatic
data = Dataset('outputs/juno_mwr-ad-main.nc', 'r')
pres = data['press'][0,:,0,0]/1.E5
tb = data['radiance'][0,:,0,0]
tb45 = tb[3::4]

tb0_ad = tb[::4]
ld0_ad = (tb0_ad - tb45)/tb0_ad*100.
temp_ad = data['temp'][0,:,0,0]

print(tb0_ad)

# super-adiabatic
data = Dataset('outputs/juno_mwr-sup-main.nc', 'r')
tb = data['radiance'][0,:,0,0]
tb45 = tb[3::4]

tb0_sup = tb[::4]
ld0_sup = (tb0_sup - tb45)/tb0_sup*100.
temp_sup = data['temp'][0,:,0,0]

print(tb0_sup)

# sub-adiabatic
data = Dataset('outputs/juno_mwr-sub-main.nc', 'r')
tb = data['radiance'][0,:,0,0]
tb45 = tb[3::4]

tb0_sub = tb[::4]
ld0_sub = (tb0_sub - tb45)/tb0_sub*100.
temp_sub = data['temp'][0,:,0,0]

print(tb0_sub)

fig, axs = subplots(1, 3, figsize = (12, 10),
  gridspec_kw = {'width_ratios': [4,1,1]})
subplots_adjust(wspace = 0.08)

ax = axs[0]
ax2 = ax.twiny()
plot_contribution_function(ax2, "juno_mwr-tau", 0.)

ax.plot(temp_ad, pres, 'k-', label='Adiabatic')
ax.plot(temp_sup, pres, 'C0-', label='Super-adiabatic')
ax.plot(temp_sub, pres, 'C1-', label='Sub-adiabatic')

ax.set_yscale('log')
ax.set_xlabel('Temperature (K)', fontsize = 14)
ax.set_ylabel('Pressure (K)', fontsize = 14)
ax.set_ylim([20., 1.])
ax.set_xlim([170., 420.])
ax.legend(loc=0, fontsize = 14)

ax = axs[1]
ax.plot([0., 0.], [0.6, 22.], 'k--')
print(tb0_sup - tb0_sub)
ax.errorbar(tb0_sup - tb0_ad, freq, xerr = 0.02*tb0_sup, 
            label='Super-adiabatic', capsize = 5, ecolor = 'C0')
ax.errorbar(tb0_sub - tb0_ad, freq, xerr = 0.02*tb0_sup, 
            label='Sub-adiabatic', capsize = 5, ecolor = 'C1')
ax.set_ylim([0.59, 22.1])
ax.set_xlim([-15., 15.])
ax.set_yscale('log')
ax.set_xlabel('Tb anomaly (K)', fontsize = 14)
ax.set_yticklabels([])

ax = axs[2]
ax.plot([0., 0.], [0.6, 22.], 'k--')
print(ld0_sup - ld0_sub)
ax.errorbar(ld0_sup - ld0_ad, freq, xerr = 0.1,
            label='Super-adiabatic', capsize = 5, ecolor = 'C0')
ax.errorbar(ld0_sub - ld0_ad, freq, xerr = 0.1, 
            label='Sub-adiabatic', capsize = 5, ecolor = 'C1')
ax.set_ylim([0.59, 22.1])
ax.set_xlim([-0.8, 0.8])
ax.yaxis.set_ticks([])
ax.set_xlabel('R45 anomaly (%)', fontsize = 14)
ax.set_yscale('log')
ax.set_yticklabels([])

labels = iter(string.ascii_uppercase)  # Generates the labels A, B, C, ...
for ax in axs:
    ax.annotate(f"({next(labels)})", xy=(0.05, 0.05), xycoords="axes fraction", fontsize=12, fontweight="bold")

#show()
savefig('three_cases_dTdz.png', bbox_inches='tight')
