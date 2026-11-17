import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import math
from pylab import *
import shutil
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

file='ar_000.txt'

with open(file, 'r') as f:
    metadata_line1 = f.readline().strip()  
#    metadata_line2 = f.readline().strip()  

metadata_values = metadata_line1.split()  

# dt=np.zeros((1))
# T=np.zeros((1))
# V0=np.zeros((1))

dt=float(metadata_values[1])
V0=float(metadata_values[4])        
k=float(metadata_values[7])    
t_max=float(metadata_values[10])

data = np.loadtxt(file, skiprows=1)
kb=1.380649e-23; e=1.60217663e-19
temp = data[:, 0]  
x = data[:, 2]/t_max
y = data[:, 3]/t_max 
ln_x = np.log(x)
ln_y = np.log(y)
temp_in=1/temp

slope1, intercept1, r_value1, p_value1, std_err1 = linregress(temp_in, ln_x)
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(temp_in, ln_y)
fitted_x = np.exp(slope1 * temp_in + intercept1)
fitted_y = np.exp(slope2 * temp_in + intercept2)
U0=-kb*slope1*1000/e
U1=-kb*slope2*1000/e

lw=1.5
fs=12
ls=12
rc('axes', linewidth=lw)
# plt.figure(figsize=(6, 5), dpi=300)
fig, axs = plt.subplots(1, 1, figsize=(5, 4), sharex=True, dpi=300)
# plt.plot(temp_in, x, 'o', color='r', alpha=1, label=f'$U_0$={U0:0.2f} meV')
plt.plot(temp_in, y, 'o', color='b', alpha=1, label=f'$U_0$= {V0:0.2f} meV')
# plt.plot(temp_in, fitted_x, 'r-')
plt.plot(temp_in, fitted_y, 'b-')
plt.text(0.010225, 4e8, f'Calculated $U$ from slope= {U1:0.2f} meV', bbox=dict(edgecolor='none',facecolor='white', alpha=1.0),
        fontsize=10)
# plt.plot(U_kbT, y, 'o', color='b')
# plt.plot(U_kbT, fitted_y, 'b-')
plt.xlabel('1/T ($K^{-1}$)', fontsize=fs)
yscale('log')
plt.ylabel('Event frequency (Hz)', fontsize=fs)
# plt.ylim(1e8,1e10)
plt.tick_params(axis='both', which='major', direction='out', width=lw, labelsize=ls)
ax = plt.gca()
secax = ax.secondary_xaxis('top', functions=(lambda temp: 1/temp, 
                                             lambda temp_in: 1/temp_in))
secax.set_xlabel("Temp (K)", fontsize=fs)
ax.tick_params(axis = 'both', which = 'major',length=4, width=lw)
ax.tick_params(axis = 'y', which = 'minor',length=3, width=1.5)
secax.tick_params(axis='x', which='major', labelsize=ls, width=lw, length=4)
# plt.suptitle("dt= %0.1e s;  $t_{max}$= %0.0f $\\mu$s;   $V_0$=%0.0f meV;   k=%0.0f meV/$rad^2$"
#     % (dt, t_max*1e6, V0, k), y=1.05, size=10, color='b')
plt.grid(True, which="both", linestyle='--')
plt.legend(fancybox=False, framealpha=1, fontsize=10, loc="upper right",edgecolor='none', bbox_to_anchor=(1, 1))
# plt.figtext(0.724,0.75, f'$U_0$:{U0:0.2f}', fontsize=8, bbox=dict(facecolor='blue', alpha=0.4))
# plt.savefig(f"ar_dt{dt}_{V0:0.0f}mV_int{k}_1.png",bbox_inches='tight',transparent=True, dpi=300)
# new_file = f"ar_dt{dt}_{V0:0.0f}mV_int{k}_1.txt"
# shutil.copy(file, new_file)
plt.show()
intercept1
