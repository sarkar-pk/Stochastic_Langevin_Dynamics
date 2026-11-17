import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.stats import linregress
from pylab import *
from scipy.ndimage import gaussian_filter1d


kb=1.380649e-23
e=1.60217663e-19

custom_colors = ['r', 'b', 'g', 'c','m', 'y', 'k','orange', 'brown', 'lime', 'purple']
results = []
# Load files
file_list = sorted(glob.glob("ar_*.txt"))
lw=2
tw=2
fs=15
ls=15
tl=6
rc('axes', linewidth=lw)
plt.figure(figsize=(8, 6), dpi=300)
fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
for idx, file in enumerate(file_list):
    color = custom_colors[idx % len(custom_colors)]  
    
    # match = re.search(r'ar1_(\d+)', file)
    # if not match:
    #     continue
    # k_val = int(match.group(1))  
    
    with open(file, 'r') as f:
        metadata_line = f.readline().strip() 
    metadata_values = metadata_line.split()
    dt=float(metadata_values[1])
    V0=float(metadata_values[4])        
    k_val=float(metadata_values[7])    
    t_max=float(metadata_values[10])
    win=float(metadata_values[13])
    # Load and filter
    data = np.loadtxt(file, skiprows=1)
    t = data[:, 0]
    c = data[:, 1]/t_max
    x = data[:, 2]/t_max
    y = data[:, 3]/t_max
    mask = c != 0
    t1 = t[mask]
    c= c[mask]
    # mask1 = c > 0
    mask2 = x > 0
    mask3 = y > 0
    # t1=t
    # filtered_data = data[data[:, 1] != 0]
    # t1 = t[mask1]
    # c = c[mask1]
    t2 = t[mask2]
    x = x[mask2]
    t3 = t[mask3]
    y = y[mask3]
    
    
    t1_in=1/t1
    t2_in=1/t2
    t3_in=1/t3
    # Fit log(x) = log(A) - B*t
    if len(c) >= 2:
        log_c = np.log(c)
        slope1, intercept1, *_ = linregress(t1_in, log_c)
        B1 = -slope1
        A1 = np.exp(intercept1)
        c_fit = A1 * np.exp(-B1 * t1_in)
        U1=kb*B1*1000/e
    else: U1=5000
        

    log_x = np.log(x)
    slope2, intercept2, *_ = linregress(t2_in, log_x)
    B2 = -slope2
    A2 = np.exp(intercept2)
    x_fit = A2 * np.exp(-B2 * t2_in)

    log_y = np.log(y)
    slope3, intercept3, *_ = linregress(t3_in, log_y)
    B3 = -slope3
    A3 = np.exp(intercept3)
    y_fit = A3 * np.exp(-B3 * t3_in)


    U2=kb*B2*1000/e
    U3=kb*B3*1000/e
    
    results.append((k_val, U1, U2, U3))
    # Plot with assigned color
    # plt.plot(t_in, x, 'o', color=color,label=f"k={k_val:.2f}, U={U:.3f}")
    # plt.plot(t_in, x_fit, '--',linewidth=lw, color=color)
    # axs[0].plot(t1_in, c, 'o', color=color, label=f"k={k_val:.2f}, U1={U1:.3f}")
    # axs[0].plot(t1_in, c_fit, '--',linewidth=lw, color=color)
    axs[0].plot(k_val, U1,'o', color='k')
    axs[1].plot(k_val, U2, 'o', color='b')
    axs[2].plot(k_val, U3, 'o', color='r')
    
    
    


  
# axs[0].set_ylim(1e1,1e9)
# axs[0].set_xlabel("1/T ($K^{-1}$)", fontsize=fs)
# axs[0].set_ylabel("Event freqency (Hz)", fontsize=fs)
# axs[0].tick_params(axis='both', which='major', direction='out', 
#                 width=tw, length=tl, labelsize=ls)
# axs[0].set_yscale('log')
# axs[0].legend(ncol=3)
# axs[1].set_ylim(60,150)
# axs[0].set_xlabel("k (meV)", fontsize=fs)
axs[0].set_ylabel("U (meV)", fontsize=fs)
axs[0].tick_params(axis='both', which='major', direction='out', 
                width=tw, labelsize=ls)
# axs[1].set_xlabel("k (meV)", fontsize=fs)
axs[1].set_ylabel("U (meV)", fontsize=fs)
axs[1].tick_params(axis='both', which='major', direction='out', 
                width=tw, labelsize=ls)

axs[2].set_xlabel("k (meV)", fontsize=fs)
axs[2].set_ylabel("U (meV)", fontsize=fs)
axs[2].tick_params(axis='both', which='major', direction='out', 
                width=tw, labelsize=ls)

# plt.grid()
plt.tight_layout()

# ax = plt.gca()
# secax = axs[0].secondary_xaxis('top', functions=(lambda t1: 1/t1, 
#                                              lambda inv_t1: 1/inv_t1))
# secax.set_xlabel("Temp (K)", fontsize=fs)
# ax.tick_params(axis = 'both', which = 'major',length=tl, width=tw)
# ax.tick_params(axis = 'y', which = 'minor',length=3, width=1.5)
# secax.tick_params(axis='x', which='major', labelsize=ls, width=tw, length=tl)
plt.suptitle("dt= %0.1e s;   t_max= %0.0f $\\mu$s;   V=%0.0f meV   win=%0.1e s"
    % (dt, t_max*1e6, V0, win), y=0.98, size=15, color='b')
plt.tight_layout()
plt.savefig('ar_both_test.png', dpi=300, bbox_inches='tight', transparent=True)
#plt.show()


# Save to file
np.savetxt("U_vs_k_win50-14_test.txt", results, header="k_val    U1    U2    U3", fmt="%.2f    %.3f    %.3f    %.3f")


################################

from scipy.interpolate import make_interp_spline
from scipy.interpolate import UnivariateSpline
file2='U_vs_k_win50-14_test.txt'
data2 = np.loadtxt(file2, skiprows=1)
k = data2[:, 0]
# u1 = data2[:, 1]
u2 = data2[:, 2]
u3 = data2[:, 3]
d_u=abs(u2-u3)
# spline1 = make_interp_spline(k, u1, k=2)  # k=3 for cubic
spline2 = make_interp_spline(k, u2, k=2)
spline3 = make_interp_spline(k, u3, k=2)
k_smooth = np.linspace(k.min(), k.max(), 300)
# c_smooth = spline1(k_smooth)
x_smooth = spline2(k_smooth)
y_smooth = spline3(k_smooth) 

coeffs1= np.polyfit(k, u2, deg=13)  # Try deg=1 to deg=5
poly1 = np.poly1d(coeffs1)
u2_fit = poly1(k)

coeffs2 = np.polyfit(k, u3, deg=13)  # Try deg=1 to deg=5
poly2 = np.poly1d(coeffs2)
u3_fit = poly2(k)

coeffs3 = np.polyfit(k, d_u, deg=10)  # Try deg=1 to deg=5
poly3 = np.poly1d(coeffs3)
du_fit = poly3(k)




lw=2
tw=2
fs=15
ls=15
tl=6
rc('axes', linewidth=lw)
# plt.figure(figsize=(8, 6), dpi=300)
fig, axs2 = plt.subplots(3, 1, figsize=(8,7), sharex=True)

axs2[0].plot(k, u2, 'o', color='b', label='$U_\\theta$', alpha=0.4)
# axs2[0].plot(k_smooth,x_smooth, color='black')
axs2[0].plot(k,u2_fit, color='black')
axs2[1].plot(k, u3, 'o', color='r', label='$U_\\phi$', alpha=0.4)
# axs2[1].plot(k_smooth,y_smooth, color='black')
axs2[1].plot(k,u3_fit, color='black')
axs2[2].plot(k, d_u, 'o', color='cyan', alpha=0.4, label='$\u0394$U=|$U_{\\theta}$ $\u2212$ $U_{\\phi}$|')
axs2[2].plot(k,du_fit, color='k')
# axs2.set_ylim(-30,500)
# axs2.set_xlim(0, 700)
axs2[2].set_xlabel("k (meV)", fontsize=fs)
axs2[0].set_ylabel("$U_{\\theta}$ (meV)", fontsize=fs)
axs2[1].set_ylabel("$U_{\\phi}$ (meV)", fontsize=fs)
axs2[2].set_ylabel("\u0394U (meV)", fontsize=fs)
# axs2[2].set_ylim(-5,50)

axs2[0].tick_params(axis='both', which='major', direction='out', 
                width=tw, length=tl, labelsize=ls)
axs2[1].tick_params(axis='both', which='major', direction='out', 
                width=tw, length=tl, labelsize=ls)
axs2[2].tick_params(axis='both', which='major', direction='out', 
                width=tw, length=tl, labelsize=ls)
axs2[0].legend(fancybox=False, framealpha=1, fontsize=12, loc="lower right", bbox_to_anchor=(1, 0))
axs2[1].legend(fancybox=False, framealpha=1, fontsize=12, loc="lower right", bbox_to_anchor=(1, 0))
axs2[2].legend(fancybox=False, framealpha=1, fontsize=12, loc="upper right", bbox_to_anchor=(1, 1))
plt.suptitle("dt= %0.1e s;   t_max= %0.0f $\\mu$s;   V=%0.0f meV   win=%0.1e s"
    % (dt, t_max*1e6, V0, win), y=0.95, size=15, color='b')
plt.savefig('U_vs_k_win_50xdt-14_test.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()
