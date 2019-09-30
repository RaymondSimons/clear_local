import matplotlib.pyplot as plt
import numpy as np
from numpy import *

z_05_10 = [-27.40, 5.02, -0.22, 8.5, 11.2 , 'blue']
z_10_15 = [-26.03, 4.62, -0.19, 9.1, 11.3,  'lightgreen']
z_15_20 = [-24.04, 4.17, -0.16, 9.2, 11.5,  'orange']
z_20_25 = [-19.99, 3.44, -0.13, 9.5, 11.5, 'red']



z_05_10 = [-27.4020,5.02150,-0.219550, 8.5, 11.2,'blue']
z_10_15 = [-26.0395,4.61789,-0.190385, 9.1, 11.3, 'lightgreen']
z_15_20 = [-24.04,4.17,-0.163760, 9.2, 11.5, 'orange']
z_20_25 = [-19.9898,3.44218,-0.129881, 9.5, 11.5,'red']



coeffs_all = [z_05_10, z_10_15, z_15_20, z_20_25]


fig, ax = plt.subplots(1,1, figsize = (5, 5.2))
for coeffs in coeffs_all:
    whit_lw = 3
    x_plot = np.linspace(coeffs[3],coeffs[4], 100)
    y_plot = coeffs[0] + coeffs[1] * x_plot + coeffs[2]*x_plot**2.
    ax.plot(x_plot, y_plot, '--', color = coeffs[-1], linewidth = whit_lw, zorder = 0)


ax.set_xlim(7.8, 12)   
ax.set_ylim(-1.4, 3.4)




fig.savefig('/Users/rsimons/Desktop/clear/figures/for_paper/m_sfr_test.png', dpi = 300)



