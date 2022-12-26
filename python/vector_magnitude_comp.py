import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

rdir = "/media/jarinf/Research1/uo2/grain_growth/cylindrical/"
data1 = np.loadtxt(rdir + "110/20degree/Basak/T2800/large_r/big/dir_1/10000_unwrappedto3000000_unwrapped_displacement_data.dat",
                   dtype = {'names': ("id", "type", "charge", "xu", "yu", "zu", "x_vec", "y_vec", "z_vec", "vec_magnitude"),
                            'formats': ('int', 'int', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float')},
                   skiprows=1)
data2 = np.loadtxt(rdir + "111/20degree/Basak/T2800/large_r/dir_1/10000_unwrappedto910000_unwrapped_displacement_data.dat",
                   dtype = {'names': ("id", "type", "charge", "xu", "yu", "zu", "x_vec", "y_vec", "z_vec", "vec_magnitude"),
                            'formats': ('int', 'int', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float')},
                   skiprows=1)

a0 = 5.454
fnn_001 = a0 / np.sqrt(2)
fnn_110 = a0 * np.sqrt(2)/2
fnn_111 = a0 * np.sqrt(3)
fnn_112 = a0 * np.sqrt(6)
select1 = (data1['vec_magnitude'] <= 3 * a0) & (data1['type'] == 1) & (data1['vec_magnitude'] >= 3.0)
select2 = (data2['vec_magnitude'] <= 3 * a0) & (data2['type'] == 1) & (data2['vec_magnitude'] >= 3.0)

fig = plt.figure(figsize=(10,8))
grid = plt.GridSpec(3,3, wspace = 0.5, hspace = 1)
ax1 = fig.add_subplot(grid[0,0])
ax2 = fig.add_subplot(grid[0,1])
ax3 = fig.add_subplot(grid[0,2])
ax4 = fig.add_subplot(grid[1:,0:])

n1, _, _ = ax1.hist(x=abs(data1['x_vec'][select1]), bins='auto', alpha = 0.5)
n2, _, _ = ax1.hist(x=abs(data2['x_vec'][select2]), bins='auto', alpha = 0.5)
ax1.axvline(np.mean(abs(data1['x_vec'][select1])), color = colors[0], linestyle = ':')
ax1.axvline(np.mean(abs(data2['x_vec'][select2])), color = colors[1], linestyle = ':')
ax1.set_title('X direction')

n1, _, _ = ax2.hist(x=abs(data1['y_vec'][select1]), bins='auto', alpha = 0.5)
n2, _, _ = ax2.hist(x=abs(data2['y_vec'][select2]), bins='auto', alpha = 0.5)
ax2.axvline(np.mean(abs(data1['y_vec'][select1])), color = colors[0], linestyle = ':')
ax2.axvline(np.mean(abs(data2['y_vec'][select2])), color = colors[1], linestyle = ':')
ax2.set_title('Y direction')

n1, _, _ = ax3.hist(x=abs(data1['z_vec'][select1]), bins='auto', alpha = 0.5)
n2, _, _ = ax3.hist(x=abs(data2['z_vec'][select2]), bins='auto', alpha = 0.5)
ax3.axvline(np.mean(abs(data1['z_vec'][select1])), color = colors[0], linestyle = ':')
ax3.axvline(np.mean(abs(data2['z_vec'][select2])), color = colors[1], linestyle = ':')
ax3.set_xlabel('Vector Magnitude')
ax3.set_title('Z direction')

n1, _, _ = ax4.hist(x=abs(data1['vec_magnitude'][select1]), bins='auto', alpha = 0.5, label = '<110> 20degree')
n2, _, _ = ax4.hist(x=abs(data2['vec_magnitude'][select2]), bins='auto', alpha = 0.5, label = '<111> 20degree')
ax4.axvline(fnn_001, 0, max(max(n1),max(n2)), color = 'k', linestyle = '--', label = '1nn distance')
ax4.axvline(np.mean(abs(data1['vec_magnitude'][select1])), color = colors[0], linestyle = ':', label = 'Average magnitude 110')
ax4.axvline(np.mean(abs(data2['vec_magnitude'][select2])), color = colors[1], linestyle = ':', label = 'Average magnitude 111')
ax4.set_xlabel('Vector Magnitude')
ax4.set_title('Overall Magnitude')
ax4.legend()
plt.show()
