import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator,FuncFormatter

x=np.arange(-5,5,.001)
eta_1 = 0.5*(1 - np.tanh(sqrt(0.5)*x))
eta_2 = 0.5*(1 + np.tanh(sqrt(0.5)*x))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-5,5)
ax.set_ylim(-0.05,1.05)

ax.xaxis.set_major_locator(MultipleLocator(5.0))
ax.yaxis.set_major_locator(MultipleLocator(1.0))
ax.tick_params(which='major', width = 1.0, length = 10.0)

ax.set_title(r'$\eta$ Profiles', fontsize=20)
ax.plot(x,eta_1, 'k', lw=2, label = r'$\eta_1$')
ax.plot(x,eta_2, 'k', lw=2, label = r'$\eta_2$')

ax.text(-4.5,1.1, r'$\eta_1$', fontsize=16, ha="right")
ax.text(4.5,1.1, r'$\eta_2$', fontsize=16, ha="left")

plt.savefig("eta_profiles.eps")
plt.show()
