from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

def RSW(x, xmin, xmax, a):
    theta = np.pi / 2 * (x - xmin) / (xmax - xmin)
    return np.sin(theta) * (1 - a * np.log(np.sin(theta)))

theta_min = 0
theta_max = 15 * np.pi / 180
a = 0.5

rswVals = []
angles = []
for i in range(1001):
    angles.append(theta_max * i / 1000*180/np.pi)
    rswVals.append(RSW(theta_max * i / 1000, theta_min, theta_max, a))

plt.plot(angles, rswVals)
plt.gca().annotate(r"$\theta_{min}, E_{min}$", xy = (0,0), xytext = (1, 0.05),
                   fontsize=22)
plt.gca().annotate(r"$\theta_{max}, E_{max}$", xy = (theta_max*180/np.pi,1), xytext = (11.5, 0.9),
                   fontsize=22)
plt.text(6, 0.6,r"$a$", fontsize=22)
plt.text(6.8, 0.6, "is the shaping parameter", fontsize=18)
plt.text(2.8, .4, r"$E_{min} + (E_{max}-E_{min})\ \mathrm{sin}\left(\frac{\pi}{2}\frac{\theta - \theta_{min}}{\theta_{max}-\theta_{min}}\right)\times$", fontsize=21)
plt.text(5.3,.25, r"$\left(1-a\ \mathrm{log}\left(\mathrm{sin}\left(\frac{\pi}{2}\frac{\theta - \theta_{min}}{\theta_{max}-\theta_{min}}\right)\right)\right)$", fontsize=21)
plt.xlabel("Angle (degrees)", fontsize=18)
plt.xlim([0,theta_max*180/np.pi + 10])
plt.ylim([0, 1.1])
plt.ylabel("RSW Value", fontsize=18)
plt.title("RSW Function", fontsize=20)
plt.savefig("../myPapers/Senior Thesis/Images/rsw.jpg")
plt.savefig("../myPapers/Senior Thesis/Images/rsw.eps")
plt.show()
