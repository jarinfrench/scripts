import matplotlib.pyplot as plt
import numpy as np
import math

data40 = []
data60 = []
data80 = []
data95 = []

i = 1
with open("post_anneal_energies_with_n.txt") as f:
   no_gb = f.readline().split()
   for line in f:
       if int(line.split()[1]) == 40:
           data40.append(line.split())
           i += 1
       elif int(line.split()[1]) == 60:
           data60.append(line.split())
           i += 1
       elif int(line.split()[1]) == 80:
           data80.append(line.split())
           i += 1
       elif int(line.split()[1]) == 95:
           data95.append(line.split())
           i += 1
       else:
           print("Unrecognized value: {}".format(line))

data40 = np.array([[int(i[0]), int(i[1]), float(i[2]), int(i[3])] for i in data40])
data60 = np.array([[int(i[0]), int(i[1]), float(i[2]), int(i[3])] for i in data60])
data80 = np.array([[int(i[0]), int(i[1]), float(i[2]), int(i[3])] for i in data80])
data95 = np.array([[int(i[0]), int(i[1]), float(i[2]), int(i[3])] for i in data95])
no_gb = [int(no_gb[0]), int(no_gb[1]), float(no_gb[2]), int(no_gb[3])]

energies40 = [(i[2] - no_gb[2]/no_gb[3]*i[3])/(2*math.pi*15*27.265000) for i in data40]
energies60 = [(i[2] - no_gb[2]/no_gb[3]*i[3])/(2*math.pi*15*27.265000) for i in data60]
energies80 = [(i[2] - no_gb[2]/no_gb[3]*i[3])/(2*math.pi*15*27.265000) for i in data80]
energies95 = [(i[2] - no_gb[2]/no_gb[3]*i[3])/(2*math.pi*15*27.265000) for i in data95]

plt.figure()
plt.plot(data40[:,0], energies40, 'r-.', label = "40%cutoff")
plt.plot(data60[:,0], energies60, 'b-', label = "60%cutoff")
plt.plot(data80[:,0], energies60, 'g.', label = "80%cutoff")
plt.plot(data95[:,0], energies60, 'k-.', label = "95%cutoff")
plt.legend()
plt.show()
