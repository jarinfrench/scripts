#! /opt/moose/miniconda/bin/python
# This script will read in the values of two output files, and then output the
# results to another file to compare the related values side by side.
# This script assumes that the output files have the same number of lines.
from sys import argv

script, file1, file2 = argv

i = 0 # An iterator
f1 = open(file1, "r")
index1 = [None]*100
index2 = [None]*100
distance1 = [None]*100
distance2 = [None]*100
ksi1 = [None]*100
ksi2 = [None]*100
eta1 = [None]*100
eta2 = [None]*100
phi1 = [None]*100
phi2 = [None]*100

# Read through the first file
while True:
    line1 = f1.readline()
    if not line1:
        break
    data1 = line1.split()
    index1[i] = float(data1[0])
    distance1[i] = float(data1[1])
    ksi1[i] = float(data1[2])
    eta1[i] = float(data1[3])
    phi1[i] = float(data1[4])
    i += 1

f1.close() # Make sure to close it

i = 0
# Open the second file
f2 = open(file2, "r")

# Read through the second file
while True:
    line2 = f2.readline()
    if not line2:
        break
    data2 = line2.split()
    index2[i] = float(data2[0])
    distance2[i] = float(data2[1])
    ksi2[i] = float(data2[2])
    eta2[i] = float(data2[3])
    phi2[i] = float(data2[4])
    i += 1
f2.close() # Close the second file

j = 0 # Another iterator
#with open('/home/frenjc/projects/marmot/examples/data_compare.tex', "w") as fp
f3 = open('/home/frenjc/projects/marmot/examples/data_compare.tex', "w")

f3.write("File1: %s\nFile2: %s\n\n"%(file1,file2))
f3.write('Index\t Distance\t Ksi\t Eta\t Phi\n')
while j < i:
    f3.write('%1.0f\t\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n'%(index1[j], distance1[j], ksi1[j], eta1[j], phi1[j]))
    f3.write('%1.0f\t\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n\n'%(index2[j], distance2[j], ksi2[j], eta2[j], phi2[j]))
    j += 1

f3.close()
