#! /usr/bin/env python3

import matplotlib.pyplot as plt
import random as rng
import numpy as np
import itertools

grid_x = 50
grid_y = 50
x_max = 5.0
y_max = 5.0
x_step = x_max / grid_x
y_step = y_max / grid_y
n_eta = 2

class Point:
    def __init__(self, x = 0, y = 0):
        self.x = float(x)
        self.y = float(y)

    def __str__(self):
        return "({},{})".format(self.x, self.y)

    def setPoint(self, list_of_two_points):
        if len(list_of_two_points) > 2:
            print("Only a list of two values may be passed in this way!")
            return
        self.x = list_of_two_points[0]
        self.y = list_of_two_points[1]

class Field:
    def __init__(self):
        self.grain_num = 0
        self.center = Point(0,0)
        self.values = [[0.0 for i in range(grid_x)] for j in range(grid_y)]

class LineInfo:
    def __init__(self, slope = 0, y_int = 0, min_x = 0, max_x = 0, is_vertical = False):
        self.slope = float(slope)
        self.y_int = float(y_int)
        self.min_x = float(min_x)
        self.max_x = float(max_x)
        self.is_vertical = is_vertical


def minPeriodicDistance(a, b):
    x = a.x - b.x
    y = a.y - b.y

    #x = x - round(x / grid_x) * grid_x
    #y = y - round(y / grid_y) * grid_y

    return x*x + y*y

def calculateLineInfo(a, b):
    top = float(a[1]) - float(b[1])
    bott = float(a[0]) - float(b[0])

    if bott == 0:
        slope = 0
        y_int = 0
        min_x = max_x = a[0]
        return LineInfo(slope, y_int, min_x, max_x, True)
    else:
        slope = top/bott
        y_int = a[1] - slope * a[0]

        min_x = min([a[0], b[0]])
        max_x = max([a[0], b[0]])
        return LineInfo(slope, y_int, min_x, max_x)



rng.seed(42)

etas = [Field() for i in range(n_eta)]
#points = [[38, 28], [14,42]]#, [7,20], [38,18], [22,10], [10,23], [35,39], [23,2], [21,1], [23,43]]

# for i in range(n_eta):
#     etas[i].center = Point(points[i][0],points[i][1])
points = []
for i in range(n_eta):
    x = rng.random() * x_max
    y = rng.random() * y_max
    etas[i].center = Point(x,y)
    points.append([x,y])
    print("{},{}".format(etas[i].center.x, etas[i].center.y))

combos = list(itertools.combinations(points,2))
lines = [LineInfo() for i in range(len(combos))]

for i in range(len(combos)):
    lines[i] = calculateLineInfo(combos[i][0], combos[i][1])

for i in range(grid_x):
    for j in range(grid_y):
        min_id = None
        _min = 1000000

        #print("Distances for {},{}:".format(i,j))
        for k in range(n_eta):
            distance = minPeriodicDistance(etas[k].center, Point(i * x_step,j * y_step))
            #print("{}: center = {}".format(k,etas[k].center))
            #print("  (index {}, {},{}) = {}".format(k,etas[k].center.x, etas[k].center.y,distance))
            if distance < _min:
                _min = distance
                min_id = k
        #print("Grid {},{} written to eta[{}] (closest point: {},{})".format(i,j,min_id,etas[min_id].center.x, etas[min_id].center.y))
        if min_id is None:
            print("ERROR!")
        etas[min_id].values[i][j] = 1.0
        etas[min_id].grain_num = min_id + 1

fout = open("field.txt", 'w')
for i in range(n_eta):
    for j in range(grid_x):
        for k in range(grid_y):
            fout.write("{} ".format(etas[i].values[j][k] * etas[i].grain_num))
        fout.write("\n")
    fout.write("\n")

fout.close()

plt.figure()
for index,value in enumerate(lines):
    if value.is_vertical:
        plt.vlines(value.min_x, min([j for k,j in combos[index]]), max([j for k,j in combos[index]]), color = 'black')
    else:
        xs = np.linspace(value.min_x, value.max_x, num=100)
        ys = [value.slope * j + value.y_int for j in xs]
        plt.plot(xs, ys, 'k')

xs,ys = zip(*points)
plt.scatter(xs,ys, c = 'green')
plt.show()
