import numpy as np
import matplotlib.pyplot as plt
from scipy import *

def scatterPlot(x,y):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    print "Scatter Plot"

    ax.scatter(x, y, color = 'b', s = 2)
    plt.show()

def imagePlot(img):    
    print "Imshow Plot"
    plt.imshow(img, cmap=plt.cm.Reds)
    plt.autoscale(True)
    plt.show()

def readFile(f):
    f = open(f, "r")
    l = []
   
    max_x = 0
    max_y = 0

    for line in f:
        if line[0] == "@":
            tokens = line.split(':')
            val1 = tokens[3]
            val2 = tokens[4]
            val3 = tokens[5]
            val4 = tokens[6].split()[0]
            if(val1 == "1" and val2 == "11101"):
                
                x = int(val3)
                y = int(val4)

                if x > max_x:
                    max_x = x

                if y > max_y:
                    max_y = y

                l.append((x, y))

    max_x = max_x/100
    max_y = max_y/100
    print max_x 
    print max_y

    print len(l)

    dims = (max_x+1,max_y+1)

    img = np.zeros(dims)
    for t in l:
        x_coord = t[0]/100
        y_coord = t[1]/100
        img[x_coord][y_coord] = 1
        
    l = None
    print img
    print "Imshow Plot"
    plt.imshow(img, cmap=plt.cm.Reds)
    plt.autoscale(True)
    plt.show()


def readFile2(f):
    f = open(f, "r")

    x_arr = []
    y_arr = []

    for line in f:
        if line[0] == "@":
            tokens = line.split(':')
            val1 = tokens[3]
            val2 = tokens[4]
            val3 = tokens[5]
            val4 = tokens[6].split()[0]
            if(val1 == "1" and val2 == "11101"):
                
                x = int(val3)
                y = int(val4)

                x_arr.append(x)
                y_arr.append(y)

    scatterPlot(x_arr,y_arr)

#readFile("./DNA8590_S47_L001_R1_001.fastq");
readFile("DNA8590_S47_L001_R1_001.fastq");
