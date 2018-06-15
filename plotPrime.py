import numpy as np
import matplotlib.pyplot as plt
from scipy import *

l = []
scale_factor = 100;


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
   
    max_x = 0
    max_y = 0
    min_x = 9999999
    min_y = 9999999
   
   
    lines = f.readlines()

    print ("length of file " + str(len(lines)))

    for i in range(0, len(lines)):
        if(i % 4 == 0):
            tokens = lines[i].split(':')
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

            if x < min_x:
                min_x = x

            if y < min_y:
                min_y = y

            l.append((x,y,lines[i+1]))


    max_x = max_x/scale_factor
    max_y = max_y/scale_factor

    min_x = min_x/scale_factor
    min_y = min_y/scale_factor
   
    #print max_x 
    #print max_y
    
    print len(l)
    createList(max_x, max_y, min_x, min_y)

def createList(max_x, max_y, min_x, min_y):

    num_of_rows = max_x - min_x
    num_of_cols = max_y - min_y

    dims = (num_of_rows + 1   , num_of_cols + 1)

    #dims = (max_x+1,max_y+1)

    img = np.zeros(dims)
    for t in l:
        x_coord = t[0]/scale_factor - num_of_rows
        y_coord = t[1]/scale_factor - num_of_cols
        img[x_coord][y_coord] = 1
        
    print img
    print "Imshow Plot"

   

    plt.imshow(crop(img, 0, 0, num_of_rows, num_of_cols), cmap=plt.cm.Reds, interpolation = "nearest")
    plt.autoscale(True)
    #ax.plot(range(11))
    plt.savefig("cluster.pdf", dpi = 300)
    plt.show()

    # plt.imshow(img, cmap=plt.cm.Reds, interpolation = "nearest")
    # extent=[min_x,max_x,min_y,max_y]
    # fig, ax = plt.subplots()
    # ax.autoscale(False)
    # ax.plot(range(11))
    # plt.savefig("cluster.pdf", dpi = 300)
    # plt.show()

def crop(img, x0, y0, height, width):
    return img[y0:y0+height , x0:x0+width][:]

def find_nearby_clusters(x_pos, y_pos, radius):
    m = []
    for t in l:
        if math.sqrt(math.pow((y_pos - t[1]), 2) + math.pow((x_pos - t[0]), 2)) < radius:
            m.append(t)

    print m
    return m

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
readFile("DNA8590_S47_L001_R1_001.fastq")
find_nearby_clusters(10000, 10000, 20000)

