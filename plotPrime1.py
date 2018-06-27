import numpy as np
import matplotlib.pyplot as plt
from scipy import *

import subprocess
import operator

scale_factor = 100;


# Read FASTQ file
# Tokenize every line that starts with '@'. if tokens[3] == '1' and tokens[4] == '11101', capture tokens[4] and tokens[5] as (x,y) coordinates
# Determine the max value of x and max value of y coordinate
# Determine the min value of x and min value of y coordinate
# Create a list of triples (x,y,z) where x is the x coordinate, y is the y coordinate and z is the sequence
# Return the list of triples l, max x coordinate, max y coordinate,  min x coordinate, min y coordinate

def readFile(f):

    max_x = 0
    max_y = 0
    min_x = 9999999
    min_y = 9999999

    f = open(f, "r")
   
    lines = f.readlines()

    print ("length of file " + str(len(lines)))

    l = []

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

    return l, max_x, max_y, min_x, min_y



# create an image matrix. Input is the list l of triples (x,y,z). Creates a 2D matrix of required size
# For every (x,y) in the list, the matrix img is set to 1 
# Returns image 2D matrix, num of rows in the matrix, num of cols in the matrix
def create_image_matrix(l, max_x, max_y, min_x, min_y):

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

    return img, num_of_rows, num_of_cols



# Display image
def display_image(img, num_of_rows, num_of_cols, color):
    plt.imshow(crop(img, 0, 0, num_of_rows, num_of_cols), cmap=color, interpolation = "nearest")
    plt.autoscale(True)
    plt.savefig("cluster.pdf", dpi = 300)
    plt.show()


def crop(img, x0, y0, height, width):
    return img[y0:y0+height , x0:x0+width][:]



# Create FASTA file 
# First input f1 is the FASTQ file, second input f2 is the FASTA file
# Read the FASTQ file 
# Tokenize every line that starts with '@'. if tokens[3] == '1' and tokens[4] == '11101', capture tokens[4] and tokens[5] as (x,y) coordinates
# Append 2 lines to the FASTA file where the first line is of the format x:y and the next line is the sequence

def create_fasta_file(f1, f2):
    f1 = open(f1, "r")
    f2 = open(f2, "w")
   
    lines = f1.readlines()

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

                f2.write(">%s:%s\n"  % (x, y) )
                f2.write("%s"  % lines[i+1] )




# Find all clusters near an x, y coordinate within a radius of r
# Append the 
def find_nearby_clusters(l, x, y, r, f3):
    f3 = open(f3, "w")
    m = []
    for t in l:
        if math.sqrt(math.pow((y - t[1]), 2) + math.pow((x - t[0]), 2)) < r:
            m.append(t)
            f3.write(">%s:%s\n"  % (t[0], t[1]) )
            f3.write("%s"  % t[2] )

    #print m
    return m



# Plot the output of bwa command run with FASTA file and FASTQ file created(Clusters.fastq)
# Iterate through the output    
# if the output starts with @ - it is a header line and therefore ignore it
# if the output does not start with @ - it is an alignment line, therefore parse it
# for each alignment line, tokenize it on \t
# select tokens[2] and tokenize on ":" to obtain x and y coordinate
# Append (x,y) to a list

def parse_bwa_output(f):

    max_x = 0
    max_y = 0
    min_x = 9999999
    min_y = 9999999

    f = open(f, "r")
   
    lines = f.readlines()

    print ("length of file bwa output" + str(len(lines)))

    l = []

    for line in lines:

        if (line[0] != "@"):
            tokens = line.split("\t")

            xy = tokens[2].split(":")

            x = int(xy[0])
            y = int(xy[1])

            if x > max_x:
                max_x = x

            if y > max_y:
                max_y = y

            if x < min_x:
                min_x = x

            if y < min_y:
                min_y = y

            l.append((x,y))

    
    max_x = max_x/scale_factor
    max_y = max_y/scale_factor

    min_x = min_x/scale_factor
    min_y = min_y/scale_factor

    print max_x 
    print max_y
    
    
    return l, max_x, max_y, min_x, min_y


# Parse BWA output to keep track of match groups
# Iterate through the output of the bwa command
# if the output starts with @ - it is a header line and therefore ignore it
# if the output does not start with @ - it is an alignment line, therefore parse it
# for each alignment line, tokenize it on \t
# if(tokens[1] == "0"), then it is a new match group
# Continue to read output lines, select tokens[2] and tokenize on ":" to obtain x and y coordinate and add to match group list

def parse_bwa_output_for_histogram(f):

    f = open(f, "r")
   
    lines = f.readlines()

    print ("length of file bwa output" + str(len(lines)))

    matchResults = []

    for i in range(0, len(lines)):

        if (lines[i][0] != "@"):
           
            tokens = lines[i].split("\t")

            if(tokens[1] == "0"):
                
                start_index = i
                i = i + 1

                while( i < len(lines) and lines[i].split("\t")[1] != "0"):
                        i = i + 1

                end_index = i - 1
                
                match_group = lines[start_index:end_index]

                #print ("match group ", match_group)

                i = i - 1

                matchResults.append(handle_match_group(match_group))

    return matchResults


def handle_match_group(lines):
    matchGroup = []

    for i in range(0, len(lines)):
        tokens = lines[i].split("\t")
        xy = tokens[2].split(":")

        x = int(xy[0])
        y = int(xy[1])

        matchGroup.append((x,y))

    return matchGroup

# Create a histogram
def create_histogram(l):
    counts = {}
    for i in range(0, len(l)):
        key = len(l[i])
        if not key in counts:
            counts[key] = 1
        else:
            counts[key] = counts[key] + 1

    print counts

    return counts

# sort the dictionary of counts and create x and y axis labels
def create_histogram_xaxis_and_yaxis_data(counts):
    # the histogram of the data
    x_label = []
    y_label = []

    for key, value in sorted(counts.iteritems()):
        print "%s: %s" % (key, value)
        x_label.append(key)
        y_label.append(value)

    return x_label, y_label

# Plot the chart
def plot_bar_x(x_label, y_label):

    index = np.arange(len(x_label))
    plt.bar(index, y_label)
    plt.xlabel('Num of matches', fontsize=5)
    plt.ylabel('Count', fontsize=5)
    plt.xticks(index, x_label, fontsize=5, rotation=30)
    plt.title('Number of sequences matching')
    plt.show()


# Read FASTQ file and plot it
l, max_x, max_y, min_x, min_y = readFile("DNA8590_S47_L001_R1_001.fastq")
img, num_of_rows, num_of_cols = create_image_matrix(l, max_x, max_y, min_x, min_y)
display_image(img, num_of_rows, num_of_cols, plt.cm.Reds)

# Create FASTA file
create_fasta_file("DNA8590_S47_L001_R1_001.fastq", "DNA8590_S47_L001_R1_001.fasta")

# Find nearby clusters(within radius 1000) to a randomly chosen point
# Create a FASTQ file out of the output
find_nearby_clusters(l, 10000, 10000, 1000, "Clusters.fastq")

# Call bwa index in python to index FASTA file
subprocess.call(['bwa', 'index',  'DNA8590_S47_L001_R1_001.fasta'])

# Call bwa in python on the FASTA file and the cluster FASTQ file
subprocess.call(['bwa', 'mem', '-a', '-o', 'outpute.txt', 'DNA8590_S47_L001_R1_001.fasta','Clusters.fastq'])

# Plot the output of bwa
l1, max_x1, max_y1, min_x1, min_y1 = parse_bwa_output("outpute.txt")
img1, num_of_rows1, num_of_cols1 = create_image_matrix(l1, max_x1, max_y1, min_x1, min_y1)
display_image(img1, num_of_rows1, num_of_cols1, plt.cm.Blues)

# Create and plot histogram
l2 = parse_bwa_output_for_histogram("outpute.txt")
counts = create_histogram(l2)
x_label, y_label = create_histogram_xaxis_and_yaxis_data(counts)
plot_bar_x(x_label, y_label)