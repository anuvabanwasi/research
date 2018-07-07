import numpy as np
import matplotlib.pyplot as plt
from scipy import *

import subprocess
import operator

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

    print  ("length of file " + str(len(lines)))

    l = []

    for i in range(0, len(lines)):
        if(i % 4 == 0):
            tokens = lines[i].split(':')
            val1 = tokens[3]
            val2 = tokens[4]
            val3 = tokens[5]
            val4 = tokens[6].split()[0]

            if(val1 == "1" and val2 == "11101"):
            #if(val1 == "1" and val2 == "1101"):
                x = int(val3)
                y = int(val4)

                l.append((x,y,lines[i+1]))

    max_x = compute_max(l, 0)
    min_x = compute_min(l, 0)

    max_y = compute_max(l, 1)
    min_y = compute_min(l, 1)

    return l, max_x, max_y, min_x, min_y



# create an image matrix. Input is the list l of triples (x,y,z). Creates a 2D matrix of required size
# For every (x,y) in the list, the matrix img is set to 1 
# Returns image 2D matrix, num of rows in the matrix, num of cols in the matrix
def create_image_matrix(l, max_x, max_y, min_x, min_y):

    num_of_rows = max_x - min_x
    num_of_cols = max_y - min_y

    dims = (num_of_rows + 1   , num_of_cols + 1)

    img = np.zeros(dims)
    for t in l:
        x_coord = t[0]/scale_factor - num_of_rows
        y_coord = t[1]/scale_factor - num_of_cols
        img[x_coord][y_coord] = 1
        
    return img, num_of_rows, num_of_cols



# Display plot
def display_image(img, num_of_rows, num_of_cols, color,title):
    global plot_id
    plt.imshow(crop(img, 0, 0, num_of_rows, num_of_cols), cmap=color, interpolation = "nearest")
    plt.autoscale(True)
    filename = "cluster"+str(plot_id)+".pdf"
    plt.title(title)
    plt.savefig(filename, dpi = 300)
    plot_id = plot_id + 1
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

    for i in range(0, len(lines)):
        if(i % 4 == 0):
            tokens = lines[i].split(':')
            val1 = tokens[3]
            val2 = tokens[4]
            val3 = tokens[5]
            val4 = tokens[6].split()[0]

            #if(val1 == "1" and val2 == "11101"):
            if(val1 == "1" and val2 == "1101"):
                x = int(val3)
                y = int(val4)

                f2.write(">%s:%s\n"  % (x, y) )
                f2.write("%s"  % lines[i+1] )




# Find all clusters near an x, y coordinate within a radius of r
# Append the 
def find_nearby_clusters(l, cluster_coords, r, f3):
    f3 = open(f3, "w")
    m = []
    for cluster_coord in cluster_coords:
        x = cluster_coord[0]
        y = cluster_coord[1]

        for t in l:
            if math.sqrt(math.pow((y - t[1]), 2) + math.pow((x - t[0]), 2)) < r:
                m.append(t)
                f3.write(">%s:%s\n"  % (t[0], t[1]) )
                f3.write("%s"  % t[2] )

    return m


# Parse BWA output to keep track of match groups
# Iterate through the output of the bwa command
# if the output starts with @ - it is a header line and therefore ignore it
# if the output does not start with @ - it is an alignment line, therefore parse it
# for each alignment line, tokenize it on \t
# if(tokens[1] == "0"), then it is a new match group
# Continue to read output lines, select tokens[2] and tokenize on ":" to obtain x and y coordinate and add to match group list

def parse_bwa_output(f):

    f = open(f, "r")
   
    lines = f.readlines()

    l = []
    plot_dict = {}
    cnt = 0
    
    for i in range(0, len(lines)):

        if (lines[i][0] != "@"):
           
            tokens = lines[i].split("\t")

            if(tokens[1] == "0"):
                
                seq_id = tokens[0]

                seq_tokens = seq_id.split(":")
                
                seq_x = int(seq_tokens[0])
                seq_y = int(seq_tokens[1])

                start_index = i
                
                i = i + 1

                while( i < len(lines) and lines[i].split("\t")[1] != "0"):
                        i = i + 1

                end_index = i - 1
                
                alignments = lines[start_index:end_index]

                i = i - 1

                l.append((seq_x, seq_y, len(alignments)))

                plot_dict[(seq_x, seq_y)] = handle_alignments(alignments)

                #print plot_dict

    return l, plot_dict


def compute_max(l, index):
    l2 = [i[index] for i in l] # change later
    l3 = np.array(l2)

    max_val = max(l3) / scale_factor

    return max_val


def compute_min(l, index):
    l2 = [i[index] for i in l]
    l3 = np.array(l2)

    min_val = min(l3) / scale_factor

    return min_val


def handle_alignments(lines):

    results = []
    for i in range(0, len(lines)):
        tokens = lines[i].split("\t")
        xy = tokens[2].split(":")

        x = int(xy[0])
        y = int(xy[1])

        results.append((x,y))

    return results


def identify_sequences(l , max_dups):
    
    # convert list l to numpy array L2
    L2 = np.array(l)

    # select all sequences where the num of alignments <= max_dups
    sequences = L2[:, 2] < max_dups

    # above statement returns array of True/False, convert to array of actual values
    valid_sequences = L2[sequences]

    # extract only the x and y coordinates
    plot_points3 = valid_sequences[:, 0:2]

    # convert np array to list of tuples
    plot_points = map(tuple,plot_points3) 

    return plot_points


def identify_alignments(plot_dict , max_val):
    
    res = []
    for key, value in plot_dict.items():
        if (len(value) <= max_val and len(value) != 0):
            res.extend(value)

    return res


#Plot a histogram on l
def plot_histogram(l, i):
    global plot_id

    L2 = np.array(l) 

    a = L2[:, 2]

    fig = plt.figure(i)

    plt.hist(a, bins = [1,2,4,8,16,32,64,128,256,512,1024]) 
    
    plt.title("histogram") 

    filename = "histogram"+str(plot_id)+".pdf"
   
    plt.savefig(filename, dpi = 300)

    plot_id = plot_id + 1

    plt.show()



# constants
scale_factor = 100
#scale_factor = 10
plot_id = 1
radius = [1000, 2000]
#radius = [200,400]
max_dups = [3, 10, 50]


# Read FASTQ file and plot it
l, max_x, max_y, min_x, min_y = readFile("DNA8590_S47_L001_R1_001.fastq")
#l, max_x, max_y, min_x, min_y = readFile("48900_S54_L001_R1_001.fastq")
#l, max_x, max_y, min_x, min_y = readFile("Nia2_ALL_R1.fastq")
img, num_of_rows, num_of_cols = create_image_matrix(l, max_x, max_y, min_x, min_y)
display_image(img, num_of_rows, num_of_cols, plt.cm.Reds,"fastq file")


# Create FASTA file
create_fasta_file("DNA8590_S47_L001_R1_001.fastq", "DNA8590_S47_L001_R1_001.fasta")
#create_fasta_file("48900_S54_L001_R1_001.fastq", "48900_S54_L001_R1_001.fasta")
#create_fasta_file("Nia2_ALL_R1.fastq", "Nia2_ALL_R1.fasta")

# Find nearby clusters(within radius 1000) to a randomly chosen point
# Create a FASTQ file out of the output

cluster_coords = [(10000,10000), (15000, 10000)]    

for r in radius:

    find_nearby_clusters(l, cluster_coords, r, "Clusters.fastq")
    #find_nearby_clusters(l, 10000, 10000, r, "Nia2_ALL_R1.fastq")

    # Call bwa index in python to index FASTA file
    subprocess.call(['bwa', 'index',  'DNA8590_S47_L001_R1_001.fasta'])
    #subprocess.call(['bwa', 'index',  '48900_S54_L001_R1_001.fasta'])
    #subprocess.call(['bwa', 'index',  'Nia2_ALL_R1.fasta'])

    # Call bwa in python on the FASTA file and the cluster FASTQ file
    subprocess.call(['bwa', 'mem', '-a', '-o', 'outpute.txt', 'DNA8590_S47_L001_R1_001.fasta','Clusters.fastq'])
    #subprocess.call(['bwa', 'mem', '-a', '-o', 'outputM.txt', '48900_S54_L001_R1_001.fasta','48900_S54_L001_R1_001.fastq'])
    #subprocess.call(['bwa', 'mem', '-a', '-o', 'outputN.txt', 'Nia2_ALL_R1.fasta','Nia2_ALL_R1.fastq'])

    # Parse the output of bwa
    l2, plot_dict  = parse_bwa_output("outpute.txt")
    
    # Plot the output of bwa
    for val in max_dups:
        plot_points = identify_sequences(l2, val)
        plot_alignments = identify_alignments(plot_dict, val)
        plot_points.extend(plot_alignments)
        img1, num_of_rows1, num_of_cols1 = create_image_matrix(plot_points, max_x, max_y, min_x, min_y)
        title = " cluster_coords : " + str(cluster_coord) + " radius : " + str(r) + " max_dups : " + str(val) 
        display_image(img1, num_of_rows1, num_of_cols1, plt.cm.Blues, title)

    # Create and plot histogram
    plot_histogram(l2, plot_id)



