import numpy as np
import matplotlib.pyplot as plt
from scipy import *

import subprocess
import operator
from random import randint

# Read FASTQ file
# Tokenize every line that starts with '@'. if tokens[3] == '1' and tokens[4] == '11101', capture tokens[4] and tokens[5] as (x,y) coordinates
# Determine the max value of x and max value of y coordinate
# Determine the min value of x and min value of y coordinate
# Create a list of triples (x,y,z) where x is the x coordinate, y is the y coordinate and z is the sequence
# Return the list of triples l, max x coordinate, max y coordinate,  min x coordinate, min y coordinate

def readFile(f, ax, ay):

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

            if(val1 == ax and val2 == ay):
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

    #num_of_rows = max_x - min_x
    #num_of_cols = max_y - min_y

    num_of_rows = max_x
    num_of_cols = max_y

    dims = (num_of_rows + 1   , num_of_cols + 1)

    img = np.zeros(dims)
    for t in l:
        #x_coord = t[0]/scale_factor - num_of_rows
        #y_coord = t[1]/scale_factor - num_of_cols
        x_coord = t[0]/scale_factor 
        y_coord = t[1]/scale_factor 
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

            if(val1 == ax and val2 == ay):
                x = int(val3)
                y = int(val4)

                f2.write(">%s:%s\n"  % (x, y) )
                f2.write("%s"  % lines[i+1] )




# Find all clusters near an x, y coordinate within a radius of r
# Append the 
def find_nearby_clusters(l, col, r, f3):
    f3 = open(f3, "w")
    m = []
    for co in col:
        x = co[0]
        y = co[1]
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
                
                start_index = i
                
                i = i + 1

                while( i < len(lines) and lines[i].split("\t")[1] != "0"):
                        i = i + 1

                end_index = i - 1
                
                alignments = lines[start_index:end_index+1]

                i = i - 1

                #print ("alignments is " , len(alignments))
                
                seq_x, seq_y, found = handle_sequences(alignments)
               
                if(found == "true"):
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

def handle_lines(lines):

    alignments = []
    seq_xy = []

    for i in range(0, len(lines)):
        tokens = lines[i].split("\t")

        if(tokens[0] == tokens[2]):
            seq_xy = tokens[0].split(":")

        align_xy = tokens[2].split(":")

        align_x = int(align_xy[0])
        align_y = int(align_xy[1])

        alignments.append((align_x,align_y))

    return seq_xy, alignments

def handle_alignments(lines):

    results = []
    for i in range(0, len(lines)):
        tokens = lines[i].split("\t")
        xy = tokens[2].split(":")

        x = int(xy[0])
        y = int(xy[1])

        results.append((x,y))

    return results

def handle_sequences(lines):

    results = []

    found = "false"

    for i in range(0, len(lines)):
        tokens = lines[i].split("\t")
        
        if (tokens[0] != tokens[2] and len(lines) == 1):
            print ("mismatch " , tokens[0] , " " , tokens[2])

        if (tokens[0] == tokens[2]):

            xy = tokens[0].split(":")

            x = int(xy[0])
            y = int(xy[1])

            found = "true"

            return x, y, found
    
    if (found == "false"):
            return 0, 0, found

        


def identify_sequences(l , max_dups):
    
    # convert list l to numpy array L2
    L2 = np.array(l)

    # select all sequences where the num of alignments <= max_dups
    sequences = L2[:, 2] <= max_dups

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

    print ("len of keys of dict ", len(plot_dict))
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


def compute_area_of_circle(r):
    area = math.pi * r * r
    return area

def compute_area_of_sector(r1,r2):
    theta = math.acos(r2/r1) * 2
    area_of_sector = (theta * r1 * r1) / 2
    return theta, area_of_sector

def compute_area_of_triangle(r1, theta):
    area_of_triangle = (r1 * r1 * math.sin(theta) ) / 2
    return area_of_triangle


# constants

plot_id = 1
max_dups = [1, 3, 10, 50]
output_file_name = 'output.txt'
output_filtered_file_name = 'output_filtered.txt'
clusters_fastq_file_name = 'Clusters.fasta'

fastq_file_name = 'DNA8590_S47_L001_R1_001.fastq'
fasta_file_name = 'DNA8590_S47_L001_R1_001.fasta'
scale_factor = 100
radius = [1000,2000]
ax = '1'
ay = '11101'

#fastq_file_name = '48900_S54_L001_R1_001.fastq'
#fasta_file_name = '48900_S54_L001_R1_001.fasta'
#scale_factor = 100
#radius = [1000,2000]
#ax = '1'
#ay = '1101'

#fastq_file_name = 'Nia2_ALL_R1.fastq'
#fasta_file_name = 'Nia2_ALL_R1.fasta'
#scale_factor = 10
#radius = [100,200]
#ax = '1'
#ay = '1101'

# Read FASTQ file and plot it
l, max_x, max_y, min_x, min_y = readFile(fastq_file_name, ax, ay)
img, num_of_rows, num_of_cols = create_image_matrix(l, max_x, max_y, min_x, min_y)
display_image(img, num_of_rows, num_of_cols, plt.cm.Reds, "fastq file")


# Create FASTA file
create_fasta_file(fastq_file_name, fasta_file_name)
    
# Find nearby clusters(within radius 1000) to a randomly chosen point
# Create a FASTQ file out of the output

#co = [(5000,5000), (7000, 7000), (9000, 9000), (10000, 10000), (10000, 12000), (12000, 12000), (15000, 15000), (18000, 18000), (18000, 19000), (16000, 15000), (15000,10000), (12000,11000), (14000, 10000), (17000, 17000), (6000, 15000), (10000, 15000), (18000, 10000), (14000, 18000), (11000, 18000), (9000, 10000), (7000,14000), (6000,11000), (11000, 6000), (14000, 12000), (14000, 14000)]

co = []
for seed_index in range(50):
    seed_x = randint(1000, 18000)
    seed_y = randint(1000, 18000)
    co.append((seed_x,seed_y))

print ("co is " , co)
#co = [(10000, 10000)]
for r in radius:

    find_nearby_clusters(l, co, r, clusters_fastq_file_name)

    # Call bwa index in python to index FASTA file
    subprocess.call(['bwa', 'index',  fasta_file_name])
    
    # Call bwa in python on the FASTA file and the cluster FASTQ file
    subprocess.call(['bwa', 'mem', '-a', '-o', output_file_name, fasta_file_name, clusters_fastq_file_name])
    
    COMMAND = "cat " +  output_file_name  + " | awk -F '\t' '{if ($1 ~ /@/) {next}; print}' | awk -F '\t' '{if ($6 ~ /H/) {next}; print}' > " + output_filtered_file_name

    print("COMMAND ", COMMAND)

    subprocess.call(COMMAND, shell=True)

    # Parse the output of bwa
    l2, plot_dict  = parse_bwa_output(output_filtered_file_name)
    
    r1 = (max_x - min_x) / 2.0
    r2 = (max_y - min_y) / 2.0

    area_circle = compute_area_of_circle(r1)
    theta, area_sector = compute_area_of_sector(r1, r2)
    area_triangle = compute_area_of_triangle(r1, theta)

    print ("area of circle ", area_circle)
    print ("area of triangle ", area_triangle)
    print ("area of sector " , area_sector)


    overall_density_denominator = area_circle - 2 * area_sector + 2 * area_triangle


    # Plot the output of bwa
    for val in max_dups:
        plot_points = identify_sequences(l2, val)

        print ("length of plot_points ", len(plot_points) )
        plot_alignments = identify_alignments(plot_dict, val)
        plot_points.extend(plot_alignments)

        print("max_x : ", max_x , " max_y : " , max_y)
        img1, num_of_rows1, num_of_cols1 = create_image_matrix(plot_points, max_x, max_y, min_x, min_y)

        # Calculate density
        numerator = len(plot_points)
        denominator = len(co) * math.pi * r * r
        density =  numerator / denominator

        overall_density_numerator = len(plot_points)
        overall_density = overall_density_numerator / overall_density_denominator

        title = "radius : " + str(r) + " max_dups : " + str(val) + " density : " + str(density) + " overall density : " + str(overall_density) 
        # print "The densities are %.3f" % (first denisty, second density)
        display_image(img1, num_of_rows1, num_of_cols1, plt.cm.Blues, title)

        #print ("plot_dict ", plot_dict)
        print( " denominator ", denominator, " numerator ", numerator , " density " , density, " overall density " , overall_density )
        

    # Create and plot histogram
    plot_histogram(l2, plot_id)




