import numpy as np
import matplotlib.pyplot as plt
from scipy import *

import subprocess
import operator
from random import randint
import os.path
import itertools
import pprint

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
    filename = image_dir + "cluster" + str(plot_id) + ".pdf"
    plt.title(title)
    plt.savefig(filename, dpi = 300)
    plot_id = plot_id + 1
    #plt.show()


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
def find_nearby_clusters(l, seeds, r, f3):
    f3 = open(f3, "w")
    m = []
    for seed in seeds:
        x = seed[0]
        y = seed[1]
        for t in l:
            if math.sqrt(math.pow((y - t[1]), 2) + math.pow((x - t[0]), 2)) <= r:
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
    seq_dict = {}
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

                #print ("length of alignments dd : " , len(alignments))
                
                seq_x, seq_y, found = handle_sequences(alignments)
               
                if(found == "true"):
                    l.append((seq_x, seq_y, len(alignments)))

                seq_dict[(seq_x, seq_y)] = handle_alignments(alignments)

    pp = pprint.PrettyPrinter(indent=4)
    #pp.pprint(seq_dict)

    return l, seq_dict


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

def handle_sequences(lines):

    results = []

    found = "false"

    for i in range(0, len(lines)):
        tokens = lines[i].split("\t")
        
        if (tokens[0] == tokens[2]):

            xy = tokens[0].split(":")

            x = int(xy[0])
            y = int(xy[1])

            found = "true"

            return x, y, found
    
    if (found == "false"):
            return 0, 0, found

        


def find_sequences(l , max_dups):
    
    # convert list l to numpy array L2
    L2 = np.array(l)

    # select all sequences where the num of alignments <= max_dups
    sequences = L2[:, 2] <= max_dups

    # above statement returns array of True/False, convert to array of actual values
    valid_sequences = L2[sequences]

    # extract only the x and y coordinates
    valid_sequences_xy = valid_sequences[:, 0:2]

    # convert np array to list of tuples
    sequences = map(tuple, valid_sequences_xy) 

    return sequences


# Identify (x,y) coordinates of alignments which have number of alignments <= specified max_dup. Create a list of (x,y) where the number of alignments for each sequence <= specified max_dup
def find_alignments(dict , max_val):
    
    alignments = []

    for key, value in dict.items():
        if (len(value) <= max_val and len(value) != 0):
            alignments.extend(value)

    #print ("len of keys of dict ", len(dict))
    return alignments


#Plot a histogram on l
def plot_histogram(l, i):
    global plot_id

    L2 = np.array(l) 

    a = L2[:, 2]

    fig = plt.figure(i)

    plt.hist(a, bins = [1,2,4,8,16,32,64,128,256,512,1024]) 
    
    plt.xlabel('No of matches')
    plt.ylabel('Frequency')

    plt.title("Histogram of number of matches") 

    filename = image_dir + "histogram" + str(plot_id) + ".pdf"
   
    plt.savefig(filename, dpi = 300)

    plot_id = plot_id + 1

    #plt.show()


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

def format(val):
    val = val * 1000000
    val = float("{0:.2f}".format(val))
    return val

def generate_random_seeds(size):
    seeds = []
    for seed_index in range(size):
        seed_x = randint(1000, 18000)
        seed_y = randint(1000, 18000)
        
        seeds.append((seed_x,seed_y))

    return seeds

def generate_evenly_spaced_seeds(size, min_x, max_x, min_y, max_y):
    seeds = []
   
    x1 = np.linspace(min_x*scale_factor, max_x*scale_factor, num=size)
    y1 = np.linspace(min_y*scale_factor, max_y*scale_factor, num=size)

    #l = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000, 15500, 16000, 16500, 17000, 17500, 18000, 18500, 19000, 195000, 20000, 20500,  21000, 21500, 22000, 22500,  23000, 23500, 24000, 24500]
    #y1 = np.array(l)

    for i in itertools.product(x1, y1):
        seeds.append(i)

    #print ("seeds is " , seeds)
    return seeds

# constants


input_data_dir = '../data/input/'
output_data_dir = '../data/output/'
image_dir = '../data/images/'

plot_id = 1
max_dups = [5]
#max_dups = [1]

'''
base_file_name = 'single_point'
fastq_file_name = input_data_dir + 'single_point.fastq'
fasta_file_name = input_data_dir + 'single_point.fasta'
scale_factor = 100
radius = [1000]
ax = '1'
ay = '11101'
'''

'''
base_file_name = 'DNA8590_S47_L001_R1_001'
fastq_file_name = input_data_dir + 'DNA8590_S47_L001_R1_001.fastq'
fasta_file_name = input_data_dir + 'DNA8590_S47_L001_R1_001.fasta'
scale_factor = 100
#radius = [1000,2000]
radius = [1000]
ax = '1'
ay = '11101'
'''


'''
base_file_name = '48900_S54_L001_R1_001'
fastq_file_name = '48900_S54_L001_R1_001.fastq'
fasta_file_name = '48900_S54_L001_R1_001.fasta'
scale_factor = 100
radius = [1000,2000]
ax = '1'
ay = '1101'
'''


base_file_name = 'Nia2_ALL_R1'
fastq_file_name = input_data_dir + 'Nia2_ALL_R1.fastq'
fasta_file_name = input_data_dir + 'Nia2_ALL_R1.fasta'
scale_factor = 10
#radius = [100,200]
radius = [50, 100]
ax = '1'
ay = '1101'



# Read FASTQ file and plot it
l, max_x, max_y, min_x, min_y = readFile(fastq_file_name, ax, ay)
img, num_of_rows, num_of_cols = create_image_matrix(l, max_x, max_y, min_x, min_y)
display_image(img, num_of_rows, num_of_cols, plt.cm.Reds, "Original fastq file")


# Create FASTA file
create_fasta_file(fastq_file_name, fasta_file_name)
   
# Call bwa index in python to index FASTA file
subprocess.call(['bwa', 'index',  fasta_file_name]) 

# Find nearby clusters(within radius 1000) to a randomly chosen point
# Create a FASTQ file out of the output

#co = [(5000,5000), (7000, 7000), (9000, 9000), (10000, 10000), (10000, 12000), (12000, 12000), (15000, 15000), (18000, 18000), (18000, 19000), (16000, 15000), (15000,10000), (12000,11000), (14000, 10000), (17000, 17000), (6000, 15000), (10000, 15000), (18000, 10000), (14000, 18000), (11000, 18000), (9000, 10000), (7000,14000), (6000,11000), (11000, 6000), (14000, 12000), (14000, 14000)]

size = [7, 10]

for N in size:

    seeds = generate_evenly_spaced_seeds(N, min_x, max_x, min_y, max_y)

    #seeds = [(1000, 1000)]

    #seeds = [(26700, 1000)]

    #seeds = [(26700, 20300)]
    
    #seeds = [(1000, 20300)]

    #seeds = [(1000, 1000), (1000, 20300), (26700, 1000), (26700, 20300)]

    #print ("seeds ", seeds)


    for r in radius:

        clusters_fasta_file_name = output_data_dir + 'clusters_' + base_file_name + "_" + str(r) + "_" + str(N) + '.fasta'

        m = find_nearby_clusters(l, seeds, r, clusters_fasta_file_name)

        print ("\n\clusters size : " , len(m))

        # Create output file names

        output_file_name = output_data_dir + 'output_'  + base_file_name + "_" + str(r) +  "_" + str(N) + '.txt'
        
        output_filtered_file_name = output_data_dir + 'output_filtered_' + base_file_name + "_" + str(r) + "_" + str(N) + '.txt'


        if not os.path.exists(output_file_name):

            # Call bwa in python on the FASTA file and the cluster FASTQ file
            subprocess.call(['bwa', 'mem', '-a', '-o', output_file_name, fasta_file_name, clusters_fasta_file_name])
        
            COMMAND = "cat " +  output_file_name  + " | awk -F '\t' '{if ($1 ~ /@/) {next}; print}' | awk -F '\t' '{if ($6 ~ /H/) {next}; print}' > " + output_filtered_file_name

            subprocess.call(COMMAND, shell=True)

        else:
            print ("\n\nAlready exists")

        # Parse the output of bwa
        sequence_list, seq_dict  = parse_bwa_output(output_filtered_file_name)
        
        r1 = (max_x * scale_factor - min_x * scale_factor) / 2.0 
        r2 = (max_y * scale_factor - min_y * scale_factor) / 2.0 

        print("max_x : ", max_x , "min_x : ", min_x , "max_y : ", max_y ," min_y : " , min_y)
        print ("r1 : " + str(r1) + " r2 : " + str(r2))

        area_circle = compute_area_of_circle(r1)
        theta, area_sector = compute_area_of_sector(r1, r2)
        area_triangle = compute_area_of_triangle(r1, theta)

        #print ("area of circle ", area_circle)
        #print ("area of triangle ", area_triangle)
        #print ("area of sector " , area_sector)


        #overall_density_denominator = area_circle - 2 * area_sector + 2 * area_triangle

        overall_density_denominator = area_circle
        
        # Plot the output of bwa
        for max_dup in max_dups:
            sequences = find_sequences(sequence_list, max_dup)
            print ("length of sequences ", len(sequences) )
            
            alignments = find_alignments(seq_dict, max_dup)
            print ("length of alignments ", len(alignments) )

            all_points = []
            all_points.extend(sequences)
            all_points.extend(alignments)
            
            all_points_set = set(all_points)
            all_points = list(all_points_set)

            print ("length of all_points ", len(all_points) )

            print ("len(seeds) : ", len(seeds) , " r : " , r)

            img1, num_of_rows1, num_of_cols1 = create_image_matrix(all_points, max_x, max_y, min_x, min_y)

            # Calculate density
            numerator = len(sequences)

            denominator = len(seeds) * math.pi * (r) * (r)
            
            density =  numerator / denominator

            #density = format(density)

            overall_density_numerator = len(all_points)
            
            overall_density = overall_density_numerator / overall_density_denominator

            #overall_density = format(overall_density)

            ratio = density / overall_density

            #ratio = float("{0:.2f}".format(ratio))

            title = "radius : " + str(r) + " max_dups : " + str(max_dup) + " density : " + str(density) + " overall density : " + str(overall_density) + " ratio: " + str(ratio)

            display_image(img1, num_of_rows1, num_of_cols1, plt.cm.Blues, title)

            print( " denominator: ", denominator, " numerator: ", numerator , " overall_density_denominator : ", overall_density_denominator,  "overall_density_numerator: ", overall_density_numerator, " density : " , density, " overall density : " , overall_density , " ratio : ", ratio )
            

        # Create and plot histogram
        plot_histogram(sequence_list, plot_id)




