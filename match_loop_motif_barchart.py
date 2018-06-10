from numpy import *
from collections import defaultdict
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# match loop anchor and motif
# for each chromose ranging from 1 .. 23
#    get row vector for motif
#    get col vector for loop anchor. Transpose operator used to create the column vector
#    compute absolute differences 
#    select only those where the distances < 20000
#    obtain corresponding i and j
def match_loop_motif(motif_dict, loop_anchor_dict):
    dict_cla_m = defaultdict(list)

    plt_x_arr = []
    for i in range(1,23):
        
        motif_vector = array([motif_dict[str(i)]])
        loop_anchor_vector = array([loop_anchor_dict[str(i)]]).T
        
        dists = absolute(loop_anchor_vector - motif_vector)

        close_ones = dists < 20000

        where_close_ones = where(close_ones)

        # where_close_ones[0] : contains indexes in motif_vector with a value of True
        # where_close_ones[1] : contains indexes in loop_anchor_vector with a value of True
        motif_index = where_close_ones[0]
        loop_anchor_index = where_close_ones[1]

        # create a dictionary
        # key : loop anchor position
        # value : motif postion
        for a,b in zip(loop_anchor_index, motif_index):
            dict_cla_m[loop_anchor_vector[b][0]].append(motif_vector[0][a])
            
        # for testing purposes
        #break

    # x axis of histogram : list of lengths of values associated with each key
    for key in dict_cla_m:
        plt_x_arr.append(len(dict_cla_m[key]))
	
    # draw the histogram
    plotHistogram(dict_cla_m, plt_x_arr)

def plotHistogram(dict_cla_m, plt_x_arr):
    num_bins = len(dict_cla_m)
    n, bins, patches = plt.hist(plt_x_arr, num_bins, facecolor='blue', alpha=0.5)
    plt.show()

# read the file f and create dictionary dict. 
# start_pos refers to the column number of the start postition and end_pos refers to the column number of the end position
def readFile(f, dict, chr_pos, start_pos, end_pos):
    f = open(f, "r")
    next(f)
    cnt = 0;
    for line in f:
        tokens = line.split('\t')
        key = tokens[chr_pos]
        if "chr" in key:
            key = key[3:]
        avg = average(tokens[start_pos], tokens[end_pos])
        #print(tokens[chr_pos] , " " , tokens[start_pos], " " , tokens[end_pos], avg)      

        dict.setdefault(key,[]).append(avg)
        cnt = cnt + 1


def get(dict, key):
    return dict[key]

    
def average(number1, number2):
  return int((int(number1) + int(number2)) / 2)

# read the motif file and create a dict
motif_dict = {}
#readFile("../data/hg19.motifs.txt", motif_dict, 1, 2, 3);
readFile("hg19.motifs.txt", motif_dict, 1, 2, 3);

# read the loop anchor file and create a dict
loop_anchor_dict = {}
#readFile("../data/GM12878_loops.txt", loop_anchor_dict, 0, 1, 2);
readFile("GM12878_loops.txt", loop_anchor_dict, 0, 1, 2);

# match loop anchors to motifs
match_loop_motif(motif_dict, loop_anchor_dict)
