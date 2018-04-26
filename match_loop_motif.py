from numpy import *

def match_loop_motif(motif_dict, loop_anchor_dict):
    for i in range(1,23):
        row_vector_motif = array([motif_dict[str(i)]])
        col_vector_loop_anchor = array([loop_anchor_dict[str(i)]]).T
        print(i, len(row_vector_motif) , row_vector_motif.shape)
        print(i, len(col_vector_loop_anchor), col_vector_loop_anchor.shape)

        dists = absolute(col_vector_loop_anchor - row_vector_motif)
        print(dists)

        close_ones = dists < 20000
        print(close_ones)

        print(where(close_ones))

        break
	
	
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

        #if cnt == 1000 :
        #    break
    


def get(dict, key):
    return dict[key]

    
def average(number1, number2):
  return int((int(number1) + int(number2)) / 2)

motif_dict = {}
readFile("hg19.motifs.txt", motif_dict, 1, 2, 3);

loop_anchor_dict = {}
readFile("GM12878_loops.txt", loop_anchor_dict, 0, 1, 2);

match_loop_motif(motif_dict, loop_anchor_dict)
