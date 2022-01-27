# given the ticket # and the fastq file lists. 
# check whether the fastqs need to be merged. 
# ASSUMMING THE THE DIRECTORY STRUCTURE AND FASTQ FILE NAMES FOR RELEASING ARE PROPERLY SET UP. 
# 1. can not have lane number larger than 10. 
# 2. for one lane it cannot have multiple flowcell. 
# it also report the space required for the merging fastqs, as we keep individual fastqs as backup
# -- written by Hongzhu Cui (IGM)  

import optparse
from os import listdir
from os.path import isfile, join, getsize
import sys

# import enum
# Enum for size units
# class SIZE_UNIT(enum.Enum):

BYTES = 1
KB = 2
MB = 3
GB = 4
TB = 5

def convert_unit(size_in_bytes, unit = TB):
    #  Convert the size from bytes to other units like KB, MB or GB
    if unit == KB:
        return size_in_bytes/1024
    elif unit == MB:
        return size_in_bytes/(1024*1024)
    elif unit == GB:
        return size_in_bytes/(1024*1024*1024)
    elif unit == TB:
        return size_in_bytes/(1024*1024*1024*1024)
    else:
        return size_in_bytes

# a function list all the fastq file for a sample. returns a list of fastqs file names. 
def list_sample_fastqs(ticket, sample):
    path = join("/nfs/tx/in", ticket, sample, "FASTQ") 
    fastqs = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith(".fastq.gz")]
    assert len(set(fastqs)) == len(fastqs) # check duplicates
    return fastqs 

# further divide the fastqs by lane number, returns a dictionary of {("L001", [f1, f2, f3, f4, ...])
def divide_fastqs_by_lanes(fastqs):
    lane_num = 1
    lane2fastqs = {}
    fastq_set = set(fastqs)
    while len(fastq_set) != 0:
        lane_str = "_L00" + str(lane_num)
        lfastqs = [f for f in fastq_set if lane_str in f]
        assert len(lfastqs) % 2 == 0 # check correct fastq number
        fastq_set.difference_update(lfastqs) # remove from fastq_set
        lane2fastqs[lane_num] = lfastqs
        lane_num += 1
        if lane_num > 50: 
            break
    return lane2fastqs

# check is the sample are sequenced on multiple flowcells.
def contain_multiple_fcs(lane2fastqs):
    for lane, fastqs in lane2fastqs.iteritems(): 
        if len(fastqs) > 2: 
            return True
    return False

# check the lane number issues. 
def lane_num_larger_than_10(lane2fastqs):
    lane_nums = [ lane_num for lane_num in lane2fastqs.keys() ]
    # return len(lane2fastqs) > 10
    return max(lane_nums) >= 10

def size_of_fastqs(ticket, sample, fastqs): 
    full_paths = [join("/nfs/tx/in", ticket, sample, "FASTQ", f) for f in fastqs]
    sizes = [getsize(p) for p in full_paths]
    return sum(sizes)

if __name__ == "__main__": 
    parser = optparse.OptionParser()
    # parser.add_option('-q', '--query',
    #     action="store", dest="query",
    #     help="query string", default="spam")

    options, args = parser.parse_args()
    ticket = args[0]
    sample_file = args[1]
    sample_file_name = sample_file.split('/')[-1]
    with open(sample_file) as f:
        samples = [line.strip() for line in f if not line.isspace() ]

    required_space = 0.0 
    problematic_samples = []
    for s in samples:
        fastqs = list_sample_fastqs(ticket, s)
        lane2fastqs = divide_fastqs_by_lanes(fastqs)

        if contain_multiple_fcs(lane2fastqs) or lane_num_larger_than_10(lane2fastqs):
            problematic_samples.append(s)
            byte_size = size_of_fastqs(ticket, s, fastqs)
            required_space += convert_unit(float(byte_size))

    outfile = join("/nfs/tx/in", ticket, ticket + "_" + sample_file_name +  "_merged_samples.txt")
    with open (outfile, 'w') as f: 
        f.write("\n".join(problematic_samples))
        f.write("\n")

    if problematic_samples: 
    	print (">>> WARNING: SAMPLES ARE NOT SAFE TO RELEASE! <<< ")
        print ("{0} out of {1} samples under ticket {2} have some issues.".format(len(problematic_samples), len(samples), ticket))
        print ("their fastqs need to be merged before lagacy release")
        print ("the problematic samples have been written to file: {0}".format(outfile))
        print ("please manually verify the fastq files to confirm the issue, then")
        print ("1. make sure you have enough free space in /nfs/tx/in, as it requires {0} TB for merging fastqs".format(required_space))
        print ("2. run the command below before extereal release: ")
        print ("	python fastq_merge.py {0} {1}".format(ticket, outfile))
    else: 
        print ("Congrats! The samples are OK and can be released without merging!")
