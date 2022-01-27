# given the ticket # and the fastq file lists. 
# check whether the fastqs files follows the NAMING CONVENTION ADN DIR SET UP INSTRUCTIONS WHEN RELEASING EXTERNAL SAMPLES. 

# 1. SAMPLE_NAME/FASTQ dir STRUCTURE 
# 2. PAIRED FASTQS.
# 3. NAMEING CONVENTION: _L001_R1_1.fastq.gz
# 4. PROPER FILE PERMISSION (READABLE) 
# 5. ENOUGH DRAGEN QUOTA? (To-Do)

# -- written by Hongzhu Cui (IGM)  

import optparse
from os import listdir
from os.path import isfile, join, getsize
import os
import sys
import re
import subprocess

RE_FASTQ_NAME = re.compile('_L00[1-9][0-9]?_(R1|R2)_(\d+)\\.fastq\\.gz$') 

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

def check_dragen_quota(): 
	HOST="dragen1"
	# Ports are handled in ~/.ssh/config since we use OpenSSH
	COMMAND="dragen_lic"

	ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                       shell=False,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
	result = ssh.stdout.readlines()

# 	if result == []:
#    		error = ssh.stderr.readlines()
#     	print >>sys.stderr, "ERROR: %s" % error
# 	else:
#     	print result        
# a function to check whether the sample have proper directory setup   
	return results
	      
def have_proper_dir_setup(ticket, s):
	path = join("/nfs/tx/in", ticket, s, "FASTQ") 
	if os.path.isdir(path):
	    return True
	return False
	
def have_proper_naming(fq_name): 
	if RE_FASTQ_NAME.search(fq_name):
	    return True
	return False

# a function list all the fastq file for a sample. returns a list of fastqs file names. 
def list_sample_fastqs(ticket, sample):
    path = join("/nfs/tx/in", ticket, sample, "FASTQ") 
    fastqs = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith(".fastq.gz")]
    # assert len(set(fastqs)) == len(fastqs) # check duplicates
    return fastqs 


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
    # sample_type = args[2] might need it to estimate the dragen quota
    
    sample_file_name = sample_file.split('/')[-1]
    with open(sample_file) as f:
        samples = [line.strip() for line in f if not line.isspace() ]
	
	# To-Do
    required_space = 0.0 
    
    samples_wihtout_proper_setup = []
    samples_with_bad_fastq_naming = []
    samples_with_unpaired_fastqs = []
    samples_without_read_permissions = []
    
    for s in samples:
    	# 1. check we have proper dir set-up for sample s
    	if not have_proper_dir_setup(ticket, s):
			samples_wihtout_proper_setup.append(s)
			continue
    	
        fastqs = list_sample_fastqs(ticket, s)

        if len(fastqs) == 0: 
        	samples_wihtout_proper_setup.append(s)
        	continue

        # check fastq naming conventions 
        proper_fqs = filter (have_proper_naming, fastqs) 
        if len(proper_fqs) != len(fastqs): 
        	samples_with_bad_fastq_naming.append(s)
        	continue
       	
        # check whether we have unpaired fastqs 
        if len(fastqs) % 2 != 0:
        	samples_with_unpaired_fastqs.append(s)
        	continue

        fq_set = set(fastqs) 

        unpair_flag = False
        for fq in fq_set: 
			if fq.find('_R1_') > -1: 
				pair_fq = fq.replace('_R1_', '_R2_', 1)
				if pair_fq not in fq_set: 
					unpair_flag = True
					break
			elif fq.find('_R2_') > -1: 
				pair_fq = fq.replace('_R2_', '_R1_', 1)
				if pair_fq not in fq_set: 
					unpair_flag = True
					break
			else: 
				raise Exception("something is wrong")  
        if unpair_flag: 
			samples_with_unpaired_fastqs.append(s)
			continue
				
        # check whether the all fastqs have proper files permissions. 
        have_not_readable = False
        for fq in fastqs: 
        	full_path_fq = join("/nfs/tx/in", ticket, s, "FASTQ", fq)
        	if not os.access(full_path_fq, os.R_OK): 
        		have_not_readable = True
        		break
        if have_not_readable:
        	samples_without_read_permissions.append(s)
        	
        # calculate the size of fastqs and check whether we have enough quota.  
        # byte_size = size_of_fastqs(ticket, s, fastqs)
        # required_space += convert_unit(float(byte_size))
        
        # TO-DO
        # estimate = estimate_dragen(required_space)		
		# result = check_dragen_quota()
		
		# TO-DO
# 	    if result == []:
#    		    # error = ssh.stderr.readlines()
#     	    print ("Failed to fetch the dragen quota information. Please check! ")
#         
#         quota = parse_dragen(result)

	

    good_samples = list( set(samples) -  set (samples_wihtout_proper_setup) - set(samples_with_bad_fastq_naming ) -set(samples_with_unpaired_fastqs) -  set( samples_without_read_permissions ))
    
    if len(good_samples) == len(samples): 
    	print ("Congrats! The samples are properly set up, and can proceed to check whether they need to be merged!")
    else: 
        print (">>> WARNING: SAMPLES ARE NOT SAFE TO RELEASE! <<< ")
    	outfile = join("/nfs/tx/in", ticket, ticket + "_" + sample_file_name +  "_qualified_samples.txt")
    	with open (outfile, 'w') as f: 
        	f.write("\n".join(good_samples))
        	f.write("\n")

    	if samples_wihtout_proper_setup: 
        	print ("{0} out of {1} samples under ticket {2} dont have proper direcotry structure setup.".format(len(samples_wihtout_proper_setup), len(samples), ticket))
        	outfile1 = join("/nfs/tx/in", ticket, ticket + "_" + sample_file_name +  "_samples_without_proper_setup.txt")
    		with open (outfile1, 'w') as f: 
        	    f.write("\n".join(samples_wihtout_proper_setup))
        	    f.write("\n")
        	print ("samples without proper setup have been written to file: {0} \n\n".format(outfile1))

    	if samples_with_bad_fastq_naming: 
        	print ("{0} out of {1} samples under ticket {2} have bad fastq file naming.".format(len(samples_with_bad_fastq_naming), len(samples), ticket))
        	outfile2 = join("/nfs/tx/in", ticket, ticket + "_" + sample_file_name +  "_samples_with_bad_naming.txt")
    		with open (outfile2, 'w') as f: 
        	    f.write("\n".join(samples_with_bad_fastq_naming))
        	    f.write("\n")
        	print ("samples with bad fastq naming have been written to file: {0} \n\n".format(outfile2))

    	if samples_with_unpaired_fastqs: 
        	print ("{0} out of {1} samples under ticket {2} have unpaired fastqs.".format(len(samples_with_unpaired_fastqs), len(samples), ticket))
        	outfile3 = join("/nfs/tx/in", ticket, ticket + "_" + sample_file_name +  "_samples_with_unpaired_fastqs.txt")
    		with open (outfile3, 'w') as f: 
        	    f.write("\n".join(samples_with_unpaired_fastqs))
        	    f.write("\n")
        	print ("samples with unpaired_fastqs have been written to file: {0} \n\n".format(outfile3))

    	if samples_without_read_permissions: 
        	print ("{0} out of {1} samples under ticket {2} dont have proper permissions.".format(len(samples_without_read_permissions), len(samples), ticket))
        	outfile4 = join("/nfs/tx/in", ticket, ticket + "_" + sample_file_name +  "_samples_without_permissions.txt")
    		with open (outfile4, 'w') as f: 
        	    f.write("\n".join(samples_without_read_permissions))
        	    f.write("\n")
        	print ("samples without permissions have been written to file: {0} \n\n".format(outfile4))

        
