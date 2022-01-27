# given the ticket # and the sample lists. 
# run fastqValidator to check the fastq files are valid. 
# -- written by Hongzhu Cui (IGM)  

import optparse
from os import listdir
from os.path import isfile, join, getsize
import sys
import subprocess

# TO-DO
# remove the ticket as argument & check samples with different file in differnt folders. 

# a function list all the fastq file for a sample. returns a list of fastqs file names. 
def list_sample_fastqs(ticket, sample):
    path = join("/nfs/tx/in", ticket, sample, "FASTQ") 
    fastqs = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith(".fastq.gz")]
    assert len(set(fastqs)) == len(fastqs) # check duplicates
    return fastqs 

# run fastq validator and return True/False based on the system call result. 
def run_validator(fastq, log_f):
    call_return = subprocess.call("/nfs/goldstein/software/fastQValidator_0.1.1a-x86_64/fastQValidator/bin/fastQValidator --file {0} --maxErrors 1".format(fastq), 
        shell=True, stdout= log_f, stderr= log_f)

#     print (call_return) 
#     sys.exit(0)
    if call_return > 0 :
        return False
    if call_return == 0: 
        return True 
    return False

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

    samples_with_fastq_missing_or_unpaired = []
    samples_with_invalid_fastqs = []

    log_file_name = ticket + "-" + sample_file_name + "-" + "VALIDATION_REPORT.txt"
    log_f = open (log_file_name, 'w')

    for s in samples:
        # log_f.write('\n' + 'CHECKING: ' +  s + '\n')
        fastqs = list_sample_fastqs(ticket, s)
        if len(fastqs) == 0 or len(fastqs) % 2 == 1 :
            # "no fastq files found!" "fastqs are not paired!"
            samples_with_fastq_missing_or_unpaired.append(s)
            continue
			
        for f in fastqs: 
            full_path_fastq = join("/nfs/tx/in", ticket, s, "FASTQ", f) 

            if run_validator(full_path_fastq, log_f):
                continue 
            else: 
                samples_with_invalid_fastqs.append(s)
		# log_f.write('\n' + 'FINISHED CHECKING: ' +  s + '\n')
    log_f.close()
    
    # valid_samples = set(samples) - set(samples_with_fastq_missing_or_unpaired.extend(samples_with_fastq_not_valid))
    valid_samples = set(samples) - set(samples_with_fastq_missing_or_unpaired) - set(samples_with_invalid_fastqs)
    if len(samples_with_fastq_missing_or_unpaired) > 0 or len(samples_with_invalid_fastqs) > 0: 
    	print (">>> WARNING: SOME SAMPLES IN THE SAMPLE LIST CONTAINS INVALID FASTQS! <<< ")
        # print ("IT MIGHT CAUSE THE PIPELIEN CRASHING! ")
        if len(samples_with_fastq_missing_or_unpaired) > 0:
            print ("1. it contains samples with missing fastqs or fastq unpaired!")
            outfile1 = ticket + sample_file_name + "_fastq_missing_or_unpaired.txt"
            with open(ticket +"/" + outfile1, 'w') as f: 
                f.write('\n'.join(samples_with_fastq_missing_or_unpaired))
            print ("samples with missing fastqs or fastq unpaired have been written to {0}".format(outfile1))
        
        if len(samples_with_invalid_fastqs) > 0:
            print ("2. it contains samples with invalid fastq files! ")
            outfile2 = ticket + sample_file_name + "_invalid_fastqs.txt"
            with open(ticket +"/" + outfile2, 'w') as f: 
                f.write('\n'.join(samples_with_invalid_fastqs))
            print ("samples with invalid fastqs have been written to {0}".format(outfile2))
        
        """
        TO-DO
        print (samples have differnt fastqs in under different ticket! )
        please rename sample folder in the ticket THAT ARE NOT THE ONE YOU INTEND TO RELEASE! 
        """

        print ("\n")
        outfile_valid = ticket + sample_file_name + "_VALID.txt"
        with open (ticket +"/" + outfile_valid, 'w' ) as f: 
            f.write("\n".join(valid_samples))
        print ("SAMPLES WITH VALID FASTQ FILES ARE WRITTEN TO {0} under".format(outfile_valid))

        print ("DETAILED INFO ABOUT FASTQ VALIDATION CAN BE FOUND AT ")

    else: 
        print ("Congrats! SAMPLES MEETS THE QUALITY CRITERIA OF FASTQ VALIDATOR! IT'S READY TO GO!")
	
	
    # print ("REMINDER:")
    # PLEASE REMOVE THE OUTPUT FILES AND VALIDATION LOG FILES AFTER YOU ARE DONE WITH THIS TICKET! 
    # KEEP THE /nfs/tx/in FOLDER CLEAN! 

