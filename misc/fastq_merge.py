# given the ticket # and the sample list file produced by fastq_check.py, 
# 1. create a backup folder under the sample directory. 
# 2. move original fastq files to the backup folder. 
# 3. merge original fastq files into two files: "R1" and "R2"
# 4. mv the two newly merged fastqs to the FASTQ FOLDER.  
import optparse
from os import listdir
from os.path import isfile, join, getsize
import sys
import subprocess


if __name__ == "__main__": 
    parser = optparse.OptionParser()
    options, args = parser.parse_args()
    ticket = args[0]
    sample_file = args[1]
    with open(sample_file) as f:
        samples = [line.strip() for line in f if not line.isspace() ]

    for s in samples: 
        fastq_dir = join("/nfs/tx/in", ticket, s, "FASTQ")
        backup_dir = join("/nfs/tx/in", ticket, s, "backup")
        subprocess.call("mkdir -p {0}".format(backup_dir), shell=True)
        subprocess.call("mv {0}/* {1}".format(fastq_dir, backup_dir), shell=True)

        # print("echo cat {0}/*_R1_*.fastq.gz >  {1}/{2}_L001_R1_1.fastq.gz".format(backup_dir, fastq_dir, s))
        subprocess.call("cat {0}/*_R1_*.fastq.gz >  {1}/{2}_L001_R1_1.fastq.gz".format(backup_dir, fastq_dir, s), shell=True)
        # print("echo cat {0}/*_R2_*.fastq.gz >  {1}/{2}_L001_R2_1.fastq.gz".format(backup_dir, fastq_dir, s))
        subprocess.call("cat {0}/*_R2_*.fastq.gz >  {1}/{2}_L001_R2_1.fastq.gz".format(backup_dir, fastq_dir, s), shell=True)

        # subprocess.call("mv %s %s" % (source_files, destination_folder), shell=True)
