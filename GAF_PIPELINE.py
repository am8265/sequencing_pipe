from CHGV_mysql import *
#from CHGV_mysql import getSequenceDB
from commands import getoutput

import getopt
import glob
import logging
import os
import sys
import subprocess
import time
from CHGV_mysql import setup_logging

def getSeqsata(sequencedb):
    seqsata_dict = {}
    for seqsata in getoutput('ls /nfs/seqsata[01]* -d').split():
        seqsata_dict[seqsata] = [os.statvfs(seqsata).f_bavail * os.statvfs(seqsata).f_frsize / 1024.0,0,0]
        #seqsata_dict[seqsata].append(getBCLjobs(seqsata))
    getPipelineExe(sequencedb,seqsata_dict)
    getBCLjobs(seqsata_dict)
    space_free = 0
    #check for seqsata drives that already have 4 or more bcl jobs running
    for key in seqsata_dict.keys():
        if seqsata_dict[key][2] > 1190000 * 1000:
            del seqsata_dict[key]

    for key in seqsata_dict.keys():
        if int(seqsata_dict[key][0]) - int(seqsata_dict[key][1]) - int(seqsata_dict[key][2]) > space_free:
            space_free = int(seqsata_dict[key][0]) - int(seqsata_dict[key][1]) - int(seqsata_dict[key][2])
            best_seqsata = key

    #print seqsata_dict
    space_free =  seqsata_dict[best_seqsata][0]

    if space_free < 600000000:
        logging.warn('No seqsata drive has greater than 600MB of space!')
        raise Exception, 'No seqsata drive has greater than 600MB of space!'
    else:
        return best_seqsata

def getPipelineExe(sequencedb,seqsata_dict):
	sequencedb.execute("SELECT CHGVID,SeqType,AlignSeqFileLoc from seqdbClone where Status = 'Pipeline executing'")
	Exe_Samples = sequencedb.fetchall()
	for sample in Exe_Samples:
		if _debug == True:
			print sample,'Pipeline executing'

		AlignSeqFileLoc = sample[2].split('/')[2]
		if sample[1] == 'Exome':
			seqsata_dict['/nfs/' + AlignSeqFileLoc][1] += 12000 * 1000
		elif sample[1] == 'Genome':
			seqsata_dict['/nfs/' + AlignSeqFileLoc][1] += 110000 * 1000
		elif sample[1] == 'RNAseq':
			seqsata_dict['/nfs/' + AlignSeqFileLoc][1] += 10000 * 1000
		elif sample[1] == 'Custom_Capture':
			seqsata_dict['/nfs/' + AlignSeqFileLoc][1] += 2000 * 1000
		elif sample[1] == 'Merged':
			seqsata_dict['/nfs/' + AlignSeqFileLoc][1] += 125000 * 1000
		else:
			logging.warn('Unknown seqtype %s for sample %s' % (sample[1],sample[0]))
			raise Exception, 'Unknown seqtype %s for sample %s' % (sample[1],sample[0])

	return seqsata_dict

def getBCLjobs(seqsata_dict):
	bcl_jobs = getoutput('qstat -r | grep "Full jobname"').split()
	for jobs in bcl_jobs:
		if 'bcl' in jobs and 'ACXX' in jobs:

			if _debug == True:
				print jobs,'bcl jobs'

			bcl_seqsata = jobs.split('_')[2]
			seqsata_dict['/nfs/' + bcl_seqsata][2] += 300000 * 1000

	return seqsata_dict


def submit(best_seqsata,bcl_drive,run_date,machine,FCID,pwd,address):
	logger = logging.getLogger('submit')
	seqsata_drive = best_seqsata.split('/')[2]

    #stage1
	os.system('python2.7 ~/github/sequencing_pipe/bcl.py --output %s' % (bcl_drive))
	logger.info('python2.7 ~/github/sequencing_pipe/bcl.py --output %s' % (bcl_drive))

	os.system('python2.7 ~/github/sequencing_pipe/bcl_mysql.py --seqsata %s' % (bcl_drive))
	logger.info('Ran python2.7 ~/github/sequencing_pipe/bcl_mysql.py --seqsata %s' % (bcl_drive))

	#stage2 wrapper
	stage2(best_seqsata,bcl_drive,run_date,machine,FCID,pwd,seqsata_drive,address)
	os.system('qsub -N %s_%s_%s_stage2 -hold_jid %s_%s_%s_bcl %s/%s_%s_%s_stage2.sh' % (machine,FCID,seqsata_drive,machine,FCID,seqsata_drive,pwd,machine,FCID,seqsata_drive))
	logger.info('Ran qsub -N %s_%s_%s_stage2 -hold_jid %s_%s_%s_bcl %s/%s_%s_%s_stage2.sh"' % (machine,FCID,seqsata_drive,machine,FCID,seqsata_drive,pwd,machine,FCID,seqsata_drive))


def header(file):
	file.write('#! /bin/bash\n')
        file.write('#\n')
        file.write('#$ -S /bin/bash -cwd\n')
        #file.write('#$ -o %s/fastqc/%s_%s_fastqc_complete.sge\n' % (seqsata,FCID,seqsata_drive))
        #file.write('#$ -e %s/fastqc/%s_%s_fastqc_out.sge\n' % (seqsata,FCID,seqsata_drive))
        file.write('#$ -V\n')
        file.write('#$ -M jb3816@cumc.columbia.edu\n')
	#file.write('#$ -M jb371@dm.duke.edu\n')
        file.write('#$ -m bea\n')
        file.write('#\n')
	file.write('\n')

def stage2(best_seqsata,bcl_drive,run_date,machine,FCID,pwd,seqsata_drive,address):
	stage2_script = open('%s_%s_%s_stage2.sh' % (machine,FCID,seqsata_drive),'w')
	header(stage2_script)

	logger = logging.getLogger('stage2')
	stage2_script.write('cd %s; python2.7 ~/github/sequencing_pipe/post_bcl.py --input %s/%s_%s_%s_Unaligned -s /nfs/%s \n' % (pwd,bcl_drive,run_date,machine,FCID,seqsata_drive))
	stage2_script.write('if [ "$(tail -1 /nfs/%s/summary/GAF_PIPELINE_LOGS/%s_%s_%s.log | grep Failure -o)" == Failure ]; then /usr/local/bin/mutt -s "Post_BCL failure: %s %s" %s < /dev/null; exit ; fi\n' % (seqsata_drive,machine,FCID,seqsata_drive,machine,FCID,address))
	stage2_script.write('sh ~/github/sequencing_pipe/storage.sh %s %s %s\n' % (FCID,bcl_drive,pwd))
	stage2_script.write('cd %s; python2.7 ~/github/sequencing_pipe/fastqc.py --input %s -f %s\n' % (pwd,bcl_drive,FCID))
	stage2_script.close()

def checkSeqsata(seqsata):
	try:
		os.path.isdir(glob.glob('/nfs/%s' % seqsata)[0])
	except:
		raise Exception, 'Seqsata input location does not exist!'


def usage():
	print '-d, --debug\t\tShows additional informaion such as Pipeline Executing and bcl jobs'
	print '-h, --help\t\tShows this help message and exit'
	print '-r, --run\t\tSubmits pipeline to the cluster'
	print '-s, --seqsata\t\tSpecify seqsata drive.  Ex. seqsata02'

	sys.exit(2)

def opts(argv):
    global run_pipeline
    run_pipeline = False
    global test
    test = True
    global _debug
    _debug = False
    global seqsata
    seqsata = 'fastq15'
    global BCL_drive
    BCL_drive = 'seqscratch09'


    try:
        opts,args = getopt.getopt(argv, "b:drs:h", ['bcl=','debug','seqsata=','help','run'])
    except getopt.GetoptError, err:
        usage()
    for o,a in opts:
        if o in ('-h','--help'):
            usage()
        elif o in ('-b','--bcl'):
            BCL_drive = a
        elif o in ('-r','--run'):
            run_pipeline = True
        elif o in ('-d','--debug'):
            _debug = True
        elif o in ('-s','--seqsata'):
            seqsata = a
        else:
            assert False, "Unhandled argument present"

def check_cwd():

        pwd = os.getcwd()
	if 'seqscratch' not in pwd or 'Runs' not in pwd or 'XX' not in pwd:
		raise Exception, 'The CWD is not within a run folder of a seqscratch drive!'
	print pwd
	return pwd

def RTA_check():
	logger = logging.getLogger('RTA_check')
	if os.path.isfile('RTAComplete.txt') == False:
		logger.warn("RTA has not completed!")
		raise Exception, "RTA has not completed!"
	else:
		logger.info('RTA has already completed')
		print "RTA has already completed"

def completeCheck(pwd):
	if os.path.isfile("%s/StorageComplete.txt" % pwd) == True:
		raise Exception, "Pipeline already completed!"

def Machine_check(sequenceDB,FCID,machine):
	sequenceDB.execute("Select Machine,Complete from Flowcell where FCillumID=%s", FCID)
	sequenceDB_machine = sequenceDB.fetchall()
	#print sequenceDB_machine
	if len(sequenceDB_machine) > 1:
		raise Exception, "Too many flowcells found during Machine_check!"
	if sequenceDB_machine == ():
		raise Exception, "Run folder's machine does not match SequenceDB machine listing!"
	if str(sequenceDB_machine[0][1]) != '1':
		raise Exception, "Flowcell has not been completed on SequenceDB!"

def print_commands(best_seqsata,bcl_drive,run_date,machine,FCID,pwd):
	print
	print "="*35+'Scripts Commands'+"="*35
	print 'python2.7 ~/github/sequencing_pipe/bcl.py --output %s' % (bcl_drive)
	print 'python2.7 ~/github/sequencing_pipe/bcl_mysql.py --seqsata %s' % (bcl_drive)
	print 'python2.7 ~/github/sequencing_pipe/post_bcl.py --input %s%s_%s_%s_Unaligned -s %s' % (bcl_drive,run_date,machine,FCID,best_seqsata)
	print 'sh ~/github/sequencing_pipe/storage.sh %s %s %s' % (FCID,bcl_drive,pwd)
	print 'python2.7 ~/github/sequencing_pipe/fastqc.py --input %s -f %s' % (bcl_drive,FCID)
	print 'python2.7 ~/github/sequencing_pipe/fastqc_mysql.py -s %s' % (bcl_drive)
	print "="*86
	print

def main():
    pwd = check_cwd()
    RTA_check()
    run_folder = pwd.split('/')[4]
    info = run_folder.split('_')
    FCID = info[3]
    run_date = info[0]
    machine = info[1]
    opts(sys.argv[1:])
    HiSeqs = {'H1A': 1,'H1B': 1,'H2A': 2,'H2B': 2,'H7A': 3,'H7B': 3,'H4A': 4,'H4B': 4,'H5A': 5,'H5B': 5,'H6A': 6,'H6B': 6,'H9A': 7,'H9B': 7,'H8A': 8,'H8B': 8}
    sequenceDB = getSequenceDB()
    address = 'jb3816@cumc.columbia.edu'

    Machine_check(sequenceDB,FCID,machine)
    best_seqsata = '/nfs/' + seqsata
    bcl_folder = '/nfs/%s/BCL/' % BCL_drive
    
    print seqsata,best_seqsata,bcl_folder
    setup_logging(machine,FCID,seqsata)
    logger = logging.getLogger(__name__)

    logger.info('Running GAF_Pipeline.py')
    if run_pipeline == True:
        completeCheck(pwd)
        logger.info('GAF_Pipeline.py in automated mode')
        submit(best_seqsata,bcl_folder,run_date,machine,FCID,pwd,address)
    else:
        print_commands(best_seqsata,bcl_folder,run_date,machine,FCID,pwd)

main()
