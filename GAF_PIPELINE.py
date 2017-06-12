#from CHGV_mysql import getSequenceDB

import getopt
import glob
import logging
import os
import sys
import subprocess
import time
from CHGV_mysql_3_6 import setup_logging
from CHGV_mysql_3_6 import getSequenceDB
from CHGV_mysql_3_6 import getTestSequenceDB
from collections import Counter

def getBestBCLDrive():
    runningSeqscratchList =[]
    bcl_jobs = getoutput('qstat -r | grep "Full jobname"').split()
    for jobs in bcl_jobs:
        if 'bcl' in jobs and 'XX' in jobs:

            if _debug == True:
                print(jobs,'bcl jobs')

            bcl_seqsata = jobs.split('_')[2]
            runningSeqscratchList.append(bcl_seqsata)
    seqscratchDrives = ['seqscratch_ssd']

    availSeqscratchDrives = set(seqscratchDrives) - set(runningSeqscratchList)

    if availSeqscratchDrives == set([]):
        seqscratchCount = Counter(runningSeqscratchList)
        #print(seqscratchCount)
        return '/nfs/%s/BCL' % (min(seqscratchCount,key=seqscratchCount.get))
    else:
        return '/nfs/' + list(availSeqscratchDrives)[0] + '/BCL'

def submit(runFolder,seqsata,run_date,machine,FCID,BCLDrive):
    logger = logging.getLogger('submit')
    address = 'jb3816@cumc.columbia.edu'
    pythonProgram = '/nfs/goldstein/software/python2.7/bin/python2.7'
    pythonProgram36 = '/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python'
    scriptLoc = '/nfs/goldstein/software/sequencing_pipe/production/'

    BCLCmd =         ("{0} {1}/bcl.py -i {2} -b {3}").format(pythonProgram,scriptLoc,runFolder,BCLDrive)
    BCLMySQLCmd =    ("{0} {1}/bcl_mysql.py -i {2}").format(pythonProgram36,scriptLoc,runFolder)
    postBCLCmd =     ("{0} {1}/post_bcl.py -i {2} -s {3} -b {4}").format(pythonProgram36,scriptLoc,runFolder,seqsata,BCLDrive)
    storageCmd =     ("{0} {1}/storage.py -i {2} -s {3} -b {4}").format(pythonProgram36,scriptLoc,runFolder,seqsata,BCLDrive)
    fastqcCmd =      ("{0} {1}/fastqc.py -i {2} -s {3}").format(pythonProgram,scriptLoc,runFolder,seqsata)
    fastqcMySQLCmd = ("{0} {1}/fastqc_mysql.py -i {2} -s {3}").format(pythonProgram,scriptLoc,runFolder,seqsata)


    if run_pipeline == False:
        #prints out all the commands to run manually
        print
        print("="*35+'Scripts Commands'+"="*35)
        print(BCLCmd)
        print(BCLMySQLCmd)
        print(postBCLCmd)
        print(storageCmd)
        print(fastqcCmd)
        #fastqcMySQLCmd is now run within the fastqcCmd
        print(fastqcMySQLCmd)
        print("="*86)
        print
        print(('Log file:  /nfs/{0}/summary/GAF_PIPELINE_LOGS/{1}_{2}_{0}.log'
            ).format(seqsata,machine,FCID))

    else:
        #submits everything to the cluster
        #stage 1
        os.system(BCLCmd)
        logger.info(BCLCmd)

        os.system(BCLMySQLCmd)
        logger.info(BCLMySQLCmd)

        #stage 2
        stage2(seqsata,runFolder,machine,FCID,address,postBCLCmd,storageCmd,fastqcCmd,fastqcMySQLCmd)
        qsubCmd =('/opt/sge6_2u5/bin/lx24-amd64/qsub -N {0}_{1}_stage2 -hold_jid {0}_{1}_{4}_bcl '
            '/nfs/seqscratch1/Runs/{3}/{0}_{1}_{2}_stage2.sh'
            ).format(machine,FCID,seqsata,runFolder,BCLDrive.split('/')[2])

        os.system(qsubCmd)
        logger.info(qsubCmd)
        print(qsubCmd)

def header(seqsata,file,FCID):
    file.write('#! /bin/bash\n')
    file.write('#\n')
    file.write('#$ -S /bin/bash -cwd\n')
    file.write('#$ -o /nfs/%s/fastqc/%s_%s_fastqc_complete.sge\n' % (seqsata,FCID,seqsata))
    file.write('#$ -e /nfs/%s/fastqc/%s_%s_fastqc_out.sge\n' % (seqsata,FCID,seqsata))
    file.write('#$ -V\n')
    file.write('#$ -M jb3816@cumc.columbia.edu\n')
    file.write('#$ -m bea\n')
    file.write('#\n')
    file.write('\n')

def stage2(seqsata,runFolder,machine,FCID,address,postBCLCmd,storageCmd,fastqCmd,fastqMySQLCmd):
    #writes out the stage 2 script
    logger = logging.getLogger('stage2')

    stage2_script = open('/nfs/seqscratch1/Runs/%s/%s_%s_%s_stage2.sh' % (runFolder,machine,FCID,seqsata),'w')
    header(seqsata,stage2_script,FCID)
    stage2_script.write("export LD_LIBRARY_PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/lib:$LD_LIBRARY_PATH ; export PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/bin:$PATH \n")
    stage2_script.write(postBCLCmd + '\n')
    stage2_script.write('if [ $? -eq 0 ] ; then echo "Post BCL completed successfully" ; else echo "Post BCL failed"\n')
    stage2_script.write('/nfs/goldstein/software/mutt-1.5.23 -s "Post_BCL failure: %s %s" %s < /dev/null; exit 1; fi\n' % (machine,FCID,address))
    stage2_script.write(storageCmd + '\n')
    stage2_script.write('cd /nfs/seqscratch1/Runs/%s \n' % (runFolder) )
    stage2_script.write(fastqCmd + '\n')
    stage2_script.close()

def checkSeqsata(seqsata):
    try:
        os.path.isdir(glob.glob('/nfs/%s' % seqsata)[0])
    except:
        raise Exception('Input location does not exist!')


def usage():
    print('-b, --bcl\t\tSpecify BCL drive.  Ex: seqscratch09')
    print('-d, --debug\t\tShows additional informaion such as Pipeline Executing and bcl jobs')
    print('-h, --help\t\tShows this help message and exit')
    print('-r, --run\t\tSubmits pipeline to the cluster')
    print('-s, --seqsata\t\tSpecify seqsata drive.  Ex: seqsata02')
    print('--FCID\t\t\tSpecify the flowcell Illumina ID')

    sys.exit(2)

def opts(argv):
    global run_pipeline
    run_pipeline = False
    global test
    test = True
    global _debug
    _debug = False
    global seqsata
    seqsata = 'igmdata01'
    global BCLDrive
    BCLDrive = ''

    global runPath
    runPath = ''

    try:
        opts,args = getopt.getopt(argv, "b:drs:h", ['bcl=','debug','FCID=','seqsata=','help','run'])
    except getopt.GetoptError:
        usage()
    for o,a in opts:
        if o in ('-h','--help'):
            usage()
        elif o in ('-b','--bcl'):
            BCLDrive = '/nfs/' + a + '/BCL'
        elif o in ('--FCID'):
            runPath = getRunPath(a)
        elif o in ('-r','--run'):
            run_pipeline = True
        elif o in ('-d','--debug'):
            _debug = True
        elif o in ('-s','--seqsata'):
            seqsata = a
        else:
            assert False, "Unhandled argument present"

def getRunPath(FCID):
    '''Grabs the run folder with the highest run number for edge cases where
       there the flowcell was rerun and there are multiple run folders.  This
       assumes the duplicate run folders have the name except for the run number
       within run folder name
    '''
    runPaths = max(glob.glob('/nfs/seqscratch1/Runs/*%s*' % FCID),key=os.path.getctime)
    runPath = runPaths
    return runPath

def check_cwd(runPath):
    print(runPath)
    if 'seqscratch' not in runPath or 'Runs' not in runPath or 'XX' not in runPath:
        raise Exception('The CWD is not within a run folder of a seqscratch drive!')

def RTA_check(runPath):
    logger = logging.getLogger('RTA_check')
    if os.path.isfile('%s/RTAComplete.txt' % runPath) == False:
        logger.warn("RTA has not completed!")
        raise Exception("RTA has not completed!")
    else:
        logger.info('RTA has already completed')
        print("RTA has already completed")

def completeCheck(runFolder):
	if os.path.isfile("/nfs/seqscratch1/Runs/%s/StorageComplete.txt" % runFolder) == True:
		raise Exception("Pipeline already completed!")
def Machine_check(sequenceDB,FCID,machine):
    query = "SELECT Machine,Complete FROM Flowcell where FCIllumID='%s'" % FCID
    #print(query)
    sequenceDB.execute(query)
    sequenceDB_machine = sequenceDB.fetchall()
    #print(sequenceDB_machine)
    if len(sequenceDB_machine) > 1:
        raise Exception("Too many flowcells found during Machine_check!")
    if sequenceDB_machine == ():
        raise Exception("Run folder's machine does not match SequenceDB machine listing!")
    if str(sequenceDB_machine[0]['Complete']) != '1':
        raise Exception("Flowcell has not been completed on SequenceDB!")


def main():
    opts(sys.argv[1:])

    global runPath
    if runPath == '':
        runFolder = os.getcwd()

    #print(runFolder)
    #check_cwd(runPath)
    RTA_check(runPath)

    runFolder = runPath.split('/')[4]
    info =  runFolder.split('_')
    FCID = info[3]
    """Sequencing output is in the following format:
    [Date]_[IGM Machine Name]_[HiSeq Run Number]_[Flowcell ID]_[Project Name]
    However in cases of Illumina maintence they will re-image the machine and
    reset the output folder to Illumina default:
    [Date]_[Illumina Machine Name]_[HiSeq Run Number]_[Machine Side][Flowcell ID]
    """
    if len(FCID) != 9:
        raise Exception("FCID %s is formatted incorrectly!" % FCID)

    run_date = info[0]
    machine = info[1]

    sequenceDB = getSequenceDB()

    Machine_check(sequenceDB,FCID,machine)
    global BCLDrive
    if BCLDrive == '':
        BCLDrive = '/nfs/seqscratch_ssd/BCL'
        #BCLDrive = getBestBCLDrive()

    setup_logging(machine,FCID,seqsata)
    logger = logging.getLogger(__name__)

    if run_pipeline == True:
        completeCheck(runFolder)
        logger.info('GAF_Pipeline.py in automated mode')

    submit(runFolder,seqsata,run_date,machine,FCID,BCLDrive)

    logger.info('Running GAF_Pipeline.py')
main()
