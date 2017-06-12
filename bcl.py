# bcl.py
# Joshua Bridgers
# 04/26/2012
# jb371@dm.duke.edu
#
# Submits finished sequenced run for BCL to fastq conversion

import csv
import glob
import logging
import os
import re
import getopt
import sys
from commands import getoutput
from CHGV_mysql import setup_logging
from CHGV_mysql import getSequenceDB
from CHGV_mysql import getTestSequenceDB
from CHGV_mysql import getUserID
from CHGV_mysql import MachineCheck
from create_sss import create_sss_bcl_2
from sss_check import sss_qc

def bcl(info,Machine,runPath,BCLDrive,seqsata,machine,sequenceDB):
    logger = logging.getLogger('bcl')
    in_dir = runPath
    script_dir = '/home/jb3816/github/sequencing_pipe'
    pwd = os.getcwd()

    #Name of folder contains flow cell info, machine info etc.  Parsing
    FCID = info[3]
    HiSeq = info[1]
    Date = info[0][3:7]
    Date_Long = info[0]

    out_dir = '{}/{}_{}_{}_Unaligned'.format(BCLDrive,Date_Long,HiSeq,FCID)
    scriptLoc = out_dir + '/' + HiSeq + '_' + Date_Long + '_BCL.sh'
    dir_check(BCLDrive + out_dir)
    base_script = '/nfs/goldstein/software/bcl2fastq2_v2.19.0/bin/bcl2fastq '

    if Machine[0] == 'A': #NovaSeq
        base_script += '--runfolder-dir %s --output-dir %s ' % (in_dir,out_dir)
    else: #HiSeq
        base_script += '--runfolder-dir %s --output-dir %s --barcode-mismatches 1 ' % (in_dir,out_dir)

    # Use if bcl or stats are missing.  They should never be missing unless 
    # there was a data transfer problem or corruption.
    #base_script += '--ignore-missing-bcl --ignore-missing-stats '

    #Depreciated.  Do not use.  If you need to re-run bcl due to failure, delete the folder and re-run
    #if forceBCL == True:
    #    base_script += ' --force'

    if sampleSheet:
        base_script += '--sample-sheet '+sampleSheet +' '
    else:
        base_script += '--sample-sheet %s/*%s*.csv ' % (runPath,FCID)

    # Tile specify what tiles on the flowcell can be converted to fastq
    # This adds tiles parameter to BCL script if specified
    if tiles != False:

        sql = ("UPDATE Flowcell "
            "SET TileCommand=CONCAT_WS(';',TileCommand,'{0}') "
            "WHERE FCillumID='{1}'"
            ).format(tiles,FCID)

        if verbose == True:
            print sql
        sequenceDB.execute(sql)

        base_script = base_script + '--tiles='+tiles+' '
        os.system('echo %s > tiles.txt' % tiles)

    # Used to mask any bases on the read
    if base_mask != False:
        base_script = base_script + '--use-bases-mask '+base_mask+' '
    #print base_script

    logger.info(base_script)
    os.mkdir(out_dir)

    #Submit bcl job to the cluster
    os.system('cp {}/SGE_header {}'.format(script_dir,scriptLoc))
    bcl_script = open(scriptLoc,'a')
    bcl_script.write(base_script + '\n')
    bcl_script.close()
    qsubLoc = '/opt/sge6_2u5/bin/lx24-amd64/qsub'
    cmd = 'cd {0} ; {1} -cwd -v PATH -N {2}_{3}_{4}_bcl {5}'.format(out_dir,qsubLoc,machine,FCID,BCLDrive.split('/')[2],scriptLoc)
    if verbose == True:
        print cmd

    logger.info(cmd)
    status = os.system(cmd)
    logger.info(status)

#checks if a bcl directory already exists
def dir_check(BCLDrive):
    logger = logging.getLogger('dir_check')
    dir_path = glob.glob(BCLDrive)
    if dir_path != []:
        logger.warn('BCL directory already exists! %s' % dir_path)
        raise Exception, 'BCL directory already exists! %s' % dir_path

def checkSataLoc(sata_loc):
	logger = logging.getLogger('checkSataLoc')
	if os.path.isdir(sata_loc) == False:
		logger.warn('Path for seqsata drive, %s is incorrect!' % sata_loc)
		raise Exception, 'Path for seqsata drive, %s is incorrect!' % sata_loc


def usage():
    print '-b, --output\t\tPath to output folder'
    print '-h, --help\t\tShows this help message and exit'
    print '-n, --noSSS\t\tForbids creation of a sample sheet'
    print '-s, --sampleSheet\tAllows user-specified sample sheet'
    print '-t, --testdb\tAll database updates occur on the test database'
    print '-v, --verbose\t\tVerbose output'
    print '--noStatus\t\tDoes not update the status of the samples on the flowcell'
    print '--tiles\t\t\tSpecifies what tiles you want BCL performed on.  Ex: s_[12]_1[123]0[1-7]'
    print '--use-bases-mask\tSpecifies what bases you want to mask.  Ex: y100n,I6n,y50n'
    sys.exit(2)

def check_sss(FCID):
    if sampleSheet:
        sample_sheets = [sampleSheet]
    else:
	sample_sheets = glob.glob('/nfs/genotyping/Sequencing_SampleSheets/*%s*.csv' % FCID)
	if len(sample_sheets) > 1:
		raise Exception, 'Multiple Sample Sheets exist with the same FCID'
	else:
		sss_qc(FCID)

def updateSamples(sequenceDB,FCID):
    logger = logging.getLogger('updateSamples')
    userID = getUserID()
    sql = ("INSERT INTO statusT "
        "(CHGVID,status_time,status,DBID,prepID,userID) "
        "SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'BCL',pt.DBID,pt.prepID,{0} "
        "FROM Flowcell f "
        "JOIN Lane l ON l.FCID=f.FCID "
        "JOIN prepT pt ON pt.prepID=l.prepID "
        "WHERE FCillumid='{1}'"
        ).format(userID,FCID)


    if verbose == True:
        print sql
    sequenceDB.execute(sql)
    logger.info(sql)

def opts(argv):
    global tiles
    tiles = False
    global base_mask
    base_mask = False
    global sata_loc
    sata_loc = ''
    global verbose
    verbose = False
    global sampleSheet
    sampleSheet = ''
    global noSSS
    noSSS = False
    global forceBCL
    forceBCL = False
    global noStatus
    noStatus = False
    global runPath
    global BCLDrive
    global testdb
    testdb = False
    try:
        opts,args = getopt.getopt(argv, "fhi:nb:vs:",
            ['input=','help','force','bcl=','tiles=','use-bases-mask=','verbose','sampleSheet=','noSSS','noStatus','testdb'])
    except getopt.GetoptError, err:
        print err
        usage()
    for o,a in opts:
        if o in ('-f','--force'):
            forceBCL = True
        elif o in ('-h','--help'):
            usage()
        elif o in ('-i','--input'):
            runPath = '/nfs/seqscratch1/Runs/' + a
        elif o in ('-n','--noSSS'):
            noSSS = True
        elif o in ('--noStatus'):
            noStatus = True
        elif o in ('-b','--bcl'):
            BCLDrive = a
        elif o in ('-s','--sampleSheet'):
            sampleSheet = a
        elif o in ('--tiles'):
            tiles = a
        elif o in ('--testdb'):
            testdb = True
        elif o in ('--use-bases-mask'):
            base_mask = a
        elif o in ('-v','--verbose'):
            verbose = True
        else:
            assert False, "Unhandled argument present"

def RTA_check(runPath):
    logger = logging.getLogger('RTA_check')
    if os.path.isfile('%s/RTAComplete.txt' % runPath) == False:
        logger.warn("RTA has not completed!")
        raise Exception, "RTA has not completed!"
    else:
        logger.info('RTA has already completed')
        print "RTA has already completed"

def main():
    opts(sys.argv[1:])
    if testdb:
        sequenceDB = getTestSequenceDB()
    else:
        sequenceDB = getSequenceDB()
    info = runPath.split('/')[4].split('_')
    Date = info[0]
    FCID = info[3]
    Machine = MachineCheck(sequenceDB,info[1],FCID)

    seqsata_drive = 'igmdata01'

    setup_logging(Machine,FCID,seqsata_drive)
    logger = logging.getLogger('main')
    logger.debug('Initializing Parameters: runPath:%s, FCID:%s, Machine:%s, seqsata_drive:%s, tiles:%s, base_mask:%s',
            (runPath,FCID,Machine,seqsata_drive,tiles,base_mask))
    print "Starting BCL job creation..."
    logger.info("Starting BCL job creation...")

    RTA_check(runPath)
    if noSSS == False:
        create_sss_bcl_2(runPath,FCID,Machine,Date,sequenceDB)
    #check_sss(FCID)
    bcl(info,Machine,runPath,BCLDrive,seqsata_drive,Machine,sequenceDB)
    if noStatus == False:
        updateSamples(sequenceDB,FCID)

    logger.info("BCL successfully started")
    print "BCL successfully started"
    sequenceDB.execute('COMMIT;')
    sequenceDB.close()

if __name__ == '__main__':
	main()
