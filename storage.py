#!/usr/bin/python
# storage.py
# Joshua Bridgers
# 03/22/2017
# jb3816@cumc.columbia.edu
#
# Moves fastq.gz files to their archival location 

import argparse
import collections
import errno
import logging
import os
import re
import subprocess
import sys
import traceback
import zipfile

from datetime import datetime
from CHGV_mysql_3_6 import getSequenceDB
from CHGV_mysql_3_6 import getTestSequenceDB
from CHGV_mysql_3_6 import setup_logging
from CHGV_mysql_3_6 import getUserID
from glob import glob

def main(noStatus,test,verbose,bclDrive,seqsataLoc,inputFolder):
    #accessing mysql gaf database
    if test:
        sequenceDB = getTestSequenceDB()
    else:
        sequenceDB = getSequenceDB()

    date,machine,runNumber,FCIllumID,*rest = inputFolder.split('_')
    logger = logging.getLogger('main')
    setup_logging(machine,FCIllumID,seqsataLoc)
    logger = logging.getLogger('main')
    logger.info('Starting Storage script')

    archiveLoc = '/nfs/' + seqsataLoc
    emailAddress = ['jb3816@cumc.columbia.edu','sif2110@cumc.columbia.edu']
    #emailAddress = ['jb3816@cumc.columbia.edu']

    logger.info('inputFolder:{}, archiveLoc:{}'.format(inputFolder,archiveLoc))

    #completeCheck(bclDrive,inputFolder)
    try:
        storage(sequenceDB,machine,FCIllumID,date,verbose,bclDrive,archiveLoc,inputFolder,noStatus)
        logger.info('Starting commit')
        sequenceDB.execute('COMMIT;')
        sequenceDB.close()
        logger.info('Closing sequenceDB cursor')
        email(emailAddress,'SUCCESS',FCIllumID)
        logger.info('Done')
        if verbose:
            print('Done')

    except:
        logger.exception('Got exception')
        traceback.print_exc()
        sequenceDB.execute('ROLLBACK;')
        sequenceDB.close()
        email(emailAddress,'FAILURE',FCIllumID)
        logger.info('Storage Failure')
        print('Storage Failure')
        sys.exit(255)

def email(emailAddress,storageStatus,FCID):
    logger = logging.getLogger('email')
    logger.info('Starting email')
    emailCmd = ('echo "BCL move {} for {}"  | mail -s "IGM:BCL move {}" {}'
        ).format(storageStatus,FCID,storageStatus,' '.join(emailAddress))
    print(emailCmd)
    logger.info(emailCmd)
    os.system(emailCmd)

def storage(sequenceDB,machine,FCIllumID,date,verbose,bclDrive,archiveLoc,inputFolder,noStatus):
    logger = logging.getLogger('storage')
    UnalignedLoc = '{}/{}_{}_{}_Unaligned'.format(bclDrive,date,machine,FCIllumID)
    folderList = glob('{}/*/*fastq.gz'.format(UnalignedLoc))
    origNumFastq = len(folderList)
    origTotalFastqSize = 0

    if len(folderList) == 0:
        raise Exception('No fastqs were found!')
    #Interate over Project*/Sample*
    for fastq in folderList:
        fastqSize = os.path.getsize(fastq)
        logger.info('{} original filesize: {}'.format(fastq,fastqSize))
        origTotalFastqSize += fastqSize
        sampleName = fastq.split('/')[6].split('_')[0]
        seqtype = getSeqtype(sequenceDB,FCIllumID,sampleName)
        GAFbin = getGAFbin(sequenceDB,sampleName)
        sampleArchiveLoc = '{}/{}/{}/{}'.format(archiveLoc,seqtype,sampleName,FCIllumID)

        logger.debug('Starting transfer of {}'.format(sampleName))
        logger.debug('mkdir -p {}'.format(sampleArchiveLoc))
        mkdir_p(sampleArchiveLoc)

        if bclDrive.split('/')[2] == sampleArchiveLoc.split('/')[2]:
            fastqMoveCmd = ['mv',fastq,sampleArchiveLoc]
        else:
            fastqMoveCmd = ['rsync','-aq',fastq,sampleArchiveLoc]
        if verbose:
            print(' '.join(fastqMoveCmd))
        logger.info(fastqMoveCmd)
        subprocess.call(fastqMoveCmd)

        reportLaneBarcodeLoc = ('{}/Reports/html/{}/{}/{}/all/laneBarcode.html'
                                ).format(UnalignedLoc,FCIllumID,GAFbin,sampleName)
        reportCpCmd = ['cp',reportLaneBarcodeLoc,sampleArchiveLoc]

        if verbose:
            print(' '.join(reportCpCmd))
        logger.info(reportCpCmd)
        subprocess.call(reportCpCmd)

    mvFolderList = glob('{}/*/*/{}/*fastq.gz'.format(archiveLoc,FCIllumID))
    mvNumFastq = len(mvFolderList)
    mvTotalFastqSize = 0
    for mvFastq in mvFolderList:
        mvFastqSize = os.path.getsize(mvFastq)
        logger.info('{} moved filesize: {}'.format(mvFolderList,mvFastqSize))
        mvTotalFastqSize += mvFastqSize

    logger.info('='*80)
    logger.info('SUM of original fastq.gz file sizes = {}'.format(origTotalFastqSize))
    logger.info('Number of original fastq.gz files = {}'.format(origNumFastq))
    logger.info('='*80)
    logger.info('SUM of moved fastq.gz file sizes = {}'.format(mvTotalFastqSize))
    logger.info('Number of moved fastq.gz files = {}'.format(mvNumFastq))
    logger.info('='*80)

    print(origNumFastq,mvNumFastq)
    print(origTotalFastqSize,mvTotalFastqSize)

    if origTotalFastqSize != mvTotalFastqSize:
        raise Exception('Total sum of files sizes after move do not match!!!')
    if origNumFastq != mvNumFastq:
        raise Exception('Total number of files after move do not match!!!')
    else:
        seqscratchBase = '/nfs/seqscratch1/Runs'
        sequenceDB = getSequenceDB()
        processSAV(verbose,FCIllumID,date,machine,archiveLoc,seqscratchBase)
        processBcl2fastqLog(verbose,FCIllumID,date,machine,archiveLoc,UnalignedLoc)
        updateStatus(verbose,sequenceDB,FCIllumID,noStatus)
        updateFlowcell(verbose,sequenceDB,FCIllumID,archiveLoc)
        createStorageCompleteFlag(verbose,seqscratchBase,FCIllumID,date,machine)
        removeFastqs(verbose,folderList)

def removeFastqs(verbose,folderList):
    logger = logging.getLogger('removeFastqs')
    for fastq in folderList:
        rmCmd = ['rm',fastq]
        if verbose:
            print(' '.join(rmCmd))
        logger.info(rmCmd)
        subprocess.call(rmCmd)

def getGAFbin(sequenceDB,sampleName):
    GAFbinQuery = ("SELECT GAFbin FROM SampleT WHERE CHGVID = '{}'").format(sampleName)
    sequenceDB.execute(GAFbinQuery)
    """ It appears that the bcl2fastq software removes whitespace from the
    inputted GAFbin which is used as the project folder name hence the
    replace """
    return sequenceDB.fetchone()['GAFbin'].replace(' ','')

def getSeqtype(sequenceDB,FCIllumID,sampleName):
    logger = logging.getLogger('getSeqtype')
    seqtypeQuery = ("SELECT DISTINCT(st.seqtype) "
            "FROM Lane l "
            "JOIN Flowcell f ON l.fcid=f.fcid "
            "JOIN SeqType st ON l.prepid=st.prepid "
            "JOIN prepT p ON l.prepid=p.prepid "
            "WHERE FCIllumID='{}' AND CHGVID='{}'"
            ).format(FCIllumID,sampleName)
    logger.info(seqtypeQuery)
    #print(seqtypeQuery)
    sequenceDB.execute(seqtypeQuery)
    seqtype = sequenceDB.fetchone()['seqtype'].upper().replace(' ','_')
    return seqtype

def createStorageCompleteFlag(verbose,seqscratchBase,FCIllumID,date,machine):
    flagloc = glob('{}/{}_{}_*_{}_*/'.format(seqscratchBase,date,machine,FCIllumID))[0]
    flagloc += 'StorageComplete'
    touchCmd = ['touch',flagloc]
    if verbose:
        print(' '.join(touchCmd))
    subprocess.call(touchCmd)

def processSAV(verbose,FCIllumID,date,machine,archiveLoc,seqscratchBase):
    logger = logging.getLogger('processSAV')
    if machine[0] == 'A':
        runParametersFile = 'RunParameters.xml'
    else:
        runParametersFile = 'runParameters.xml'
    RunInfoLoc = glob('{}/Runs/{}_{}_*_{}*/RunInfo.xml'.format('/nfs/seqscratch1',date,machine,FCIllumID))[0]
    runParaLoc= glob('{}/Runs/{}_{}_*_{}*/{}'.format('/nfs/seqscratch1',date,machine,FCIllumID,runParametersFile))[0]
    InteropsLoc = glob('{}/Runs/{}_{}_*_{}*/InterOp'.format('/nfs/seqscratch1',date,machine,FCIllumID))[0]
    SAVLoc = '{}/summary/SAV/{}_{}_{}_SAV.tar.gz'.format(archiveLoc,FCIllumID,date,machine)
    tarCmd = ['tar','czf',SAVLoc,RunInfoLoc,runParaLoc,InteropsLoc]

    if verbose:
        print(' '.join(tarCmd))
    logger.info(tarCmd)
    returnCode = subprocess.call(tarCmd)
    logger.info(returnCode)

def processBcl2fastqLog(verbose,FCIllumID,date,machine,archiveLoc,UnalignedLoc):
    logger = logging.getLogger('processBcl2fastqLog')
    bcl2fastqLogLoc = '{}/nohup.sge'.format(UnalignedLoc)
    zipLoc = '{}/summary/bcl_nohup/{}_{}_{}_bcl2fastqLog.zip'.format(archiveLoc,FCIllumID,date,machine)
    zipCmd = ['zip',zipLoc,bcl2fastqLogLoc]

    if verbose:
        print(' '.join(zipCmd))
    logger.info(zipCmd)
    returnCode = subprocess.call(zipCmd)
    logger.info(returnCode)

def updateStatus(verbose,sequenceDB,FCID,noStatus):
    """sequenceDB Sample Update"""
    logger = logging.getLogger('updateStatus')
    userID = getUserID(sequenceDB)
    #Status update for entire flowcell
    if noStatus == False:
        query = ("INSERT INTO statusT "
                "(CHGVID,status_time,status,DBID,prepID,userID) "
                "SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'Storage',pt.DBID,pt.prepID,'{}' "
                "FROM Flowcell f "
                "JOIN Lane l on l.FCID=f.FCID "
                "JOIN prepT pt on pt.prepID=l.prepID "
                "WHERE FCillumid='{}'").format(userID,FCID)

        if verbose == True:
            print(query)
        logger.info(query)
        sequenceDB.execute(query)

def updateFlowcell(verbose,sequenceDB,FCID,archiveLoc):
    logger = logging.getLogger('updateFlowcell')
    query = ("UPDATE Flowcell "
             "SET DateStor=CURRENT_TIMESTAMP(), "
             "SeqsataLoc={} "
             "WHERE FCIllumid='{}'").format(archiveLoc,FCID)
    if verbose:
        print(query)
    logger.info(query)
    exception_count=0
    try:
        sequenceDB.execute(query)
    except():
        if exception_count <=3:
            sequenceDB.execute(query)
            exception_count += 1
        else:
            raise Exception('Could not complete Flowcell MySQL update')

def mkdir_p(path):
    try:
        print(path)
        os.makedirs(path)
        #os.chmod(path,0775)

    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def completeCheck(bclDrive,inputFolder):
    pass

if __name__ == "__main__":
    desc = ("Program that takes the fastqs from a completed bcl2fastq run and"
            "processes and moves them to the appropriate archival location")
    parser = argparse.ArgumentParser(
        description=desc, epilog=None,
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=40))
    parser.add_argument('-i',
        dest='inputFolder', help='input folder',
        required=True)
    parser.add_argument('-s',
        dest='seqsataLoc', help='Archival Location',
        required=True)
    parser.add_argument('-b',
        dest='bclDrive', help='BCL root dir',
        required=True)
    parser.add_argument('--test', action='store_true',
        help='update test SequenceDB instead')
    parser.add_argument('--noStatus', action='store_true',
        help='do not update sample statuses')
    parser.add_argument('--version', action='version',
        version='%(prog)s v1.3')
    parser.add_argument('--verbose', action='store_true',
        help='print verbose output')
    args=parser.parse_args()

    main(args.noStatus,args.test,args.verbose,args.bclDrive,args.seqsataLoc,args.inputFolder)
