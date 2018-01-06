#!/usr/bin/python
# storage.py
"""Moves fastq.gz files to their archival location"""

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
from glob import glob
from utilities import *

def main():
    config = get_config()
    args = parse_arguments()

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    run_info_dict = parse_run_parameters_xml(args.fcillumid,database)
    setup_logging(run_info_dict['machine'],args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    logger.info('Starting storage script')
    fcillumid=run_info_dict['FCIllumID']
    machine=run_info_dict['machine']
    unaligned_dir = '{}/{}_{}_{}_Unaligned'.format(config.get('locs','bcl2fastq_scratch_dir'),
                                             run_info_dict['runDate'],
                                             machine,fcillumid)

    logger.info('Starting Storage script')
    # Josh Bridgers, Sophia Frantz, Brett Copeland
    emailAddress = ['jb3816@cumc.columbia.edu','sif2110@cumc.columbia.edu','bc2675@cumc.columbia.edu']
    emailAddressFail = ['jb3816@cumc.columbia.edu']

    #completeCheck(bcl_drive,inputFolder)
    try:
        storage(machine,fcillumid,args,config,run_info_dict,database)
        email(emailAddress,'SUCCESS',fcillumid)
        logger.info('Done')
        if args.verbose:
            print('Done')

    except Exception:
        logger.exception('Got exception')
        traceback.print_exc()
        email(emailAddressFail,'FAILURE',fcillumid)
        logger.info('Storage Failure')
        print('Storage Failure')
        sys.exit(255)


def email(emailAddress,storageStatus,FCID):
    logger = logging.getLogger(__name__)
    logger.info('Starting email')
    emailCmd = ('echo "BCL move {} for {}"  | mail -s "IGM:BCL move {}" {}'
        ).format(storageStatus,FCID,storageStatus,' '.join(emailAddress))
    print(emailCmd)
    logger.info(emailCmd)
    os.system(emailCmd)

def regenerate_archive_tuples_list(config,fcillumid,database):
    archive_tuples_list = []
    query = """SELECT CHGVID,FCILLUMID,LANENUM,SAMPLE_TYPE
               FROM prepT p
               JOIN Lane l ON l.PREPID=p.PREPID
               JOIN Flowcell f ON l.FCID=f.FCID
               WHERE f.fcillumid='{}'
            """.format(fcillumid)
    sample_info_on_flowcell = run_query(query,database)
    for sample in sample_info_on_flowcell:
        sampleArchiveLoc = '/nfs/{}/{}/{}/{}'.format(config.get('locs','fastq_archive_drive'),
                                                     sample['SAMPLE_TYPE'].upper(),
                                                     sample['CHGVID'],fcillumid)
        for read in range(1,3):
            scratchFastq = ('/nfs/{archive_drive}/{sample_type}/{sample_name}/{fcillumid}/{sample_name}_*L00{lanenum}_R{read}_*.fastq.gz'
                           ).format(archive_drive=config.get('locs','bcl2fastq_scratch_drive'),
                                    sample_type=sample['SAMPLE_TYPE'].upper(),
                                    sample_name=sample['CHGVID'],
                                    lanenum=sample['LANENUM'],
                                    fcillumid=fcillumid,
                                    read=read)
            if glob(scratchFastq) == []:
                print(scratchFastq)
                raise Exception('Fastq not found for regenerated tuple list')
            else:
                if len(glob(scratchFastq)) > 1:
                    raise Exception('Too many fastqs found for regenerated tuple list')
                archive_tuples_list.append((glob(scratchFastq)[0],sampleArchiveLoc))
    return archive_tuples_list

def archiveFastqs(config,archive_tuples_list,fcillumid,verbose,database):
    if archive_tuples_list == []:
        archive_tuples_list = regenerate_archive_tuples_list(config,fcillumid,database)
    logger = logging.getLogger(__name__)
    origNumFastq = 0
    origTotalFastqSize = 0
    for orig_dest_pair in archive_tuples_list:
        origNumFastq += 1

        fastqSize = os.path.getsize(orig_dest_pair[0])
        logger.info('{} original filesize: {}'.format(orig_dest_pair[0],fastqSize))
        origTotalFastqSize += fastqSize

        subprocess.check_call(['mkdir','-p',orig_dest_pair[1]])
        fastqMoveCmd = ['rsync','--inplace','-aq',orig_dest_pair[0],orig_dest_pair[1]]
        if verbose:
            print(' '.join(fastqMoveCmd))
        logger.info(fastqMoveCmd)
        subprocess.check_call(fastqMoveCmd)
    if verbose:
        print('Counting archived fastqs...')
    dest_drive = config.get('locs','fastq_archive_drive')
    fastqList = glob('/nfs/{}/*/*/{}/*.fastq.gz'.format(dest_drive,fcillumid))
    mvTotalFastqSize = 0
    for archive_fastq in fastqList:
        mvTotalFastqSize += os.path.getsize(archive_fastq)
    mvNumFastq = len(fastqList)
    check_orig_dest_transfer(dest_drive,origNumFastq,origTotalFastqSize,mvNumFastq,mvTotalFastqSize)

def get_mv_fastq_size_sum(drive,fcillumid):
    logger = logging.getLogger(__name__)
    mvFastqList = glob('/nfs/{}/*/*/{}/*fastq.gz'.format(drive,fcillumid))
    mvNumFastq = len(mvFastqList)
    mvTotalFastqSize = 0
    for mvFastq in mvFastqList:
        mvFastqSize = os.path.getsize(mvFastq)
        logger.info('{} moved filesize: {}'.format(mvFastqList,mvFastqSize))
        mvTotalFastqSize += mvFastqSize
    return mvTotalFastqSize,mvNumFastq

def check_orig_dest_transfer(dest_drive,origNumFastq,origTotalFastqSize,mvNumFastq,mvTotalFastqSize):
    logger = logging.getLogger(__name__)
    logger.info('{} check'.format(dest_drive))
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
        raise Exception('Total sum of files sizes after {} move do not match!!!'.format(dest_drive))
    if origNumFastq != mvNumFastq:
        raise Exception('Total number of files after {} move do not match!!!'.format(dest_drive))

def updateLaneStatus(verbose,fcillumid,noStatus,database):
    logger = logging.getLogger(__name__)
    if noStatus == False:

        update_lane_step_status_query = ("UPDATE Lane l "
                                         "JOIN Flowcell f ON l.fcid=f.fcid "
                                         "SET step_status='in storage' "
                                         "WHERE FCILLUMID = '{}'"
                                        ).format(fcillumid)
        if verbose:
            print(update_lane_step_status_query)
        logger.info(update_lane_step_status_query)
        run_query(update_lane_step_status_query,database)

def storage(machine,fcillumid,args,config,run_info_dict,database):
    logger = logging.getLogger(__name__)
    bcl_drive = args.bcl_drive
    archive_drive = args.archive_dir
    noStatus = args.noStatus
    verbose = args.verbose
    UnalignedLoc = '/nfs/{}/BCL/{}_{}_{}_Unaligned'.format(bcl_drive,run_info_dict['runDate'],
                                                  machine,fcillumid)
    fastq_list = glob('{}/*/*fastq.gz'.format(UnalignedLoc))
    origNumFastq = len(fastq_list)
    origTotalFastqSize = 0
    archive_tuples_list = []
    if args.rerun == False:
        if len(fastq_list) == 0:
            raise Exception('No fastqs were found!')

        for fastq in fastq_list:
            fastqSize = os.path.getsize(fastq)
            logger.info('{} original filesize: {}'.format(fastq,fastqSize))
            origTotalFastqSize += fastqSize
            sampleName = '_'.join(fastq.split('/')[6].split('_')[0:-4])
            seqtype = getSeqtype(fcillumid,sampleName,database)
            dest_drive = config.get('locs','bcl2fastq_scratch_drive')
            gaf_bin = run_query(GET_GAFBIN_FROM_SAMPLE_NAME.format(CHGVID=sampleName),database)
            sampleArchiveLoc = '/nfs/{}/{}/{}/{}'.format(archive_drive,seqtype,sampleName,fcillumid)
            scratchArchiveLoc = '/nfs/{}/{}/{}/{}'.format(dest_drive,seqtype,
                                                          sampleName,fcillumid)
            archive_tuples_list.append((scratchArchiveLoc,sampleArchiveLoc))
            logger.debug('Starting transfer of {}'.format(sampleName))
            logger.debug('mkdir -p {}'.format(scratchArchiveLoc))
            mkdir_p(scratchArchiveLoc,verbose)
            fastqMoveCmd = ['mv',fastq,scratchArchiveLoc]
            if verbose:
                print(' '.join(fastqMoveCmd))
            logger.info(fastqMoveCmd)
            subprocess.check_call(fastqMoveCmd)

            reportLaneBarcodeLoc = ('{}/Reports/html/{}/{}/{}/all/laneBarcode.html'
                                    ).format(UnalignedLoc,fcillumid,gaf_bin[0]['GAFBIN'],sampleName)
            reportCpCmd = ['cp',reportLaneBarcodeLoc,scratchArchiveLoc]

            if verbose:
                print(' '.join(reportCpCmd))
            logger.info(reportCpCmd)
            subprocess.call(reportCpCmd)

        mvTotalFastqSize,mvNumFastq = get_mv_fastq_size_sum(dest_drive,fcillumid)
        check_orig_dest_transfer(dest_drive,origNumFastq,origTotalFastqSize,mvNumFastq,mvTotalFastqSize)

    seqscratchBase = config.get('locs','bcl_dir')
    processSAV(verbose,fcillumid,machine,run_info_dict,config,archive_drive,seqscratchBase)
    processBcl2fastqLog(verbose,fcillumid,machine,archive_drive,UnalignedLoc)
    archiveFastqs(config,archive_tuples_list,fcillumid,verbose,database)
    updateLaneStatus(verbose,fcillumid,noStatus,database)
    updateStatus(verbose,fcillumid,noStatus,database)
    updateFlowcell(verbose,fcillumid,archive_drive,database)
    createStorageCompleteFlag(verbose,config,run_info_dict)

def getSeqtype(fcillumid,sampleName,database):
    logger = logging.getLogger(__name__)
    seqtypeQuery = ("SELECT DISTINCT(SAMPLE_TYPE) AS SAMPLE_TYPE "
            "FROM Lane l "
            "JOIN Flowcell f ON l.fcid=f.fcid "
            "JOIN prepT p ON l.prepid=p.prepid "
            "WHERE FCIllumID='{}' AND CHGVID='{}'"
            ).format(fcillumid,sampleName)
    logger.info(seqtypeQuery)
    #print(seqtypeQuery)
    seqtype = run_query(seqtypeQuery,database)
    seqtype = seqtype[0]['SAMPLE_TYPE'].upper().replace(' ','_')
    return seqtype

def createStorageCompleteFlag(verbose,config,run_info_dict):
    raw_sequencing_folder_loc = '{}/{}'.format(config.get('locs','bcl_dir'),run_info_dict['runFolder'])
    flagloc = '{}/StorageComplete'.format(raw_sequencing_folder_loc)

    touchCmd = ['touch',flagloc]
    if verbose:
        print(' '.join(touchCmd))
    subprocess.check_call(touchCmd)

def processSAV(verbose,fcillumid,machine,run_info_dict,config,archive_drive,seqscratchBase):
    logger = logging.getLogger(__name__)
    if machine[0] == 'A':
        runParametersFile = 'RunParameters.xml'
    else:
        runParametersFile = 'runParameters.xml'
    raw_sequencing_folder_loc = '{}/{}'.format(config.get('locs','bcl_dir'),run_info_dict['runFolder'])
    RunInfoLoc = glob('{}/RunInfo.xml'.format(raw_sequencing_folder_loc))[0]
    runParaLoc= glob('{}/{}'.format(raw_sequencing_folder_loc,runParametersFile))[0]
    InteropsLoc = glob('{}/InterOp'.format(raw_sequencing_folder_loc))[0]
    SAVLoc = '/nfs/{}/summary/SAV/{}_{}_SAV.tar.gz'.format(archive_drive,fcillumid,machine)
    tarCmd = ['tar','czf',SAVLoc,RunInfoLoc,runParaLoc,InteropsLoc]

    if verbose:
        print(' '.join(tarCmd))
    logger.info(tarCmd)
    returnCode = subprocess.call(tarCmd)
    logger.info(returnCode)

def processBcl2fastqLog(verbose,fcillumid,machine,archive_drive,UnalignedLoc):
    logger = logging.getLogger(__name__)
    bcl2fastqLogLoc = '{}/nohup.sge'.format(UnalignedLoc)
    zipLoc = '/nfs/{}/summary/bcl_nohup/{}_{}_bcl2fastqLog.zip'.format(archive_drive,fcillumid,machine)
    zipCmd = ['zip',zipLoc,bcl2fastqLogLoc]

    if verbose:
        print(' '.join(zipCmd))
    logger.info(zipCmd)
    returnCode = subprocess.call(zipCmd)
    logger.info(returnCode)

def updateStatus(verbose,FCID,noStatus,database):
    """sequenceDB Sample Update"""
    logger = logging.getLogger(__name__)
    userID = get_user_id(database)
    #Status update for entire flowcell
    if noStatus == False:
        query = ("INSERT INTO statusT "
                "(CHGVID,status_time,status,DBID,prepID,userID) "
                "SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'Storage',pt.DBID,pt.prepID,'{}' "
                "FROM Flowcell f "
                "JOIN Lane l on l.FCID=f.FCID "
                "JOIN prepT pt on pt.prepID=l.prepID "
                "WHERE FCillumid='{}'").format(userID,FCID)
        prepT_status_update_query = """UPDATE prepT p
                                       JOIN Lane l ON p.PREPID=l.PREPID
                                       JOIN Flowcell f ON f.FCID=l.FCID
                                       SET STATUS='Storage'
                                       WHERE FCILLUMID='{}'
                                    """.format(fcillumid)

        if verbose == True:
            print(query)
            print(prepT_status_update_query)
        run_query(prepT_status_update_query,database)
        logger.info(prepT_status_update_query)
        run_query(query,database)
        logger.info(query)

def updateFlowcell(verbose,fcillumid,archive_drive,database):
    logger = logging.getLogger(__name__)
    query = ("UPDATE Flowcell "
             "SET DateStor=CURRENT_TIMESTAMP(), "
             "SeqsataLoc='{}',"
             "PIPELINECOMPLETE=1 "
             "WHERE FCIllumid='{}'").format(archive_drive,fcillumid)
    if verbose:
        print(query)
    logger.info(query)
    run_query(query,database)

def mkdir_p(path,verbose):
    try:
        if verbose:
            print("Creating dir:{}".format(path))
        os.makedirs(path)
        #os.chmod(path,0775)

    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def completeCheck(bcl_drive):
    pass

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", dest='fcillumid', required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-b','--bcl_drive', default='seqscratch_ssd', dest='bcl_drive',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('-a','--archive_dir', default='igmdata01', dest='archive_dir',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('--rerun',default=False, action='store_true',
                        help="""Allows for re-archival of sample even after its been
                              moved to its scratch location""")
    parser.add_argument("--test", default=False, action="store_true",
                        help="""Query and updates to the database occur on the
                             test server""")
    parser.add_argument('--noStatus', action='store_true',
                        help='do not update sample statuses')
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    return args

if __name__ == "__main__":
    main()
