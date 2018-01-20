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
    check_bcl_complete(unaligned_dir)
    logger.info('Starting Storage script')
    # Josh Bridgers, Sophia Frantz, Brett Copeland
    emailAddress = config.get('emails','storage_success')
    emailAddressFail = config.get('emails','storage_failure')

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
        ).format(storageStatus,FCID,storageStatus,emailAddress)
    print(emailCmd)
    logger.info(emailCmd)
    os.system(emailCmd)

def generate_archive_tuples_list(rerun,config,run_info_dict,fcillumid,archive_drive,database):
    bcl_drive = config.get('locs','bcl2fastq_scratch_drive')
    UnalignedLoc = '/nfs/{}/BCL/{}_{}_{}_Unaligned'.format(bcl_drive,run_info_dict['runDate'],
                                                           machine,fcillumid)
    archive_tuples_list = []
    query = """SELECT CHGVID,FCILLUMID,LANENUM,SAMPLE_TYPE
               FROM prepT p
               JOIN Lane l ON l.PREPID=p.PREPID
               JOIN Flowcell f ON l.FCID=f.FCID
               WHERE f.fcillumid='{}'
            """.format(fcillumid)
    sample_info_on_flowcell = run_query(query,database)
    for sample in sample_info_on_flowcell:
        if rerun == True:
            bcl_loc = ''
        else:
            gaf_bin_query = GET_GAFBIN_FROM_SAMPLE_NAME.format(CHGVID=sample_info_on_flowcell['CHGVID'])
            gaf_bin = run_query(gaf_bin_query,database)
            bcl_loc = UnalignedLoc + '/' + gaf_bin_query[0]['GAFBIN'].upper()
        scratchLoc = '/nfs/{}/{}/{}/{}'.format(config.get('locs','bcl2fastq_scratch_drive'),
                                                     sample['SAMPLE_TYPE'].upper(),
                                                     sample['CHGVID'],fcillumid)
        sampleArchiveLoc = '/nfs/{}/{}/{}/{}'.format(config.get('locs','fastq_archive_drive'),
                                                     sample['SAMPLE_TYPE'].upper(),
                                                     sample['CHGVID'],fcillumid)
        archive_tuples_list.append(bcl_loc,scratchLoc,sampleArchiveLoc,sample['CHGVID'])
    return archive_tuples_list

def get_distinct_seqtype_on_flowcell(fcillumid,database):
    seqtypes_on_flowcell_query = ("SELECT DISTINCT(SAMPLE_TYPE) AS SAMPLE_TYPE "
                                  "FROM Lane l "
                                  "JOIN Flowcell f ON l.fcid=f.fcid "
                                  "JOIN prepT p ON l.prepid=p.prepid "
                                  "WHERE FCIllumID = '{}'").format(fcillumid)
    seqtypes_on_flowcell = run_query(seqtypes_on_flowcell_query,database)

    distinct_seqtype = []
    for entry in seqtypes_on_flowcell:
        distinct_seqtype.append(entry['SAMPLE_TYPE'].upper())
    distinct_seqtype = list(set(distinct_seqtype))
    if len(distinct_seqtype) == 1:
        return ','.join(distinct_seqtype)
    elif len(distinct_seqtype) > 1:
        return  '{' + ','.join(distinct_seqtype) + '}'
    else:
        raise Exception('No seqtypes were found!')

def archiveFastqs(config,archive_tuples_list,fcillumid,verbose,database):
    logger = logging.getLogger(__name__)
    bcl_fastqs = glob('{}/*/*fastq.gz'.format(orig_dest_triplet[0])
    bcl_total_num_fastq = len(bcl_fastqs)
    bcl_total_fastq_size = 0
    for bcl_fastq in bcl_fastqs:
        bcl_total_fastq_size += os.path.getsize(bcl_fastq)
    distinct_seqtype = get_distinct_seqtype_on_flowcell(fcillumid,database)
    scratch_check_drive = orig_dest_triplet[offset+1].split('/')[2]
    logger.info('checking {}'.format(scratch_check_drive))
    scratch_fastqs = glob('/nfs/{}/{}/*/{}/*.fastq.gz'
                      ).format(scratch_check_drive,distinct_seqtype,fcillumid)
    scratch_total_num_fastq = len(scratch_fastqs)
    scratch_total_fastq_size = 0
    for scratch_fastq in scratch_fastqs:
        scratch_total_fastq_size += os.path.getsize(scratch_fastq)
    check_orig_dest_transfer('scratch',bcl_total_num_fastq,bcl_total_fastq_size,
        scratch_total_num_fastq,scratch_total_fastq_size)

    if rerun is False:
        offset = 0
        for orig_dest_triplet in archive_tuples_list:
            mv_rsync_fastq(orig_dest_triplet,offset,verbose)

    offset = 1

    for orig_dest_triplet in archive_tuples_list:
        mv_rsync_fastq(orig_dest_triplet,offset,verbose)
    logger.info('checking {}'.format(archive_check_drive))
    archive_check_drive = orig_dest_triplet[offset+1].split('/')[2]
    logger.info('checking {}'.format(dest_check_drive))
    archive_fastqs = glob('/nfs/{}/{}/*/{}/*.fastq.gz'
                      ).format(archive_check_drive,distinct_seqtype,fcillumid)
    archive_total_num_fastq = len(archive_fastqs)
    archive_total_fastq_size = 0
    for archive_fastq in archive_fastqs:
        archive_total_fastq_size += os.path.getsize(archive_fastq)
    check_orig_dest_transfer('archive',bcl_total_num_fastq,bcl_total_fastq_size,
        archive_total_num_fastq,archive_total_fastq_size)


def mv_rsync_fastq(orig_dest_triplet,offset,verbose):
    fastqs = glob('{}/*fastq.gz'.format(orig_dest_triplet[offset])
    logger.debug('Starting transfer of {}'.format(orig_dest_triplet[offset]))
    mkdir_p(orig_dest_triplet[offset + 1],verbose)
    for fastq in fastqs:
        fastq_filename=os.path.basename(fastq)
        fastqSize = os.path.getsize(fastq)
        origTotalFastqSize += fastqSize
        logger.info('{} original filesize: {}'.format(fastq,fastqSize))

        if offset = 0:
            fastqMoveCmd = ['mv', fastq ,orig_dest_pair[offset + 1]]
        else:
            fastqMoveCmd = ['rsync','--inplace', '-aq', fastq ,orig_dest_pair[offset + 1]]

            gaf_bin_query = GET_GAFBIN_FROM_SAMPLE_NAME.format(CHGVID=orig_dest_triplet[3])
            gaf_bin = run_query(gaf_bin_query,database)
            reportLaneBarcodeLoc = ('{}/Reports/html/{}/{}/{}/all/laneBarcode.html'
                               ).format(UnalignedLoc,fcillumid,gaf_bin[0]['GAFBIN'],sampleName)
            if verbose:
                print(' '.join(reportLaneBarcodeLoc)
            logger.info(reportLaneBarcodeLoc)
            subprocess.check_call(reportLaneBarcodeLoc)

        if verbose:
            print(' '.join(fastqMoveCmd)
        logger.info(fastqMoveCmd)
        subprocess.check_call(fastqMoveCmd)a

def get_mv_fastq_size_sum(drive,fcillumid,database):
    logger = logging.getLogger(__name__)
    distinct_seqtype = get_distinct_seqtype_on_flowcell(fcillumid,database)
    dest_check_loc ='/nfs/{}/{}/*/{}/*.fastq.gz'.format(drive,distinct_seqtype,fcillumid)
    logger.info('checking {}'.format(dest_check_loc))
    mvFastqList = glob(dest_check_loc)
    mvNumFastq = len(mvFastqList)
    mvTotalFastqSize = 0
    for mvFastq in mvFastqList:
        mvFastqSize = os.path.getsize(mvFastq)
        logger.info('{} moved filesize: {}'.format(mvFastq,mvFastqSize))
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
    archive_drive = config.get('locs','fastq_archive_drive')
    noStatus = args.noStatus
    verbose = args.verbose

    fastq_list = glob('{}/*/*fastq.gz'.format(UnalignedLoc))
    origNumFastq = len(fastq_list)
    origTotalFastqSize = 0
    archive_tuples_list = generate_archive_tuples_list(args.rerun,config,
            run_info_dict,fcillumid,archive_drive,database)

    archiveFastqs(config,archive_tuples_list,fcillumid,verbose,database)
    processSAV(verbose,fcillumid,machine,run_info_dict,config,archive_drive)
    processBcl2fastqLog(verbose,fcillumid,machine,archive_drive,UnalignedLoc)
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

def processSAV(verbose,fcillumid,machine,run_info_dict,config,archive_drive):
    logger = logging.getLogger(__name__)
    if machine[0] == 'A':
        runParametersFile = 'RunParameters.xml'
    else:
        runParametersFile = 'runParameters.xml'
    seqscratchBase = config.get('locs','bcl_dir')
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
    zipCmd = ['/usr/bin/zip',zipLoc,bcl2fastqLogLoc]

    if verbose:
        print(' '.join(zipCmd))
    logger.info(zipCmd)
    returnCode = subprocess.call(zipCmd)
    logger.info(returnCode)

def updateStatus(verbose,fcillumid,noStatus,database):
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
                "WHERE FCillumid='{}'").format(userID,fcillumid)
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
            msg = "Creating dir:{}".format(path)
            print(msg)
            logger.debug(msg)
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
