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
    ######## check for touchfile
    check_bcl_complete(unaligned_dir)

    logger.info('Starting Storage script')

    # Josh Bridgers, Sophia Frantz, Brett Copeland
    emailAddress = config.get('emails','storage_success')
    emailAddressFail = config.get('emails','storage_failure')

    #completeCheck(bcl_drive,inputFolder)
    try:
        storage(machine,fcillumid,args,config,run_info_dict,database)
        email(emailAddress,'SUCCESS',fcillumid,args.test)
        logger.info('Done')
        if args.verbose:
            print('Done')

    except Exception:
        logger.exception('Got exception')
        traceback.print_exc()
        email(emailAddressFail,'FAILURE',fcillumid,args.test)
        logger.info('Storage Failure')
        print('Storage Failure')
        sys.exit(255)


def email(emailAddress,storageStatus,fcillumid,test):
    logger = logging.getLogger(__name__)
    logger.info('Starting email')
    subject ='"IGM:BCL move {}'.format(storageStatus)
    if test == True:
        subject += ' TEST"'
    else:
        subject += '"'
    emailCmd = ('echo "BCL move {} for {}"  | mail -s {} {}'
        ).format(storageStatus,fcillumid,subject,emailAddress)
    print(emailCmd)
    logger.info(emailCmd)
    op = os.system(emailCmd)
    if op != 0:
        raise Exception("{} incorrectly returned {}".format(emailCmd,op))

def generate_archive_tuples_list(rerun,config,run_info_dict,fcillumid,
        UnalignedLoc,archive_drive,database):
    machine = run_info_dict['machine']
    archive_tuples_list = []
    query = """SELECT DISTINCT CHGVID,FCILLUMID,SAMPLE_TYPE,EXOMEKIT
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
            gaf_bin_query = (GET_GAFBIN_FROM_SAMPLE_NAME
                            ).format(CHGVID=sample['CHGVID'])
            gaf_bin = run_query(gaf_bin_query,database)
            bcl_loc = UnalignedLoc + '/' + gaf_bin[0]['GAFBIN']
        scratchLoc = '/nfs/{}/{}/{}/{}'.format(config.get('locs','bcl2fastq_scratch_drive'),
                                                     sample['SAMPLE_TYPE'].upper(),
                                                     sample['CHGVID'],fcillumid)
        sampleArchiveLoc = '/nfs/{}/{}/{}/{}'.format(config.get('locs','fastq_archive_drive'),
                                                     sample['SAMPLE_TYPE'].upper(),
                                                     sample['CHGVID'],fcillumid)
        archive_tuples_list.append((bcl_loc,scratchLoc,sampleArchiveLoc,sample['CHGVID']))
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
    if len(distinct_seqtype) < 1:
        raise Exception('No seqtypes were found!')
    else:
        return distinct_seqtype

def archiveFastqs(offset,rerun,config,archive_tuples_list,UnalignedLoc,fcillumid,verbose,database):

    logger = logging.getLogger(__name__)
    distinct_seqtypes = get_distinct_seqtype_on_flowcell(fcillumid,database)
    if offset==0: # offset = 0
        if rerun is False:
            bcl_fastqs = glob('{}/*/*fastq.gz'.format(UnalignedLoc))
            bcl_total_num_fastq = len(bcl_fastqs)
            if bcl_total_num_fastq == 0:
                raise Exception("no fastqs to move from {0}".format(UnalignedLoc))
            bcl_total_fastq_size = 0
            for bcl_fastq in bcl_fastqs:
                bcl_total_fastq_size += os.path.getsize(bcl_fastq)
            if bcl_total_fastq_size == 0:
                raise Exception("total fastq size = 0 in {0}".format(UnalignedLoc))

            for archive_locs_list in archive_tuples_list:
                mv_rsync_fastq(
                archive_locs_list,
                UnalignedLoc,
                fcillumid,
                offset,
                verbose,
                database
                )
    elif offset==1:
        for archive_locs_list in archive_tuples_list:
            mv_rsync_fastq(
            archive_locs_list,
            UnalignedLoc,
            fcillumid,
            offset,
            verbose,
            database
            )
        archive_check_drive = config.get('locs','fastq_archive_drive')
        logger.info('checking {}'.format(archive_check_drive))
        archive_total_fastq_size = 0
        archive_total_num_fastq = 0
        for seqtype in distinct_seqtypes:
            archive_fastq_loc = ('/nfs/{}/{}/*/{}/*.fastq.gz'
                                ).format(archive_check_drive,seqtype,fcillumid)
            logger.info('checking {} for {} fastqs'.format(archive_fastq_loc,fcillumid))
            archive_fastqs = glob(archive_fastq_loc)
            if len(archive_fastqs) == 0:
                raise Exception("no fastqs of flowcell {} of type {} in {}".format(fcillumid,seqtype,archive_fastq_loc))
            archive_total_num_fastq += len(archive_fastqs)
            for archive_fastq in archive_fastqs:
                sz_archive_fastq = os.path.getsize(archive_fastq)
                if sz_archive_fastq == 0:
                    raise Exception("{} is empty. What's happening here?".format(archive_fastq))
                archive_total_fastq_size += sz_archive_fastq
    else:
        raise ValueError("Like, WTF.?!?")

    scratch_check_drive = config.get('locs','bcl2fastq_scratch_drive')
    logger.info('checking {}'.format(scratch_check_drive))
    scratch_total_fastq_size = 0
    scratch_total_num_fastq = 0
    for seqtype in distinct_seqtypes:
            scratch_fastq_loc = ('/nfs/{}/{}/*/{}/*.fastq.gz').format(scratch_check_drive,seqtype,fcillumid)
            logger.info('checking {} for {} fastqs'.format(scratch_fastq_loc,fcillumid))
            scratch_fastqs = glob(scratch_fastq_loc)
            scratch_total_num_fastq += len(scratch_fastqs)
            for scratch_fastq in scratch_fastqs:
                sz_scratch_fastq = os.path.getsize(scratch_fastq)
                if sz_scratch_fastq == 0:
                    raise Exception("{} is empty. What's happening here?".format(scratch_fastq))
                scratch_total_fastq_size += sz_scratch_fastq
    if offset == 0 and rerun is False:
        check_orig_dest_transfer('scratch',bcl_total_num_fastq,bcl_total_fastq_size,scratch_total_num_fastq,scratch_total_fastq_size)
    if offset == 1:
        check_orig_dest_transfer('archive',scratch_total_num_fastq,scratch_total_fastq_size,archive_total_num_fastq,archive_total_fastq_size)



def mv_rsync_fastq( archive_locs_list, UnalignedLoc, fcillumid, offset, verbose, database ):

    logger = logging.getLogger(__name__)

    if offset == 0:
        fastqs = glob('{}/{}_*fastq.gz'.format(archive_locs_list[offset],
                                             archive_locs_list[-1]))
    else:
        fastqs = glob('{}/*fastq.gz'.format(archive_locs_list[offset]))

    msg = 'Starting transfer of {}'.format(archive_locs_list[offset])
    logger.debug(msg)
    if verbose:
        print(msg)
    mkdir_p(archive_locs_list[offset + 1],verbose)
    gaf_bin_query = GET_GAFBIN_FROM_SAMPLE_NAME.format(CHGVID=archive_locs_list[3])
    gaf_bin = run_query(gaf_bin_query,database)
    reportLaneBarcodeLoc = ('{}/Reports/html/{}/{}/{}/all/laneBarcode.html'
                       ).format(UnalignedLoc,fcillumid,gaf_bin[0]['GAFBIN'],
                               archive_locs_list[3])

    reportCpCmd = ['cp',reportLaneBarcodeLoc,archive_locs_list[offset+1]]

    if verbose:
        print(' '.join(reportCpCmd))
    logger.info(reportCpCmd)
    op = subprocess.call(reportCpCmd)
    if op != 0:
        raise Exception("{} exited incorrectly: {}".format(reportCpCmd,op))
    for fastq in fastqs:
        fastq_filename=os.path.basename(fastq)
        fastqSize = os.path.getsize(fastq)
        logger.info('{} original filesize: {}'.format(fastq,fastqSize))

        if offset == 0:
            fastqMoveCmd = ['mv', fastq ,archive_locs_list[offset + 1]]
        else:
            fastqMoveCmd = ['rsync','--inplace', '-aq', fastq ,archive_locs_list[offset + 1]]

        if verbose:
            print(' '.join(fastqMoveCmd))
        logger.info(fastqMoveCmd)
        subprocess.check_call(fastqMoveCmd)

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

def check_orig_dest_transfer(dest_drive,origNumFastq,origTotalFastqSize,
        mvNumFastq,mvTotalFastqSize):
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
    if 0 in (origNumFastq,mvNumFastq,origTotalFastqSize,mvTotalFastqSize):
        raise Exception("number or size of fastq is 0: {0} {1} {2} {3}".format(origNumFastq,mvNumFastq,origTotalFastqSize,mvTotalFastqSize))
    if origTotalFastqSize != mvTotalFastqSize:
        raise Exception('Total sum of files sizes after {} move do not match!!!'.format(dest_drive))
    if origNumFastq != mvNumFastq:
        raise Exception('Total number of files after {} move do not match!!!'.format(dest_drive))

def update_rg_aka_lane_table_to_in_storage(WHAT,WHAT_rg,verbose,fcillumid,noStatus,database):
    logger = logging.getLogger(__name__)
    if noStatus == False:

        update_lane_step_status_query = ("UPDATE Lane l "
                                         "JOIN Flowcell f ON l.fcid=f.fcid "
                                         "SET step_status='{}', rg_status='{}' "
                                         "WHERE FCILLUMID = '{}'"
                                        ).format(WHAT,WHAT_rg,fcillumid)
        if verbose:
            print(update_lane_step_status_query)
        logger.info(update_lane_step_status_query)
        run_query(update_lane_step_status_query,database)

def storage(machine,fcillumid,args,config,run_info_dict,database):
    logger = logging.getLogger(__name__)
    archive_drive = config.get('locs','fastq_archive_drive')
    noStatus = args.noStatus
    verbose = args.verbose
    rerun = args.rerun
    bcl_drive = config.get('locs','bcl2fastq_scratch_drive')
    UnalignedLoc = ('/nfs/{}/BCL/{}_{}_{}_Unaligned'
                   ).format(bcl_drive,run_info_dict['runDate'],
                            machine,fcillumid)

    archive_tuples_list = generate_archive_tuples_list(args.rerun,config,
            run_info_dict,fcillumid,UnalignedLoc,archive_drive,database)

    ####################################################################################
    archiveFastqs(0, rerun, config, archive_tuples_list, UnalignedLoc, fcillumid, verbose, database )

    query = ("UPDATE Flowcell SET DateStor=CURRENT_TIMESTAMP(), SeqsataLoc='{}', PIPELINECOMPLETE=1 WHERE FCIllumid='{}'").format(archive_drive,fcillumid)
    logger.info(query)
    run_query(query,database)

    update_rg_aka_lane_table_to_in_storage('in storage','fastq_ready',verbose,fcillumid,noStatus,database)
    update_status_to_storage('Archiving',verbose,fcillumid,noStatus,database)

    ####################################################################################
    archiveFastqs(1, rerun, config, archive_tuples_list, UnalignedLoc, fcillumid, verbose, database )

    processSAV(verbose,fcillumid,machine,run_info_dict,config,archive_drive)
    processBcl2fastqLog(verbose,fcillumid,machine,archive_drive,UnalignedLoc)

    update_rg_aka_lane_table_to_in_storage('in storage','fastq_copied',verbose,fcillumid,noStatus,database)
    update_status_to_storage('Storage',verbose,fcillumid,noStatus,database)
    create_checkpoint_file(verbose,config,run_info_dict,UnalignedLoc)

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

def create_checkpoint_file(verbose,config,run_info_dict,UnalignedLoc):
    raw_sequencing_folder_loc = '{}/{}'.format(config.get('locs','bcl_dir'),
                                              run_info_dict['runFolder'])
    flagloc = '{}/StorageComplete'.format(raw_sequencing_folder_loc)
    flagloc_ar = '{}/StorageComplete'.format(UnalignedLoc)
    touchCmd = ['touch',flagloc]
    if verbose:
        print(' '.join(touchCmd))
    subprocess.check_call(touchCmd)
    subprocess.check_call(['touch',flagloc_ar])

def processSAV(verbose,fcillumid,machine,run_info_dict,config,archive_drive):
    logger = logging.getLogger(__name__)
    if machine[0] == 'A':
        runParametersFile = 'RunParameters.xml'
    else:
        runParametersFile = 'runParameters.xml'
    seqscratchBase = config.get('locs','bcl_dir')
    raw_sequencing_folder_loc = '{}/{}'.format(config.get('locs','bcl_dir'),
                                              run_info_dict['runFolder'])
    RunInfoLoc = glob('{}/RunInfo.xml'.format(raw_sequencing_folder_loc))[0]
    runParaLoc= glob('{}/{}'.format(raw_sequencing_folder_loc,runParametersFile))[0]
    InteropsLoc = glob('{}/InterOp'.format(raw_sequencing_folder_loc))[0]
    SAVLoc = '/nfs/{}/summary/SAV/{}_{}_SAV.tar.gz'.format(archive_drive,fcillumid,machine)
    if not os.path.isdir('/nfs/{}/summary/SAV'.format(archive_drive)):
        raise Exception("/nfs/{}/summary/SAV does not exist. Create it before tarring SAV".format(archive_drive))
    tarCmd = ['tar','czf',SAVLoc,RunInfoLoc,runParaLoc,InteropsLoc]
    if verbose:
        print(' '.join(tarCmd))
    logger.info(tarCmd)
    try:
        returnCode = subprocess.check_output(tarCmd)
        logger.info(returnCode)
    except subprocess.CalledProcessError as e:
        sys.exit("'%s' failed, returned code %d" % (tarCmd,e.returncode))
    except OSError as e:
        sys.exit("failed to execute program '%s': '%s'" % (tarCmd, str(e)))

def processBcl2fastqLog(verbose,fcillumid,machine,archive_drive,UnalignedLoc):
    logger = logging.getLogger(__name__)
    bcl2fastqLogLoc = '{}/nohup.sge'.format(UnalignedLoc)
    zipLoc = '/nfs/{}/summary/bcl_nohup/{}_{}_bcl2fastqLog.zip'.format(archive_drive,
                                                                      fcillumid,machine)
    if not os.path.isdir('/nfs/{}/summary/bcl_nohup'.format(archive_drive)):
        raise Exception("/nfs/{}/summary/bcl_nohup does not exist. Create it before zipping log".format(archive_drive))

    zipCmd = ['/usr/bin/zip',zipLoc,bcl2fastqLogLoc]

    if verbose:
        print(' '.join(zipCmd))
    logger.info(zipCmd)
    try:
        returnCode = subprocess.check_output(zipCmd)
        logger.info(returnCode)
    except subprocess.CalledProcessError as e:
        sys.exit("'%s' failed, returned code %d" % (zipCmd,e.returncode))
    except OSError as e:
        sys.exit("failed to execute program '%s': '%s'" % (zipCmd, str(e)))


def update_status_to_storage(WHAT,verbose,fcillumid,noStatus,database):
    """sequenceDB Sample Update"""
    logger = logging.getLogger(__name__)
    userID = get_user_id(database)
    #Status update for entire flowcell
    if noStatus == False:
        query = ("INSERT INTO statusT "
                "(CHGVID,STATUS_TIME,STATUS,sample_id,PREPID,USERID,POOLID,SEQID,PLATENAME) "
                "SELECT DISTINCT(pt.CHGVID),UNIX_TIMESTAMP(),'{}',pt.sample_id,"
                "pt.PREPID,'{}',0,0,' ' "
                "FROM Flowcell f "
                "JOIN Lane l ON l.FCID=f.FCID "
                "JOIN prepT pt ON pt.prepID=l.prepID "
                "WHERE FCillumid='{}'").format(WHAT,userID,fcillumid)
        prepT_status_update_query = """UPDATE prepT p
                                       JOIN Lane l ON p.PREPID=l.PREPID
                                       JOIN Flowcell f ON f.FCID=l.FCID
                                       SET STATUS='{}', status_time=unix_timestamp()
                                       WHERE FCILLUMID='{}'
                                    """.format(WHAT,fcillumid)

        if verbose == True:
            print(query)
            print(prepT_status_update_query)
        run_query(prepT_status_update_query,database)
        logger.info(prepT_status_update_query)
        run_query(query,database)
        logger.info(query)

def mkdir_p(path,verbose):
    logger = logging.getLogger(__name__)
    try:
        if verbose:
            msg = "Creating dir:{}".format(path)
            print(msg)
            logger.debug(msg)
        os.makedirs(path,0o770)
        #os.chmod(path,0775)

    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise Exception

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
