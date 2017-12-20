import argparse
import configparser
import logging
import os
from copy import copy
from dragen_sample import dragen_sample
from glob import glob
from utilities import *

def main():
    config = get_config()
    args = parse_arguments()

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    auto_release_flag = False

    rejected_samples = []
    if args.fcillumid:
        get_sample_info_from_flowcell = """SELECT DISTINCT CHVID,p.PREPID,SEQTYPE,
                                         p.EXOMEKIT,s.PRIORITY
                                         FROM Lane l
                                         JOIN SampleT ON p.dbid=s.dbid
                                         JOIN Flowcell f ON f.fcid=l.fcid
                                         JOIN SeqType st on st.prepid=l.prepid
                                         JOIN prepT p ON p.prepid=l.prepid
                                         WHERE FCILLUMID='{}'
                                         """.format(args.fcillumid)
        print("Starting sample release for Flowcell: {}".format(args.fcillumid))
        for sample_info in run_query(get_sample_info_from_flowcell):
            seqtype = sample_info['SeqType']
            sample_name = sample_info['CHVID']
            priority=sample_info['PRIORITY']
            if seqtype.lower() == 'exome' or seqtype.lower() == 'custom_capture':
                rejected_samples = run_sample(sample_name,seqtype,
                                              sample_info['EXOMEKIT'],priority,
                                              auto_release_flag,rejected_samples,
                                              database)
        email(fcillumid,rejected_samples,config)

    else:
        run_sample(args,config,auto_release_flag,rejected_samples,database)

def email(fcillumid,rejected_samples,config):
    email_program=config.get('emails','email_program')
    failure_email_addy = config.get('emails','release_failure')


    """
    today = datetime.today().strftime("%y.%m.%d")
    if SeqType == 'exome':
        Subject = 'Exome Release '+today+' ('+release_loc+')'
    elif SeqType == 'genome':
        Subject = 'Whole Genome Release '+today+' ('+release_loc+')'
    elif SeqType == 'Custom Capture':
        Subject = 'Custom Capture Release '+today+' ('+release_loc+')'
    elif SeqType == 'RNASeq':
        Subject = 'RNASeq Release '+today+' ('+release_loc+')'
    return Subject
    """
    sampleNumber = len(rejected_samples)


    #name = get nameAasdl;fkjal
    release_email = open('release_email.txt','w')
    release_email.write('The following {} sample(s) were not release for Flowcell: {}'.format(sampleNumber,fcillumid))
    release_email.write('\n')

    for samp in rejected_samples:
        release_email.write("%s\t%s\n" % (samp))
    release_email.write('\nThanks,\n')
    release_email.write('\n')
    release_email.write('%s\n' % name)

    release_email.close()
    email_cmd = "{} -s 'Sample Release' {} < email.tmp".format(email_program,failure_email_addy)

    #os.system(email_cmd)
    #os.remove('release_email.txt')


def update_queues(sample_name,sample_type,capture_kit,ppid,priority,database):
    table = 'tmp_dragen'
    queue_query = """SELECT *
                     FROM {table}
                     WHERE PSEUDO_PREPID = {ppid}
                  """
    tmp_dragen_query = queue_query.format(table=table,ppid=ppid)
    in_tmp_dragen = run_query(tmp_dragen_query,database)
    if in_tmp_dragen:
        print("Moving sample to dragen_queue")
        insert_query = "INSERT INTO dragen_queue " + tmp_dragen_query
        run_query(insert_query,database)
        ppid = int(ppid)
        rm_query = """DELETE FROM {table}
                      WHERE PSEUDO_PREPID = {ppid}
                   """.format(table=table,ppid=ppid)
        table = 'dragen_queue'
        dragen_queue_query = queue_query.format(table=table,ppid=ppid)
        in_dragen_queue = run_query(dragen_queue_query,database)
        if in_dragen_queue:
            print("Removing samples from tmp_dragen")
            run_query(rm_query,database)

    else:
        table = 'dragen_queue'
        dragen_queue_query = queue_query.format(table=table,ppid=ppid)
        in_dragen_queue = run_query(dragen_queue_query,database)
        if not in_dragen_queue:
            print("Inserting sample into dragen_queue")
            insert_query = """INSERT INTO dragen_queue
                              (PSEUDO_PREPID,SAMPLE_NAME,SAMPLE_TYPE,
                              CAPTURE_KIT,PRIORITY)
                              VALUES ({ppid},'{sample_name}','{sample_type}',
                              '{capture_kit}',{priority})
                           """.format(sample_name=sample_name,
                                      sample_type=sample_type,
                                      capture_kit=capture_kit,
                                      ppid=ppid,
                                      priority=priority)
            run_query(insert_query,database)

    print("Updating sample {} priority to {}".format(sample_name,priority))
    update_dragen_queue = """UPDATE dragen_queue
                             SET PRIORITY={priority}
                             WHERE PSEUDO_PREPID = {ppid}
                          """.format(priority=priority,ppid=ppid)
    run_query(update_dragen_queue,database)



def check_db_check_seqscratch(config,sample_name,sample_type,ppid,database):
    ssd_dir = '{}/{}/{}.{}'.format(config.get('locs','dragen_aligment_dir'),
                                         sample_type.upper(),sample_name,ppid)
    ft_dir = '{}/{}/{}.{}'.format('/nfs/fastq_temp/ALIGNMENT/',sample_type.upper(),sample_name,ppid)
    ft2_dir = '{}/{}/{}.{}'.format('/nfs/fastq_temp2',sample_type.upper(),sample_name,ppid)
    alignment_dirs = [ssd_dir,ft_dir,ft2_dir]
    for alignment_dir in alignment_dirs:
        sample_alignment = (glob(alignment_dir))
        if sample_alignment != []:
            raise ValueError("Scratch alignment director found.  Cleanup required first!")

    sample_status_query = """SELECT MAX(PIPELINE_STEP_ID) as PIPELINE_STEP_ID
                             FROM dragen_pipeline_step
                             WHERE STEP_STATUS='completed' AND
                             PSEUDO_PREPID={ppid}
                          """.format(ppid=ppid)
    sample_status = run_query(sample_status_query,database)

    if sample_status[0]['PIPELINE_STEP_ID'] != None:
        raise ValueError("Sample {} has already been run in GATK pipe.  Cleanup required first".format(sample_name))

    dsm_query = """SELECT * FROM dragen_sample_metadata
                   WHERE PSEUDO_PREPID={ppid}
                """.format(ppid=ppid)
    dsm_status = run_query(dsm_query,database)
    if dsm_status:
        if dsm_status[0]['is_merged'] != 0:
            raise ValueError("Sample {} has an attempted merge already! Cleanup required first".format(sample_name))
        if dsm_status[0]['component_bams']:
            print(dsm_status[0]['component_bams'])
            raise ValueError("Sample {} already has component bams!".format(sample_name))


def run_sample(args,config,auto_release_flag,rejected_samples,database):
    sample_name = args.sample_name
    sample_type = args.sample_type
    capture_kit = args.capture_kit
    ppid,pids = update_ppid(sample_name,sample_type,capture_kit,database)
    sample = dragen_sample(sample_name,sample_type,ppid,
                           capture_kit,database)
    if args.priority is None:
        priority = sample.metadata['priority']
    else:
        priority = args.priority

    check_db_check_seqscratch(config,sample_name,sample_type,ppid,database)

    print("Starting sample release for sample: {}, {}, {}".format(sample_name,
                                                                  sample_type,
                                                                  capture_kit))

    #every external sample should only have one pids
    if is_external_or_legacy_sample(pids[0],database) == True:
        pass
    else:
        rejected_samples = check_yield(ppid,sample_type,
                                       auto_release_flag,
                                       rejected_samples,database)
    update_queues(sample_name,sample_type,capture_kit,ppid,priority,database)
    ##update statusT
    ##update prepT.status
    return rejected_samples

def check_yield(ppid,sample_type,auto_release_flag,rejected_samples,database):
   query = GET_YIELD_FROM_PPID.format(ppid=ppid)
   lane_yield_sum = run_query(query,database)[0]['LANE_YIELD_SUM']
   lane_yield_sum = int(lane_yield_sum)

   if sample_type.lower() == 'exome' and lane_yield_sum > 7000:
       return rejected_samples
   elif sample_type.lower() == 'genome' and lane_yield_sum > 100000:
       return rejected_samples
   elif sample_type.lower() == 'custom_capture' and lane_yield_sum > 150:
       return rejected_samples
   else:
       if auto_release_flag == False:
           pass
       else:

           rejected_samples.append(ppid)
       return rejected_samples

def update_prepT_ppid(ppid,pid,database):
    INSERT_ppid_on_pid = """UPDATE prepT
                            SET P_PREPID = {}
                            WHERE PREPID = {}
                         """.format(ppid,pid)
    print("Updating prepT with ppid: {}".format(ppid))
    run_query(INSERT_ppid_on_pid,database)


def update_ppid(sample_name,sample_type,capture_kit,database):
    GET_PID_PPID_FROM_TRIPLET= """
        SELECT P_PREPID,p.PREPID
        FROM prepT p
        JOIN SeqType st on p.prepid=st.prepid
        WHERE CHGVID='{chgvid}'
        AND SEQTYPE='{seqtype}'
        AND FAILEDPREP = 0
        """.format(chgvid=sample_name,
                   seqtype=sample_type)

    if sample_type != 'genome':
        GET_PID_PPID_FROM_TRIPLET += ("AND p.EXOMEKIT='{exomekit}'"
                                     ).format(exomekit=capture_kit)
    pid_and_ppid = run_query(GET_PID_PPID_FROM_TRIPLET,database)

    if pid_and_ppid == ():
        raise Exception('No preps were found!')
    else:
        ppid_from_triplet = [x['P_PREPID'] for x in pid_and_ppid]
        prepid_from_triplet = [x['PREPID'] for x in pid_and_ppid]
        if len(set(ppid_from_triplet)) != 1:
            raise Exception('Too many ppids found from triplet:{}'.format(set(ppid_from_triplet)))

        return list(set(ppid_from_triplet))[0],prepid_from_triplet


        """  ---->  check difference from DSM <------"""

def remove_pid_from_ppid(pid,ppid,database):
    check_fail_query = """SELECT FAILPREP
                          FROM PREPT
                          WHERE PREPID={pid}
                       """.format(pid=pid)
    fail_flag = run_query(check_fail_query,database)[0]['FAILPREP']
    if fail_flag != '1':
        raise ValueError('prepid flagged for removal from ppid but its not failed: {}'.format(pid))

    rm_query = """DELETE FROM PSEUDO_PREPID
                  WHERE PREPID={pid}
               """.format(pid=pid)
    run_query(rm_query,database)

def add_pid_to_ppid(pid,ppid,database):
    query = """INSERT INTO PSEUDO_PREPID
               (PSEUDO_PREPID,PREPID)
               VALUES ({ppid},{pid})
            """.format(ppid=ppid,pid=pid)

    run_query(query,database)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fcillumid",
                        help="Specify Illumina's Flowcell ID")
    group.add_argument('-s', '--sample_name',
                        help="Specify sample name")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-t','--sample_type',
                        help="Specify sample type (exome,genome,etc)")
    parser.add_argument('-k','--capture_kit',
                        help="Specify capture kit ")
    parser.add_argument("-p", "--priority", type=int)
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    if args.sample_name and (bool(args.sample_type) == False or  bool(args.capture_kit) == False) :
        parser.error('--sample_name and --sample_type and --capture_kit must be given together')
    if args.fcillumid and (bool(args.priority)) == True:
        parser.error('--priority is not compatible with --flowcell')
    return args
if __name__ == '__main__':
    main()
