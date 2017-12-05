import argparse
import configparser
import logging
import os
from copy import copy
from dragen_sample import dragen_sample
from utilities import *

def main():
    config = get_config()
    args = parse_arguments()

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    auto_release_flag = True

    rejected_samples = []
    if args.fcillumid:
        auto_release_flag = False
        get_sample_info_from_flowcell = """SELECT DISTINCT CHVID,p.PREPID,SEQTYPE,
                                         p.EXOMEKIT,s.PRIORITY
                                         FROM Lane l
                                         JOIN SampleT ON p.dbid=s.dbid
                                         JOIN Flowcell f ON f.fcid=l.fcid
                                         JOIN SeqType st on st.prepid=l.prepid
                                         JOIN prepT p ON p.prepid=l.prepid
                                         WHERE FCILLUMID='{}'
                                         """.format(args.fcillumid)
        print(get_sample_info_from_flowcell)
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
        run_sample(args.sample_name,args.sample_type,args.capture_kit,
                   args.priority,auto_release_flag,rejected_samples,database)

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
                                      priority=40)
            run_query(insert_query,database)
    if priority is not None:
        print("Updating sample {} priority to {}".format(sample_name,priority))
        update_dragen_queue = """UPDATE dragen_queue
                                 SET PRIORITY={priority}
                                 WHERE PSEUDO_PREPID = {ppid}
                              """.format(priority=priority,ppid=ppid)
        run_query(update_dragen_queue,database)

def run_sample(sample_name,sample_type,capture_kit,priority,
               auto_release_flag,rejected_samples,database):
    ppid,pids = update_ppid(sample_name,sample_type,capture_kit,database)
    sample_status_query = """ SELECT MAX(PIPELINE_STEP_ID) as PIPELINE_STEP_ID
                              FROM dragen_pipeline_step
                              WHERE STEP_STATUS='completed' AND
                              PSEUDO_PREPID={ppid}
                          """.format(ppid=ppid)
    sample_status = run_query(sample_status_query,database)[0]['PIPELINE_STEP_ID']
    if sample_status:
        raise ValueError("Sample {} has already been run.  Cleanup required first".format(sample_name))
    else:
        print("Starting sample release for sample: {}, {}, {}".format(sample_name,
                                                                      sample_type,
                                                                      capture_kit))
        sample = dragen_sample(sample_name,sample_type,ppid,
                               capture_kit,database)
        #every external sample should only have one pids
        is_external = run_query(IS_SAMPLE_EXTERNAL_FROM_PREPID.format(prepid=pids[0]),database)
        if is_external == True or int(pids[0]) < 20000:
            print("Sample is external or legacy")
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
   lane_yield_sum = run_query(GET_YIELD_FROM_PPID.format(ppid=ppid),database)[0]['LANE_YIELD_SUM']
   lane_yield_sum = int(lane_yield_sum)

   if sample_type.lower() == 'exome' and lane_yield_sum > 7500:
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

def insert_pid_into_ppid(sample_name,sample_type,capture_kit,database):
    GET_PIDS_QUERY = """SELECT PREPID
                        FROM prepT p
                        WHERE CHGVID='{chgvid}'
                        AND SAMPLE_TYPE='{seqtype}'
                     """.format(chgvid=sample_name,seqtype=sample_type)
    if sample_type.lower() != 'genome':
        GET_PIDS_QUERY += ("AND EXOMEKIT='{exomekit}'"
                          ).format(exomekit=capture_kit)

    pids = run_query(GET_PIDS_QUERY,database)
    pid_counter = 0
    if pids:
        print(pids)
    else:
        raise ValueError('No sample found!')
    for pid in pids:

        INSERT_NEW_PID_INTO_PPID = """
            INSERT INTO pseudo_prepid
            (pseudo_prepid,prepid)
            VALUES (NULL,"{pid}")
            """
        query = INSERT_NEW_PID_INTO_PPID.format(pid=pid['PREPID'])
        run_query(query,database)
        if pid_counter == 0:
            pseudo_prepid_query = "SELECT LAST_INSERT_ID()"
            ppid = run_query(pseudo_prepid_query,database)
            pid_counter += 1
        else:
            INSERT_PID_INTO_PPID = """
                INSERT INTO pseudo_prepid
                ("{ppid"},"{pid}")
                VALUES (pseudo_prepid,prepid)
                """
            INSERT_PID_INTO_PPID.format(pid=pid['PREPID'],
                                        ppid=ppid[0]['PSEUDO_PREPID'])
            run_query(query,database)
    return ppid

def update_ppid(sample_name,sample_type,capture_kit,database):
    GET_PID_PPID_FROM_TRIPLET= """
        SELECT P_PREPID,p.PREPID
        FROM prepT p
        JOIN pseudo_prepid pp on p.prepid=pp.prepid
        WHERE CHGVID='{chgvid}'
        AND SAMPLE_TYPE='{seqtype}'
        """.format(chgvid=sample_name,
                   seqtype=sample_type)

    if sample_type != 'genome':
        GET_PID_PPID_FROM_TRIPLET += ("AND EXOMEKIT='{exomekit}'"
                                     ).format(exomekit=capture_kit)
    pid_and_ppid = run_query(GET_PID_PPID_FROM_TRIPLET,database)

    if pid_and_ppid == ():
        ppid = insert_pid_into_ppid(sample_name,sample_type,capture_kit,database)
        return ppid
    else:
        ppid_from_triplet = [x['P_PREPID'] for x in pid_and_ppid]
        prepid_from_triplet = [x['PREPID'] for x in pid_and_ppid]
        final_pids = copy(prepid_from_triplet)
        if len(set(ppid_from_triplet)) != 1:
            raise Exception('Too many ppids found from triplet:{}'.format(set(ppid_from_triplet)))
        query = GET_PREPID_FROM_PSEUDO_PREPID.format(list(set(ppid_from_triplet))[0])
        prepid_from_ppid = run_query(query,database)
        prepid_from_ppid = [x['PREPID'] for x in prepid_from_ppid]

        while len(prepid_from_triplet) > 0:
            pid = prepid_from_triplet.pop()
            try:
                prepid_from_ppid.remove(pid)
            except ValueError:
                raise Exception("Really add {} to pid?".format(pid))
                #add_pid_to_ppid(prepid,database)
        while len(prepid_from_ppid) > 0:
            pid2 = prepid_from_ppid.pop()
            raise Exception("Really remove {} from pid?".format(pid2))
            #remove_pid_from_ppid(prepid,database)

        return list(set(ppid_from_triplet))[0],final_pids


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
    run_query(INSERT_PID_INTO_PPID.format(ppid=ppid,pid=pid),database)
    run_query(UPDATE_PSEUDO_PREPID_IN_PREPT.format(ppid=ppid,pid=pid),database)

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
