import argparse
import configparser
import logging
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

    if args.fcillumid:
        auto_release_flag = False
        get_sample_info_from_flowcell = ("SELECT DISTINCT CHVID,p.PREPID,SEQTYPE,p.EXOMEKIT"
                                         "FROM Lane l "
                                         "JOIN Flowcell f ON f.fcid=l.fcid "
                                         "JOIN SeqType st on st.prepid=l.prepid "
                                         "JOIN prepT p ON p.prepid=l.prepid "
                                         "WHERE FCILLUMID='{}'"
                                        ).format(args.fcillumid)

        print("Starting sample release for Flowcell: {}".format(args.fcillumid))
        rejected_samples = []
        for sample_info in run_query(get_sample_info_from_flowcell):
            seqtype = sample_info['SeqType']
            sample_name = sample_info['CHVID']
            if seqtype.lower() == 'exome' or seqtype.lower() == 'custom_capture':
                check_yield(sample_name,fcillumid)
                rejected_samples = run_sample(sample_name,seqtype,
                                              sample_info['EXOMEKIT'],
                                              auto_release_flag,rejected_samples,
                                              database)
        email(fcillumid,rejected_samples,config)

    else:
        run_sample(args.sample_name,args.sample_type,args.capture_kit,auto_release_flag,database)

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
    os.system("{} -s 'Sample Release' {} < email.tmp".format(email_program,failure_email_addy)

    os.system('rm release_email.txt')


def update_dragen_queue(ppid,database):
    table = 'tmp_dragen'
    queue_query = """SELECT PSEUDO_PREPID
                     FROM {table}
                     WHERE PSEUDO_PREPID = {ppid}
                  """
    in_tmp_dragen = run_query(queue_query.format(table=table,ppid=ppid),database)
    if in_tmp_dragen:
        insert_query = "INSERT INTO dragen_queue " + queue_query
    else:
        table = 'dragen_queue'
        in_dragen_queue = run_query(queue_query.format(table=table,ppid=ppid),database)
        if not in_dragen_queue:
            insert_query = """INSERT INTO TABLE dragen_queue
                              (PSEUDO_PREPID,SAMPLE_NAME,SAMPLE_TYPE,
                              CAPTURE_KIT,PRIORITY)
                              VALUES ({ppid},'{sample_name}','{sample_type}',
                              '{capture_kit}',{priority})
                           """.format(sample_name=sample_name,
                                      sample_type=sample_type,
                                      capture_kit=capture_kit,
                                      priority=40)
            run_query(insert_query,database)

def run_sample(sample_name,sample_type,capture_kit,auto_release_flag,rejected_samples,database)
    ppid,pids = update_ppid(sample_name,sample_type,capture_kit,database)
    try:
        print("Starting sample release for Flowcell: {},{},{}".format(sample_name,
                                                                      sample_type,
                                                                      capture_kit))
        sample = dragen_sample(sample_name,sample_type,ppid
                               capture_kit,database)
        sample_status_query = """ SELECT MAX(PIPELINE_STEP_ID) as PIPELINE_STEP_ID
                                  FROM DRAGEN_PIPELINE_STEP
                                  WHERE STEP_STATUS='complete' AND
                                  PSEUDO_PREPID={ppid}
                              """.format(ppid=ppid)
        sample_status = run_query(sample_status_query,database)[0]['PIPELINE_STEP_ID']
        if sample_status:
            raise ValueError("Sample has already been run.  Cleanup required first")
        else:
            #every external sample should only have one pids
            is_external = run_query(IS_SAMPLE_EXTERNAL_FROM_PREPID.format(prepid=pids[0])
            if is_external == True and pids[0] < 20000:
                pass
            else:
                rejected_samples = check_yield(ppid,sample_type,
                                               auto_release_flag,
                                               rejected_samples,database)
            update_dragen_queue(ppid,database)
            
            ##update statusT
            ##update prepT.status
            
            return rejected_samples

    except FastqError:
        pass

def check_yield(ppid,sample_type,auto_release_flag,rejected_samples,database):
   lane_yield_sum = run_query(GET_YIELD_FROM_PPID.format(ppid=ppid),database)
   lane_yield_sum = int(lane_yield_sum)

   if sample_type.lower() == 'exome' and lane_yield_sum > 7500:
       return rejected_samples
   elif sample_type.lower() == 'genome' and lane_yield_sum > 100000:
       return rejected_samples
   elif sample_type.lower() == 'custom_capture' and lane_yield_sum > 150:
       return rejected_samples
   else:
       if auto_release_flag == True:
           pass
       else:

       rejected_samples.append(ppid)
       return rejected_samples

def insert_pid_into_ppid(sample_name,sample_type,capture_kit)
    GET_PIDS_QUERY = """SELECT PSEUDO_PREPID,PREPID
                        FROM prepT p
                        JOIN pseudo_prepid pp on p.prepid=pp.prepid
                        WHERE CHGVID='{chgvid}'
                        AND SEQTYPE='{seqtype}'
                        AND EXOMEKIT='{exomekit}'
                     """
    pids = run_query(GET_PIDS_QUERY.format(chgvid=sample_name,
                                           seqtype=sample_type,
                                           exomekit=capture_kit),
                     database)
    pid_counter = 0
    for pid in pids:
        query = INSERT_PID_INTO_PPID.format(pid[0]['PREPID']
        run_query(query,database)
        if pid_counter = 0:
            pseudo_prepid_query = "SELECT LAST_INSERT_ID()"
            ppid = run_query(pseudo_prepid_query,database)
    return ppid

def update_ppid(sample_name,sample_type,capture_kit,database):
    pid_and_ppid = run_query(GET_PID_PPID_FROM_TRIPLET.format(chgvid=sample_name,
                                                              seqtype=sample_type,
                                                              exomekit=capture_kit),
                             database)
    if pid_and_ppid == ():
        ppid = insert_pid_into_ppid(sample_name,sample_type,capture_kit)
        return ppid
    else:
        ppid_from_triplet = [x['PSEUDO_PREPID'] from x in pid_and_ppid]
        prepid_from_triplet = [x['PREPID'] from x in pid_and_ppid]
        final_pids = copy.copy(prepid_from_triplet)
        if len(set(ppid_from_triplet)) != 1:
            raise Exception('Too many ppids found from triplet:{}'.format(set(ppid_from_triplet))
        prepid_from_ppid = run_query(GET_PREPID_FROM_PSEUDO_PREPID.format(set(ppid_from_triplet))[0])
        prepid_from_ppid = [x['PREPID'] from x in prepid_from_ppid]


        while len(prepid_from_triplet) > 0:
            pid = prepid_from_triplet.pop()
            try:
                prepid_from_ppid.remove(pid)
            except ValueError:
                raise Exception("Really add {} to ppid?".format(prepid)
                #add_pid_to_ppid(prepid,database)
        while len(prepid_from_ppid) > 0:
            pid = prepid_from_ppid.pop()
            raise Exception("Really remove {} from ppid?".format(prepid2)
            #remove_pid_from_ppid(prepid,database)

        return set(ppid_from_triplet),final_pids

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
    group.add_argument("-f", "--fcillumid", dest='fcillumid',
                        help="Specify Illumina's Flowcell ID")
    group.add_argument('-s', '--sample_name',dest='sample_name',
                        help="Specify sample name")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-t','--sample_type', dest='sample_type',
                        help="Specify sample type (exome,genome,etc)")
    parser.add_argument('-k','--capture_kit', dest='capture_kit',
                        help="Specify capture kit ")
    parser.add_argument("-p", "--priority", default=99)
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    if args.sample_name and (bool(args.sample_type) == False or  bool(args.capture_kit) == False) :
        parser.error('--sample_name and --sample_type and --capture_kit must be given together')
    return args

if __name__ == '__main__':
    main()
