#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Create a configuration file used by the Dragen system to processed IGM samples
based on their priority level in the Dragen pipeline queue
"""

import argparse
import os
import sys
import subprocess
import traceback
from create_align_config import create_align_config
from datetime import datetime
from db_statements import *
from dragen_sample import dragen_sample
from glob import glob
from utilities import *

def main(run_type_flag, debug, dontexecute, database, seqscratch_drive):
    config = get_config()
    if isinstance(run_type_flag,str) == True: #pseudo_prepid present
        pseudo_prepid = run_type_flag
        sample_name, sample_type, pseudo_prepid, capture_kit = get_next_sample(pseudo_prepid,database,debug)
        sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,database)
        run_sample(sample,dontexecute,config,seqscratch_drive,database,debug)

    elif run_type_flag == True: #Automated run
        pseudo_prepid = 0
        sample_name, sample_type, pseudo_prepid, capture_kit = get_next_sample(pseudo_prepid,database,debug)
        while sample_name is not None:
            sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,database)
            run_sample(sample,dontexecute,config,seqscratch_drive,database,debug)
            """ These two functions will need to be moved to the dragen_release
                script when available.  For now will run on every bam in scratch
                folder """

            try:
                sample_name, sample_type, pseudo_prepid, capture_kit  = get_next_sample(0,database,debug)
            except:
                print("No more samples in the queue")
                sys.exit()

def update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid,database):
    #read_group_insert = ("UPDATE Lane l "
    #                     "JOIN Flowcell f on l.FCID=f.FCID "
    #                     "SET tep_status='completed' "
    #                     "WHERE FCILLUMD='{} AND "
    #                     "LANENUM={} and "
    #                     "PREPID={} "
    #                    ).format(rg_fcillumid,rg_lane_num,rg_prepid)
    # Auto adds can_release = 1 until HTS starts approving samples
    read_group_insert = ("UPDATE Lane l "
                         "JOIN Flowcell f on l.FCID=f.FCID "
                         "SET released=1,step_status='completed' "
                         "WHERE FCILLUMID='{}' AND "
                         "LANENUM={} AND "
                         "failr1 IS NULL AND "
                         "PREPID={} "
                        ).format(rg_fcillumid,rg_lane_num,rg_prepid)
    run_query(read_group_insert,database)

def check_bam_found_vs_bam_db(sample,qualified_bams_found):
    max_prepid = max(map(lambda x:int(x),sample.metadata['prepid']))
    if is_external_or_legacy_sample(max_prepid,database) == True:
       print("Checking bams found vs RGs form fastqs")
       for laneFCID in sample.metadata['lane'][0]: #loop over read groups
            rg_lane_num,rg_fcillumid,rg_prepid = laneFCID

            bam_loc = ("{output_dir}/{sample_name}.{pseudo_prepid}.{rg_fcillumid}.{rg_lane_num}.*bam"
                      ).format(rg_fcillumid=rg_fcillumid,
                               rg_lane_num=rg_lane_num,
                               **sample.metadata)
            if glob(bam_loc) == []:
                raise ValueError("Bam {} not found!".format(bam_loc))

    else:
        print("Checking bams found vs RGs in db")
        query = GET_QUALIFIED_BAMS.format(**sample.metadata)
        qualified_bams_in_db = run_query(query,database)
        if len(qualified_bams_in_db) != len(qualified_bams_found):
            print(qualified_bams_found,qualified_bams_in_db)
            raise ValueError("Number of bams found != number of RGs in database!")
        for db_bam in qualified_bams_in_db:
            db_bam = ("{output_dir}/{sample_name}.{pseudo_prepid}.{fcillumid}.{lane_num}.bam"
                     ).format(fcillumid=db_bam['fcillumid'],lane_num=db_bam['laneNum'],**sample.metadata)
            if glob(db_bam) == []:
                raise ValueError("Bam {} not found!".format(db_bam))

def get_component_bams(sample,debug):
    qualified_bams_found = glob('{output_dir}/*bam'.format(**sample.metadata))
    check_bam_found_vs_bam_db(sample,qualified_bams_found)
    if len(qualified_bams_found) < 1:
        raise Exception("No qualified bams were found!")
    tmp_bam_list = []
    for bam in qualified_bams_found:
        bam_name = bam.split('/')[-1]
        tmp_bam_list.append(bam_name)
    component_bams = ','.join(sorted(tmp_bam_list))
    return component_bams

def write_sge_header(sample,step,script_loc):
    script = open(script_loc,'w')
    script.write("#! /bin/bash\n")
    script.write("#\n")
    script.write("#$ -S /bin/bash -cwd\n")
    script.write("#$ -j y\n")
    script.write("#$ -l mem_free=8G\n")
    script.write("#$ -l h_vmem=8G\n")
    script.write("#$ -o {log_dir}/{sample_name}.{pseudo_prepid}.{step}.out\n".format(step=step,**sample.metadata))
    script.write("#$ -e {log_dir}/{sample_name}.{pseudo_prepid}.{step}.err\n".format(step=step,**sample.metadata))
    script.write("#$ -V\n")
    script.write("#$ -N {sample_name}.{pseudo_prepid}.{step}\n".format(step=step,**sample.metadata))
    script.write("#$ -M jb3816@cumc.columbia.edu\n")
    script.write("#$ -m ea\n")
    return script


def run_sample(sample,dontexecute,config,seqscratch_drive,database,debug):
    check_Fastq_Total_Size(sample,debug)
    setup_dir(seqscratch_drive,sample,debug)
    existing_bams_check = True
    output_dir = sample.metadata['output_dir']
    pseudo_prepid = sample.metadata['pseudo_prepid']
    submit_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    existing_bams = glob("{}/*.bam".format(output_dir))
    if existing_bams == [] or existing_bams_check == False:
        update_queue(pseudo_prepid,database,debug)
        for laneFCID in sample.metadata['lane'][0]: #loop over read groups
            rg_lane_num,rg_fcillumid,rg_prepid = laneFCID
            if debug:
                print("RG info: lane {}, fcillumd {}, prepid {}".format(rg_lane_num,rg_fcillumid,rg_prepid))
            setup_first_read_RG(sample,rg_lane_num,rg_fcillumid,rg_prepid,debug)
            set_seqtime(rg_fcillumid,sample,database)
            create_align_config(sample,rg_lane_num,rg_fcillumid,rg_prepid)
            if dontexecute == False:
               run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug)
            update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid,database)

        component_bams = get_component_bams(sample,debug)
        update_dragen_metadata(sample,component_bams,database,debug)
        rm_query = "DELETE FROM {0} WHERE pseudo_prepid={1}".format("tmp_dragen",pseudo_prepid)
        if debug:
            print(rm_query)
        run_query(rm_query,database)

    else:
        print("Sample with bam files already exists!")
        sys.exit(1)

def run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug):
    output_dir = sample.metadata['output_dir']
    dragen_cmd = ['dragen', '-f', '-v', '-c', sample.metadata['conf_file'], '--watchdog-active-timeout', '600']
    if debug:
        print(' '.join(dragen_cmd))
    stderr_file_loc = ('{}/{}.{}.{}.{}.dragen.err'
                      ).format(sample.metadata['log_dir'],sample.metadata['sample_name'],
                               sample.metadata['pseudo_prepid'],rg_fcillumid,rg_lane_num)

    stdout_file_loc = ('{}/{}.{}.{}.{}.dragen.out'
                      ).format(sample.metadata['log_dir'],sample.metadata['sample_name'],
                               sample.metadata['pseudo_prepid'],rg_fcillumid,rg_lane_num)

    dragen_stderr = open(stderr_file_loc,'a')
    with open(stdout_file_loc,'a') as dragen_stdout:
        process = subprocess.Popen(dragen_cmd,stdout=subprocess.PIPE,stderr=dragen_stderr)
        while True:
            output =  process.stdout.readline().decode()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output)
                dragen_stdout.write(output)
        rc = process.poll()
    dragen_stdout.close()
    dragen_stderr.close()
    if rc != 0:
        raise Exception("Dragen alignment did not complete successfully!")
    try:
        subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
    except:
        pass

def update_dragen_metadata(sample,component_bams,database,debug):
    seqscratch_drive = sample.metadata['output_dir'].split('/')[2]
    get_dragen_metadata_query = ("SELECT * from dragen_sample_metadata "
                                 "WHERE pseudo_prepid = {pseudo_prepid} "
                                ).format(**sample.metadata)
    dbInfo = run_query(get_dragen_metadata_query,database)
    if dbInfo:
        print("Performing dragen_sample_metadata update")
        query = ("UPDATE dragen_sample_metadata "
                 "set sample_name='{sample_name}',pseudo_prepid={pseudo_prepid},"
                 "sample_type='{sample_type}',capture_kit='{capture_kit}',"
                 "priority={priority},seqscratch_drive='{seqscratch_drive}',"
                 "is_merged=0,component_bams='{component_bams}' "
                 "WHERE pseudo_prepid={pseudo_prepid}"
                ).format(seqscratch_drive=seqscratch_drive,
                         component_bams=component_bams,**sample.metadata)
    else:
        print("Performing dragen_sample_metadata insert")
        query = ("INSERT INTO dragen_sample_metadata "
                       "(sample_name,pseudo_prepid,sample_type,"
                       "capture_kit,priority,seqscratch_drive,"
                       "is_merged,component_bams) "
                       "VALUES ('{sample_name}',{pseudo_prepid},"
                       "'{sample_type}','{capture_kit}',{priority},"
                       "'{seqscratch_drive}',0,'{component_bams}'"
                       ") "
                      ).format(seqscratch_drive=seqscratch_drive,component_bams=component_bams,**sample.metadata)
    if debug:
        print(query)
    run_query(query,database)


def update_queue(pseudo_prepid,database,debug):
    insert_query = ("INSERT INTO tmp_dragen "
                    "SELECT * FROM dragen_queue WHERE pseudo_prepid={0}"
                   ).format(pseudo_prepid)
    rm_query = ("DELETE FROM {0} WHERE pseudo_prepid={1}").format("dragen_queue",pseudo_prepid)
    run_query(insert_query,database)
    run_query(rm_query,database)
    if debug:
        print(insert_query)
        print(rm_query)

def set_seqtime(rg_fcillumid,sample,database):
    query = ("SELECT FROM_UNIXTIME(Seqtime) AS SEQTIME "
            "FROM Flowcell WHERE FCIllumID='{}'").format(rg_fcillumid)
    seqtime = run_query(query,database)
    if seqtime:
        seqtime = seqtime[0]['SEQTIME'].date().isoformat() #ISO8601 format for RGDT field
    else: #In case there is not flowcell information (Ex. old and external samples)
        seqtime = '1970-1-1'
    sample.set('seqtime',seqtime)


def setup_dir(seqscratch_drive,sample,debug):
    output_dir = ('/nfs/{}/ALIGNMENT/BUILD37/DRAGEN/{sample_type}/{sample_name}.{pseudo_prepid}/'
                                    ).format(seqscratch_drive,
                                             sample_type=sample.metadata['sample_type'].upper(),
                                             sample_name=sample.metadata['sample_name'],
                                             pseudo_prepid=sample.metadata['pseudo_prepid'])
    sample.metadata['output_dir'] = output_dir
    sample.metadata['script_dir'] = sample.metadata['output_dir']+'scripts'
    sample.metadata['log_dir'] = sample.metadata['output_dir']+'logs'
    dragen_stdout = "{log_dir}/{sample_name}.{pseudo_prepid}.dragen.out".format(**sample.metadata)
    sample.metadata['dragen_stdout'] = dragen_stdout
    dragen_stderr = "{log_dir}/{sample_name}.{pseudo_prepid}.dragen.err".format(**sample.metadata)
    sample.metadata['dragen_stderr'] = dragen_stderr

    mkdir_cmd = ['mkdir','-p',sample.metadata['script_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['log_dir']]
    subprocess.call(mkdir_cmd)

def setup_first_read_RG(sample,rg_lane_num,rg_fcillumid,rg_prepid,debug):
    first_read1,first_read2 = get_first_read(sample,rg_lane_num,rg_fcillumid,debug)
    sample.set('first_fastq1',first_read1)
    sample.set('first_fastq2',first_read2)

def get_first_read(sample,rg_lane_num,rg_fcillumid,debug):
    #Using first fastq as template for all fastq.gz
    for fastqLoc in sample.metadata['fastq_loc']:
        if rg_fcillumid in fastqLoc.split('/')[-1]:
            RGfastqStr ='{0}/*L00{1}_R1_*.fastq.gz'.format(fastqLoc,rg_lane_num)
            if debug:
                print(RGfastqStr)
            RGfastqs = glob(RGfastqStr)

    #The fastqs might still be on the quantum tapes
    first_read = sorted(RGfastqs)[0]
    second_read = first_read.replace('_R1_','_R2_')
    return first_read,second_read

def user_input_fastq_size(lg_or_gt):
    if lg_or_gt == 'lt':
        userInput = input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
        if userInput == 'n':
            raise Exception("Sum of fastq files sizes is too small for a {} sample!".format(sample_type))
    elif lg_or_gt == 'gt':
        userInput = input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
        if userInput == 'n':
             raise Exception("Sum of fastq files sizes is too big for a {} sample!".format(sample_type))
    else:
        raise Exception('Unhandled lg_or_gt variable')

def check_Fastq_Total_Size(sample,debug):
    sample_type = sample.metadata['sample_type']
    fastq_filesize_sum = 0
    for fastq_loc in sample.metadata['fastq_loc']:
        for fastq in glob(fastq_loc + '/*gz'):
            fastq_filesize = os.stat(os.path.realpath(fastq)).st_size
            fastq_filesize_sum += fastq_filesize
    print("Sum of Fastq size: {}".format(fastq_filesize_sum))
    if sample_type == 'genome':
        if fastq_filesize_sum < 42949672960: # < 40GB
            user_input_fastq_size('lt')
    elif sample_type == 'exome':
        if fastq_filesize_sum > 32212254720: # > 30GB
            user_input_fastq_size('gt')
        elif fastq_filesize_sum < 1073741824: # < 1GB
            user_input_fastq_size('lt')
    elif sample_type == 'rnaseq':
        if fastq_filesize_sum > 32212254720: # > 30GB
            user_input_fastq_size('gt')
        elif fastq_filesize_sum < 1073741824: # < 1GB
            user_input_fastq_size('lt')
    elif sample_type == 'custom_capture':
        if fastq_filesize_sum > 10737418240: # > 10GB
            user_input_fastq_size('gt')
    else:
        raise Exception('Unhandled sample_type found: {}!'.format(sample_type))
    if debug:
        print('Fastq size',fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024)

def get_reads(sample,read_number,debug):
    if debug:
        print(sample.metadata['fastq_loc'])
    fastq_loc = sample.metadata['fastq_loc'][0]
    if debug:
        print('{0}/*L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number))
    read = glob('{0}/*L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number))
    #The fastqs might still be on the quantum tapes
    if read == []:
        raise Exception("Fastq file not found!")
    else:
        return sorted(read)[0]

def get_next_sample(pseudo_prepid,database,debug):
    if pseudo_prepid == 0:
        query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit "
            "FROM dragen_queue "
            "WHERE PRIORITY < 99 "
            "ORDER BY PRIORITY ASC LIMIT 1 ")
        sample_info = run_query(query,database)
    else:
        query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit "
            "FROM dragen_queue "
            "WHERE pseudo_prepid={}"
            ).format(pseudo_prepid)
        sample_info = run_query(query,database)
    if debug:
        print(query)
        print('Dragen_queue info: {0}'.format(sample_info))
    print(sample_info)
    return sample_info[0]['sample_name'],sample_info[0]['sample_type'],sample_info[0]['pseudo_prepid'],sample_info[0]['capture_kit']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--pseudo_prepid", type=str,
                        help="Run dragen in single sample mode with provided pseudo_prepid")
    group.add_argument("-a", "--auto", default=False, action="store_true",
                        help="Run pipeline in automated mode")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    parser.add_argument("--seqscratch_drive", default='seqscratch_ssd',
                        action='store', help="Set output destination")
    parser.add_argument("--dontexecute", default=False, action="store_true",
                        help="Perform setup but do not start a Dragen run")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', '-v', action='version',
                        version='%(prog)s v2.0')
    args=parser.parse_args()

    if args.pseudo_prepid:
        run_type_flag = args.pseudo_prepid
    else:
        run_type_flag = args.auto
    if args.test:
        database="testDB"
    else:
        database="sequenceDB"

    try:
        main(run_type_flag, args.debug, args.dontexecute, database, args.seqscratch_drive)
    except Exception:
        tb = traceback.format_exc()
        print(tb)
        emailAddresses = ['jb3816@cumc.columbia.edu']
        emailCmd = ('echo "Dragen Pipe failure" | mail -s "Dragen Pipe Failure" {}'
                   ).format(' '.join(emailAddresses))
        print(emailCmd)
        os.system(emailCmd)
