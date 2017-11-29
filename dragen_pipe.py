#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Create a configuration file used by the Dragen system to processed IGM samples
based on their priority level in the Dragen pipeline queue
"""

import argparse
import MySQLdb
import os
import sys
import subprocess
import traceback
from create_align_config import create_align_config
from ConfigParser import SafeConfigParser
from datetime import datetime
from db_statements import *
from dragen_sample import dragen_sample
from glob import glob


def main(run_type_flag, debug, dontexecute, seqscratch_drive):
    CNF = "/nfs/goldstein/software/dragen_pipe/resources/dragen.cnf" # defaults file for pipeline
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    if isinstance(run_type_flag,str) == True: #pseudo_prepid present
        pseudo_prepid = run_type_flag
        sample_name, sample_type, pseudo_prepid, capture_kit = get_next_sample(pseudo_prepid,debug)
        sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,get_curs())
        setup_dir(sample,debug)
        run_sample(sample,dontexecute,config_parser,debug)

    elif run_type_flag == True: #Automated run
        pseudo_prepid = 0
        sample_name, sample_type, pseudo_prepid, capture_kit = get_next_sample(pseudo_prepid,debug)
        while sample_name is not None:
            sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,get_curs())
            sample.metadata['output_dir'] = ('/nfs/{}/ALIGNMENT/BUILD37/DRAGEN/{}/{}.{}/'
                                            ).format(seqscratch_drive,sample_type.upper(),
                                                     sample_name,pseudo_prepid)
            setup_dir(sample,debug)
            run_sample(sample,dontexecute,config_parser,debug)
            """ These two functions will need to be moved to the dragen_release
                script when available.  For now will run on every bam in scratch
                folder """

            try:
                sample_name, sample_type, pseudo_prepid, capture_kit  = get_next_sample(0,debug)
            except:
                print "No more samples in the queue"
                sys.exit()

def get_curs():
    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
        read_default_file='~/.my.cnf')

    curs = db.cursor()
    return curs

def update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid):
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
                         "LANENUM={} and "
                         "PREPID={} "
                        ).format(rg_fcillumid,rg_lane_num,rg_prepid)
    #print read_group_insert
    run_query(read_group_insert)

def check_bam_found_vs_bam_db(sample,qualified_bams_found):
    if is_external_or_legacy_sample(sample) == True:
        #skips db vs found bam sanity check
        pass
    else:
        qualified_bams_in_db = run_query(GET_QUALIFIED_BAMS.format(**sample.metadata))
        if len(qualified_bams_in_db) != len(qualified_bams_found):
            #print qualified_bams_found,qualified_bams_in_db
            raise ValueError, "Number of bams found != number of RGs in database!"
        for db_bam in qualified_bams_in_db:
            db_bam = ("{output_dir}/{sample_name}.{pseudo_prepid}.{fcillumid}.{lane_num}.bam"
                     ).format(fcillumid=db_bam[0],lane_num=db_bam[1],**sample.metadata)
            if glob(db_bam) == []:
                raise ValueError, "Bam {} not found!".format(db_bam)

def get_component_bams(sample,debug):
    qualified_bams_found = glob('{output_dir}/*bam'.format(**sample.metadata))
    check_bam_found_vs_bam_db(sample,qualified_bams_found)
    if len(qualified_bams_found) < 1:
        raise Exception, "No qualified bams were found!"
    component_bams = ','.join(sorted(qualified_bams_found))
    return component_bams

def is_external_or_legacy_sample(sample):
    """In cases of legecy samples (prepid < 20000) or external samples. There are no cases of legecy samples with
    multiple preps"""
    if max(map(lambda x:int(x),sample.metadata['prepid'])) < 20000:
        print "Sample has a prepID < 20000"
        return True

    else:
        is_external_sample = int(run_query(IS_SAMPLE_EXTERNAL_FROM_PREPID.format(prepid=sample.metadata['prepid'][0]))[0][0])
        if is_external_sample == 1:
            print "Sample is an external sample!"
            return True
        elif is_external_sample == 0:
            print "Sample is not an external sample"
            return False
        else: #Returns None
            raise ValueError, "No value found.  Does prepID exist?"


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

def run_query(query):
    database='sequenceDB'
    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
        read_default_file='~/.my.cnf')
    curs = db.cursor()
    curs.execute(query)
    results = curs.fetchall()
    db.commit()
    db.close()
    return results

def run_sample(sample,dontexecute,config_parser,debug):
    output_dir = sample.metadata['output_dir']
    pseudo_prepid = sample.metadata['pseudo_prepid']
    submit_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    existing_bams = glob("{}/*.bam".format(output_dir))
    if existing_bams == []:
        update_queue(pseudo_prepid,debug)
        for laneFCID in sample.metadata['lane'][0]: #loop over read groups
            rg_lane_num,rg_fcillumid,rg_prepid = laneFCID
            if debug:
                print "RG info: lane {}, fcillumd {}, prepid {}".format(rg_lane_num,rg_fcillumid,rg_prepid)
            setup_first_read_RG(sample,rg_lane_num,rg_fcillumid,rg_prepid,debug)
            set_seqtime(rg_fcillumid,sample)
            create_align_config(sample,rg_lane_num,rg_fcillumid,rg_prepid)
            if dontexecute == False:
               run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug)
            update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid)
        rm_query = "DELETE FROM {0} WHERE pseudo_prepid={1}".format("tmp_dragen",pseudo_prepid)
        if debug:
            print rm_query
        run_query(rm_query)

        component_bams = get_component_bams(sample,debug)
        update_dragen_metadata(sample,component_bams,debug)

    else:
        print ("Sample {sample_name} bam file already exists!"
              ).format(sample_name=sample.metadata['sample_name'])
        sys.exit(1)

def run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug):
    output_dir = sample.metadata['output_dir']
    dragen_cmd = ['dragen', '-f', '-v', '-c', sample.metadata['conf_file'], '--watchdog-active-timeout', '30000']
    if debug:
        print ' '.join(dragen_cmd)
    stderr_file_loc = ('{}/{}.{}.{}.{}.dragen.err'
                      ).format(sample.metadata['log_dir'],sample.metadata['sample_name'],
                               sample.metadata['pseudo_prepid'],rg_fcillumid,rg_lane_num)

    stdout_file_loc = ('{}/{}.{}.{}.{}.dragen.out'
                      ).format(sample.metadata['log_dir'],sample.metadata['sample_name'],
                               sample.metadata['pseudo_prepid'],rg_fcillumid,rg_lane_num)

    dragen_stderr = open(stdout_file_loc,'a')
    with open(stdout_file_loc,'a') as dragen_stdout:
        process = subprocess.Popen(dragen_cmd,stdout=subprocess.PIPE,stderr=dragen_stderr)
        for line in iter(process.stdout.readline, ''):
            if debug:
                sys.stdout.write(line)
            dragen_stdout.write(line)
        process.communicate()
        error_code = process.wait()
        subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
    if error_code != 0:
        raise Exception, "Dragen alignment did not complete successfully!"
    dragen_stdout.close()
    dragen_stderr.close()

def update_dragen_metadata(sample,component_bams,debug):
    seqscratch_drive = sample.metadata['output_dir'].split('/')[2]
    get_dragen_metadata_query = ("SELECT * from dragen_sample_metadata "
                                 "WHERE pseudo_prepid = {pseudo_prepid} "
                                ).format(**sample.metadata)
    dbInfo = run_query(get_dragen_metadata_query)
    if dbInfo:
        print "Performing dragen_sample_metadata update"
        query = ("UPDATE dragen_sample_metadata "
                 "set sample_name='{sample_name}',pseudo_prepid={pseudo_prepid},"
                 "sample_type='{sample_type}',capture_kit='{capture_kit}',"
                 "priority={priority},seqscratch_drive='{seqscratch_drive}',"
                 "is_merged=0,component_bams='{component_bams}' "
                 "WHERE pseudo_prepid={pseudo_prepid}"
                ).format(seqscratch_drive=seqscratch_drive,
                         component_bams=component_bams,**sample.metadata)
    else:
        print "Performing dragen_sample_metadata insert"
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
        print query
    run_query(query)

def get_pipeline_version():
    version = subprocess.check_output(
        ["/nfs/goldstein/software/git-2.5.0/bin/git", "describe", "--always"]).strip()
    if version:
        return version
    else:
            raise ValueError("Could not get the version # of the pipeline; "
                             "Check cwd?")

def update_queue(pseudo_prepid,debug):
    insert_query = ("INSERT INTO tmp_dragen "
                    "SELECT * FROM dragen_queue WHERE pseudo_prepid={0}"
                   ).format(pseudo_prepid)
    rm_query = ("DELETE FROM {0} WHERE pseudo_prepid={1}").format("dragen_queue",pseudo_prepid)
    run_query(insert_query)
    run_query(rm_query)
    if debug:
        print insert_query
        print rm_query

def set_seqtime(rg_fcillumid,sample):
    query = ("SELECT FROM_UNIXTIME(Seqtime) "
            "FROM Flowcell WHERE FCIllumID='{}'").format(rg_fcillumid)
    seqtime = run_query(query)

    if seqtime:
        seqtime = seqtime[0][0].date().isoformat() #ISO8601 format for RGDT field
    else: #In case there is not flowcell information (Ex. old and external samples)
        seqtime = '1970-1-1'
    sample.set('seqtime',seqtime)

def check_dir_contents(sample,debug):
    bam_bai_files = glob('{output_dir}/*.ba[im]'.format(**sample.metadata))
    if len(bam_bai_files) > 1:
        raise ValueError, "Previous bam/bai files were found!"

def setup_dir(sample,debug):
    mkdir_cmd = ['mkdir','-p',sample.metadata['script_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['log_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['fastq_dir']]
    subprocess.call(mkdir_cmd)
    #check_dir_contents(sample,debug)
    """Removes any fastq.gz files in the fastq folder. Symlink-only fastqs
        should be in this folder.  This is just in case the sample was run
        through once"""
    symLinkedFastqFiles = glob(str(sample.metadata['fastq_dir']) + '/*fastq.gz')
    for file in symLinkedFastqFiles:
        #remove any symlinked fastq.gz files
        try:
            os.remove(file)
        except OSError:
            pass

    first_read1 = get_reads(sample,1,debug)
    first_read2 = get_reads(sample,2,debug)
    fastq_counter=0
    for fastq_loc in sample.metadata['fastq_loc']:
        fastqs = glob(fastq_loc + '/*_R1_*fastq.gz')
        for fastq in fastqs:
            if 'fastq16-rsync' in fastq:
                pass
            else:
                fastq_counter+=1
                new_fastq_read1 = first_read1.split('/')[-1].replace(
                        '001.fastq.gz','{0:03d}.fastq.gz'.format(fastq_counter))
                new_fastq_read2 = first_read2.split('/')[-1].replace(
                        '001.fastq.gz','{0:03d}.fastq.gz'.format(fastq_counter))
                ln_cmd1 = ['ln','-s',fastq,sample.metadata['fastq_dir']+'/'+new_fastq_read1]
                ln_cmd2 = ['ln','-s',fastq.replace('_R1_','_R2_'),sample.metadata['fastq_dir']+'/'+new_fastq_read2]

                if fastq_counter == 1:
                    sample.set('first_fastq1',sample.metadata['fastq_dir']+'/'+new_fastq_read1)
                    sample.set('first_fastq2',sample.metadata['fastq_dir']+'/'+new_fastq_read2)
                    sample.set('first_lane',sample.metadata['lane'][0][0][0])
                    sample.set('first_flowcell',sample.metadata['lane'][0][0][1])

                subprocess.call(ln_cmd1)
                subprocess.call(ln_cmd2)
    check_Fastq_Total_Size(sample,debug)

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
                print RGfastqStr
            RGfastqs = glob(RGfastqStr)

    #The fastqs might still be on the quantum tapes
    first_read = sorted(RGfastqs)[0]
    second_read = first_read.replace('_R1_','_R2_')
    return first_read,second_read


def check_Fastq_Total_Size(sample,debug):
    fastq_dir = sample.metadata['fastq_dir']
    sample_type = sample.metadata['sample_type']
    fastq_filesize_sum = 0
    for fastq in glob(fastq_dir + '/*gz'):
        fastq_filesize = os.stat(os.path.realpath(fastq)).st_size
        fastq_filesize_sum += fastq_filesize
    print "Sum of Fastq size: {}".format(fastq_filesize_sum)
    if sample_type == 'genome':
        if fastq_filesize_sum < 42949672960: # < 40GB
            print fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024
            userInput = raw_input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
            if userInput == 'n':
                raise Exception, "Sum of fastq files sizes is too small for a {} sample!".format(sample_type)
    elif sample_type == 'exome':
        if fastq_filesize_sum > 32212254720: # > 30GB
            print fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024
            userInput = raw_input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
            if userInput == 'n':
                 raise Exception, "Sum of fastq files is too big for a {} sample!".format(sample_type)
        elif fastq_filesize_sum < 1073741824: # < 1GB
            userInput = raw_input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
            if userInput == 'n':
                 raise Exception, "Sum of fastq files is too small for a {} sample!".format(sample_type)
    elif sample_type == 'rnaseq':
        if fastq_filesize_sum > 32212254720: # > 30GB
            print fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024
            userInput = raw_input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
            if userInput == 'n':
                 raise Exception, "Sum of fastq files sizes is too big for a {} sample!".format(sample_type)
        elif fastq_filesize_sum < 1073741824: # < 1GB
            userInput = raw_input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
            if userInput == 'n':
                 raise Exception, "Sum of fastq files is too small for a {} sample!".format(sample_type)
    elif sample_type == 'custom_capture':
        if fastq_filesize_sum > 10737418240: # > 10GB
            print fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024
            userInput = raw_input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
            if userInput == 'n':
                 raise Exception, "Sum of fastq files sizes is too big for a {} sample!".format(sample_type)
    else:
        raise Exception, 'Unhandled sample_type found: {}!'.format(sample_type)
    if debug:
        print 'Fastq size',fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024

def get_reads(sample,read_number,debug):
    if debug:
        print sample.metadata['fastq_loc']
    fastq_loc = sample.metadata['fastq_loc'][0]
    if debug:
        print '{0}/*L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number)
    read = glob('{0}/*L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number))
    #The fastqs might still be on the quantum tapes
    if read == []:
        raise Exception, "Fastq file not found!"
    else:
        return sorted(read)[0]

def get_next_sample(pseudo_prepid,debug):
    if pseudo_prepid == 0:
        query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit "
            "FROM dragen_queue "
            "WHERE PRIORITY < 99 "
            "ORDER BY PRIORITY ASC LIMIT 1 ")
        sample_info = run_query(query)
    else:
        query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit "
            "FROM dragen_queue "
            "WHERE pseudo_prepid={}"
            ).format(pseudo_prepid)
        sample_info = run_query(query)
    if debug:
        print query
        print 'Dragen_queue info: {0}'.format(sample_info)
    print sample_info
    return sample_info[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--pseudo_prepid", type=str, 
                        help="Run dragen in single sample mode with provided pseudo_prepid")
    group.add_argument("-a", "--auto", default=False, action="store_true",
                        help="Run pipeline in automated mode")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    parser.add_argument("--seqscratch_drive", default='seqscratch_ssd', action="store_true",
                        help="Set output destination")
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
        args.database="testDB"
    else:
        args.database="sequenceDB"
    global database
    database = args.database

    try:
        main(run_type_flag, args.debug, args.dontexecute, args.seqscratch_drive)
    except Exception:
        tb = traceback.format_exc()
        print tb
        emailAddresses = ['jb3816@cumc.columbia.edu']
        emailCmd = ('echo "Dragen Pipe failure" | mail -s "Dragen Pipe Failure" {}'
                   ).format(' '.join(emailAddresses))
        print emailCmd
        os.system(emailCmd)
