#!/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python
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
import socket
import smtplib,email,email.encoders,email.mime.text,email.mime.base
from email.mime.multipart import MIMEMultipart as MM
import pprint 

def emailit(rarp,arse):
    smtpserver = 'localhost'
    fromAddr = 'dh2880@cumc.columbia.edu'
    emailMsg = MM('alternative') #email.MIMEMultipart.MIMEMultipart('alternative')
    emailMsg['Subject'] = rarp
    emailMsg['From'] = fromAddr
    to=['dsth@cantab.net','dh2880@cumc.columbia.edu']
    emailMsg['To'] = ', '.join(to)
    body='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '
    body +='"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"><html xmlns="http://www.w3.org/1999/xhtml">'
    body +='<body style="font-size:12px;font-family:Verdana"><PRE>'
    body += arse
    body += '</PRE></body></html>'
    emailMsg.attach(email.mime.text.MIMEText(body,'html'))
    server = smtplib.SMTP(smtpserver)
    server.sendmail(fromAddr,to,emailMsg.as_string())
    server.quit()

def main(run_type_flag, debug, dontexecute, database, seqscratch_drive):

    os.system("dragen_reset")

    config = get_config()
    if isinstance(run_type_flag,str) == True: #pseudo_prepid present
        pseudo_prepid = run_type_flag
        sample_name, sample_type, pseudo_prepid, capture_kit = get_next_sample(pseudo_prepid,database,debug)
        sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,database)
        try:
            run_sample(sample,dontexecute,config,seqscratch_drive,database,debug)
        except:
            tb = traceback.format_exc()
            print(tb)
            status='Dragen Alignment Failure'
            update_status(sample,status,database)
            emailit('alignment issue:2 ' + sample_name,pprint.pformat(vars(sample)))
            # sick of all the emails!?!
            # email_failure()

#############################################################################
    elif run_type_flag == True: #Automated run
#############################################################################

        pseudo_prepid = 0
        sample_name, sample_type, pseudo_prepid, capture_kit = get_next_sample(pseudo_prepid,database,debug)

###### USE STADNARD PYTHON EMAIL!?!?!
###### USE STADNARD PYTHON EMAIL!?!?!
###### USE STADNARD PYTHON EMAIL!?!?!

        # should set up seqscratch but atm fastq_temp & fastq_temp2 are set up
        # if sample_type == "Genome":
            # print("UPDATING DIR FOR GENOMES!?!")
            # seqscratch_drive = "fastq_temp"
            # seqscratch_drive = "fastq_temp2"

        while sample_name is not None:

            print("running pp= {}".format(pseudo_prepid))
            ##### this is the silly, ott bit so lock first and update to avoid cycling forever with probs...?!?
            ##### - in the end locked in get_next_sample

            sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,database)

            try:
                run_sample(sample,dontexecute,config,seqscratch_drive,database,debug)
            except:
                tb = traceback.format_exc()
                print(tb)
                status='Dragen Alignment Failure'
                update_status(sample,status,database)
                # sick of all the emails!?!
                emailit('alignment issue:1 ' + sample_name,pprint.pformat(vars(sample)))
                # email_failure()
                sys.exit()

            try:
                sample_name, sample_type, pseudo_prepid, capture_kit  = get_next_sample(0,database,debug)
            except:
                print("No more samples in the queue")
                sys.exit()
#############################################################################

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
       print("Checking bams found vs RGs from fastqs")
       for laneFCID in sample.metadata['lane'][0]: #loop over read groups
            rg_lane_num,rg_fcillumid,rg_prepid = laneFCID

            bam_loc = ("{output_dir}/{sample_name}.{pseudo_prepid}.{rg_fcillumid}.{rg_lane_num}.*bam"
                      ).format(rg_fcillumid=rg_fcillumid,
                               rg_lane_num=rg_lane_num,
                               **sample.metadata)
            if glob(bam_loc) == []:
                raise Exception("Bam {} not found!".format(bam_loc))

    else:
        print("Checking bams found vs RGs in db")
        query = GET_QUALIFIED_BAMS.format(**sample.metadata)
        qualified_bams_in_db = run_query(query,database)
        if len(qualified_bams_in_db) != len(qualified_bams_found):
            print(qualified_bams_found,qualified_bams_in_db)
            raise Exception("Number of bams found != number of RGs in database!")
        for db_bam in qualified_bams_in_db:
            db_bam = ("{output_dir}/{sample_name}.{pseudo_prepid}.{fcillumid}.{lane_num}.bam"
                     ).format(fcillumid=db_bam['fcillumid'],lane_num=db_bam['laneNum'],**sample.metadata)
            if glob(db_bam) == []:
                raise Exception("Bam {} not found!".format(db_bam))

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
    script.write("#$ -M {}@cumc.columbia.edu\n".format(get_user()))
    script.write("#$ -m ea\n")
    return script


def run_sample(sample,dontexecute,config,seqscratch_drive,database,debug):

    ###### should really update these to fastq error!?!
    check_Fastq_Total_Size(sample,debug)
    setup_dir(seqscratch_drive,sample,debug)
    existing_bams_check = False
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

            if dontexecute == False: #For test purposes
               run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug)

            update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid,database)

        component_bams = get_component_bams(sample,debug)
        update_dragen_metadata_prepT_status(sample,component_bams,database,pseudo_prepid,debug)
        # rm_query = "DELETE FROM {0} WHERE pseudo_prepid={1}".format("tmp_dragen",pseudo_prepid)
        # if debug:
        #    print(rm_query)
        # run_query(rm_query,database)

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
        ##### this we do want an email about!?!
        email_failure(dragen_cmd)
        emailit('alignment issue ' + stderr_file_loc,dragen_cmd)
        raise Exception("Dragen alignment did not complete successfully!")
    try:
        subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
    except:
        pass

def update_dragen_metadata_prepT_status(sample,component_bams,database,pseudo_prepid,debug):

    seqscratch_drive = sample.metadata['output_dir'].split('/')[2]
    get_dragen_metadata_query = ("SELECT * from dragen_sample_metadata WHERE pseudo_prepid = {} ").format(pseudo_prepid)
    dbInfo = run_query(get_dragen_metadata_query,database)

    ###### since we don't actually check for consistency shouldn't it just use replace into?!?
    ###### but more fundamentally, they really shouldn't already exist as ths 'should' be the only thing populating dsm?!?
    if dbInfo:
        print("Performing dragen_sample_metadata update")
        query = ("UPDATE dragen_sample_metadata set "
                 " sample_name='{sample_name}', "
                 " seqscratch_drive='{seqscratch_drive}',"
                 " is_merged=0,component_bams='{component_bams}' "
                 " WHERE pseudo_prepid={pseudo_prepid}"
                ).format(seqscratch_drive=seqscratch_drive,component_bams=component_bams,**sample.metadata)
    else:
        raise ValueError("pseudo_prepid {} is missing from dragen_sample_metadata".format(pseudo_prepid));

    print(query)
    run_query(query,database)
    status='Dragen Alignment Completed'
    update_status(sample,status,database)

def update_status(sample,status,database):

    prepT_query = """UPDATE prepT
                     SET status='{status}',status_time=unix_timestamp()
                     WHERE p_prepid={pseudo_prepid}
                  """.format(status=status,**sample.metadata)
    run_query(prepT_query,database)

    prepT_public_query = """UPDATE prepT_public
                            SET prepT_status='{status}',
                            prepT_status_time=unix_timestamp()
                            WHERE p_prepid={pseudo_prepid}
                         """.format(status=status,**sample.metadata)
    run_query(prepT_public_query,database)

    userID = get_user_id(database)
    statusT_insert = """INSERT INTO statusT
                        (CHGVID,STATUS,STATUS_TIME,DBID,PREPID,USERID,POOLID,SEQID,PLATENAME)
                        VALUES ('{sample_name}','{status}',unix_timestamp(),
                                {dbid},{prepid},{userID},0,0,'')
                     """.format(userID=userID,status=status,
                                prepid=sample.metadata['prepid'][0],
                                dbid=sample.metadata['dbid'],
                                sample_name=sample.metadata['sample_name'])
    run_query(statusT_insert,database)

########## should attempt to lock sample immediately at selection?!?
def update_queue(pseudo_prepid,database,debug):

####################### need to make it change the value so that we don't end up in a permanent loop - i.e. in def get_next_sample(pid,database,debug):
####################### need to make it change the value so that we don't end up in a permanent loop - i.e. in def get_next_sample(pid,database,debug):
####################### need to make it change the value so that we don't end up in a permanent loop - i.e. in def get_next_sample(pid,database,debug):

    who=socket.gethostname()
    state=0;
    if who == "dragen1.igm.cumc.columbia.edu":
        state=80011
    elif who == "dragen2.igm.cumc.columbia.edu":
        state=80012
    else:
        emailit('alignment issue','wrong server')
        raise ValueError("{} is not allowed to run this".format(who))

    connection = get_connection(database)
    print(" db= {} and pp= {}".format(database,pseudo_prepid))

    ##### make sure dragen isn't actually running too!?!?
    with connection.cursor() as cur:

        ###### clearly anything with this state is in an error state if we're running!?!
        ###### also check for conf files etc. - i.e. really make sure not getting races?!?
        q="update dragen_sample_metadata set is_merged = {} WHERE is_merged = {}".format( state+10, state )
        print(q)
        cur.execute(q)
        print("we marked {} problem samples".format(cur.rowcount))

        q="update dragen_sample_metadata set is_merged = {} WHERE pseudo_prepid = {} and is_merged = 80010".format(state,pseudo_prepid) 
        print(q)
        cur.execute(q)
        if cur.rowcount != 1:
            raise ValueError("seems we couldn't get a lock on sample {}".format(pseudo_prepid));
        else:
            print("we likely have a lock")

    connection.commit()
    connection.close()

    # rm_query = ("DELETE FROM dragen_queue WHERE pseudo_prepid={}").format(pseudo_prepid)
    # print(rm_query)
    # run_query(rm_query,database)

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

def user_input_fastq_size(lg_or_gt,sample_type):
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

    gb=(1024*1024*1024);
    sample_type = sample.metadata['sample_type']
    print("> have '{}'".format(sample_type))
    fastq_filesize_sum = 0
    for fastq_loc in sample.metadata['fastq_loc']:
        print(" > checking location '{}'".format(fastq_loc))
        for fastq in glob(fastq_loc + '/*gz'):
            fastq_filesize = os.stat(os.path.realpath(fastq)).st_size
            fastq_filesize_sum += fastq_filesize
            print("  > checking fastq '{}' ({}GB)".format(fastq,'%.2f'%(fastq_filesize/gb)))
    print("Sum of Fastq size: {}".format(fastq_filesize_sum))
    if sample_type == 'genome':
        if fastq_filesize_sum < 42949672960: # < 40GB
            raise ValueError('\n\nnot sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('lt',sample.metadata['sample_type'])
    elif sample_type == 'exome':
        #### why would it stall?!? at most it should register and warning and move on not block entire procedure for user input!?!
        # if fastq_filesize_sum > 32212254720: # > 30GB
        if fastq_filesize_sum > (150*gb): # > 30GB
            # user_input_fastq_size('gt',sample.metadata['sample_type'])
            raise ValueError('\n\nnot sure this is necessary but it is definitely triggered too easily for external data - are they genomes?!? ({}/{})'.format(90*gb,fastq_filesize_sum))
        elif fastq_filesize_sum < 1073741824: # < 1GB
            raise ValueError('\n\nnot sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('lt',sample.metadata['sample_type'])
    elif sample_type == 'rnaseq':
        if fastq_filesize_sum > 32212254720: # > 30GB
            raise ValueError('\n\nnot sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('gt',sample.metadata['sample_type'])
        elif fastq_filesize_sum < 1073741824: # < 1GB
            raise ValueError('\n\nnot sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('lt',sample.metadata['sample_type'])
    elif sample_type == 'custom_capture':
        if fastq_filesize_sum > 10737418240: # > 10GB
            raise ValueError('\n\nnot sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('gt',sample.metadata['sample_type'])
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

######## need to make it change the value so that we don't end up in a permanent loop
######## need to make it change the value so that we don't end up in a permanent loop
######## need to make it change the value so that we don't end up in a permanent loop

##################### is this doing anything but loop on the same samples forever?!?
##################### is this doing anything but loop on the same samples forever?!?
##################### is this doing anything but loop on the same samples forever?!?

def get_next_sample(pid,database,debug):

    q="SELECT d.sample_name,d.sample_type,d.capture_kit,d.pseudo_prepid FROM dragen_sample_metadata d "
    q+=" join prepT p on p.p_prepid=d.pseudo_prepid "

    if pid == 0:
        # q=q+"WHERE is_merged = 80000 ORDER BY PSEUDO_PREPID asc LIMIT 1 "
        # q=q+"WHERE is_merged = 80000 ORDER BY PSEUDO_PREPID desc LIMIT 1 "
        ###### need to pick up old als samples
        # q=q+"WHERE is_merged = 80000 ORDER BY priority asc, sample_type desc LIMIT 1 "
        ####### really old, so far not run samples are actually rather high in terms of ppid...?!?
        # q=q+"WHERE is_merged = 80000 and d.sample_type != 'Genome' ORDER BY p.prepid asc LIMIT 1 "
        who=socket.gethostname()
        if who == "dragen2.igm.cumc.columbia.edu":
            q=q+"WHERE is_merged = 80000 ORDER BY p.prepid desc LIMIT 1 "
        else:
            q=q+"WHERE is_merged = 80000 ORDER BY p.prepid desc LIMIT 1 "
            # q=q+"WHERE is_merged = 80000 and d.sample_type != 'Genome' ORDER BY p.prepid desc LIMIT 1 "
        # q=q+"WHERE is_merged = 80000 and d.sample_type != 'Genome' ORDER BY p.prepid desc LIMIT 1 "
        # q=q+"WHERE is_merged = 80000 and sample_type != 'Genome' ORDER BY pseudo_prepid desc LIMIT 1 "
        # q=q+"WHERE is_merged = 80000 and sample_type != 'Genome' ORDER BY priority asc, sample_type desc LIMIT 1 "
        # q=q+"WHERE is_merged = 80000 and sample_type = 'Exome' order by pseudo_prepid desc LIMIT 1 "
    else:
        q=q+("WHERE P_PREPID={}".format(pid))

    connection = get_connection(database)

    with connection.cursor() as cur:

        cur.execute(q)
        if cur.rowcount != 1:
            print("there's nothing to do")
            exit(0)
            # raise ValueError("couldn't get a sample")

        sample=cur.fetchone()

        cur.execute("update dragen_sample_metadata set is_merged = 80010 WHERE pseudo_prepid = {} and is_merged = 80000".format(sample['pseudo_prepid']) )
        if cur.rowcount != 1:
            raise ValueError("seems we couldn't get a lock on sample {}".format(pseudo_prepid));

        connection.commit()
        connection.close()

        return sample['sample_name'], sample['sample_type'], sample['pseudo_prepid'], sample['capture_kit']



def email_failure(stuff):
    emailAddresses = '{}@cumc.columbia.edu'.format(get_user())
    emailCmd = ('echo "Dragen Pipe failure" | mail -s "Dragen Pipe Failure" {} <<< {}'
               ).format(emailAddresses,stuff)
    print(emailCmd)
    os.system(emailCmd)

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
                        version='%(prog)s v3.0')
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
        email_failure(tb)
        emailit('alignment issue',tb)

