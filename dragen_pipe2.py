#!/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python -u
"""
Create a configuration file used by the Dragen system to processed IGM samples
based on their priority level in the Dragen pipeline queue
"""

import argparse
import time
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
import pymysql
import json
from dragen_sample import get_bed_file_loc

twat_global=0
twat_global_restage_list_no_mapping=False
twat_global_no_cram_conversion=False

# -input} = {file}
raw_config="""#================================================================================
# Dragen 2.5 Configuration File
#================================================================================
# SAMPLE SETUP
intermediate-results-dir = /staging/tmp
ref-dir = /staging/REF/b37_decoy/
{type}
# fastq-offset = 33 		
enable-auto-multifile = true
#================================================================================
# OUTPUT
#================================================================================
output-file-prefix = {sample_name}.{pseudo_prepid}.prerel # bypassing merge...
output-format = BAM
output-directory = {output_dir}
#================================================================================
# ALIGNMENT
#================================================================================

enable-map-align=true
enable-map-align-output=true      
enable-bam-indexing = true
enable-sort = true
enable-duplicate-marking = true
remove-duplicates = false
enable-sampling = true 			# automatically detect paired-end parameters with aligner test
enable-deterministic-sort = true 	# ensure sort order is completely repeatable at cost of a small decrease in speed

[Aligner]

match-score = 1 	# Score increment for matching reference nucleotide
mismatch-pen = 4 	# Score penalty for a mismatch
gap-open-pen = 6 	# Score penalty for opening a gap (insertion or deletion)
gap-ext-pen = 1 	# Score penalty for gap extension
unclip-score = 5 	# Score bonus for reaching each edge of the read
global = 0 		# If alignment is global (N-W) rather than local (S-W)
pe-orientation = 0 	# Expected paired-end orientation: 0=FR, 1=RF, 2=FF
pe-max-penalty = 60 	# Maximum pairing score penalty, for unpaired or distant ends
mapq-max = 60 		# Ceiling on reported MAPQ
supp-aligns = 3 	# Maximum supplimentary (chimeric) alignments to report per read
sec-aligns = 0 		# Maximum secondary (suboptimal) alignments to report per read
supp-as-sec = 0 	# If supplementary alignments should be reported with secondary flag
hard-clips = 6 		# Flags for hard clipping: 0 primary, 1 supplementary, 2 secondary
unpaired-pen = 80 	# Penalty for unpaired alignments in Phred scale
dedup-min-qual = 15 	# Minimum base quality for calculating read quality metric for deduplication
no-unpaired = 0 	# If only properly paired alignments should be reported for paired reads

#================================================================================
# MAPPER
#================================================================================

[Mapper]

seed-density = 0.5 	# Requested density of seeds from reads queried in the hash table
edit-mode = 0 		# 0 = No edits, 1 = Chain len test, 2 = Paired chain len test, 3 = Edit all std seeds
edit-seed-num = 6 	# For edit-mode 1 or 2: Requested number of seeds per read to allow editing on 
edit-read-len = 100 	# For edit-mode 1 or 2: Read length in which to try edit-seed-num edited seeds
edit-chain-limit = 29 	# For edit-mode 1 or 2: Maximum seed chain length in a read to qualify for seed editing
map-orientations = 0 	# 0=Normal, 1=No Rev Comp, 2=No Forward  (paired end requires Normal)
#================================================================================
# FIN
#================================================================================
"""

def emailit(rarp,arse):

    if twat_global!=0:

        print("not sending email")
        import os
        os._exit(1)

    smtpserver = 'localhost'
    fromAddr = 'igm-bioinfo@columbia.edu'
    emailMsg = MM('alternative') #email.MIMEMultipart.MIMEMultipart('alternative')
    emailMsg['Subject'] = rarp
    emailMsg['From'] = fromAddr
    # to=['dh2880@cumc.columbia.edu']
    # to=['dh2880@cumc.columbia.edu','ce2373@cumc.columbia.edu','nb2975@cumc.columbia.edu'] # 'mml2204@cumc.columbia.edu']
    to=['igm-bioinfo@columbia.edu'] # '880@cumc.columbia.edu','mml2204@cumc.columbia.edu']
    # to=['dsth@cantab.net','dh2880@cumc.columbia.edu']
    # to=['igm-bioinfo@columbia.edu','nb2975@cumc.columbia.edu'] # 'mml2204@cumc.columbia.edu']
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

def check_space(vol):
    df = subprocess.Popen(["df",vol], stdout=subprocess.PIPE)
    # print(df.stdout.readline()) # print(df.stdout.readline())
    # split eithout arg for arbitrary whitespace  
    # only want first of stdout, stderr # output = df.communicate()[0] 
    # output = df.communicate()[0].decode("utf-8").split("\n")[2].split()
    # all in 1K-blocks by default
    # size, used, available, percent, mount = df.communicate()[0].decode("utf-8").split("\n")[2].split()
    filesystem, size, used, available, percent, mount = df.communicate()[0].decode("utf-8").split("\n")[1].split()
    available_t=float(int(available)/(1024*1024*1024))
    # print('type={}, output="{}"'.format(type(output),output))
    # [1].split()
    print("size={}, used={}, available={}, percent={}, mountpoint={}".format(size,used,available,percent,mount));
    return available_t
    # return rg[0]["count"]==0 and dsm[0]["count"]==0

def no_work(database):
    # from pprint import pprint as pp
    rg = run_query("select count(1) count from Lane where rg_status = 'fastq_ready'",database)
    # pp(rg)
    dsm = run_query("select count(1) count from dragen_sample_metadata where is_merged = 80000",database)
    # pp(dsm)
    return rg[0]["count"]==0 and dsm[0]["count"]==0

def main(reset_dragen,no_prerelease_align,experiment_id,no_gvcf):

    if "dragen2.igm.cumc.columbia.edu" == socket.gethostname():
        subprocess.run(["echo -e '\e[38;5;196m'"],shell=True)

    ###### some of these should 'perhaps' be arguments?!?
    seqscratch_drive = 'seqscratch_ssd'
    database="sequenceDB"
    dontexecute = False
    run_type_flag = True
    debug = True

    print('reset_dragen={}\nno_prerelease_align={}\nexperiment_id={}\ntwat_global={}'.format(reset_dragen,no_prerelease_align,experiment_id,twat_global))

    config = get_config()

    # saga sample errors causing system to hang on reset with restarts...
    if reset_dragen:
        print("reseting dragen")
        os.system("dragen_reset")
    else:
        print("not reseting system")

    ######### need to force this when space goes < 5T to make sure ssd drains off - start sending emails and force this
    # if check_space("/nfs/seqscratch_ssd")<5.0 or (no_work(database) and no_gvcf==False):
    ###### consider having a max number in queue too?!?
    if check_space("/nfs/seqscratch_ssd")<12.0:
    # if check_space("/nfs/seqscratch_ssd")<12.0:
        print("space is low - need to clear any backlog and continue on to alignment...")
        # time.sleep(30)
        os.system("/nfs/seqscratch09/informatics/tmp/RunAux.pl")
        # os.system("/nfs/seqscratch_ssd/informatics/tmp/RunAux.pl")
    elif no_work(database) and no_gvcf==False:
        print("nothing to align - run variant calling")
        os.system("/nfs/seqscratch09/informatics/tmp/RunAux.pl")
        # os.system("/nfs/seqscratch_ssd/informatics/tmp/RunAux.pl")
        exit(1)
    else:
        print("not running gvcf generation")

    if check_space("/nfs/seqscratch_ssd")<2.5:
    # if check_space("/nfs/seqscratch_ssd")<3.0:
        print("there's insufficient space to continue")
        time.sleep(30)
        exit(1)
    else:
        print("will run")

    try:

        pseudo_prepid = 0
        sample_name, sample_type, pseudo_prepid, capture_kit, tx, mi = get_next_sample(pseudo_prepid,database,debug,no_prerelease_align,experiment_id)

        while sample_name is not None:

            print(" > running experiment_id = {}".format(pseudo_prepid))

            if tx>1:

                j=json.loads(mi)

                if len(j)!=1:
                    raise ValueError("not doing this yet!?!")

                filey=''

                if j[0]['format']=='bam':
                    filey='bam-input = {}/{}.{}'.format(j[0]['path'],j[0]['name'],j[0]['format'])
                elif j[0]['format']=='cram':

                    scratch='/nfs/{}/ALIGNMENT/BUILD37/DRAGEN/{}/{}.{}/'.format(seqscratch_drive,sample_type.upper(),sample_name,pseudo_prepid)

                    if not os.path.isdir(scratch):
                        os.mkdir(scratch)

                    cmd='time dragen -f --enable-map-align=false --file-conversion=true --cram-input {}/{}.{} \
                      --output-format=bam -r {} --output-directory {} --output-file-prefix {}.{}.cram2bam'.format(
                      j[0]['path'],     j[0]['name'], j[0]['format'],       j[0]['ref'],        scratch,      sample_name,pseudo_prepid
                    )

                    print('will run\n\n\t{}\n'.format(cmd))
                    if twat_global_no_cram_conversion==False and os.system(cmd):
                        raise ValueError("error in dragen cram2bam conversion.")

                    filey='bam-input = {}/{}.{}.cram2bam.bam'.format(scratch,sample_name,pseudo_prepid)

                else:
                    raise ValueError("not doing this yet!?!")

                prepid = run_query("select prepid from prepT where p_prepid = {} and failedprep = 0 or failedprep = 11 or failedprep >= 100 ".format(pseudo_prepid),database)
                print('prepid= {}'.format(prepid))
                run_sample_external(config,database,seqscratch_drive,sample_type,capture_kit,sample_name,pseudo_prepid,filey,prepid[0]['prepid'])

            else:

                try:
                    sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,database)
                except Exception as e:
                    raise Exception("unable to construct dragen sample for {} : {}".format(sample_name,e))

                try:
                    run_sample(sample,dontexecute,config,seqscratch_drive,database,debug)
                except Exception as e:
                    status='Dragen Alignment Failure'
                    update_status(sample,status,database)
                    raise Exception("unable to run dragen for {} : {}".format(sample_name,e))
                    # email_failure()
                    # sys.exit(1)

            try:
                if experiment_id!=0:
                    print("gonna stop")
                    os._exit(0)
                sample_name, sample_type, pseudo_prepid, capture_kit, tx, mi = get_next_sample(0,database,debug,no_prerelease_align,experiment_id)
            except:
                print("No more samples in the queue")
                sys.exit(0)

    except Exception as e:

        msg="Fastq error ({})".format(str(e).replace('\n',' ').replace('  ',' ').lstrip(' '))
        tb = traceback.format_exc()
        print("sending email with '{}'\n\nand '{}'".format(msg,tb))
        emailit(msg,tb)
        print(sample_name+", "+sample_type+", "+str(pseudo_prepid)+", "+capture_kit);

        connection = get_connection(database)
        with connection.cursor() as cursor:
            # if wipe:
                # cursor.execute("update prepT set is_released = 0, failedprep = 11, status = %s where p_prepid = %s and sample_internal_name = %s",(msg,pseudo_prepid,sample_name))
                # cursor.execute("delete from dragen_sample_metadata where pseudo_prepid = %s and sample_name = %s",(pseudo_prepid,sample_name))
            # else:
            cursor.execute("update prepT set status = %s where p_prepid = %s and sample_internal_name = %s",(msg,pseudo_prepid,sample_name))
            cursor.execute("update dragen_sample_metadata set is_merged = 80100 where pseudo_prepid = %s and sample_name = %s",(pseudo_prepid,sample_name))

            connection.commit()

        raise Exception("\n\n> "+msg)

        exit(1)

def update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid,database):
    read_group_insert = ("UPDATE Lane l "
                         "JOIN Flowcell f on l.FCID=f.FCID "
                         "SET released=1,step_status='completed' "
                         "WHERE FCILLUMID='{}' AND "
                         "LANENUM={} AND "
                         "failr1 IS NULL AND "
                         "PREPID={} "
                        ).format(rg_fcillumid,rg_lane_num,rg_prepid)
    run_query(read_group_insert,database)

def check_bam_found_vs_bam_db(sample,qualified_bams_found,database):

    max_prepid = max(map(lambda x:int(x),sample.metadata['prepid']))

    if is_external_or_legacy_sample(max_prepid,database) == True:
       print("(A) Checking bams found vs RGs from fastqs")
       for laneFCID in sample.metadata['lane']: #[0]: #loop over read groups
            rg_lane_num,rg_fcillumid,rg_prepid = laneFCID

            bam_loc = ("{output_dir}/{sample_name}.{pseudo_prepid}.{rg_fcillumid}.{rg_lane_num}.*bam"
                      ).format(rg_fcillumid=rg_fcillumid,
                               rg_lane_num=rg_lane_num,
                               **sample.metadata)
            if glob(bam_loc) == []:
                raise Exception("Bam {} not found!".format(bam_loc))

    else:

        query = "SELECT fcillumid,laneNum FROM prepT p JOIN Lane l on p.PREPID=l.PREPID JOIN Flowcell f on l.fcid=f.fcid WHERE f.fail=0 and p.experiment_id = {pseudo_prepid} AND FCILLUMID NOT LIKE 'X%' AND FAILR1 IS NULL AND FAILR2 IS NULL AND (FAILEDPREP = 0 or FAILEDPREP>=100)".format(**sample.metadata)

        qualified_bams_in_db = run_query(query,database)

        if len(qualified_bams_in_db) != len(qualified_bams_found):

            from pprint import pprint as pp
            pp(qualified_bams_found)
            pp(qualified_bams_in_db)
            raise Exception("Number of bams found ({}) != number of RGs in database ({})!".format(len(qualified_bams_found),len(qualified_bams_in_db)))

        for db_bam in qualified_bams_in_db:
            db_bam = ("{output_dir}/{sample_name}.{pseudo_prepid}.{fcillumid}.{lane_num}.bam"
                     ).format(fcillumid=db_bam['fcillumid'],lane_num=db_bam['laneNum'],**sample.metadata)
            if glob(db_bam) == []:
                raise Exception("Bam {} not found!".format(db_bam))

def get_component_bams(sample,debug,database):

    qualified_bams_found = glob('{output_dir}/*bam'.format(**sample.metadata))

    check_bam_found_vs_bam_db(sample,qualified_bams_found,database)
    if len(qualified_bams_found) < 1:
        raise Exception("No qualified bams were found!")
    tmp_bam_list = []
    for bam in qualified_bams_found:
        bam_name = bam.split('/')[-1]
        tmp_bam_list.append(bam_name)
    component_bams = ','.join(sorted(tmp_bam_list))
    return component_bams

def run_sample(sample,dontexecute,config,seqscratch_drive,database,debug):

    from pprint import pprint as pp
    pp(vars(sample))
    if 'bams' in sample.metadata:
        if twat_global_restage_list_no_mapping:
            return
            # os._exit(1)
    else:
        check_Fastq_Total_Size(sample,debug)

    setup_dir(seqscratch_drive,sample,debug)

    print("using scratch dir {}".format(sample.metadata['output_dir']))
    if not os.path.isdir(sample.metadata['output_dir']):
        print("this is wrong!")
        os._exit(1)

    # rfo
    r_list='{}/restage.txt'.format(sample.metadata['output_dir'])
    if twat_global_restage_list_no_mapping:
        print("only generating restaging list {}".format(r_list))
        if(os.path.isfile(r_list)):
            rfo=open(r_list,"w+")
            rfo.truncate()
            rfo.close()
        # else:
            # rfo=open(list,"w+")
        # rfo.close()

    existing_bams_check = False
    output_dir = sample.metadata['output_dir']
    pseudo_prepid = sample.metadata['pseudo_prepid']
    submit_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    existing_bams = glob("{}/*.bam".format(output_dir))

    arse1='{}/{}.{}.bam'.format(sample.metadata['output_dir'],sample.metadata['sample_name'],sample.metadata['pseudo_prepid'])
    arse2='{}/{}.{}.merge.bam'.format(sample.metadata['output_dir'],sample.metadata['sample_name'],sample.metadata['pseudo_prepid'])

    if os.path.exists(arse1):
        os.remove(arse1)
    if os.path.exists(arse2):
        os.remove(arse2)

    arse3='{}/{}.{}.realn.recal.bam'.format(sample.metadata['output_dir'],sample.metadata['sample_name'],sample.metadata['pseudo_prepid'])
    if os.path.exists(arse3):
        os.remove(arse3)
    arse4='{}/{}.{}.realn.bam'.format(sample.metadata['output_dir'],sample.metadata['sample_name'],sample.metadata['pseudo_prepid'])
    if os.path.exists(arse4):
        os.remove(arse4)

    qualified_bams_found = glob('{output_dir}/*bam'.format(**sample.metadata))

    if existing_bams == [] or existing_bams_check == False:

        update_queue(pseudo_prepid,database)

        if 'bams' in sample.metadata:

            print("\n---------------------------------------------------------------");
            print("GO THROUGH MAPPED RG COMPONENTS DIRECT FROM FC - i.e. we ONLY link files over!?!")
            print("---------------------------------------------------------------\n");

            bc=0
            for b in sample.metadata['bams']:

                bc=bc+1

                bam='{}/{}'.format(output_dir,os.path.basename(b[2]))
                if not os.path.exists(bam):
                    os.symlink(b[2],bam)

                index='{}/{}.bai'.format(output_dir,os.path.basename(b[2]))
                if not os.path.exists(index):
                    os.symlink(b[2]+'.bai',index)

        else:

            ### this could never have worked?!?
            # for laneFCID in sample.metadata['lane'][0]: #loop over read groups
            for laneFCID in sample.metadata['lane']: 

                rg_lane_num,rg_fcillumid,rg_prepid = laneFCID

                lane_table = run_query("select count(1) count, data from Lane l join Flowcell f on l.fcid=f.fcid where prepid={} and fcillumid='{}' and lanenum='{}'".format(rg_prepid,rg_fcillumid,rg_lane_num),database)[0]
                if lane_table["count"]==1:
                    if lane_table["data"] == None:
                        print("legacy internal - just map it...")
                    else:
                        print("This seems to have a data field and yet we're running in legacy mode - make this an error after the hybrids are finished!?!")
                elif lane_table["count"]==0:
                    print("legacy external procedure")
                else:
                    raise ValueError("what the heck is going on")

                setup_first_read_RG(sample,rg_lane_num,rg_fcillumid,rg_prepid,debug)
                set_seqtime(rg_fcillumid,sample,database)
                create_align_config(sample,rg_lane_num,rg_fcillumid,rg_prepid)

                run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug,r_list)

                #### utterly pointless
                update_lane_metrics(sample,rg_lane_num,rg_fcillumid,rg_prepid,database)

        ####### f' me - just globbing...
        if twat_global_restage_list_no_mapping==False:
            component_bams = get_component_bams(sample,debug,database)
            update_dragen_metadata_prepT_status(sample,component_bams,database,pseudo_prepid,debug)
        else:
            print("simply build restage list {}".format(r_list));
            run_query("update prepT p join dragen_sample_metadata d on p.experiment_id=d.experiment_id set status = 'Marked_For_Restaging', is_merged = 81000 where failedprep = 0  and d.experiment_id = {}".format(pseudo_prepid),database)
            # os._exit(0)

        # rm_query = "DELETE FROM {0} WHERE pseudo_prepid={1}".format("tmp_dragen",pseudo_prepid)
        # if debug:
        #    print(rm_query)
        # run_query(rm_query,database)

    else:
        raise Exception("[legacy] Sample with bam files already exists!")

def run_dragen_on_read_group(sample,rg_fcillumid,rg_lane_num,debug,r_list):

    output_dir = sample.metadata['output_dir']

    dragen_cmd = ['dragen', '-f', '-v', '-c', sample.metadata['conf_file'], '--qc-cross-cont-vcf', '/opt/edico/config/sample_cross_contamination_resource_GRCh37.vcf.gz', '--watchdog-active-timeout', '600']

    if debug:
        print(' '.join(dragen_cmd))
    stderr_file_loc = ('{}/{}.{}.{}.{}.dragen.err'
                      ).format(sample.metadata['log_dir'],sample.metadata['sample_name'],
                               sample.metadata['pseudo_prepid'],rg_fcillumid,rg_lane_num)

    stdout_file_loc = ('{}/{}.{}.{}.{}.dragen.out'
                      ).format(sample.metadata['log_dir'],sample.metadata['sample_name'],
                               sample.metadata['pseudo_prepid'],rg_fcillumid,rg_lane_num)

    if twat_global_restage_list_no_mapping==False:

        bored = ( sample.metadata['output_dir'], sample.metadata['sample_name'], sample.metadata['pseudo_prepid'], rg_fcillumid, rg_lane_num )
        cp1='{}/{}.{}.{}.{}.bam'.format(*bored)
        cp2='{}/{}.{}.{}.{}.bam.bai'.format(*bored)
        cp3='{}/{}.{}.{}.{}.bam.md5sum'.format(*bored)
        cp4='{}/{}.{}.{}.{}.time_metrics.csv'.format(*bored)
        cp5='{}/{}.{}.{}.{}-replay.json'.format(*bored)

        if os.path.isfile(cp1) and os.path.isfile(cp2) and os.path.isfile(cp3) and os.path.isfile(cp4) and os.path.isfile(cp5) \
          and os.path.getmtime(cp5)>= os.path.getmtime(cp4) \
          and os.path.getmtime(cp4)>= os.path.getmtime(cp3) \
          and os.path.getmtime(cp3)>= os.path.getmtime(cp2) \
          and os.path.getmtime(cp2)>= os.path.getmtime(cp1):

            print(">>>>> appears to be done <<<<<\n");
            return
            # pass
        else:
            print(">>>>> need to run it\n");

        sys.stdout.flush()
        sys.stderr.flush()

        # os._exit(1)

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
            # email_failure(dragen_cmd)
            # emailit('alignment issue ' + stderr_file_loc,dragen_cmd)
            raise Exception("Dragen alignment did not complete successfully (check log)")
    else:
        rfo=open(r_list,"a+")
        rfo.write("{}\n{}\n".format(sample.metadata['first_fastq1'],sample.metadata['first_fastq2']))
        rfo.close()
        # os._exit(0)

    try:
        subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
    except:
        pass

def update_dragen_metadata_prepT_status(sample,component_bams,database,pseudo_prepid,debug):

    seqscratch_drive = sample.metadata['output_dir'].split('/')[2]
    get_dragen_metadata_query = ("SELECT * from dragen_sample_metadata WHERE pseudo_prepid = {} ").format(pseudo_prepid)
    dbInfo = run_query(get_dragen_metadata_query,database)

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
                        (STATUS,STATUS_TIME,PREPID,USERID,POOLID,SEQID)
                        VALUES ('{status}',unix_timestamp(),{prepid},{userID},0,0)
                     """.format(userID=userID,status=status,
                                prepid=sample.metadata['prepid'][0],
                                sample_name=sample.metadata['sample_name'])
    run_query(statusT_insert,database)

def update_queue(pseudo_prepid,database):

    who=socket.gethostname()
    state=0;
    # if who == "dragen1.igm.cumc.columbia.edu":
    state=80011
    # elif who == "dragen2.igm.cumc.columbia.edu":
        # state=80012
    # elif twat_global_restage_list_no_mapping:
        # print("we don't care")
    # else:
        # raise ValueError("{} is not allowed to run this".format(who))

    connection = get_connection(database)
    print(" db= {} and pp= {}".format(database,pseudo_prepid))

    with connection.cursor() as cur:

        ########### clear previous errors - shouldn't be there!?!
        # q="update dragen_sample_metadata set is_merged = {} WHERE is_merged = {}".format( state+10, state )
        ###### usual gtac ball dropping situation - need so isolate propblem samples asap.
        q="update dragen_sample_metadata set is_merged = {} WHERE is_merged = {}".format( state+50, state )
        print(q)
        cur.execute(q)

        q="update dragen_sample_metadata set is_merged = {} WHERE pseudo_prepid = {} and is_merged = 80010".format(state,pseudo_prepid) 
        print(q)
        cur.execute(q)
        if cur.rowcount != 1 and twat_global==0:
            raise ValueError("seems we couldn't get a lock on sample {} : {} : {}".format(pseudo_prepid,cur.rowcount,q));
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
        ######## this isn't correct but it doesn't matter - i.e. it's giving the date only and not the time/timezone but merging puts that in on it's own
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
    first_read1,first_read2 = get_first_read(sample,rg_lane_num,rg_fcillumid)
    sample.set('first_fastq1',first_read1)
    sample.set('first_fastq2',first_read2)

def get_first_read(sample,rg_lane_num,rg_fcillumid):
    #Using first fastq as template for all fastq.gz
    for fastqLoc in sample.metadata['fastq_loc']:
        if rg_fcillumid in fastqLoc.split('/')[-1]:
            ######## AZIPF have FASTQ filenames of \w+L\d{3} names so just use standard '_' delimiter to avoid - this is terrible stuff!?!
            # RGfastqStr ='{0}/*L00{1}_R1_*.fastq.gz'.format(fastqLoc,rg_lane_num)
            RGfastqStr ='{0}/*_L00{1}_R1_*.fastq.gz'.format(fastqLoc,rg_lane_num)
            print(RGfastqStr)
            RGfastqs = glob(RGfastqStr)
            if len(RGfastqs)==0:
                raise Exception("read1 of pair missing '{}'".format(RGfastqStr))
            else:
                print("read1 of pair found '{}'".format(RGfastqStr))
            print("have RGfastqs = '{}'".format(RGfastqs))

    first_read = sorted(RGfastqs)[0]
    second_read = first_read.replace('_R1_','_R2_')
    if not os.path.isfile(second_read): 
        print("read2 of pair missing '{}'".format(RGfastqStr))
        raise Exception("read2 of pair missing '{}'".format(RGfastqStr))
    else:
        print("read2 of pair found '{}'".format(RGfastqStr))
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
            if not os.path.isfile(fastq):
                raise Exception("file seems to have gone AWOL : '{}'".format(fastq))
            fastq_filesize = os.stat(os.path.realpath(fastq)).st_size
            fastq_filesize_sum += fastq_filesize
            print("  > checking fastq '{}' ({}GB)".format(fastq,'%.2f'%(fastq_filesize/gb)))
    print("Sum of Fastq size: {}".format(fastq_filesize_sum))

    print("THIS IS REALLY A BITR POINTLESS NOW - let them go through and get rejected as absurd smaples!?!")
    ##################### cannot recall where the update is done?!?
    if sample_type == 'genome' or sample_type == 'genome_as_fake_exome':

        min = 2594 # min = 25949672960
        if fastq_filesize_sum < min:
            # if fastq_filesize_sum < 42949672960: # < 40GB
            gb=min/(1024*1024*1024)
            # run_query("UPDATE dragen_sample_metadata SET is_merged = 80113 where pseudo_prepid = {}".format(pseudo_prepid),database)
            raise ValueError('Release Error; Genome FASTQ size too low ({0:.2f} GB).'.format(gb))
            # user_input_fastq_size('lt',sample.metadata['sample_type'])

    elif sample_type == 'exome':

        if fastq_filesize_sum > (150*gb): # > 30GB
            # user_input_fastq_size('gt',sample.metadata['sample_type'])
            raise ValueError('not sure this is necessary but it is definitely triggered too easily for external data - are they genomes?!? ({}/{})'.format(90*gb,fastq_filesize_sum))
        # elif fastq_filesize_sum < 573741824: # < 1GB
        elif fastq_filesize_sum < 773741824: 
        # elif fastq_filesize_sum < 1073741824: # < 1GB
            raise ValueError('exome fastq file sum is rather small : {}'.format(fastq_filesize_sum) )
            user_input_fastq_size('lt',sample.metadata['sample_type'])
    elif sample_type == 'rnaseq':
        if fastq_filesize_sum > 32212254720: # > 30GB
            raise ValueError('not sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('gt',sample.metadata['sample_type'])
        elif fastq_filesize_sum < 1073741824: # < 1GB
            raise ValueError('not sure this is necessary - changing to raise to get email and let it move on...')
            user_input_fastq_size('lt',sample.metadata['sample_type'])
    elif sample_type == 'custom_capture':
        if fastq_filesize_sum > 10737418240: # > 10GB
            raise ValueError('not sure this is necessary - changing to raise to get email and let it move on...')
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
        print('{0}/*_L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number))
    read = glob('{0}/*_L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number))
    #The fastqs might still be on the quantum tapes
    if read == []:
        raise Exception("Fastq file not found!")
    else:
        return sorted(read)[0]

def get_next_sample(pid,database,debug,no_prerelease_align,experiment_id):

    if no_prerelease_align==False:

        if False:
            if os.system("/nfs/seqscratch09/informatics/bin/GoButtonAlign align")!=0:
                raise Exception("\n\nproblem with se alignment process!\n\n")
        else:
            os.system("/nfs/seqscratch09/informatics/bin/GoButtonAlign align")
            # os.system("/nfs/seqscratch_ssd/informatics/bin/GoButton align")


    else:
        print("no prerelease align")

    q="SELECT d.sample_name,d.sample_type,d.experiment_id,d.capture_kit,d.pseudo_prepid,d.is_external ticket_num,d.mapping_input FROM dragen_sample_metadata d "
    q+=" join prepT p on p.experiment_id=d.experiment_id where (failedprep=0 or failedprep=11 or failedprep>=100) and "

    if experiment_id == 0:

        who=socket.gethostname()
        if who == "dragen1.igm.cumc.columbia.edu":
        # if who != "dragen2.igm.cumc.columbia.edu":
        ###### just get the diagseq done before the AZ...?!?
            # q=q+"is_merged = 80000 ORDER BY d.experiment_id desc LIMIT 1 "
            q=q+"is_merged = 80000 ORDER BY d.sample_name desc LIMIT 1 "
            # q=q+"is_merged = 80000 and externaldata is null ORDER BY d.experiment_id desc LIMIT 1 "
        else:
            # q=q+"is_merged = 80000 ORDER BY d.experiment_id desc LIMIT 1 "
            q=q+"is_merged = 80000 ORDER BY d.experiment_id asc LIMIT 1 "

    else:
        q=q+("pseudo_prepid={} group by experiment_id".format(experiment_id))

    connection = get_connection(database)

    with connection.cursor() as cur:

        cur.execute(q)
        if cur.rowcount != 1:
            print("[dp] didn't retrive a sample to map from dsm ({})".format(cur.rowcount))
            sys.exit(0)
            # raise ValueError("couldn't get a sample")

        sample=cur.fetchone()

        if sample['pseudo_prepid']!=sample['experiment_id']:
            raise ValueError("experiment_id does not match legacy pseudo_prepid")

        cur.execute("update dragen_sample_metadata set is_merged = 80010 WHERE pseudo_prepid = {} and is_merged = 80000".format(sample['pseudo_prepid']) )

        if cur.rowcount != 1 and experiment_id==0:
            raise ValueError("seems we couldn't get a lock on sample {}".format(sample['pseudo_prepid']));

        connection.commit()
        connection.close()

        from pprint import pprint as pp
        pp(q)
        print("we have {}".format(no_prerelease_align))
        pp(sample)
        # exit(1)
        # wtf?!?
        return sample['sample_name'], sample['sample_type'], sample['pseudo_prepid'], sample['capture_kit'], int(sample['ticket_num']),sample['mapping_input']

def run_sample_external(config,database,seqscratch_drive,sample_type,capture_kit,sample_name,pseudo_prepid,filey,prepid):

    output_dir      = '/nfs/{}/ALIGNMENT/BUILD37/DRAGEN/{}/{}.{}/'.format(seqscratch_drive,sample_type.upper(),sample_name,pseudo_prepid)
    script_dir      = output_dir+'/scripts'
    log_dir         = output_dir+'/logs'
    dragen_stdout   = "{}/{}.{}.dragen.out".format(log_dir,sample_name,pseudo_prepid)
    dragen_stderr   = "{}/{}.{}.dragen.err".format(log_dir,sample_name,pseudo_prepid)

    subprocess.call(['mkdir','-p',script_dir])
    subprocess.call(['mkdir','-p',log_dir])

    x = glob("{}/*.bam".format(output_dir))
    if len(x)!=0:

        if True:
            x = glob("{}/*.bam".format(output_dir))
            for i in x:
                print( "TEMPORARY NASTY HACK : {}".format(i) );
                if i.find('cram2bam')==-1:
                    print( "wiping bam : {}".format(i) );
                    os.system( "rm -f {}".format(i) )
                else:
                    print( "ignoring cram2bam intermediate : {}".format(i) );
                # os.system("rm -f {}/*.bam".format(output_dir))
            # exit(1)
            # os.system("rm -f {}/*.bam".format(output_dir))
        else:
            run_query("UPDATE dragen_sample_metadata SET is_merged = 80113 where pseudo_prepid = {}".format(pseudo_prepid),database)
            exit(1)
            raise Exception("EXTERNAL_BAM2BAM : Sample with bam files already exists!")

    update_queue(pseudo_prepid,database)

    conf_file = ("{}/{}.{}.DragenAlignment.conf").format(script_dir,sample_name,pseudo_prepid)

    final_config_cont = raw_config.format(
       # file=bam_file,
       sample_name=sample_name,
       pseudo_prepid=pseudo_prepid,
       output_dir=output_dir,
       type=filey
    )

    with open(conf_file,'w') as cf:
        cf.write(final_config_cont)

    dragen_cmd = ['dragen', '-f', '-v', '-c', conf_file, '--qc-cross-cont-vcf', '/opt/edico/config/sample_cross_contamination_resource_GRCh37.vcf.gz', '--watchdog-active-timeout', '600']

    stderr_file_loc = ('{}/{}.{}.dragen.err').format(log_dir,sample_name,pseudo_prepid)
    stdout_file_loc = ('{}/{}.{}.dragen.out').format(log_dir,sample_name,pseudo_prepid)

    # makes no sense here?!?
    if twat_global_restage_list_no_mapping:
        print("this makes no sense at all!?!")
        os._exit(1)
    else:
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

        # with?!?
        dragen_stderr.close()

        if rc != 0:
            run_query("UPDATE dragen_sample_metadata SET is_merged = 80112 where pseudo_prepid = {}".format(pseudo_prepid),database)
            exit(1)
            raise Exception("EXTERNAL_BAM2BAM : Dragen alignment did not complete successfully : {} ".format(dragen_cmd))

    try:
        subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
    except:
        pass

######## f'ing nasty!?!
    qualified_bams_found = glob('{}/*prerel*bam'.format(output_dir))
    if len(qualified_bams_found) != 1:
        raise Exception("EXTERNAL: this is wrong : " + qualified_bams_found)

    print("Performing dragen_sample_metadata update")
    run_query("UPDATE dragen_sample_metadata set seqscratch_drive='{}',is_merged=0,component_bams='{}' WHERE pseudo_prepid={}".format(
      seqscratch_drive,qualified_bams_found[0].split('/')[-1],pseudo_prepid
    ),database)

    status='Dragen Alignment Completed'

    run_query("UPDATE prepT SET status='{}',status_time=unix_timestamp() WHERE p_prepid={}".format(status,pseudo_prepid),database)
    run_query("UPDATE prepT_public SET prepT_status='{}', prepT_status_time=unix_timestamp() WHERE p_prepid={}".format(status,pseudo_prepid),database)

    userID = get_user_id(database)
    # ['prepid'][0],
    run_query("INSERT INTO statusT (STATUS,STATUS_TIME,PREPID,USERID,POOLID,SEQID) VALUES ('{status}',unix_timestamp(),{prepid},{userID},0,0)".format(
      userID=userID,
      status=status,
      prepid=prepid,
      sample_name=sample_name
    ),database)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)

    ###### get rid of some of these options!?!
    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument("-e", "--experiment_id", type=int)
    # group.add_argument("-a", "--auto",          default=False,              action="store_true",    help="Run pipeline in automated mode")
    # group.add_argument("-p", "--pseudo_prepid", type=str,                                           help="Run dragen in single sample mode with provided pseudo_prepid")
    # parser.add_argument("-d", "--debug",        default=False,              action="store_true",    help="Verbose output for debugging")
    # parser.add_argument("--seqscratch_drive",   default='seqscratch_ssd',   action='store',         help="Set output destination")

    parser.add_argument("-r", "--reset", default=False, action="store_true") # make it an option flag with store_true?!?
    parser.add_argument("--no_prerelease_align", default=False, action="store_true") # make it an option flag with store_true?!?
    parser.add_argument("--no_gvcf", default=False, action="store_true") # make it an option flag with store_true?!?
    parser.add_argument("-e", "--experiment_id", default=0, type=int)

    parser.add_argument("--restage_list", default=False, action="store_true") # make it an option flag with store_true?!?

    parser.add_argument("--no_cram_conversion", default=False, action="store_true") # make it an option flag with store_true?!?

    arg_list=parser.parse_args()

    if arg_list.experiment_id!=0:
        twat_global=arg_list.experiment_id
        arg_list.no_prerelease_align=False

    if arg_list.restage_list==True:
        twat_global_restage_list_no_mapping=arg_list.restage_list
    
    if arg_list.no_cram_conversion==True:
        twat_global_no_cram_conversion=arg_list.no_cram_conversion

    main(arg_list.reset,arg_list.no_prerelease_align,arg_list.experiment_id,arg_list.no_gvcf)
    # main(run_type_flag, args.debug, args.dontexecute, database, args.seqscratch_drive)


