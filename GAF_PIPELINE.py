#from CHGV_mysql import getSequenceDB

import argparse
import getopt
import glob
import logging
import os
import sys
import subprocess
import time
from collections import Counter
from operator import itemgetter
from utilities import *

def submit(config,args,run_info_dict,database):
    bcl2fastq_scratch_drive = config.get('locs','bcl2fastq_scratch_drive')
    archive_dir = config.get('locs','fastq_archive_drive')
    logger = logging.getLogger(__name__)
    address = '{}@cumc.columbia.edu'.format(get_user_id(database))
    python36_program = config.get('programs','python36_program')
    scriptLoc = '/nfs/goldstein/software/sequencing_pipe/dev/'
    #scriptLoc = '/nfs/goldstein/software/sequencing_pipe/production/'
    runFolder = run_info_dict['runFolder']
    machine = run_info_dict['machine']
    fcillumid=args.fcillumid

    BCLCmd =         ("{} {}/bcl.py -f {}").format(python36_program,scriptLoc,fcillumid)
    BCLMySQLCmd =    ("{} {}/bcl_mysql.py -f {}").format(python36_program,scriptLoc,fcillumid)
    postBCLCmd =     ("{} {}/post_bcl.py -f {}").format(python36_program,scriptLoc,fcillumid)
    storageCmd =     ("{} {}/storage.py -f {}").format(python36_program,scriptLoc,fcillumid)
    fastqcCmd =      ("{} {}/fastqc.py -f {}").format(python36_program,scriptLoc,fcillumid)
    fastqcMySQLCmd = ("{} {}/fastqc_mysql.py -f {}").format(python36_program,scriptLoc,fcillumid)

    if args.test == True:
        BCLCmd += ' --test'
        BCLMySQLCmd += ' --test'
        postBCLCmd += ' --test'
        storageCmd += ' --test'
        fastqcCmd += ' --test'
        fastqcMySQLCmd += ' --test'

    if args.run == False:
        #prints out all the commands to run manually
        print
        print("="*35+'Scripts Commands'+"="*35)
        print(BCLCmd)
        print(BCLMySQLCmd)
        print(postBCLCmd)
        print(storageCmd)
        print(fastqcCmd)
        #fastqcMySQLCmd is now run within the fastqcCmd
        print(fastqcMySQLCmd)
        print("="*86)
        print
        print(('Log file: {}' 
            ).format(logging.getLoggerClass().root.handlers[0].baseFilename))

    else:
        #run bcl2fastq
        os.system(BCLCmd) #bcl2fastq auto submits to cluster
        logger.info(BCLCmd)
        os.system(BCLMySQLCmd) #quick msyql updates.  
        logger.info(BCLMySQLCmd)

        #run post_bcl2fastq
        create_post_bcl_script(config,archive_dir,runFolder,machine,fcillumid,address,postBCLCmd)
        create_run_qsub_command(sample,'post_bcl','bcl',run_info_dict)
        create_storage_script(config,archive_dir,runFolder,machine,fcillumid,address,storageCmd)
        create_run_qsub_command(sample,'storage','post_bcl',run_info_dict)
        #create_fastqc_script(config,archive_dir,runFolder,machine,fcillumid,address,fastqcCmd)
        #create_run_qsub_command(sample,'fastqc','storage',run_info_dict)

def create_run_qsub_command(sample,step,fcillumid,run_info_dict):
    logger = logging.getLogger(__name__)
    qsub_cmd = ('{qsub_program} -N {machine}_{fcillumid}_{bcl2fastq_scratch_drive}_{step} '
                '/nfs/seqscratch1/Runs/{runFolder}/{machine}_{fcillumid}_{archive_dir}_{step}.sh'
               ).format(qsub_program=config.get('programs','qsub_program'),
                        machine=run_info_dict['machine'],
                        fcillumid=fcillumid,
                        archive_dir=run_info_dict['marchive_dir'],
                        runFolder=run_info_dict['runFolder'],
                        bcl2fastq_scratch_drive=config.get('locs','bcl2fastq_scratch_drive'))

    logger.info(qsub_cmd)
    os.system(qsub_cmd)

def add_sge_header(config,file,fcillumid,step):
    log_dir = config.get('locs','logs_dir')
    file.write('#! /bin/bash\n')
    file.write('#\n')
    file.write('#$ -S /bin/bash -cwd\n')
    file.write('#$ -o {}/{}_{}.out\n'.format(log_dir,fcillumid,step))
    file.write('#$ -e {}/{}_{}.err\n'.format(log_dir,fcillumid,step))
    file.write('#$ -V\n')
    file.write('#$ -M jb3816@cumc.columbia.edu\n')
    file.write('#$ -m bea\n')
    file.write('#\n')
    file.write('\n')

def create_post_bcl_script(config,archive_dir,runFolder,machine,fcillumid,address,postBCLCmd):
    #writes out the stage 2 script
    logger = logging.getLogger(__name__)

    post_bcl_script = open('/nfs/seqscratch1/Runs/%s/%s_%s_%s_post_bcl.sh' % (runFolder,machine,fcillumid,archive_dir),'w')
    add_sge_header(config,post_bcl_script,fcillumid,'post_bcl')
    post_bcl_script.write("export LD_LIBRARY_PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/lib:$LD_LIBRARY_PATH ; export PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/bin:$PATH \n")
    post_bcl_script.write(postBCLCmd + '\n')
    post_bcl_script.write('if [ $? -eq 0 ] ; then echo "Post BCL completed successfully" ; else echo "Post BCL failed"\n')
    post_bcl_script.write('/nfs/goldstein/software/mutt-1.5.23 -s "Post_BCL failure: %s %s" %s < /dev/null; exit 1; fi\n' % (machine,fcillumid,address))

def create_storage_script(config,archive_dir,runFolder,machine,fcillumid,address,storageCmd):
    storage_script = open('/nfs/seqscratch1/Runs/%s/%s_%s_%s_storage.sh' % (runFolder,machine,fcillumid,archive_dir),'w')
    add_sge_header(config,storage_script,fcillumid,'storage')
    storage_script.write(storageCmd + '\n')

def create_fastqc_script(config,archive_dir,runFolder,machine,fcillumid,address,fastqcCmd):
    add_sge_header(config,fastqc_script,fcillumid,'fastqc')
    fastqc_script.write(fastqCmd + '\n')
    fastqc_script.close()

def checkSeqsata(archive_dir):
    try:
        os.path.isdir(glob.glob('/nfs/%s' % archive_dir)[0])
    except:
        raise Exception('Input location does not exist!')

def arg_parser(config):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fcillumid", dest='fcillumid', required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-b','--bcl_drive', default=config.get('locs','bcl2fastq_scratch_drive'),
                        dest='bcl_drive',help="Specify scratch dir for bcl2fastq")
    parser.add_argument('-a','--archive_dir', default=config.get('locs','fastq_archive_drive'),
                        dest='archive_dir',help="Specify scratch dir for bcl2fastq")
    parser.add_argument("-r","--run", default=False, action="store_true",
                        help="Run bcl2fastq pipeline")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    return args


def main():
    config = get_config()
    args = arg_parser(config)
    run_info_dict = parse_run_parameters_xml(args.fcillumid)
    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    # Ex A00123 + B = A00123B 
    machine = run_info_dict['machine']
    setup_logging(machine,args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    logger.info('Running GAF_Pipeline.py')

    check_flowcell_complete(config.get('locs','bcl_dir'),run_info_dict['runFolder'])
    check_machine(machine,args.fcillumid,database)
    """Sequencing output is in the following format:
    [Date]_[IGM Machine Name]_[HiSeq Run Number]_[Flowcell ID]_[Project Name]
    However in cases of Illumina maintence they will re-image the machine and
    reset the output folder to Illumina default:
    [Date]_[Illumina Machine Name]_[HiSeq Run Number]_[Machine Side][Flowcell ID]
    """

    if args.run == True:
        logger.info('GAF_Pipeline.py in automated mode')
    submit(config,args,run_info_dict,database)

    #submit(run_info_dict['runFolder'],args.archive_dir,run_info_dict['runDate'],
    #       machine,args.fcillumid,config.get('locs','bcl2fastq_scratch_drive'))

if __name__ == '__main__':
    main()
