#!/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python
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
    address = '{}@cumc.columbia.edu'.format(get_user())

    runFolder = run_info_dict['runFolder']
    machine = run_info_dict['machine']
    fcillumid=args.fcillumid

    ##### please get ride of me!
    python36_program = config.get('programs','python36_program')

    #scriptLoc = '/nfs/goldstein/software/sequencing_pipe/dev/' # config.get('locs','scr_dir')
    rarp = (python36_program, os.path.dirname(os.path.realpath(__file__)), fcillumid)

    update_flowcell_metrics = "{} {}/1_update_flowcell_metrics.py -f {}".format( *rarp )

    BCLCmd                  = "{} {}/2_bcl.py -f {}".format(      *rarp )
    update_yields           = "{} {}/3a_update_yields_and_fractions.py -f {}".format( *rarp )
    storageCmd              = "{} {}/3b_storage.py -f {}".format(  *rarp )

    if args.test == True:
        BCLCmd += ' --test'; update_flowcell_metrics += ' --test'; update_yields += ' --test'; storageCmd += ' --test'

    if args.run == False:
        #prints out all the commands to run manually
        print
        print("="*35+'Scripts Commands'+"="*35)
        print(update_flowcell_metrics)

        print(BCLCmd)
        print(update_yields)
        print(storageCmd)

        print("="*86)
        print
        print(('Log file: {}'
            ).format(logging.getLoggerClass().root.handlers[0].baseFilename))

    else:
        ##### run BCLCmd, update_flowcell_metrics locally. submit 
        #run bcl2fastq
        out_dir = run_info_dict['out_dir']

        check_exist_bcl_dir(fcillumid,out_dir,database)
        if os.system(BCLCmd) != 0: #bcl2fastq auto submits to cluster
            raise Exception("\n\nreally? : '{}'".format(BCLCmd))
        logger.info(BCLCmd)

        x=os.system(update_flowcell_metrics)
        if x != 0: #quick msyql updates.  I
            raise Exception("\n\nwooooooow, look, we have a return value '{}' from '{}'".format(x,update_flowcell_metrics))
        logger.info(update_flowcell_metrics)

        #run post_bcl2fastq
        create_post_bcl_script(config,archive_dir,runFolder,machine,fcillumid,address,update_yields)
        run_qsub_command(config,fcillumid,'post_bcl',run_info_dict)

        create_storage_script(config,archive_dir,runFolder,machine,fcillumid,address,storageCmd)
        run_qsub_command(config,fcillumid,'storage',run_info_dict)

def run_qsub_command(config,fcillumid,step,run_info_dict):
    logger = logging.getLogger(__name__)
    qsub_cmd = ('{qsub_program} {bcl_dir}/{runFolder}/{machine}_{fcillumid}_{step}.sh'
               ).format(qsub_program=config.get('programs','qsub_program'),
                        machine=run_info_dict['machine'],
                        fcillumid=fcillumid,
                        bcl_dir=config.get('locs','bcl_dir'),
                        runFolder=run_info_dict['runFolder'],
                        step=step)
    print(qsub_cmd)
    logger.info(qsub_cmd)
    op = os.system(qsub_cmd)
    if op != 0:
        raise Exception("{} didn't execute".format(qsub_cmd))

def add_sge_header(config,file,fcillumid,step,hold):
    log_dir = config.get('locs','logs_dir')
    file.write('#! /bin/bash\n')
    file.write('#\n')
    file.write('#$ -S /bin/bash -cwd\n')
    file.write('#$ -o {}/{}_{}.out\n'.format(log_dir,fcillumid,step))
    file.write('#$ -e {}/{}_{}.err\n'.format(log_dir,fcillumid,step))
    file.write('#$ -V\n')
    file.write('#$ -N {}_{}\n'.format(step,fcillumid))
    file.write('#$ -M {}@cumc.columbia.edu\n'.format(get_user()))
    file.write('#$ -m bea\n')
    file.write('#\n')

    if hold == True:
        file.write('#$ -hold_jid bcl_{}\n'.format(fcillumid))
    file.write('\n')

def create_post_bcl_script(config,archive_dir,runFolder,machine,fcillumid,address,update_yields):
    #writes out the stage 2 script
    logger = logging.getLogger(__name__)
    script_loc = ('/{bcl_dir}/{runFolder}/{machine}_{fcillumid}_post_bcl.sh'
                 ).format(bcl_dir=config.get('locs','bcl_dir'),
                          runFolder=runFolder,
                          machine=machine,
                          fcillumid=fcillumid)
    post_bcl_script = open(script_loc,'w')
    add_sge_header(config,post_bcl_script,fcillumid,'post_bcl',True)
    post_bcl_script.write("export LD_LIBRARY_PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/lib:$LD_LIBRARY_PATH\n")
    post_bcl_script.write("export PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/bin:$PATH \n")
    post_bcl_script.write(update_yields + ' -v\n')
    post_bcl_script.write('if [ $? -eq 0 ] \n')
    post_bcl_script.write('then echo "Post BCL completed successfully"\n')
    post_bcl_script.write('else echo "Post BCL failed"\n')
    post_bcl_script.write('/nfs/goldstein/software/mutt-1.5.23/bin/mutt -s "Post_BCL failure: %s %s" %s < /dev/null; exit 1; fi\n' % (machine,fcillumid,address))

def create_storage_script(config,archive_dir,runFolder,machine,fcillumid,address,storageCmd):
    script_loc = ('/{bcl_dir}/{runFolder}/{machine}_{fcillumid}_storage.sh'
                 ).format(bcl_dir=config.get('locs','bcl_dir'),
                          runFolder=runFolder,
                          machine=machine,
                          fcillumid=fcillumid)
    storage_script = open(script_loc,'w')
    add_sge_header(config,storage_script,fcillumid,'storage',True)
    storage_script.write('export LD_LIBRARY_PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/lib:$LD_LIBRARY_PATH\n')
    storage_script.write(storageCmd + ' -v \n')

def checkSeqsata(archive_dir):
    try:
        os.path.isdir(glob.glob('/nfs/%s' % archive_dir)[0])
    except:
        raise Exception('Input location does not exist!')

def arg_parser(config):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", dest='fcillumid', required=True,
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
    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    run_info_dict = parse_run_parameters_xml(args.fcillumid,database)
    machine = run_info_dict['machine'] # Ex A00123 + B = A00123B 
    setup_logging(machine,args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    logger.info('Running bcl2fastq_pipeline_wrapper.py')

    check_flowcell_complete(args.fcillumid,config.get('locs','bcl_dir'),
                            run_info_dict['runFolder'],run_info_dict['type'],database)
    check_machine(machine,args.fcillumid,database)
    """Sequencing output is in the following format:
    [Date]_[IGM Machine Name]_[HiSeq Run Number]_[Flowcell ID]_[Project Name]
    However in cases of Illumina maintence they will re-image the machine and
    reset the output folder to Illumina default:
    [Date]_[Illumina Machine Name]_[HiSeq Run Number]_[Machine Side][Flowcell ID]
    """

    if args.run == True:
        logger.info('bcl2fastq_pipeline_wrapper.py in automated mode')
    submit(config,args,run_info_dict,database)

    #submit(run_info_dict['runFolder'],args.archive_dir,run_info_dict['runDate'],
    #       machine,args.fcillumid,config.get('locs','bcl2fastq_scratch_drive'))

if __name__ == '__main__':
    os.umask(int("002",8))
    main()
