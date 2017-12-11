# bcl.py
# Submits finished sequenced run for BCL to fastq conversion

import argparse
import logging
import os
from datetime import datetime
from glob import glob
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary
from utilities import *

def run_bcl2fastq(args,run_info_dict,config,sss_loc,database,verbose):
    fcillumid = args.fcillumid
    logger = logging.getLogger(__name__)
    script_dir = config.get('locs','scr_dir')
    bcl2fastq_loc = config.get('programs','bcl2fastq_program')
    bcl_dir = '{}/{}'.format(config.get('locs','bcl_dir'),run_info_dict['runFolder'])
    out_dir = '{}/{}_{}_{}_Unaligned'.format(config.get('locs','bcl2fastq_scratch_dir'),
                                             run_info_dict['runDate'],
                                             run_info_dict['machine'],
                                             run_info_dict['FCIllumID'])
    script_loc = '{}/{}_{}_BCL.sh'.format(out_dir,run_info_dict['FCIllumID'],run_info_dict['machine'])

    check_exist_bcl_dir(out_dir)
    bcl2fastq_cmd = build_bcl2fastq_cmd(args,fcillumid,bcl2fastq_loc,sss_loc,bcl_dir,out_dir,database)
    print(bcl2fastq_cmd)
    logger.info(bcl2fastq_cmd)
    os.mkdir(out_dir)

    #Submit bcl job to the cluster
    with open(script_loc,'w') as bcl_script:
        add_sge_header_to_script(bcl_script,fcillumid)
        add_libraries_to_script(bcl_script)
        bcl_script.write(bcl2fastq_cmd + '\n')

    qsub_loc = config.get('programs','qsub_program')
    qsub_cmd = ("cd {} ; {} -v PATH  {}"
               ).format(out_dir,qsub_loc,script_loc)
    if verbose == True:
        print(qsub_cmd)
    logger.info(qsub_cmd)
    pid = os.system(qsub_cmd)
    logger.info(pid)

def build_bcl2fastq_cmd(args,fcillumid,bcl2fastq_loc,sss_loc,bcl_dir,out_dir,database):
    logger = logging.getLogger(__name__)
    bcl2fastq_cmd = ("{} --runfolder-dir {} --output-dir {} --barcode-mismatches 1 "
                    ).format(bcl2fastq_loc,bcl_dir,out_dir)
    """ Use if bcl or stats are missing.  They should never be missing unless
        there was a data transfer problem or corruption."""
    #bcl2fastq_cmd += '--ignore-missing-bcl --ignore-missing-stats '

    bcl2fastq_cmd += '--sample-sheet '+ sss_loc +' '
    """ Tile specify what tiles on the flowcell can be converted to fastq
        This adds tiles parameter to BCL script if specified.  Multiple
        bcl2fastq runs might be required"""
    if args.tiles:
        bcl2fastq_cmd += '--tiles={} '.format(args.tiles)
        update_flowcell_tile_query = ("UPDATE Flowcell "
                                      "SET TileCommand=CONCAT_WS(';',TileCommand,'{0}') "
                                      "WHERE FCillumID='{1}'"
                                     ).format(args.tiles,fcillumid)
        logger.info(update_flowcell_tile_query)
        run_query(update_flowcell_tile_query,database)

    # Used to mask any bases on the read
    if args.use_bases_mask != False:
        bcl2fastq_cmd += '--use-bases-mask {} '.format(args.use_bases_mask)
    return bcl2fastq_cmd

def add_sge_header_to_script(bcl_script,fcillumid):
    bcl_script.write("#! /bin/bash\n")
    bcl_script.write("#\n")
    bcl_script.write("#$ -S /bin/bash -cwd\n")
    bcl_script.write("#$ -j y\n")
    bcl_script.write("#$ -o nohup.sge\n")
    bcl_script.write("#$ -e error.sge\n")
    bcl_script.write("#$ -V\n")
    bcl_script.write("#$ -M jb3816@cumc.columbia.edu\n")
    bcl_script.write("#$ -m bea\n")
    bcl_script.write("#$ -N bcl_{}\n".format(fcillumid))

def add_libraries_to_script(bcl_script):
    bcl_script.write("export PATH=/nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/bin:"
                     "/nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/lib:"
                     "/nfs/goldstein/software/libxml2-2.9.4-x86_64/bin:"
                     "/nfs/goldstein/software/libxml2-2.9.4-x86_64/lib:"
                     "/nfs/goldstein/software/binutils-2.26.1/bin:"
                     "/nfs/goldstein/software/binutils-2.26.1/lib:"
                     "/nfs/goldstein/software/boost_1_54_0/include:"
                     "/nfs/goldstein/software/boost_1_54_0/lib:"
                     "/nfs/goldstein/software/gcc-4.9.3/bin:"
                     "/nfs/goldstein/software/gcc-4.9.3/lib64:"
                     "$PATH\n")
    bcl_script.write("export LD_LIBRARY_PATH="
                     "/nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/lib:"
                     "/nfs/goldstein/software/libxml2-2.9.4-x86_64/lib:"
                     "/nfs/goldstein/software/binutils-2.26.1/lib:"
                     "/nfs/goldstein/software/boost_1_54_0/lib:"
                     "/nfs/goldstein/software/gcc-4.9.3/lib64:"
                     "$LD_LIBRARY_PATH\n")

def update_sample_status(database,fcillumid,verbose):
    logger = logging.getLogger(__name__)
    userID = get_user_id(database)
    sample_status_update_query = ("INSERT INTO statusT "
                                  "(CHGVID,status_time,status,DBID,prepID,userID) "
                                  "SELECT DISTINCT(pt.CHGVID),unix_timestamp(),"
                                  "'BCL',pt.DBID,pt.prepID,{0} "
                                  "FROM Flowcell f "
                                  "JOIN Lane l ON l.FCID=f.FCID "
                                  "JOIN prepT pt ON pt.prepID=l.prepID "
                                  "WHERE FCillumid='{1}'"
                                 ).format(userID,fcillumid)

    if verbose == True:
        print(sample_status_update_query)
    run_query(sample_status_update_query,database)
    logger.info(sample_status_update_query)

def check_fcillumid(inputted_fcillumid,xml_fcillumid):
    if inputted_fcillumid != xml_fcillumid:
        raise ValueError("FCIllumIDs don't match! {} != {}".format(inputted_fcillumid,xml_fcillumid))

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-b','--bcl_drive', default='seqscratch_ssd',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('-a','--archive_dir', default='igmdata01',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('--tiles',
                        help="Specify which tiles of the flowcell to convert")
    parser.add_argument('--use_bases_mask', action='store_true', default=False,
                        help="Specify any base positions to mask")
    parser.add_argument('--sss',
                        help="Specify your own sequencing sample sheet")
    parser.add_argument('--noStatus', action='store_true', default=False,
                        help="Do not update status")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    return args

def main():
    config = get_config()
    args = parse_arguments()
    run_info_dict = parse_run_parameters_xml(args.fcillumid)

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    setup_logging(run_info_dict['machine'],args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    logger.info("Starting bcl2fastq job creation for Flowcell: {}".format(args.fcillumid))
    print("Starting bcl2fastq job creation for Flowcell: {}".format(args.fcillumid))

    check_fcillumid(args.fcillumid,run_info_dict['FCIllumID'])
    check_machine(run_info_dict['machine'],args.fcillumid,database)
    check_flowcell_complete(config.get('locs','bcl_dir'),run_info_dict['runFolder'])

    if args.sss_loc is None:
        sss_loc = create_sss_from_database(args.fcillumid,run_info_dict['machine'],run_info_dict,config,database)
        cp_cmd= ('cp {} /nfs/igmdata01/Sequencing_SampleSheets/'
             ).format(sss_loc)
        os.system(cp_cmd)

    else:
        sss_loc = args.sss_loc
        print("Using SSS: {}".format(sss_loc))
    #check_sss(sss_loc)
    run_bcl2fastq(args,run_info_dict,config,sss_loc,database,args.verbose)
    if args.noStatus == False:
        update_sample_status(database,args.fcillumid,args.verbose)
    logger.info("bcl2fastq successfully started")
    print("bcl2fastq successfully started")

if __name__ == '__main__':
	main()
