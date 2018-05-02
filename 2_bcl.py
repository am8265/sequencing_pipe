# bcl.py
# Submits finished sequenced run for BCL to fastq conversion

import argparse
import logging
import os
from datetime import datetime
from utilities import *
import warnings


def run_bcl2fastq(args,run_info_dict,config,database,verbose):

    fcillumid = args.fcillumid
    logger = logging.getLogger(__name__)
    #script_dir = config.get('locs','scr_dir')
    bcl2fastq_loc = config.get('programs','bcl2fastq_program')
    bcl_dir = '{}/{}'.format(config.get('locs','bcl_dir'),run_info_dict['runFolder'])
    out_dir = run_info_dict['out_dir'] 
    script_loc = '{}/{}_{}_BCL.sh'.format(out_dir,run_info_dict['FCIllumID'],run_info_dict['machine'])

    check_exist_bcl_dir(fcillumid,out_dir,database)
    sss_loc = create_sss_from_database(args.fcillumid,run_info_dict['machine'],
                                       run_info_dict,config,database)
    if not os.path.exists(sss_loc):
        raise("yo' punk. where's my sample sheet? : '{}".format(sss_loc))

    print("we should have an sample sheet now? : '{}'".format(sss_loc))
    
    
    cp_cmd= ( 'cp {} /nfs/{}/Sequencing_SampleSheets/' ).format(
      sss_loc,config.get('locs','fastq_archive_drive') 
    )

    if os.system(cp_cmd) != 0:
        print("i beg your pardon '{}'".format(cp_cmd))

    logger.info("Using SSS: {}".format(sss_loc))

    bcl2fastq_cmd = build_bcl2fastq_cmd(args,fcillumid,bcl2fastq_loc,sss_loc,bcl_dir,out_dir,database)
    print(bcl2fastq_cmd)
    logger.info(bcl2fastq_cmd)
    if os.path.exists(out_dir):
        warnings.warn("directory {} already exists".format(out_dir))
    else:    
        try:
            os.mkdir(out_dir,0o770)
        except:
            print("error creating {}".format(out_dir))

    #Submit bcl job to the cluster
    with open(script_loc,'w') as bcl_script:
        add_sge_header_to_script(bcl_script,fcillumid)
        add_libraries_to_script(bcl_script)
        bcl_script.write(bcl2fastq_cmd + '\n')
        bcl_script.write(' if [ $? -eq 0 ] ; then touch {}/bcl_complete ; fi\n'.format(out_dir))
    qsub_loc = config.get('programs','qsub_program')
    qsub_cmd = ("cd {} ; {} -v PATH  {}"
               ).format(out_dir,qsub_loc,script_loc)
    if verbose == True:
        print(qsub_cmd)
    logger.info(qsub_cmd)
    pid = os.system(qsub_cmd)
    if pid != 0:
        raise Exception("{} didnt execute correctly. return value is {}".format(qsub_cmd,pid))
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
    if args.use_bases_mask is not None:
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
    bcl_script.write("#$ -M {}@cumc.columbia.edu".format(get_user()))
    bcl_script.write("#$ -m bea\n")
    ##### perhaps we should mention that this name is rather important for when we get to -hold_jid bcl_{} in the wrapper...?!?
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

def set_flowcell_samples_and_flowcell_status_to_bclstarted(database,fcillumid,verbose):
    logger = logging.getLogger(__name__)
    userID = get_user_id(database)
    sample_status_insert_query = ("INSERT INTO statusT "
                                  "(STATUS_TIME,STATUS,PREPID,USERID,POOLID,SEQID) "
                                  "SELECT DISTINCT unix_timestamp(),"
                                  "'BCL Started',pt.prepID,{},0,0 "
                                  "FROM Flowcell f "
                                  "JOIN Lane l ON l.FCID=f.FCID "
                                  "JOIN prepT pt ON pt.prepID=l.prepID "
                                  "WHERE FCillumid='{}'"
                                 ).format(userID,fcillumid)
    prepT_status_update_query = """UPDATE prepT p
                                   JOIN Lane l ON p.PREPID=l.PREPID
                                   JOIN Flowcell f ON f.FCID=l.FCID
                                   SET STATUS='BCL Started', status_time=unix_timestamp()
                                   WHERE FCILLUMID='{}'
                                """.format(fcillumid)
    if verbose == True:
        print(sample_status_insert_query)
        print(prepT_status_update_query)
    run_query(sample_status_insert_query,database)
    logger.info(sample_status_insert_query)
    run_query(prepT_status_update_query,database)
    logger.info(prepT_status_update_query)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('--tiles',
                        help="Specify which tiles of the flowcell to convert")
    parser.add_argument('--use_bases_mask', dest='use_bases_mask',
                        help="""Specify any base positions to mask"
                             Ex: y150n,I8n,y150n""")
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

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'

    run_info_dict = parse_run_parameters_xml(args.fcillumid,database)
    setup_logging(run_info_dict['machine'],args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    logger.info("Starting bcl2fastq job creation for Flowcell: {}".format(args.fcillumid))
    print("Starting bcl2fastq job creation for Flowcell: {}".format(args.fcillumid))

    # check_fcillumid(args.fcillumid,run_info_dict['FCIllumID'])
    # def check_fcillumid(inputted_fcillumid,xml_fcillumid):

    ######## lock flowcell supplied doesn't match that of run dir flowcell
    if args.fcillumid != run_info_dict['FCIllumID']:
        update_pipelinecomplete( args.fcillumid, '-2', database )
        raise Exception("FCIllumIDs don't match! {} != {}".format(args.fcillumid,run_info_dict['FCIllumID']))

    check_machine(run_info_dict['machine'],args.fcillumid,database)

    check_flowcell_complete(args.fcillumid,config.get('locs','bcl_dir'),run_info_dict['runFolder'],run_info_dict['type'],database)

    ######## lock flowcell
    mofo=update_pipelinecomplete( args.fcillumid, '-1', database )

    print(mofo)

    #check_sss(sss_loc)

    ######## finally, we do something...
    run_bcl2fastq(args,run_info_dict,config,database,args.verbose)

    if args.noStatus == False:
        set_flowcell_samples_and_flowcell_status_to_bclstarted(database,args.fcillumid,args.verbose)

    logger.info("bcl2fastq successfully started")
    print("bcl2fastq successfully started")

if __name__ == '__main__':
    os.umask(int("002",8))
    main()
