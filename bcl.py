# bcl.py
# Submits finished sequenced run for BCL to fastq conversion

import logging
import os
import utilities
from CHGV_mysql import (create_sss_bcl_2)
from glob import glob
from xml.etree import ElementTree

config = get_config()
args = parse_arguments()
run_info_dict = parse_run_parameters_xml(args.fcillumid)

def run_bcl2fastq(machine,database)
    logger = logging.getLogger(__name__)
    script_dir = config['scr_dir']
    bcl2fastq_loc = config['bcl2fastq_program']
    bcl_dir = '{}/{}'.format(config['bcl_dir'],run_info_dict['runFolder'])
    out_dir = '{}/{}_{}_{}_Unaligned'.format(config['bcl2fastq_scratch_dir'],
                                             run_info_dict['runDate'],
                                             machine,
                                             run_info_dict['FCIllumID'])
    script_loc = '{}/{}_{}_BCL.sh'.format(out_dir,run_info_dict['FCIllumID'],machine)

    check_exist_bcl_dir(out_dir)
    bcl2fastq_cmd = build_bcl2fastq_cmd(bcl2fastq_loc,bcl_dir,out_dir,database)
    print bcl2fastq_cmd
    logger.info(bcl2fastq_cmd)
    os.mkdir(out_dir)

    #Submit bcl job to the cluster
    os.system('cp {}/SGE_header {}'.format(script_dir,script_loc))
    with open(script_loc,'w') as bcl_script
        add_libraries_to_script(bcl_script)
        bcl_script.write(bcl2fastq_cmd + '\n')

    qsub_loc = config['qsub_program']
    qsub_cmd = ("cd {0} ;"
                "{} -cwd -v PATH -N {}_{}_{}_bcl {}"
               ).format(out_dir,
                        qsub_loc,
                        machine,
                        FCID,
                        bcl_dir.split('/')[2],
                        script_loc)

    if verbose == True:
        print cmd
    logger.info(cmd)
    status = os.system(cmd)
    logger.info(status)

def build_bcl2fastq_cmd(bcl2fastq_loc,bcl_dir,out_dir,database):
    bcl2fastq_cmd = ("{} --runfolder-dir {} --output-dir {} --barcode-mismatches 1' "
                    ).format(bcl2fastq_loc,bcl_dir,out_dir)
    """ Use if bcl or stats are missing.  They should never be missing unless
        there was a data transfer problem or corruption."""
    #bcl2fastq_cmd += '--ignore-missing-bcl --ignore-missing-stats '

    if sample_sheet: #for custom sample_sheets
        sss_qc(sample_sheet)
        bcl2fastq_cmd += '--sample-sheet '+sample_sheet +' '
    else:
        sample_sheets = glob('%s/*%s*.csv '.format(runPath,FCID))
        if len(sample_sheets) != 1:
            raise ValueError('Incorrect number of sequencing '
                             'sample sheets found in run folder!')
        else:
            sss_qc(sample_sheet)
            bcl2fastq_cmd += '--sample-sheet %s/*%s*.csv ' % (runPath,FCID)

    """ Tile specify what tiles on the flowcell can be converted to fastq
        This adds tiles parameter to BCL script if specified.  Multiple
        bcl2fastq runs might be required"""
    if tiles != False:
        bcl2fastq_cmd += '--tiles={}'.format(args.tiles)
        update_flowcell_tile_query = ("UPDATE Flowcell "
                                      "SET TileCommand=CONCAT_WS(';',TileCommand,'{0}') "
                                      "WHERE FCillumID='{1}'"
                                     ).format(tiles,FCID)
        logger.info(update_flowcell_tile_query)
        sequenceDB.execute(update_flowcell_tile_query)

    # Used to mask any bases on the read
    if base_mask != False:
        bcl2fastq_cmd += '--use-bases-mask {} '.format(base_mask)


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
    return bcl_script


def update_sample_status(database,FCID):
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
                                 ).format(userID,FCID)

    if verbose == True:
        print sample_status_update_query
    sequenceDB.execute(sample_status_update_query)
    logger.info(sample_status_update_query)

def check_fcillumid(inputtedFCID,xmlFCID):
    if inputtedFCID != xmlFCID:
        raise ValueError("FCIllumIDs don't match! {} != {}".format(inputtedFCID,xmlFCID))

def parse_run_parameters_xml(FCID):
    logger = logging.getLogger(__name__)
    runParametersXML = glob('/nfs/seqscratch1/Runs/*{}*'.format(FCID))
    if len(runParametersXML) != 1:
        raise ValueError("Too many run folders found!")
    else:
        runParametersXML = runParametersXML[0]

    with open(runParametersXML,'r') as xml
        run_info_dict = {}
        tree = ElementTree.parse(xml)
        run_info_dict['rtaVersion'] = tree.findall('RtaVersion')[0].text
        run_info_dict['experimentName'] = tree.findall('ExperimentName')[0].text
        run_info_dict['runFolder'] = tree.findall('RunId')[0].text
        run_info_dict['runDate'] = runFolder.split('_')[0]
        run_info_dict['machineName'] = tree.findall('InstrumentName')[0].text
        run_info_dict['ControlSoftwareVer'] = tree.findall('ApplicationVersion')[0].text
        for FCInfo in tree.findall('RfidsInfo'):
            run_info_dict['FCIllumID'] = FCInfo.find('FlowCellSerialBarcode')
            run_info_dict['machineSide'] = FCInfo.find('FlowCellPartNumber')

    return run_info_dict

def create_sss(noSSS):
    if noSSS == False:
        create_sss_bcl_ver2(machine,run_info_dict,database)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    group.add_argument("-f", "--fcillumid", type='str', required=True
                        help="Specify Illumina's Flowcell ID")
    group.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-b','--bcl_drive', default='seqscratch_ssd', type='str',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('-a','--archive_dir', default='igmdata01', type='str',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('--tiles', type='str',
                        help="Specify which tiles of the flowcel to convert")
    parser.add_argument('--use_bases_mask',  type='str',
                        help="Specify any base positions to mask")
    parser.add_argument('-s','--sss', type='str',
                        help="Specify your own sequencing sample sheet")
    parser.add_argument('--noStatus', type='str',
                        help="Do not update status")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v2.0')
    args=parser.parse_args()

def main():
    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    machine = machineName + machineSide # Ex A00123 + B = A00123B
    setup_logging(machine,args.fcillumid,args.bcl_drive)
    logger = logging.getLogger(__name__)
    logger.info("Starting bcl2fastq job creation...")
    print "Starting bcl2fastq job creation..."
    check_fcillumid(args.fcillumid,run_info_dict['FCIllumID'])
    check_machine(machine,args.fcillumid,database)
    check_rta_complete(run_folder_path)
    create_sss(args.noSSS)
    run_bcl2fastq(machine,database)
    if noStatus == False:
        update_sample_status(database,FCID)
    logger.info("bcl2fastq successfully started")
    print "bcl2fastq successfully started"

if __name__ == '__main__':
	main()
