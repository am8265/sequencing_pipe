import argparse
import csv
import getopt
import logging
import math
import subprocess
import sys
import traceback
from datetime import datetime
from glob import glob
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary
from utilities import *

def get_cur_datetime():
    tDateBCL = datetime.now()
    tDateBCL = str(tDateBCL).split('.')
    DateBCL = tDateBCL[0]
    return DateBCL

def get_rta_date(run_folder):
    cmd = ['ls','-l','--full-time','{}/RTAComplete.txt'.format(run_folder)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    DateRTA = proc.stdout.read().decode().split(' ')
    date = DateRTA[5]
    time = DateRTA[6].split('.')[0]
    return '%s %s' % (date,time)

def get_read1_date(run_folder):
    cmd = ['ls','{}/RunInfo.xml'.format(run_folder),'-l','--time-style=full-iso']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    tmp = proc.stdout.read().decode().split()
    Read1Date = tmp[5] + ' ' + tmp[6].split('.')[0]
    return Read1Date

def fill_in_some_flowcell_table_info(database,fcillumid,run_info_dict,config,run_folder,verbose):
    logger = logging.getLogger(__name__)
    fastq_archive_loc = config.get('locs','fastq_archive_drive')
    rta_ver = run_info_dict['RtaVersion']
    hsc_ver = run_info_dict['ControlSoftwareVer']
    DateBcl = get_cur_datetime()
    DateRTA = get_rta_date(run_folder)
    DateRead1 = get_read1_date(run_folder)

    sql = ("UPDATE Flowcell "
        "SET RTAVer='{}', HCSVer='{}', DateRead1='{}', DateRTA='{}', "
        "DateBcl='{}', SeqsataLoc='{}' "
        "WHERE FCillumID='{}'"
        ).format(rta_ver,hsc_ver,DateRead1,DateRTA,
                 DateBcl,fastq_archive_loc,fcillumid)
    if verbose == True:
        print(sql)

    run_query(sql,database)
    logger.info(sql)


def update_per_sample_illumina_metrics(database,fcillumid,Machine,run_folder,run_summary,verbose):
    """Updates ClusterDensity for Lane Entries"""
    logger = logging.getLogger(__name__)
    logger.debug('Running updateLane')

    sss_lanes = run_summary.lane_count()

    for LaneNum in list(range(0,sss_lanes)):
        ClustDen = run_summary.at(0).at(LaneNum).density().mean() / 1000
        ClustDenStDev = run_summary.at(0).at(LaneNum).density().stddev() / 1000
        ClustPF = run_summary.at(0).at(LaneNum).density_pf().mean() / 1000
        ClustPFStDev = run_summary.at(0).at(LaneNum).density_pf().stddev() /1000


        #print(ClustPF,ClustPFStDev)
        #interOP LaneNum is zero-indexed while sql LaneNum is not.
        LaneNum += 1
        sql = ("UPDATE Lane l "
            "JOIN Flowcell f on l.FCID=f.FCID "
            "SET ClustDen='{0}',ClustDenStDev='{1}',"
            "ClusterPF='{2}',ClusterPFStDev=round('{3}',3) "
            "WHERE LaneNum={4} AND f.FCillumID='{5}'"
            ).format(ClustDen,ClustDenStDev,ClustPF,ClustPFStDev,str(LaneNum),fcillumid)

        if verbose == True:
            print(sql)
        run_query(sql,database)
        logger.info(sql)

    totalNumLanes = totalLanesCheck(sss_lanes,fcillumid)
    qmets(database,run_folder,totalNumLanes,Machine,fcillumid,run_summary,verbose)
    update_per_sample_illumina_metrics_fraction(database,fcillumid,run_folder,verbose)

def totalLanesCheck(sss_lanes,fcillumid):
    """Check if # of lanes in generated sequencing sample sheet matches actual
    number"""
    if fcillumid[-3] == 'N': #Normal Flowcells have 8 lanes
        actualNumLanes = [8]
    elif fcillumid[0] == 'H': #Rapid Runs, NovaSeq S1-4 can have 2 or 4 lanes
        actualNumLanes = [2,4]
    else:
        raise Exception('Unhandled fcillumid for flowcell %s' % fcillumid)
    for lanes in actualNumLanes:
        match = 0
        if sss_lanes in actualNumLanes:
            match = 1
            return sss_lanes
    if match == 0 :
        raise Exception( 'Number of lanes in SSS is incorrect!')

def update_per_sample_illumina_metrics_fraction(database,fcillumid,run_folder,verbose):
    logger = logging.getLogger(__name__)
    logger.debug('Running update_per_sample_illumina_metrics_fraction')
    sql = ("SELECT l.DBID,CHGVID,l.FCID,l.lanenum,FCillumID "
        "FROM Lane l "
        "JOIN Flowcell f on l.FCID=f.FCID "
        "JOIN prepT p on p.prepID=l.prepID "
        "WHERE FCillumID='{0}'"
        ).format(fcillumid)
    info = run_query(sql,database)
    #print(info)

    sampleSheet = glob('{}/*{}*.csv'.format(run_folder,fcillumid))
    if len(sampleSheet) == 0:
        raise IOError("No sample sheet found!")
    elif len(sampleSheet) > 1:
        raise IOError("too many sheets found!")

    sss_samples_info_dict = {}
    #print(sampleSheet)
    with open(sampleSheet[0]) as csvfile:
        skipped_ahead_sample_sheet = csvfile.readlines()[17:]
        reader = csv.DictReader(skipped_ahead_sample_sheet)
        for row in reader:
            key_from_sss = '{}_{}'.format(row['Sample_ID'],row['Lane'])
            sss_samples_info_dict[key_from_sss] = row
    if sss_samples_info_dict == {}:
        raise ValueError("No samples were found in the sample sheet!")
    for samp in info:
        key_from_db = '{}_{}'.format(samp['CHGVID'],samp['lanenum'])
        LnFraction = sss_samples_info_dict[key_from_db]['Description'].split('_')[0]
        #LnFraction = sss_samples_info_dict[key_from_db]['Description'].split("'")[1].split('_')[0]

        sql = ("UPDATE Lane l "
            "SET LnFraction={0} "
            "WHERE FCID={1} and LaneNum={2} and DBID={3}"
            ).format(LnFraction,samp['FCID'],samp['lanenum'],samp['DBID'])

        logger.info(sql)
        if verbose == True:
            print(sql)
        run_query(sql,database)

def qmets(database,run_folder,total_lanes,machine,fcillumid,run_summary,verbose):
    logger = logging.getLogger(__name__)
    logger.info('Running loadQualityMetrics.pl')

    #checkes number of index reads on flowcell
    indexLengths = run_query(GET_FLOWCELL_RECIPE.format(fcillumid=fcillumid),database)
    indexOneLength = indexLengths[0]['LenI1']
    indexTwoLength = indexLengths[0]['LenI2']

    for LaneNum in list(range(0,total_lanes)):
        readNumberOffset = 0

        if indexTwoLength != 0: # no read2 index (standard FC run)
            readNumberOffset += 2 # 4 total reads
        elif indexOneLength != 0:
            readNumberOffset += 1 # 3 total reads
        elif indexOneLength == 0:
            readNumberOffset = 0 # 2 total reads 
        else:
            raise Exception("Unhandled index reads")

        perQ30R1 = run_summary.at(0).at(LaneNum).percent_gt_q30()
        errorR1 = run_summary.at(0).at(LaneNum).error_rate().mean()
        percentAlignR1 = run_summary.at(0).at(LaneNum).percent_aligned().mean()
        perQ30R2 = run_summary.at(1 + readNumberOffset).at(LaneNum).percent_gt_q30()
        errorR2 = run_summary.at(1 + readNumberOffset).at(LaneNum).error_rate().mean()
        percentAlignR2 = run_summary.at(1 + readNumberOffset).at(LaneNum).percent_aligned().mean()

        if readNumberOffset == 0:
            perQ30I1 = '0'
            perQ30I2 = '0'
        elif readNumberOffset == 1:
            perQ30I1 = run_summary.at(1).at(LaneNum).percent_gt_q30()
        elif readNumberOffset == 2:
            perQ30I2 = run_summary.at(3).at(LaneNum).percent_gt_q30()
        else:
            raise Exception("Unhandled readNumberOffset value")

        # when read2 fails, perQ30R2 == 'nan' 
        if math.isnan(perQ30R2):
            perQ30R2 = '0'

        #SAV's interOP LaneNum is zero-indexed while sql LaneNum is not.
        LaneNum += 1
        sql = ("UPDATE Lane l "
            "join Flowcell f on l.FCID=f.FCID "
            "SET perQ30R1={},perQ30R2={},perQ30I1={},"
            "errorR1={},errorR2={},percentAlignR1={},percentAlignR2={} "
            "WHERE LaneNum={} AND f.FCillumID='{}'"
            ).format(perQ30R1,perQ30R2,perQ30I1,errorR1,errorR2,percentAlignR1,percentAlignR2,LaneNum,fcillumid)

        if verbose == True:
            print(sql)
        run_query(sql,database)
        logger.info(sql)

def parse_illumina_metrics_from_run_dir(run_folder):
    #run_folder is suppose to be a raw string literal
    run_folder = run_folder.replace("\\", "\\\\")
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_folder = run_metrics.read(run_folder, valid_to_load)
    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)
    return summary


def usage():
	print('-h, --help\t\tShows this help message and exit')
	print('-v, --verbose\t\tPrints out MYSQL injection commands')
	sys.exit(2)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", dest='fcillumid', required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
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
    # Ex A00123 + B = A00123B 
    setup_logging(run_info_dict['machine'],args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    print('BCL MySQL updates started')
    logger.info('BCL MySQL updates started')
    try:

        bcl_dir = config.get('locs','bcl_dir')

        run_folder = bcl_dir + run_info_dict['runFolder']

        fill_in_some_flowcell_table_info(database,args.fcillumid,run_info_dict, config,run_folder,args.verbose)

        run_summary = parse_illumina_metrics_from_run_dir(run_folder)

        update_per_sample_illumina_metrics(database,args.fcillumid,run_info_dict['machine'], run_folder,run_summary,args.verbose)

        #mysql updates happen now on BCL
        logger.info('BCL MySQL updates completed')
        print('BCL MySQL updates completed')
    except:
        traceback.print_exc()
        logger.info('BCL MySQL updates failure')
        print('BCL MySQL update failure')


main()
