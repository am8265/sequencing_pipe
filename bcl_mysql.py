#!/usr/bin/python
# bcl_mysql.py
# Joshua Bridgers
# 04/26/2012
# jb371@dm.duke.edu
#
# Submits finished sequenced run for BCL to fastq conversion

import csv
import getopt
import glob
import logging
import math
import subprocess
import sys
import traceback
from CHGV_mysql_3_6 import getTestSequenceDB, getSequenceDB, setup_logging
from datetime import datetime
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary

#gets current time and date

def cur_datetime():
    tDateBCL = datetime.now()
    tDateBCL = str(tDateBCL).split('.')
    DateBCL = tDateBCL[0]
    return DateBCL

def getRTAdate(pwd):
    cmd = ['ls','-l','--full-time','{}/RTAComplete.txt'.format(pwd)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    DateRTA = proc.stdout.read().decode().split(' ')
    date = DateRTA[5]
    time = DateRTA[6].split('.')[0]
    return '%s %s' % (date,time)

def getRead1Date(pwd):
    cmd = ['ls','{}/RunInfo.xml'.format(pwd),'-l','--time-style=full-iso']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    tmp = proc.stdout.read().decode().split()
    Read1Date = tmp[5] + ' ' + tmp[6].split('.')[0]
    return Read1Date

def ver_num(pwd,Machine):
    if Machine[0] == 'A':
        runParametersFile = 'RunParameters.xml'
    else:
        runParametersFile = 'runParameters.xml'
    RTAcmd = ['grep','-i','-m1','rta','{}/{}'.format(pwd,runParametersFile)]
    proc = subprocess.Popen(RTAcmd, stdout=subprocess.PIPE)
    tmp = proc.stdout.read().decode().split()
    RTAVer = tmp[0].split('>')[1].split('<')[0]
    HCScmd = ['grep','ApplicationVersion','{}/{}'.format(pwd,runParametersFile),'-m1']
    proc2 = subprocess.Popen(HCScmd, stdout=subprocess.PIPE)
    tmp2 = proc2.stdout.read().decode().split()
    HCSVer = tmp2[0].split('>')[1].split('<')[0]
    return RTAVer,HCSVer


def updateFC(sequenceDB,FCID,Machine,pwd,seqLoc):
    logger = logging.getLogger('updateFC')
    DateBcl = cur_datetime()
    DateRTA = getRTAdate(pwd)
    DateRead1 = getRead1Date(pwd)
    RTAVer,HCSver = ver_num(pwd,Machine)
    hasIndexRead = 1

    #print(DateRTA)

    sql = ("UPDATE Flowcell "
        "SET RTAVer='{}', HCSVer='{}', DateRead1='{}', DateRTA='{}', "
        "DateBcl='{}', SeqsataLoc='{}' "
        "WHERE FCillumID='{}'"
        ).format(RTAVer,HCSver,DateRead1,DateRTA,DateBcl,seqLoc,FCID)
    if verbose == True:
        print(sql)

    sequenceDB.execute(sql)
    logger.info(sql)


#Updates ClusterDensity for Lane Entries
def updateLane(sequenceDB,FCID,Machine,pwd,run_summary):
    logger = logging.getLogger('updateLane')
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
            ).format(ClustDen,ClustDenStDev,ClustPF,ClustPFStDev,str(LaneNum),FCID)

        if verbose == True:
            print(sql)
        sequenceDB.execute(sql)
        logger.info(sql)

    totalNumLanes = totalLanesCheck(sss_lanes,FCID)
    qmets(sequenceDB,pwd,totalNumLanes,Machine,FCID,run_summary)
    updateLnFraction(sequenceDB,FCID,pwd)

def totalLanesCheck(sss_lanes,FCID):
    """Check if # of lanes in generated sequencing sample sheet matches actual
    number"""
    if FCID[-3] == 'N': #Normal Flowcells have 8 lanes
        actualNumLanes = [8]
    elif FCID[0] == 'H': #Rapid Runs, NovaSeq S1-4 can have 2 or 4 lanes
        actualNumLanes = [2,4]
    else:
        raise Exception('Unhandled FCID for flowcell %s' % FCID)
    for lanes in actualNumLanes:
        match = 0
        if sss_lanes in actualNumLanes:
            match = 1
            return sss_lanes
    if match == 0 :
        raise Exception( 'Number of lanes in SSS is incorrect!')

def updateLnFraction(sequenceDB,FCID,pwd):
    logger = logging.getLogger('updateLnFraction')
    logger.debug('Running updateLnFraction')
    sql = ("SELECT l.DBID,CHGVID,l.FCID,l.lanenum,FCillumID "
        "FROM Lane l "
        "JOIN Flowcell f on l.FCID=f.FCID "
        "JOIN prepT p on p.prepID=l.prepID "
        "WHERE FCillumID='{0}'"
        ).format(FCID)
    sequenceDB.execute(sql)
    info = sequenceDB.fetchall()
    #print(info)

    sampleSheet = glob.glob('{}/*{}*.csv'.format(pwd,FCID))[0]
    sampleDict = {}
    with open(sampleSheet) as csvfile:
        skipAheadCSV = csvfile.readlines()[17:]
        reader = csv.DictReader(skipAheadCSV)
        for row in reader:
            #print(row)
            key = '{}_{}'.format(row['Sample_ID'],row['Lane'])
            sampleDict[key] = row
    for samp in info:
        key = '{}_{}'.format(samp['CHGVID'],samp['lanenum'])
        LnFraction = sampleDict[key]['Description'].split('_')[0]
        sql = ("UPDATE Lane l "
            "SET LnFraction={0} "
            "WHERE FCID={1} and LaneNum={2} and DBID={3}"
            ).format(LnFraction,samp['FCID'],samp['lanenum'],samp['DBID'])

        logger.info(sql)
        if verbose == True:
            print(sql)
        sequenceDB.execute(sql)

def qmets(sequenceDB,run_folder,total_lanes,machine,FCID,run_summary):
    logger = logging.getLogger('qmets')
    logger.info('Running loadQualityMetrics.pl')

    #checkes number of index reads on flowcell 
    indexCheckSQL = ("SELECT LenI1,LenI2 "
                    "from Flowcell "
                    "WHERE FCILLUMID='{}'"
                    ).format(FCID)
    sequenceDB.execute(indexCheckSQL)
    indexLengths = sequenceDB.fetchone()
    indexOneLength = indexLengths['LenI1']
    indexTwoLength = indexLengths['LenI2']

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
            perQ30I1 = run_summary.at(1).at(lanenum).percent_gt_q30()
        elif readNumberOffset == 2:
            perQ30I2 = run_summary.at(3).at(lanenum).percent_gt_q30()
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
            ).format(perQ30R1,perQ30R2,perQ30I1,errorR1,errorR2,percentAlignR1,percentAlignR2,LaneNum,FCID)

        if verbose == True:
            print(sql)
        sequenceDB.execute(sql)
        logger.info(sql)

def getMetricsSummary(run_folder):
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


def opts(argv):
    global verbose
    verbose = False
    global runPath
    runPath = ''

    try:
        opts,args = getopt.getopt(argv, "hvi:", ['help','verbose','input='])
    except getopt.GetoptError(err):
        print(str(err))
        usage()
    for o,a in opts:
        if o in ('-v','--verbose'):
            verbose = True
        elif o in ('-h','--help'):
            usage()
        elif o in ('-i','--input'):
            runPath = a
        else:
            assert False, "Unhandled argument present"

def main():
    #accessing mysql gaf database
    sequenceDB = getSequenceDB()
    opts(sys.argv[1:])

    pwd = '/nfs/seqscratch1/Runs/' + runPath
    info = pwd.split('/')[4].split('_')
    #print(info)
    FCID = info[3]
    Machine = info[1]
    seqloc = 'igmdata01'
    setup_logging(Machine,FCID,seqloc)
    print('BCL MySQL updates started')
    logger = logging.getLogger('main')
    logger.info('BCL MySQL updates started')
    logger.debug('Initializing Parameters: pwd:%s, FCID:%s, Machine:%s, seqsata_drive:%s', (pwd,FCID,Machine,seqloc))
    try:
        updateFC(sequenceDB,FCID,Machine,pwd,seqloc)
        run_summary = getMetricsSummary(pwd)
        updateLane(sequenceDB,FCID,Machine,pwd,run_summary)

        #mysql updates happen now on BCL
        logger.info('BCL MySQL updates completed')
        print('BCL MySQL updates completed')
        sequenceDB.execute('COMMIT')
        sequenceDB.close()
    except:
        traceback.print_exc()
        sequenceDB.execute('ROLLBACK')
        sequenceDB.close()
        logger.info('BCL MySQL updates failure')
        print('BCL MySQL update failure')


main()
