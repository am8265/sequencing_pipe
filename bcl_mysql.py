#!/usr/bin/python
# bcl_mysql.py
# Joshua Bridgers
# 04/26/2012
# jb371@dm.duke.edu
#
# Submits finished sequenced run for BCL to fastq conversion

import glob
import os
import getopt
import logging
import sys
import traceback
from CHGV_mysql import *
from commands import getoutput
from MySQLdb import Connect
from datetime import datetime
from CHGV_mysql import getUserID

#gets current time and date


def cur_datetime():
        tDateBCL = datetime.now()
        tDateBCL = str(tDateBCL).split('.')
        DateBCL = tDateBCL[0]
        return DateBCL

def getRTAdate(pwd):
	DateRTA = getoutput('ls -l --full-time %s/RTAComplete.txt' % pwd).split(' ')
	date = DateRTA[5]
	time = DateRTA[6].split('.')[0]
	return '%s %s' % (date,time)

def getRead1Date(FCID):
	Read1Date = getoutput('ls -l /nfs/seqscratch1/Runs/*'+FCID+'*/RunInfo.xml --time-style=full-iso 2>/dev/null | cut -d\  -f6,7 | sort -rk6 | head -1 | cut -d. -f1')
	return Read1Date

def ver_num(pwd):
    tRTA = getoutput('grep RTA %s/ImageAnalysis_Netcopy_complete.txt' % pwd)
    tRTA2 = tRTA.split(' ')
    RTAVer = tRTA2[2].strip()

    HCSver = getoutput('grep ApplicationVersion %s/runParameters.xml -m1 | cut -d\> -f2 | cut -d\< -f1' % pwd)
    #print RTAVer,HCSver
    return RTAVer,HCSver

def getReads(sequenceDB,FCID):
    sequenceDB.execute("SELECT recipe from Flowcell where FCillumID=%s", FCID)
    Recipe = sequenceDB.fetchone()
    Recipe = Recipe[0]

    #generated str to be inserted into mysql statement.
    if Recipe == '1':
        return "LenR1='101', LenR2='101', LenI1='7'" #HiSeq 2000 V3 and Rapid V1 chemistry
    elif Recipe =='2':
        return "LenR1='101', LenR2='101'," #No index
    elif Recipe =='3':
        return "LenR1='101', LenR2='101', LenI1='7', LenI2='7'" #Dual Indexed v3, Nova S3,S4
    elif Recipe =='4':
        return "LenR1='100', LenR2='100', LenI1='9'" #Extended Index read
    elif Recipe =='5':
        return "LenR1='126', LenR2='126', LenI1='7'" #HiSeq 2500
    elif Recipe =='6':
        return "LenR1='151', LenR2='151', LenI1='7'" #HiSeq X,Nova S1,S2,S3,S4
    elif Recipe =='8':
        return "LenR1='251', LenR2='251', LenI1='7'" #Rapid V2
    elif Recipe =='9':
        return "LenR1='50', LenR2='50', LenI1='7'" #RNASeq
    elif Recipe =='10':
        return "LenR1='50', LenR2='50', LenI1='7', LenI2='7'" #Highly multiplexed RNASeq
    elif Recipe =='11':
        return "LenR1='151', LenR2='151', LenI1='7', LenI2='7'" #Highly multiplexed Nova S3,S4
    else:
        raise Exception, "Unhandled Recipe code:{}".format(Recipe)

def updateFC(sequenceDB,FCID,Machine,pwd):
    logger = logging.getLogger('updateFC')
    seqLoc = 'fastq18'
    DateBcl = cur_datetime()
    DateRTA = getRTAdate(pwd)
    DateRead1 = getRead1Date(FCID)
    RTAVer,HCSver = ver_num(pwd)
    Read_SQL_Command = getReads(sequenceDB,FCID)
    #print DateRTA

    sql = ("UPDATE Flowcell "
        "SET {0}, RTAVer='{1}', HCSVer='{2}', DateRead1='{3}', DateRTA='{4}', "
        "DateBcl='{5}', SeqsataLoc='{6}' "
        "WHERE FCillumID='{7}'"
        ).format(Read_SQL_Command,RTAVer,HCSver,DateRead1,DateRTA,DateBcl,seqLoc,FCID)
    if verbose == True:
        print sql

    sequenceDB.execute(sql)
    logger.info(sql)

#Updates ClusterDensity for Lane Entries
def updateLane(sequenceDB,FCID,Machine,pwd):
    logger = logging.getLogger('updateLane')
    logger.debug('Running updateLane')

    #tr command: remove everything but printable ascii characters
    tLaneKDen = getoutput("cat %s/First_Base_Report.htm | /usr/bin/w3m -T text/html | tr -cd '\11\12\15\40-\176' | sed 's,(k/mm2),,g' | grep -o Density.*" % pwd)
    LaneKDen = tLaneKDen.split('\n')
    Read1KDen = LaneKDen[0].split()
    Read2KDen = LaneKDen[1].split()

    #need to remove blank lines from sss_lanes
    sss_lanes = getoutput("cut -d, -f2 /nfs/genotyping/Sequencing_SampleSheets/*%s* | sort -u | grep -v Lane" % FCID).split('\n')
    #uses script that parses binary file for cluster density information
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadTileMetrics.pl %s/InterOp/TileMetricsOut.bin | grep 100[[:space:]] > %s/TileClusterDensity.txt' % (pwd,pwd))
    # %ClusterPF metrics are determined at a tile level and then the %Cluster PF is the average of the %Cluster PF of all tiles.
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadTileMetrics.pl %s/InterOp/TileMetricsOut.bin | grep [[:space:]]102[[:space:]] | sort -k2n > %s/TotalCluster.txt' % (pwd,pwd))
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadTileMetrics.pl %s/InterOp/TileMetricsOut.bin | grep [[:space:]]103[[:space:]] | sort -k2n > %s/ClusterPF.txt' % (pwd,pwd))
    for LaneNum in range(1,len(sss_lanes)+1):
        ClustDen,ClustDenStDev = getClustDen(LaneNum,pwd)
        ClustPF,ClustPFStDev = getClustPF(LaneNum,pwd)
        AvgKDen = (float(Read1KDen[LaneNum])+float(Read2KDen[LaneNum]))/2
        ApproxClustDen = str(1.1296*AvgKDen+237.4)
        #print ClustPF,ClustPFStDev
        sql = ("UPDATE Lane l "
            "JOIN Flowcell f on l.FCID=f.FCID "
            "SET ClustDen='{0}',ClustDenStDev='{1}',"
            "ClusterPF='{2}',ClusterPFStDev=round('{3}',3) "
            "WHERE LaneNum={4} AND f.FCillumID='{5}'" 
            ).format(ClustDen,ClustDenStDev,ClustPF,ClustPFStDev,str(LaneNum),FCID)

        if verbose == True:
            print sql
        sequenceDB.execute(sql)
        logger.info(sql)

    totalNumLanes = totalLanesCheck(sss_lanes,FCID)
    qmets(sequenceDB,pwd,totalNumLanes,Machine,FCID)
    updateLnFraction(sequenceDB,FCID,pwd)

def totalLanesCheck(sss_lanes,FCID):
    """Check if # of lanes in generated sequencing sample sheet matches actual
    number"""
    if FCID[-3] == 'N': #Normal Flowcells
        actualNumLanes = 8
    elif FCID[0] == 'H': #Rapid Runs
        actualNumLanes = 2
    else:
        raise Exception, 'Unhandled FCID for flowcell %s' % FCID

    if len(sss_lanes) != actualNumLanes:
        raise Exception, 'Number of lanes in SSS is incorrect!'
    else:
        return len(sss_lanes)

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
    #print info
    #print info
    for samp in info:
        LnFraction = getoutput("grep ,%s, %s/*%s*.csv | grep ,%s, | cut -d, -f6 | cut -d_ -f1" % (samp[1],pwd,samp[4],samp[3]))
        #print LnFraction
        sql = ("UPDATE Lane l "
            "SET LnFraction={0} "
            "WHERE FCID={1} and Lanenum={2} and DBID={3}" 
            ).format(LnFraction,samp[2],samp[3],samp[0])
        logger.info(sql)

        if verbose == True:
            print "grep ,%s, %s/*%s*.csv | grep ,%s, | cut -d, -f6 | cut -d_ -f1" % (samp[1],pwd,samp[4],samp[3])
            print sql
        sequenceDB.execute(sql)

#gets average of cluster density values from all tiles and also calcs standard deviation 
def getClustPF(LaneNum,pwd):
	info = getoutput("paste %s/ClusterPF.txt %s/TotalCluster.txt | grep ^%s | awk '{sum +=$4/$8; sumsq+=($4/$8)**2} END { avg = sum/NR*100; stdev=sqrt(sumsq/NR - (sum/NR)**2)*100; print avg\"\t\"stdev}'" % (pwd,pwd,LaneNum))
	ClustPF = info.split()[0]
	ClustPFStDev = info.split()[1]
	return ClustPF,ClustPFStDev

#gets average of cluster density values from all tiles and also calcs standard deviation 
def getClustDen(LaneNum,pwd):
	info = getoutput("egrep ^%s %s/TileClusterDensity.txt | cut -f4 | awk '{SUM+=$1; SUMSQ+=$1*$1} END {avg = SUM/NR/1000; stdev =sqrt(SUMSQ/NR - (SUM/NR)**2)/1000; print avg\"\t\"stdev}'" % (LaneNum,pwd))
	ClustDen = info.split()[0]
	ClustDenStDev = info.split()[1]
	return ClustDen,ClustDenStDev


def qmets(sequenceDB,run_folder,total_lanes,machine,FCID):
    logger = logging.getLogger('qmets')
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadQualityMetrics.pl %s/InterOp/QMetricsOut.bin > %s/QMetrics.txt' % (run_folder,run_folder))
    logger.info('Running loadQualityMetrics.pl')
    if os.path.isfile('%s/QMetrics.txt' % run_folder) == False:
        logger.info('No QMetrics file found!  Skipping Q30 info insertion')
        return 0
    QMetrics = open('%s/QMetrics.txt' % (run_folder),'r')
    met_dictR1 = {}
    met_dictI1 = {}
    met_dictR2 = {}

    for lane in range(1,total_lanes+1):
        met_dictR1[lane] = [0]*50
        met_dictI1[lane] = [0]*50
        met_dictR2[lane] = [0]*50
    for line in QMetrics.readlines()[1:]:
        info = line.split()
        lane = info[0]
        cycle = int(info[2])
        qscores = [int(y) for y in info[3:]]
        #print qscores
        #V4 Flowcell
        if FCID[-3] == 'N':
            if cycle > 134:
                met_dictR2[int(lane)] = [sum(x) for x in zip(met_dictR2[int(lane)],qscores)]
            elif cycle > 126:
                met_dictI1[int(lane)] = [sum(x) for x in zip(met_dictI1[int(lane)],qscores)]
            elif cycle > 0:
                met_dictR1[int(lane)] = [sum(x) for x in zip(met_dictR1[int(lane)],qscores)]

        else:
            if cycle > 108:
                met_dictR2[int(lane)] = [sum(x) for x in zip(met_dictR2[int(lane)],qscores)]
            elif cycle > 101:
                met_dictI1[int(lane)] = [sum(x) for x in zip(met_dictI1[int(lane)],qscores)]
            elif cycle > 0:
                met_dictR1[int(lane)] = [sum(x) for x in zip(met_dictR1[int(lane)],qscores)]
        #print met_dictR1
    for lane in range(1,total_lanes+1):
        totalR1 = 0
        qsumR1 = 0
        q30R1 = 0
        totalI1 = 0
        qsumI1 = 0
        q30I1 = 0
        totalR2 = 0
        qsumR2 = 0
        q30R2 = 0
        #print met_dictR1
        for pos in range(0,50):
            totalR1 += met_dictR1[int(lane)][pos]
            qsumR1 += met_dictR1[int(lane)][pos]*pos
            if pos >= 29:
                q30R1 += met_dictR1[int(lane)][pos]
                
            totalI1 += met_dictI1[int(lane)][pos]
            qsumI1 += met_dictI1[int(lane)][pos]*pos
            if pos >= 29:
                q30I1 += met_dictI1[int(lane)][pos]
                                
            totalR2 += met_dictR2[int(lane)][pos]
            qsumR2 += met_dictR2[int(lane)][pos]*pos
            if pos >= 29:
                q30R2 += met_dictR2[int(lane)][pos]
        if totalR1 !=0 or totalI1 != 0 or totalR2 != 0:
            #print totalR1,lane,met_dictR1
            avgQ30R1 = float(qsumR1)/totalR1
            perQ30R1 = float(q30R1)/totalR1
            avgQ30I1 = float(qsumI1)/totalI1
            perQ30I1 = float(q30I1)/totalI1
            avgQ30R2 = float(qsumR2)/totalR2
            perQ30R2 = float(q30R2)/totalR2
        else:
            avgQ30R1 = avgQ30I1 = avgQ30I2 = perQ30R1 = perQ30R2 = perQ30I1 = 0
        sql = ("UPDATE Lane l "
            "join Flowcell f on l.FCID=f.FCID "
            "SET mnQscR1='{0}',mnQscI1='{1}',mnQscR2='{2}',"
            "perQ30R1='{3}',perQ30I1='{4}',perQ30R2='{5}' "
            "WHERE LaneNum='{6}' AND f.FCillumID='{7}'"
            ).format(avgQ30R1,avgQ30I1,avgQ30R2,perQ30R1,perQ30I1,perQ30R2,lane,FCID)

        if verbose == True:

            print(sql)
        sequenceDB.execute(sql)
        logger.info(sql)

def usage():
	print '-h, --help\t\tShows this help message and exit'
	print '-v, --verbose\t\tPrints out MYSQL injection commands'
	sys.exit(2)


def opts(argv):
	global verbose
	verbose = False
        global runPath
	runPath = ''

	try:
		opts,args = getopt.getopt(argv, "hvi:", ['help','verbose','input='])
	except getopt.GetoptError, err:
		print str(err)
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
    #print info
    FCID = info[3]
    Machine = info[1]

    setup_logging(Machine,FCID,'fastq18')
    print 'BCL MySQL updates started'
    logger = logging.getLogger('main')
    logger.info('BCL MySQL updates started')
    logger.debug('Initializing Parameters: pwd:%s, FCID:%s, Machine:%s, seqsata_drive:%s', (pwd,FCID,Machine,'fastq18'))
    try:
        updateFC(sequenceDB,FCID,Machine,pwd)
        updateLane(sequenceDB,FCID,Machine,pwd)

        #mysql updates happen now on BCL
        logger.info('BCL MySQL updates completed')
        print 'BCL MySQL updates completed'
        sequenceDB.execute('COMMIT;')
        sequenceDB.close()
    except:
        traceback.print_exc()
        sequenceDB.execute('ROLLBACK;')
        sequenceDB.close()
        logger.info('BCL MySQL updates failure')
        print 'BCL MySQL update failure'


main()
