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

def getRTAdate():
	DateRTA = getoutput('ls -l --full-time RTAComplete.txt').split(' ')
	date = DateRTA[5]
	time = DateRTA[6].split('.')[0]
	return '%s %s' % (date,time)

def getRead1Date(FCID):
	Read1Date = getoutput('ls -l /nfs/seqscratch1/Runs/*'+FCID+'*/RunInfo.xml --time-style=full-iso 2>/dev/null | cut -d\  -f6,7 | sort -rk6 | head -1 | cut -d. -f1')
	return Read1Date

"""
def getChemVer(Machine):
	LaneCount = getoutput('grep FlowcellLayout RunInfo.xml | cut -d= -f2').split('"')[1]

	if LaneCount == "2" and Machine[0:2] == "H8":
		return 'v3r'
	elif LaneCount == '8':
		return 'v3'
	else:
		raise Exception, "LaneCount incorrect!"
"""
def ver_num():
	tRTA = getoutput('grep RTA ImageAnalysis_Netcopy_complete.txt')
	tRTA2 = tRTA.split(' ')
	RTAVer = tRTA2[2].strip()
	
	HCSver = getoutput('grep ApplicationVersion runParameters.xml -m1 | cut -d\> -f2 | cut -d\< -f1')
	return RTAVer,HCSver

def getReads(sequenceDB,FCID):
    sequenceDB.execute("SELECT recipe from Flowcell where FCillumID=%s", FCID)
    Recipe = sequenceDB.fetchone()
    Recipe = Recipe[0]

    if Recipe == '1':
        return "LenR1='101', LenR2='101', LenI1='7',"
    elif Recipe =='2':
        return "LenR1='101', LenR2='101',"
    elif Recipe =='3':
        return "LenR1='101', LenR2='101', LenI1='7', LenI2='7',"
    elif Recipe =='4':
        return "LenR1='100', LenR2='100', LenI1='9',"
    elif Recipe =='5':
        return "LenR1='126', LenR2='126', LenI1='7',"
    elif Recipe =='6':
        return "LenR1='101', LenR2='101', LenI1='8',"
    elif Recipe =='6':
        return "LenR1='151', LenR2='151', LenI1='7',"


def updateFC(sequenceDB,sataloc,FCID,Machine):
    logger = logging.getLogger('updateFC')
    sataloc = sataloc.split('/')[2]
    DateBcl = cur_datetime()
    DateRTA = getRTAdate()
    DateRead1 = getRead1Date(FCID) 
    #ChemVer = getChemVer(Machine)
    RTAVer,HCSver = ver_num()
    Read_SQL_Command = getReads(sequenceDB,FCID)
    if verbose == True:
        print "UPDATE Flowcell SET %s RTAVer='%s', HCSVer='%s', DateRead1='%s', DateRTA='%s', DateBcl='%s', SeqsataLoc='%s' WHERE FCillumID='%s'" % (Read_SQL_Command,RTAVer,HCSver,DateRead1,DateRTA,DateBcl,sataloc,FCID)
    sequenceDB.execute("UPDATE Flowcell SET "+Read_SQL_Command+" RTAVer=%s, HCSVer=%s, DateRead1=%s, DateRTA=%s, DateBcl=%s, SeqsataLoc=%s WHERE FCillumID=%s", (RTAVer,HCSver,DateRead1,DateRTA,DateBcl,sataloc,FCID))
    logger.info("UPDATE Flowcell SET %s RTAVer='%s', HCSVer='%s', DateRead1='%s', DateRTA='%s', DateBcl='%s', SeqsataLoc='%s' WHERE FCillumID='%s'" % (Read_SQL_Command,RTAVer,HCSver,DateRead1,DateRTA,DateBcl,sataloc,FCID))

#Updates ClusterDensity for Lane Entries
def updateLane(sequenceDB,FCID,Machine,pwd):
    logger = logging.getLogger('updateLane')
    logger.debug('Running updateLane')

    #tr command: remove everything but printable ascii characters
    tLaneKDen = getoutput("cat First_Base_Report.htm | w3m -T text/html | tr -cd '\11\12\15\40-\176' | sed 's,(k/mm2),,g' | grep -o Density.*")
    LaneKDen = tLaneKDen.split('\n')
    Read1KDen = LaneKDen[0].split()
    Read2KDen = LaneKDen[1].split()

    #need to remove blank lines from sss_lanes
    sss_lanes = getoutput("cut -d, -f2 /nfs/genotyping/Sequencing_SampleSheets/*%s* | sort -u | grep -v Lane" % FCID).split('\n')
    #uses script that parses binary file for cluster density information
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadTileMetrics.pl InterOp/TileMetricsOut.bin | grep 100[[:space:]] > TileClusterDensity.txt')
    # %ClusterPF metrics are determined at a tile level and then the %Cluster PF is the average of the %Cluster PF of all tiles.
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadTileMetrics.pl InterOp/TileMetricsOut.bin | grep [[:space:]]102[[:space:]] | sort -k2n > TotalCluster.txt')
    os.system('perl /nfs/goldstein/goldsteinlab/GAF/scripts/SAV_parse/loadTileMetrics.pl InterOp/TileMetricsOut.bin | grep [[:space:]]103[[:space:]] | sort -k2n > ClusterPF.txt')
    for LaneNum in range(1,len(sss_lanes)+1):
        ClustDen,ClustDenStDev = getClustDen(LaneNum)
        ClustPF,ClustPFStDev = getClustPF(LaneNum)
        AvgKDen = (float(Read1KDen[LaneNum])+float(Read2KDen[LaneNum]))/2
        ApproxClustDen = str(1.1296*AvgKDen+237.4)
        #print ClustPF,ClustPFStDev
        if verbose == True:
            print "UPDATE Lane l join Flowcell f on l.FCID=f.FCID SET ClustDen='%s',ClustDenStDev='%s',ClusterPF='%s',ClusterPFStDev=round('%s',3) WHERE LaneNum='%s' AND f.FCillumID='%s'" % (ClustDen,ClustDenStDev,ClustPF,ClustPFStDev,str(LaneNum),FCID)
        sequenceDB.execute("UPDATE Lane l join Flowcell f on l.FCID=f.FCID SET ClustDen=%s,ClustDenStDev=%s,ClusterPF=%s,ClusterPFStDev=round(%s,3) WHERE LaneNum=%s AND f.FCillumID=%s", (ClustDen,ClustDenStDev,ClustPF,ClustPFStDev,str(LaneNum),FCID))
        logger.info("UPDATE Lane l join Flowcell f on l.FCID=f.FCID SET ClustDen='%s',ClustDenStDev='%s',ClusterPF='%s',ClusterPFStDev=round('%s',3) WHERE LaneNum='%s' AND f.FCillumID='%s'" % (ClustDen,ClustDenStDev,ClustPF,ClustPFStDev,str(LaneNum),FCID))
    qmets(sequenceDB,pwd,len(sss_lanes),Machine,FCID)
    updateLnFraction(sequenceDB,FCID)

def updateLnFraction(sequenceDB,FCID):
	logger = logging.getLogger('updateLnFraction')
	logger.debug('Running updateLnFraction')
	sequenceDB.execute("SELECT l.DBID,s.CHGVID,l.FCID,l.lanenum,FCillumID from Lane l join Flowcell f on l.FCID=f.FCID join SampleT s on s.DBID=l.DBID where FCillumID=%s", FCID)
	info = sequenceDB.fetchall()
	#print info
	for samp in info:
		LnFraction = getoutput("grep %s *%s*.csv | grep ,%s, | cut -d, -f6 | cut -d_ -f1" % (samp[1],samp[4],samp[3]))
		logger.info("UPDATE Lane l SET LnFraction='%s' WHERE FCID='%s' and Lanenum='%s' and DBID='%s'" % (LnFraction,samp[2],samp[3],samp[0]))
		if verbose == True:
			print "UPDATE Lane l SET LnFraction='%s' WHERE FCID='%s' and Lanenum='%s' and DBID='%s'" % (LnFraction,samp[2],samp[3],samp[0])
		sequenceDB.execute("UPDATE Lane l SET LnFraction=%s WHERE FCID=%s and Lanenum=%s and DBID=%s", (LnFraction,samp[2],samp[3],samp[0]))

#gets average of cluster density values from all tiles and also calcs standard deviation 
def getClustPF(LaneNum):
	info = getoutput("paste ClusterPF.txt TotalCluster.txt | grep ^%s | awk '{sum +=$4/$8; sumsq+=($4/$8)**2} END { avg = sum/NR*100; stdev=sqrt(sumsq/NR - (sum/NR)**2)*100; print avg\"\t\"stdev}'" % (LaneNum))
	ClustPF = info.split()[0]
	ClustPFStDev = info.split()[1]
	return ClustPF,ClustPFStDev


#gets average of cluster density values from all tiles and also calcs standard deviation 
def getClustDen(LaneNum):
	info = getoutput("egrep ^%s TileClusterDensity.txt | cut -f4 | awk '{SUM+=$1; SUMSQ+=$1*$1} END {avg = SUM/NR/1000; stdev =sqrt(SUMSQ/NR - (SUM/NR)**2)/1000; print avg\"\t\"stdev}'" % (LaneNum))
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
        
        if verbose == True:
            print "UPDATE Lane l join Flowcell f on l.FCID=f.FCID SET mnQscR1='%s',mnQscI1='%s',mnQscR2='%s',perQ30R1='%s',perQ30I1='%s',perQ30R2='%s' WHERE LaneNum='%s' AND f.FCillumID='%s'" % (avgQ30R1,avgQ30I1,avgQ30R2,perQ30R1,perQ30I1,perQ30R2,lane,FCID)
        sequenceDB.execute("UPDATE Lane l join Flowcell f on l.FCID=f.FCID SET mnQscR1=%s,mnQscI1=%s,mnQscR2=%s,perQ30R1=%s,perQ30I1=%s,perQ30R2=%s WHERE LaneNum=%s AND f.FCillumID=%s", (avgQ30R1,avgQ30I1,avgQ30R2,perQ30R1,perQ30I1,perQ30R2,lane,FCID))
        logger.info("UPDATE Lane l join Flowcell f on l.FCID=f.FCID SET mnQscR1='%s',mnQscI1='%s',mnQscR2='%s',perQ30R1='%s',perQ30I1='%s',perQ30R2='%s' WHERE LaneNum='%s' AND f.FCillumID='%s'" % (avgQ30R1,avgQ30I1,avgQ30R2,perQ30R1,perQ30I1,perQ30R2,lane,FCID))

def usage():
	print '-h, --help\t\tShows this help message and exit'
	print '-v, --verbose\t\tPrints out MYSQL injection commands'
	sys.exit(2)


def opts(argv):
	global verbose 
	verbose = False
        global sata_loc
	sata_loc = ''

	try:
		opts,args = getopt.getopt(argv, "hvs:", ['help','verbose','seqsata='])
	except getopt.GetoptError, err:
		print str(err)
		usage()
	for o,a in opts:
		if o in ('-v','--verbose'):
			verbose = True
		elif o in ('-h','--help'):
			usage()
		elif o in ('-s','--seqsata'):
			sata_loc = a
		else:
			assert False, "Unhandled argument present"

def main():
    #accessing mysql gaf database
    sequenceDB = getSequenceDB()

    pwd = os.getcwd()
    Info = pwd.split('/')[4].split('_')
    #print Info
    FCID = Info[3]
    Machine = Info[1]
    opts(sys.argv[1:])
    seqsata_drive = sata_loc.split('/')[2]
    setup_logging(Machine,FCID,seqsata_drive)
    print 'BCL MySQL updates started'
    logger = logging.getLogger('main')
    logger.info('BCL MySQL updates started')
    logger.debug('Initializing Parameters: pwd:%s, FCID:%s, Machine:%s, seqsata_drive:%s', (pwd,FCID,Machine,seqsata_drive))

    try:
        updateFC(sequenceDB,sata_loc,FCID,Machine)
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
