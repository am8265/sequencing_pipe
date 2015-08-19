# bcl.py
# Joshua Bridgers
# 04/26/2012
# jb371@dm.duke.edu
#
# Submits finished sequenced run for BCL to fastq conversion

import csv
import glob
import logging
import os
import re
import getopt
import sys
from commands import getoutput
from CHGV_mysql import getGAFdb
from CHGV_mysql import setup_logging
from CHGV_mysql import getSequenceDB
from CHGV_mysql import getUserID
from sss_check import sss_qc

def bcl(info,sata_loc,seqsata,machine,sequenceDB):
    logger = logging.getLogger('bcl')
    in_dir = 'Data/Intensities/BaseCalls/'
    script_dir = '/home/jb3816/github/sequencing_pipe'
    pwd = os.getcwd()

    #Name of folder contains flow cell info, machine info etc.  Parsing
    FCID = info[3]
    HiSeq = info[1]
    Date = info[0][3:7]
    Date_Long = info[0]
    Script = HiSeq+'_'+Date+'_BCL.sh'

    dir_check(sata_loc,FCID)
    out_dir = sata_loc + '/' + Date_Long + '_' + HiSeq + '_' + FCID + '_Unaligned'
    base_script =  '/nfs/goldstein/software/bcl2fastq_v1.8.4/bin/configureBclToFastq.pl --input-dir %s --output-dir %s --mismatches 1 ' % (in_dir,out_dir)

    #adds tiles parameter to BCL script if specified
    if sampleSheet:
        base_script += '--sample-sheet '+sampleSheet +' '
    else:
        base_script += '--sample-sheet *%s*.csv ' % FCID

    if tiles != False:
        GAFdb = getGAFdb()

        print "UPDATE Flowcell SET TileCommand=CONCAT_WS(';',TileCommand,'%s') WHERE FCID='%s'" % (tiles,FCID)
        GAFdb.execute("UPDATE Flowcell SET TileCommand=CONCAT_WS(';',TileCommand,%s) WHERE FCID=%s", (tiles,FCID))
        print "UPDATE Flowcell SET TileCommand=CONCAT_WS(';',TileCommand,'%s') WHERE FCillumID='%s'" % (tiles,FCID)
        sequenceDB.execute("UPDATE Flowcell SET TileCommand=CONCAT_WS(';',TileCommand,%s) WHERE FCillumID=%s", (tiles,FCID))

        GAFdb.execute('COMMIT;')
        sequenceDB.execute('COMMIT;')
        GAFdb.close()
        sequenceDB.close()

        base_script = base_script + '--tiles='+tiles+' '
        os.system('echo %s > tiles.txt' % tiles)

    if base_mask != False:
        base_script = base_script + '--use-bases-mask '+base_mask+' '
    print base_script
    #raw_input(base_script + ' 2>> %s/summary/GAF_PIPELINE_LOGS/%s_%s_%s.log' % (sata_loc,machine,FCID,seqsata))
    logger.info(base_script)
    os.system(base_script + ' 2>> %s/summary/GAF_PIPELINE_LOGS/%s_%s_%s.log' % (sata_loc,machine,FCID,seqsata))
    #Submit bcl job to the cluster
    os.system('cp '+script_dir+'/run_C1.8.sh '+out_dir+'/'+Script)
    #print 'echo "cd %s ; qsub -N %s_bcl -pe make 32 %s/%s" | ssh solexa3.lsrc.duke.edu ' % (out_dir,FCID,out_dir,Script)
    logger.info('cd %s ; qsub -cwd -v PATH -N %s_%s_%s_bcl %s/%s' % (out_dir,machine,FCID,seqsata,out_dir,Script))
    os.system('cd %s ; qsub -cwd -v PATH -N %s_%s_%s_bcl %s/%s' % (out_dir,machine,FCID,seqsata,out_dir,Script))


#checks if a bcl directory already exists
def dir_check(sata_loc,FCID):
    logger = logging.getLogger('dir_check')
    dir_path = glob.glob('/nfs/%s/*%s*Unaligned' % (sata_loc,FCID))
    if dir_path != []:
        logger.warn('BCL directory already exists! %s' % dir_path)
        raise Exception, 'BCL directory already exists! %s' % dir_path

def storage(info,sata_loc,seqsata,machine,sequenceDB):
    if sata_loc == '':
        raise exception, 'Storage Location not specified!'
    checkSataLoc(sata_loc)
    bcl(info,sata_loc,seqsata,machine,sequenceDB)


def checkSataLoc(sata_loc):
	logger = logging.getLogger('checkSataLoc')
	if os.path.isdir(sata_loc) == False:
		logger.warn('Path for seqsata drive, %s is incorrect!' % sata_loc)
		raise Exception, 'Path for seqsata drive, %s is incorrect!' % sata_loc


def usage():
    print '-h, --help\t\tShows this help message and exit'
    print '-o, --output\t\tPath to output folder'
    print '-s, --sampleSheet\tAllows user-specified sample sheet'
    print '--noSSS\t\t\tForbids creation of a sample sheet'
    print '--tiles\t\t\tSpecifies what tiles you want BCL performed on.  Ex: s_[12]_1[123]0[1-7]'
    print '--use-bases-mask\tSpecifies what bases you want to mask.  Ex: y100n,I6n,y50n'
    sys.exit(2)

def check_sss(FCID):
    if sampleSheet:
        sample_sheets = [sampleSheet]
    else:
	sample_sheets = glob.glob('/nfs/genotyping/Sequencing_SampleSheets/*%s*.csv' % FCID)
	if len(sample_sheets) > 1:
		raise Exception, 'Multiple Sample Sheets exist with the same FCID'
	else:
		sss_qc(FCID)


def getSSSLaneFraction(DBID,FCID,LaneNum,sequenceDB):
    sequenceDB.execute("SELECT SeqType FROM Lane l  JOIN SeqType st ON st.prepID=l.prepID JOIN Flowcell f on l.FCID=f.FCID WHERE l.DBID=%s AND LaneNum=%s AND FCillumID=%s", (DBID,LaneNum,FCID))
    seqtype = sequenceDB.fetchone()[0]
    sequenceDB.execute("SELECT laneFraction FROM Lane l JOIN Flowcell f ON f.FCID=l.FCID JOIN samplesTOrun s2r ON l.seqID = s2r.seqID WHERE l.DBID=%s AND FCIllumID=%s AND laneNum=%s", (DBID,FCID,LaneNum))
    laneFraction = sequenceDB.fetchone()[0]

    sequenceDB.execute("SELECT COUNT(DISTINCT(PoolID)) FROM Lane l JOIN Flowcell f ON f.FCID=l.FCID WHERE FCIllumID=%s AND laneNum=%s AND poolID!=0", (FCID,LaneNum))
    numPools = sequenceDB.fetchone()[0]

    sequenceDB.execute("SELECT COUNT(PoolID) FROM Lane l JOIN Flowcell f ON f.FCID=l.FCID WHERE FCIllumID=%s AND laneNum=%s AND poolID=0", (FCID,LaneNum))
    NumOtherSamples = sequenceDB.fetchone()[0]

    sequenceDB.execute("SELECT COUNT(*) FROM Lane l join Flowcell f ON l.FCID=f.FCID WHERE l.prepID=(SELECT prepID FROM Lane l JOIN Flowcell f ON l.FCID=f.FCID WHERE l.LaneNum=%s AND f.FCillumID=%s AND DBID=%s) AND FCillumID=%s", (LaneNum,FCID,DBID,FCID))
    NumLanesSampleOn = sequenceDB.fetchone()[0] #The number of times a sample is on a flowcell (for genomes and RNASeq)

    """
    sequenceDB.execute("SELECT poolID FROM Lane l JOIN Flowcell f on l.fcid=f.fcid WHERE LaneNum=%s AND FCIllumID=%s AND DBID=%s", (LaneNum,FCID,DBID))
    poolID = sequenceDB.fetchone()[0]
    sequenceDB.execute("select count(distinct LaneNum) from Lane l JOIN Flowcell f on l.fcid=f.fcid where FCillumID=%s and poolID=%s", (FCID,poolID))
    NumLanesExomeSampleOn = sequenceDB.fetchone()[0] #The number of times a sample is on a flowcell for Exomes and Custom Capture samples
    """
    #            JOIN poolMembers pM ON l.dbid=pM.dbid \

    #print DBID,FCID,LaneNum,DBID,FCID,LaneNum
    if numPools > 0:
        sequenceDB.execute("SELECT \
            (CASE \
                WHEN l.poolID=0 THEN 1 \
                WHEN l.poolID!=0 THEN \
                    (SELECT COUNT(DISTINCT pM.prepID) FROM poolMembers pM WHERE pM.poolid=\
                        (SELECT DISTINCT l.poolID FROM Lane l JOIN Flowcell f ON f.FCID=l.FCID WHERE l.DBID=%s AND FCillumID=%s AND LaneNum=%s)) \
            END) \
            FROM Lane l \
                JOIN Flowcell f ON f.FCID=l.FCID \
            WHERE \
                l.DBID=%s AND \
                FCIllumid=%s AND \
                LaneNum=%s", \
                (DBID,FCID,LaneNum,DBID,FCID,LaneNum))
        NumPoolSamples = sequenceDB.fetchone()[0]
    else:
        numPools = 0
        NumPoolSamples = 1
    """
    if LaneNum==3:
        print DBID,seqtype,LaneNum,float(laneFraction),numPools,NumPoolSamples,NumOtherSamples,NumLanesSampleOn,seqtype
        print laneFraction,NumPoolSamples,NumOtherSamples+1
    #SampleLaneFraction = 1.0/numPools/NumPoolSamples
    """

    if seqtype == 'Genome':
        SampleLaneFraction = float(laneFraction)/(numPools+1)/NumPoolSamples/NumLanesSampleOn
        #SampleLaneFraction = float(laneFraction)
    elif seqtype == 'RNAseq':
        SampleLaneFraction = float(laneFraction)/(numPools+1)/NumPoolSamples/NumLanesSampleOn
    elif seqtype == 'Exome':
        SampleLaneFraction = float(laneFraction)/NumLanesSampleOn/NumPoolSamples/(NumOtherSamples+1)
    elif seqtype == 'Custom Capture':
        SampleLaneFraction = float(laneFraction)/NumPoolSamples/(NumOtherSamples+1+numPools-1)
        #SampleLaneFraction = float(laneFraction)/NumPoolSamples

    return SampleLaneFraction

def create_sss(FCID,Machine,date,sequenceDB):
    logger = logging.getLogger('create_sss')
    sample_sheet = []
    outfile=open('%s_%s_%s.csv' % (Machine,date,FCID),'w')
    outfile.write("FCID,Lane,SampleID,SampleReference,Index,Description,Control,Recipe,Operator,Project\n")
    lane_num = 8
    if FCID[0] == 'H':
        lane_num = 2
    for LaneNum in range(1,lane_num+1):
        sequenceDB.execute("SELECT DBID FROM Lane l join Flowcell f ON l.FCID=f.FCID where FCillumID=%s AND f.Machine=%s AND l.laneNum=%s ORDER BY LaneNum", (FCID,Machine,LaneNum))
        ss_samples = sequenceDB.fetchall()
        if len(ss_samples) == 0:
            print FCID,Machine,LaneNum
            raise Exception, "No Samples were found with FCillumID: %s" % FCID
        #print ss_samples
        for DBID in ss_samples:
            DBID = str(DBID[0])
            laneFraction = getSSSLaneFraction(DBID,FCID,LaneNum,sequenceDB)
            logger.info("SELECT f.FCillumID FCID,l.LaneNum Lane,pt.CHGVID SampleID,'human' SampleReference,(case when f.recipe=2 then '' when s.SeqType='Exome' then pm.adapterlet when s.Seqtype='RNAseq' then s2r.adapterlet when s.Seqtype='Genome' then s2r.adapterlet when s.SeqType='Custom Capture' then pm.adapterlet END) 'Index',CONCAT('%s','_',round(s2r.picomoles,1),'pM) Description,'N' Control,(case when f.recipe=1 then '101x7x101' when f.recipe=2 then '101x101' when f.recipe=3 then '101x7x7x101' END) Recipe,replace(u.name,' ','') Operator,replace(s.GAFbin,' ','') Project FROM Lane l join Flowcell f ON f.FCID=l.FCID join prepT pt ON l.prepID=pt.prepID join samplesTOrun s2r ON s2r.seqID=l.seqID join users u ON u.userID=f.userID join SampleT s ON s.DBID=pt.DBID LEFT JOIN poolMembers pm ON (pm.DBID=pt.DBID AND pm.poolID=l.poolID) where l.dbid='%s' AND f.FCillumID='%s' AND LaneNum='%s'" % (laneFraction,DBID,FCID,LaneNum))

            sequenceDB.execute("SELECT f.FCillumID FCID,\
                l.LaneNum Lane,\
                pt.CHGVID SampleID,\
                'human' SampleReference,\
                (case \
                    when f.recipe=2 then '' \
                    when st.SeqType='Exome' then pm.adapterlet \
                    when st.Seqtype='RNAseq' then s2r.adapterlet \
                    when st.Seqtype='Genome' then s2r.adapterlet \
                    when st.SeqType='Custom Capture' then pm.adapterlet \
                END) 'Index',\
                CONCAT(round(%s,4),'_',round(s2r.picomoles,1),'pM') Description, \
                'N' Control,\
                (case \
                    when f.recipe=1 then '101x7x101'\
                    when f.recipe=2 then '101x101'\
                    when f.recipe=3 then '101x7x7x101'\
                    when f.recipe=4 then '100x9x100'\
                    when f.recipe=5 then '126x7x126'\
                    when f.recipe=6 then '101x8x101'\
                    when f.recipe=7 then '151x7x151'\
                END) Recipe,\
                replace(u.name,' ','') Operator,\
                replace(s.GAFbin,' ','') Project \
                FROM Lane l \
                    join Flowcell f ON f.FCID=l.FCID \
                    join prepT pt ON l.prepID=pt.prepID \
                    join samplesTOrun s2r ON s2r.seqID=l.seqID \
                    join users u ON u.userID=f.userID \
                    join SeqType st ON l.prepID=st.prepID \
                    join SampleT s ON s.DBID=pt.DBID \
                    LEFT JOIN poolMembers pm ON \
                        (pm.DBID=pt.DBID AND pm.poolID=l.poolID) \
                where \
                    l.dbid=%s AND \
                    f.FCillumID=%s AND \
                    LaneNum=%s", \
            (laneFraction,DBID,FCID,LaneNum))
            ss_line = sequenceDB.fetchone()

            sample_sheet.append(ss_line)
            #print ss_line,DBID,FCID,LaneNum
            outfile.write(",".join(map(str,ss_line))+'\n')
    outfile.close()
    #os.chmod('%s_%s_%s.csv' % (Machine,date,FCID),0775) #changes permissions of *csv file
    #copies sequencing sample sheet to genotyping location
    os.system('cp %s_%s_%s.csv /nfs/genotyping/Sequencing_SampleSheets/' % (Machine,date,FCID))
    logger.info('cp %s_%s_%s.csv /nfs/genotyping/Sequencing_SampleSheets/' % (Machine,date,FCID))

def updateSamples(sequenceDB,FCID):
    logger = logging.getLogger('updateSamples')
    userID = getUserID()
    if verbose == True:
        print "INSERT INTO statusT (CHGVID,status_time,status,DBID,prepID,userID) SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'BCL',pt.DBID,pt.prepID,'%s' FROM Flowcell f join Lane l ON l.FCID=f.FCID join prepT pt ON pt.prepID=l.prepID where FCillumid='%s'" % (userID,FCID)
    sequenceDB.execute("INSERT INTO statusT (CHGVID,status_time,status,DBID,prepID,userID) SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'BCL',pt.DBID,pt.prepID,%s FROM Flowcell f join Lane l ON l.FCID=f.FCID join prepT pt ON pt.prepID=l.prepID where FCillumid=%s", (userID,FCID))
    logger.info("INSERT INTO statusT (CHGVID,status_time,status,DBID,prepID,userID) SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'BCL',pt.DBID,pt.prepID,'%s' FROM Flowcell f join Lane l ON l.FCID=f.FCID join prepT pt ON pt.prepID=l.prepID where FCillumid='%s'" % (userID,FCID))


def opts(argv):
    global tiles
    tiles = False
    global base_mask
    base_mask = False
    global sata_loc
    sata_loc = ''
    global verbose
    verbose = False
    global sampleSheet
    sampleSheet = ''
    global noSSS
    noSSS = False

    try:
        opts,args = getopt.getopt(argv, "bhno:s:v", ['help','output=','tiles=','use-bases-mask=','verbose','sampleSheet=','noSSS'])
    except getopt.GetoptError, err:
        usage()
    for o,a in opts:
        if o in ('-h','--help'):
            usage()
        elif o in ('--tiles'):
            tiles = a
        elif o in ('--use-bases-mask'):
            base_mask = a
        elif o in ('-v','--verbose'):
            verbose = True
        elif o in ('-s','--sampleSheet'):
            sampleSheet = a
        elif o in ('-n','--noSSS'):
            noSSS = True
        elif o in ('-o','--output'):
            if 'fastq' not in sata_loc:
                sata_loc = a.rstrip('/')
            else:
                raise Exception, 'Output location %s is not in the whole_genome folder within a seqsata drive!' % a
        else:
            assert False, "Unhandled argument present"

def RTA_check():
	logger = logging.getLogger('RTA_check')
	if os.path.isfile('RTAComplete.txt') == False:
		logger.warn("RTA has not completed!")
		raise Exception, "RTA has not completed!"
	else:
		logger.info('RTA has already completed')
		print "RTA has already completed"

def main():
    pwd = os.getcwd()
    opts(sys.argv[1:])
    info = pwd.split('/')[4].split('_')
    Date = info[0]
    FCID = info[3]
    Machine = info[1]
    seqsata_drive = sata_loc.split('/')[2]

    setup_logging(Machine,FCID,seqsata_drive)
    logger = logging.getLogger('main')
    logger.debug('Initializing Parameters: pwd:%s, FCID:%s, Machine:%s, seqsata_drive:%s, tiles:%s, base_mask:%s', (pwd,FCID,Machine,seqsata_drive,tiles,base_mask))
    print "Starting BCL job creation..."
    logger.info("Starting BCL job creation...")

    RTA_check()
    sequenceDB = getSequenceDB()
    if noSSS == False:
        create_sss(FCID,Machine,Date,sequenceDB)
    #check_sss(FCID)
    storage(info,sata_loc,seqsata_drive,Machine,sequenceDB)
    updateSamples(sequenceDB,FCID)

    logger.info("BCL successfully started")
    print "BCL successfully started"
    sequenceDB.execute('COMMIT;')
    sequenceDB.close()

if __name__ == '__main__':
	main()
