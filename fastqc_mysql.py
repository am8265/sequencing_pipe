#!/bin/python2.7/python
#fastqc_mysql.py
#updates Lane table of sequenceDB with information from fastqc
#Run from run folder

import getopt
import glob
import logging
import os
import sys
from CHGV_mysql import getGAFdb
from CHGV_mysql import getSequenceDB
from CHGV_mysql import setup_logging
from commands import getoutput

def fastqc_mysql(sequenceDB,samp,seqsata,read,FCID):
    logger = logging.getLogger('fastqc_mysql')
    SampleID = samp[0]
    LaneNum = samp[1]
    Lane = samp[1]
    r = 'R'+str(read)

    sequenceDB.execute("SELECT s.Seqtype from Lane l join SampleT s on l.dbid=s.dbid join Flowcell f on f.fcid=l.fcid where s.CHGVID=%s and f.FCillumID=%s and LaneNum=%s", (SampleID,FCID,LaneNum))
    Seqtype = sequenceDB.fetchone()[0]

    summaryFile = '%s/%s/%s/%s/%s_*L00%s_R%s_fastqc/summary.txt'  % (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read)
    fastqcDataFile = '%s/%s/%s/%s/%s_*L00%s_R%s_fastqc/fastqc_data.txt' % (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read)

    #print fastqcDataFile
    PCRDup = getoutput('grep "Total Duplicate Percentage" %s | cut -f2' % fastqcDataFile)
    BasicStat = getoutput('grep "Basic Statistics" %s | cut -f1' % summaryFile)
    PerBaseQual = getoutput('grep "Per base sequence quality" %s | cut -f1' % summaryFile)
    PerSeqQual = getoutput('grep "Per sequence quality scores" %s | cut -f1' % summaryFile)
    PerBaseContent = getoutput('grep "Per base sequence content" %s | cut -f1' % summaryFile)
    PerBaseGCContent = getoutput('grep "Per base GC content" %s | cut -f1' % summaryFile)
    PerSeqGCContent = getoutput('grep "Per sequence GC content" %s | cut -f1' % summaryFile)
    PerBaseNContent = getoutput('grep "Per base N content" %s | cut -f1' % summaryFile)
    SeqLenDist = getoutput('grep "Sequence Length Distribution" %s | cut -f1' % summaryFile)
    OverRepSeq = getoutput('grep "Overrepresented sequences" %s | cut -f1' % summaryFile)
    KmerContent = getoutput('grep "Kmer Content" %s | cut -f1' % summaryFile)

    #print '%s/%s/%s/%s_*L00%s_R%s_fastqc/fastqc_data.txt' % (seqsata,SampleID,FCID,SampleID,Lane,read)	
    if glob.glob('%s/%s/%s/%s/%s_*L00%s_R%s*_fastqc/fastqc_data.txt' %
            (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read)) != []:
        if r == 'R1':
            sql= ("UPDATE Lane l "
                "join prepT pt on l.prepID=pt.prepID "
                "join Flowcell f on l.FCID=f.FCID "
                "SET l.BasicStatR1='{0}',l.PerBaseQualR1='{1}',"
                "l.PerSeqQualR1='{2}',l.PerBaseContentR1='{3}',"
                "l.PerBaseGCContentR1='{4}',l.PerSeqGCContentR1='{5}',"
                "l.SeqDupR1='{6}',l.OverRepSeqR1='{7}',l.KmerContentR1='{8}' "
                "where pt.chgvid='{9}' and f.FCillumID='{10}' and "
                "l.LaneNum='{11}'"
                ).format(BasicStat,PerBaseQual,PerSeqQual,PerBaseContent,PerBaseGCContent,PerSeqGCContent,PCRDup,OverRepSeq,KmerContent,SampleID,FCID,Lane)

            if verbose == True:
                print sql
            sequenceDB.execute(sql)
            logger.info(sql)
        elif r == 'R2':
            sql = ("UPDATE Lane l "
                "join prepT pt on l.prepID=pt.prepID "
                "join Flowcell f on l.FCID=f.FCID "
                "SET l.BasicStatR2='{0}',l.PerBaseQualR2='{1}',"
                "l.PerSeqQualR2='{2}',l.PerBaseContentR2='{3}',"
                "l.PerBaseGCContentR2='{4}',l.PerSeqGCContentR2='{5}',"
                "l.SeqDupR2='{6}',l.OverRepSeqR2='{7}',l.KmerContentR2='{8}' "
                "where pt.chgvid='{9}' and f.FCillumID='{10}' and "
                "l.LaneNum='{11}'"
                ).format(BasicStat,PerBaseQual,PerSeqQual,PerBaseContent,PerBaseGCContent,PerSeqGCContent,PCRDup,OverRepSeq,KmerContent,SampleID,FCID,Lane)

            if verbose == True:
                print sql
            sequenceDB.execute(sql)
            logger.info(sql)

        else:
            raise Exception, "Could not find fastqc_data in %s!" % seqsata


    else:
        print os.path.isfile('%s/%s/%s/%s/%s_*L00%s_R%s*_fastqc/fastqc_data.txt' % (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read))
        print '%s/%s/%s/%s/%s_*L00%s_R%s*_fastqc/fastqc_data.txt' % (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read)
        print seqsata,SampleID,FCID,SampleID,Lane,read,Seqtype
        raise Exception
    #print SampleID,Seqtype
    #alerts for fastqc checks
    #print '%s/%s/%s/%s/%s_*L00%s_R%s_fastqc/fastqc_data.txt' % (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read)
    if glob.glob('%s/%s/%s/%s/%s_*L00%s_R%s_fastqc/fastqc_data.txt' % (seqsata,Seqtype.upper(),SampleID,FCID,SampleID,Lane,read)) != []:
        if PerBaseQual == 'FAIL':
            report(PerBaseQual,r,LaneNum,SampleID,FCID,Seqtype,'Per Base Quality Fail')
        if PerBaseNContent == 'FAIL':
            report(PerBaseNContent,r,LaneNum,SampleID,FCID,Seqtype,'Per Base N Content Fail')
        if OverRepSeq == 'FAIL':
            report(OverRepSeq,r,LaneNum,SampleID,FCID,Seqtype,'Over Represented Sequences Fail')
        if float(PCRDup) > 25 and Seqtype == 'Exome':
            report(PCRDup,r,LaneNum,SampleID,FCID,Seqtype,'Percent PCR duplicates Fail')
        if float(PCRDup) > 7.5 and Seqtype == 'Genome':
            report(PCRDup,r,LaneNum,SampleID,FCID,Seqtype,'Percent PCR duplicates Fail')
        if float(PCRDup) > 60 and Seqtype == 'RNAseq':
            report(PCRDup,r,LaneNum,SampleID,FCID,Seqtype,'Percent PCR duplicates Fail')

def report(Value,Read,LaneNum,SampleID,FCID,Seqtype,failure):
	os.system('echo "%s: %s %s %s Ln%s %s %s" >> fastqc_failures.txt ' % (failure,SampleID,Seqtype,FCID,LaneNum,Read,Value))


def fastqc(sequenceDB,seqsata,FCID):
	sequenceDB.execute("SELECT s.CHGVID,l.LaneNum from Lane l join SampleT s on l.dbid=s.dbid join Flowcell f on f.fcid=l.fcid where f.FCillumID=%s and (l.FailR1 IS NULL and l.FailR2 IS NULL)", FCID) 
	Samples = sequenceDB.fetchall()
	os.system('rm fastqc_failures.txt')
	#iterate over every sample, read 1 and 2 in a FCID in every lane.  
	for samp in Samples:
		for read in range(1,3):
			fastqc_mysql(sequenceDB,samp,seqsata,read,FCID)


def xenlims_cp(FCID,seqsata):
	logger = logging.getLogger('xenlims_cp')
	#print 'ssh apache@10.73.50.38 mkdir /home/dev/public_html/FastQCfiles/%s' % (FCID)
	#print 'scp %s/*/%s/*fastqc apache@10.73.50.38:/home/dev/public_html/FastQCfiles/%s' % (seqsata,FCID,FCID)
	os.system('ssh apache@10.73.50.38 mkdir /home/dev/public_html/FastQCfiles/%s' % (FCID))
	os.system('scp -r %s/*/*/%s/*fastqc apache@10.73.50.38:/home/dev/public_html/FastQCfiles/%s' % (seqsata,FCID,FCID)) 
	logger.info('scp -r %s/*/*/%s/*fastqc apache@10.73.50.38:/home/dev/public_html/FastQCfiles/%s' % (seqsata,FCID,FCID))

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
    pwd = os.getcwd()
    info = pwd.split('_')
    Machine = info[1]
    FCID = info[3]
    sequenceDB = getSequenceDB()
    email = 'jb3816@.columbia.edu'
    opts(sys.argv[1:])
    
    global sata_loc
    if 'fastq' not in sata_loc:
        seqsata_drive = sata_loc.split('/')[2]
    else:
        seqsata_drive = 'fastq15'
        sata_loc = '/nfs/fastq15'

    
    seqsata_drive = sata_loc.split('/')[2]

    setup_logging(Machine,FCID,seqsata_drive)
    logger = logging.getLogger('main')
    logger.info("Starting Mysql injection of Fastqc results...")
    fastqc(sequenceDB,sata_loc,FCID)
    #xenlims_cp(FCID,sata_loc)

    if os.path.isfile('fastqc_failures.txt'):
        os.system('cat fastqc_failures.txt | mail -s "FASTQC failure `pwd`" %s' % email)

    sequenceDB.execute('COMMIT;')
    sequenceDB.close()
    logger.info("Fastqc result upload completed")
main()
