#Runs fastqc on the first last and middle fastq.gz file for each Sample, for each Lane, for each Read, for a flowcell.  Fastqc is run on the first last and middle fastq.gz files to compensate for biases due to read positions on the flowcell.  Reads from the first fastq.gz file should be from the first upper sections (artifically making the run look better) while reads fromthe last fastq.gz will be from the lower sections.
#
#joshbridgers@gmail.com
#01.05.2013

from commands import getoutput
from CHGV_mysql import *
from MySQLdb import Connect

import getopt
import itertools
import os
import re
import sys

def fastqc(script,sample,seqsata,FCID,seqtype):
    seqsata_drive = 'fastq16'
    seqsata = '/nfs/fastq16/'

    #print seqsata,sample,FCID,seqtype
    #print 'cd %s/%s/%s/%s' % (seqsata,seqtype,sample,FCID)
    os.chdir('%s/%s/%s/%s' % (seqsata,seqtype,sample,FCID))

    #print script,sample,seqsata,FCID
    Lanes = getoutput("ls *.fastq.gz | grep -o L00[0-9] | sed 's/L00//g' | sort -u").split('\n')
    script.write('cd ' + getoutput('pwd')+'\n')
    for lane in Lanes:
        for read in range(1,3):
            numFile = getoutput('ls %s*L00%s_R%s_*.fastq.gz | wc -l' % (sample,lane,read))
            first = '001'

            #print 'ls %s*L00%s_R%s_*.fastq.gz | wc -l' % (sample,lane,read)

            #gets penultimate fastq.gz in the chance that the last fastq.gz has few reads 
            last = str(int(numFile)-1).zfill(3)

            if last == first:
                last = numFile.zfill(3)

            #fastqc on first, last fastq.gz of the read 
            script.write('/nfs/goldstein/software/FastQC-0.10.1/fastqc -t 20 --casava %s*L00%s_R%s_%s.fastq.gz %s*L00%s_R%s_%s.fastq.gz\n' % (sample,lane,read,first,sample,lane,read,last))


def script_header(seqsata,FCID):

    script = open("/nfs/%s/fastqc/%s_%s_fastqc_script.sh" % (seqsata,FCID,seqsata),'w')
    script.write('#! /bin/bash\n')
    script.write('#\n')
    script.write('#$ -S /bin/bash -cwd\n')
    script.write('#$ -o /nfs/%s/fastqc/%s_%s_fastqc_complete.sge\n' % (seqsata,FCID,seqsata))
    script.write('#$ -e /nfs/%s/fastqc/%s_%s_fastqc_out.sge\n' % (seqsata,FCID,seqsata))
    script.write('#$ -V\n')
    script.write('#$ -M jb3816@cumc.columbia.edu\n')
    script.write('#$ -m bea\n')
    script.write('#\n')
    script.write('\n')

    return script

def submit(seqsata,FCID):
    #print '/opt/sge6_2u5/bin/lx24-amd64/qsub -N %s_%s_fastqc /nfs/%s/fastqc/%s_%s_fastqc_script.sh' % (FCID,seqsata,seqsata,FCID,seqsata)
    os.system('/opt/sge6_2u5/bin/lx24-amd64/qsub -N %s_%s_fastqc /nfs/%s/fastqc/%s_%s_fastqc_script.sh' % (FCID,seqsata,seqsata,FCID,seqsata))

def opts(argv):
    global seqsata
    seqsata = ''
    global runFolder
    try:
        opts,args = getopt.getopt(argv, "s:i:", ['input=','seqsata=',])
    except getopt.GetoptError, err:
                print str(err)
                usage()
    for o,a in opts:
        if o in ('-i','--input'):
            runFolder = a
        elif o in ('-s','--seqsata'):
            seqsata = a
        else:
            assert False, "Unhandled argument present"

def main():
    opts(sys.argv[1:])
    sequenceDB = getSequenceDB()
    info = runFolder.split('_')
    FCID = info[3]


    #print 'ls %s/*%s*/Project_*/Sample_*/ -d | cut -d/ -f6 | cut -d_ -f2' % (seqsata,FCID)
    Samples = getoutput('cat /nfs/seqscratch1/Runs/%s/*csv | cut -d, -f3 | sort -u | grep -v SampleID' % runFolder).split('\n')
    script = script_header(seqsata,FCID)

    sql = ("SELECT p.CHGVID from Lane l "
            "join prepT p on l.prepid=p.prepid "
            "join Flowcell f on l.fcid = f.fcid "
            "where failr1 = 1 and failr2 = 1 and FCIllumID = '{0}'"
            ).format(FCID)
    sequenceDB.execute(sql)
    failedSamples = sequenceDB.fetchall()
    failedSamples = list(itertools.chain(*failedSamples))
    for failSample in failedSamples:
        Samples.remove(failSample)

    for samp in Samples:
        #print samp
        query = ' '.join((
            "SELECT DISTINCT(REPLACE(UPPER(st.seqtype),' ','_')) FROM Lane l JOIN Flowcell",
            "f ON l.fcid=f.fcid JOIN SeqType st ON l.prepid=st.prepid JOIN",
            "prepT p ON l.prepid=p.prepid WHERE FCILLUMID='"+FCID+"' AND",
            "CHGVID='"+samp+"'"))
        #print query
        sequenceDB.execute(query)
        seqtype=sequenceDB.fetchall()

        if len(seqtype) > 1:
            print seqtype
            raise excecption, 'Too many seqtypes returned for sample %s' % samp
        #print seqtype,samp,Samples

        fastqc(script,samp,seqsata,FCID,seqtype[0][0])
    script.write('cd /nfs/seqscratch1/Runs/%s\n' % runFolder)
    script.write('/nfs/goldstein/software/python2.7/bin/python2.7 ~/github/sequencing_pipe/fastqc_mysql.py -i %s -s %s\n' % (runFolder,seqsata))
    script.close()
    submit(seqsata,FCID)
main()
