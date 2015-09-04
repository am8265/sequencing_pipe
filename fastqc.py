#Runs fastqc on the first last and middle fastq.gz file for each Sample, for each Lane, for each Read, for a flowcell.  Fastqc is run on the first last and middle fastq.gz files to compensate for biases due to read positions on the flowcell.  Reads from the first fastq.gz file should be from the first upper sections (artifically making the run look better) while reads fromthe last fastq.gz will be from the lower sections.
#
#joshbridgers@gmail.com
#01.05.2013

from commands import getoutput
from CHGV_mysql import *
from MySQLdb import Connect

import getopt
import os
import re
import sys

def fastqc(script,sample,seqsata,FCID,seqtype):
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
	seqsata_drive = seqsata.split('/')[2]
	script = open("%s/fastqc/%s_%s_fastqc_script.sh" % (seqsata,FCID,seqsata_drive),'w')
	script.write('#! /bin/bash\n')
	script.write('#\n')
	script.write('#$ -S /bin/bash -cwd\n')
	script.write('#$ -o %s/fastqc/%s_%s_fastqc_complete.sge\n' % (seqsata,FCID,seqsata_drive))
	script.write('#$ -e %s/fastqc/%s_%s_fastqc_out.sge\n' % (seqsata,FCID,seqsata_drive))
	script.write('#$ -V\n')
	script.write('#$ -M jb3816@cumc.columbia.edu\n')
	script.write('#$ -m bea\n')
	script.write('#\n')
	script.write('\n')
	
	return script,seqsata_drive

def submit(seqsata,FCID,seqsata_drive):
	print 'qsub -N %s_%s_fastqc %s/fastqc/%s_%s_fastqc_script.sh' % (FCID,seqsata_drive,seqsata,FCID,seqsata_drive)
	os.system('qsub -N %s_%s_fastqc %s/fastqc/%s_%s_fastqc_script.sh' % (FCID,seqsata_drive,seqsata,FCID,seqsata_drive))

def opts(argv):
	global FCID
	FCID = ''
	global seqsata
	seqsata = ''

	try:
		opts,args = getopt.getopt(argv, "f:i:", ['input=','FCID=',])
        except getopt.GetoptError, err:
                print str(err)
                usage()
	for o,a in opts:
		if o in ('-i','--input'): 
			seqsata = a
		elif o in ('-f','--FCID'):
			FCID = a
		else:
			assert False, "Unhandled argument present"

def main():
    opts(sys.argv[1:])
    pwd = os.getcwd()
    sequenceDB = getSequenceDB()

    #print 'ls %s/*%s*/Project_*/Sample_*/ -d | cut -d/ -f6 | cut -d_ -f2' % (seqsata,FCID)
    Samples = getoutput('ls %s/*%s*/Project_*/Sample_*/ -d | cut -d/ -f6 | cut -d_ -f2-' % (seqsata,FCID)).split('\n')

    script,seqsata_drive = script_header(seqsata,FCID)
    for samp in Samples:

        query = ' '.join((
            "SELECT distinct(UPPER(st.seqtype)) FROM Lane l JOIN Flowcell",
            "f ON l.fcid=f.fcid JOIN SeqType st ON l.prepid=st.prepid JOIN",
            "prepT p ON l.prepid=p.prepid WHERE FCILLUMID='"+FCID+"' AND",
            "CHGVID='"+samp+"'"))

        sequenceDB.execute(query)
        seqtype=sequenceDB.fetchall()
        if len(seqtype) > 1:
            print seqtype
            raise excecption, 'Too many seqtypes returned for sample %s' % samp
        fastqc(script,samp,seqsata,FCID,seqtype[0][0])

    script.write('cd %s; python2.7 ~/github/sequencing_pipe/fastqc_mysql.py -s %s\n' % (pwd,seqsata))
    script.close()
    submit(seqsata,FCID,seqsata_drive)
main()
