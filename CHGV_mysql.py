#CHGV_gaf.py module
from MySQLdb import Connect
from datetime import datetime
from commands import getoutput
import logging
import os
#mysql cursors

def getseqdb():
    seqdb=Connect(read_default_group='seqdb')
    seq_curs=seqdb.cursor()
    seq_curs.execute('BEGIN;')
    return seq_curs

def getGAFdb():
    gafdb=Connect(read_default_group="gaf")
    gaf_curs=gafdb.cursor()
    gaf_curs.execute('BEGIN;')
    return gaf_curs

def getTestSequenceDB():
    sequencedb=Connect(read_default_group='sequenceDBtest')
    sequence_curs=sequencedb.cursor()
    sequence_curs.execute('BEGIN;')
    return sequence_curs

def getSequenceDB():
    """Get the SequenceDB Cursor"""
    seqdb=Connect(read_default_group="sequenceDB")
    seq_curs=seqdb.cursor()
    seq_curs.execute('BEGIN;')
    return seq_curs

def getSampleID(gaf_curs,idnum):
        gaf_curs.execute("SELECT SampleID from Sample WHERE idnum=%s", (idnum))
        SampleID = gaf_curs.fetchone()
        return SampleID[0]

def getStatus(gaf_curs,idnum):
        gaf_curs.execute("SELECT Status from Sample WHERE idnum=%s", (idnum))
        Status = gaf_curs.fetchone()
        return Status[0]

def getSeqType(gaf_curs,idnum):
        gaf_curs.execute("SELECT SeqType from Sample WHERE idnum=%s", (idnum))
        Seqtype = gaf_curs.fetchone()
        return Seqtype[0]

def getSamYield(gaf_curs,idnum):
    gaf_curs.execute("SELECT SampleYield FROM Sample WHERE idnum=%s", (idnum))
    SamYield = gaf_curs.fetchone()
    if SamYield == None:
        SamYield = '0'
        return str(SamYield[0])

def getUserID():
    sequenceDB = getSequenceDB()
    userName = getoutput('echo $LOGNAME')
    if userName == 'solexa':
        userName = 'jb3816'
    sequenceDB.execute("SELECT userID FROM users WHERE netid='%s'" % (userName))
    userID = sequenceDB.fetchone()
    return userID[0]


def MachineCheck(sequenceDB,Machine,FCID):
    sql = """SELECT MACHINE
            FROM Flowcell
            WHERE FCIllumID like'%{0}%'
            """.format(FCID)

    sequenceDB.execute(sql)
    MachineFromDB = sequenceDB.fetchall()
    
    if len(MachineFromDB) != 1:
        raise Exception, "Incorrect number of FCIDs found!"
    else:
        if MachineFromDB[0][0] != Machine:

            #Updates database entry if incorrect.  Need to add logging
            sql = """UPDATE Flowcell
                    SET Machine='{0}' 
                    WHERE FCIllumID = '{1}'
                    """.format(Machine,FCID)

            sequenceDB.execute(sql)
            sequenceDB.execute('COMMIT;')

    return Machine

def getPredYield(gaf_curs,idnum):
        gaf_curs.execute("SELECT PredYield FROM Sample WHERE idnum=%s", (idnum))
        PredYield = gaf_curs.fetchone()
        if PredYield == None:
                PredYield = '0'
        return str(PredYield[0])

def setup_logging(machine,FCID,seqsata_drive):
	logging.basicConfig(level=logging.INFO,format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',filename='/nfs/%s/summary/GAF_PIPELINE_LOGS/%s_%s_%s.log' % (seqsata_drive,machine,FCID,seqsata_drive))
	logger = logging.getLogger(__name__)
	#create a file handler
	#handler = logging.FileHandler('/nfs/%s/summary/GAF_PIPELINE_LOGS/%s_%s_%s.log' % (seqsata_drive,machine,FCID,seqsata_drive))
	#handler.setLevel(logging.INFO)

	#create a logging format
	#formatter = logging.Formatter('%(asctime)s: %(name)s: %(levelname)s - %(message)s')

	#logger.addHandler(handler)

def getSamMultiplex(curs,Idnum):
	multiplexSamples = []
	curs.execute("SELECT LaneNum,FCID from Lane where samidnum='"+Idnum+"';")
	info = curs.fetchall()
	for lanes in info:
		curs.execute("SELECT SampleID from Lane where LaneNum='"+lanes[0]+"' and FCID='"+lanes[1]+"';")
		allSamp = curs.fetchall()
		for samp in allSamp:
			multiplexSamples.append(samp[0])

	return ' '.join([str(x) for x in set(multiplexSamples)])

def calcPredYield(gaf_curs,idnum):
        gaf_curs.execute("SELECT Sum(LnFraction) FROM Lane WHERE samidnum=%s and Fail='0'", (idnum))
        PredYield = gaf_curs.fetchone()
        gaf_curs.execute("SELECT LnNeedMan From S2R WHERE SampleID=%s", (getSampleID(gaf_curs,idnum)))
	S2Rlane = gaf_curs.fetchone()
	if PredYield == None or PredYield[0] == None or PredYield == ():
		PredYield = 0
	else:
		PredYield = PredYield[0]

        if S2Rlane == None or S2Rlane[0] == None or S2Rlane == ():
                S2Rlane = 0
        else:
                S2Rlane = float(S2Rlane[0])
        NewPredYield = (PredYield * 35000) + (S2Rlane * 35000)

	#print idnum,PredYield,S2Rlane,NewPredYield
        return str(int(round(float(NewPredYield))))

def getClusterDensity(gaf_curs,FCID,LaneNum,SampleID):
	gaf_curs.execute("SELECT ClustDen FROM Lane WHERE FCID=%s and LaneNum=%s and SampleID=%s", (FCID,LaneNum,SampleID))
	ClusterDensity = gaf_curs.fetchone()
	return str(ClusterDensity[0])

def getInsertSize(gaf_curs,idumn):
	GAFdb.execute("SELECT ROUND(AVG(MedInsSize)) from Lane where samidnum=%s", (idnum))
	avgInsertSize = gaf_curs.fetchone()
	return avgInsertSize

def getDatePool(gaf_curs,idnum):
        gaf_curs.execute("SELECT DatePool from Sample WHERE idnum=%s", (idnum))
        DatePool = gaf_curs.fetchone()
        return DatePool[0]

def getIdnum(gaf_curs,SampleID,SeqType):
	gaf_curs.execute("SELECT Idnum FROM Sample WHERE SampleID=%s and SeqType=%s", (SampleID,SeqType))
	idnum = gaf_curs.fetchone()
	#print SampleID,SeqType
	return str(int(idnum[0]))

def getIndexBaseFromIdnum(gaf_curs,idnum):
	#print idnum
	gaf_curs.execute("SELECT IndexBase FROM Sample WHERE idnum=%s", (idnum))
	IndexBase = gaf_curs.fetchone()
	return IndexBase[0]

def getPL(seq_curs,SampleID):
	seq_curs.execute("SELECT CurrProjLeader from build37 where CHGVID=%s", (SampleID))
	PL = seq_curs.fetchone()
	if PL == None:
		raise Exception, "%s or PL does not exist in SeqDB build37" % (SampleID)
	else:
		return PL[0]

def getIndexBase(IndexNum,SampleID):
        if IndexNum == 'na' or IndexNum == None or IndexNum == '0' or IndexNum == 'none':
                IndexNum = ''
                IndexBase = ''
        elif IndexNum == '':
                IndexBase = ''
        elif int(IndexNum) < 28 :
                IndexBase = getoutput('egrep ^'+IndexNum+ ', /home/solexa/indices.csv | cut -d, -f2')
        elif int(IndexNum) > 27:
                raise Exception, 'Index Number '+IndexNum+' does not exist!'
        else:
                raise Exception, 'Check Index for '+SampleID+'!'
        return IndexBase

def getNM(gaf_curs,idnum):
	gaf_curs.execute("SELECT KapaNM from Sample WHERE idnum=%s", (idnum))
	NM = gaf_curs.fetchone()
	if NM == None:
		return 'NULL'
	else:
		return NM[0]

def cur_datetime():
        tDate = datetime.now()
        Date = str(tDateBCL).split('.')[0]
	return Date

def getPM(gaf_curs,idnum):
        gaf_curs.execute("SELECT KapaPM from Sample WHERE idnum=%s", (idnum))
        PM = gaf_curs.fetchone()
        if PM == None:
                return 'NULL'
        else:
                return PM[0]

def Idnum_checker(gaf_curs,SampleID):
        entry_num = gaf_curs.execute("SELECT Idnum FROM Sample WHERE SampleID=%s", (SampleID))
	idnum = gaf_curs.fetchone()
        if entry_num > 2:
                raise Exception, 'Error with %s: 2 entries for %s already exist!' % (SampleID,SampleID)
        elif entry_num == 2:
                genome_count = gaf_curs.execute("SELECT Idnum FROM Sample Where SampleID=%s and SeqType='genome'", (SampleID))
		idnum = gaf_curs.fetchone()
                if genome_count > 1:
                        raise Exception, 'More than one genome entry exists for '+SampleID+'!'
                else:
			#code for cases when exome and genome entry exists
        		raise Exception, 'Exome and genome exist, code needed for this situation'
	elif entry_num == 1:
                return str(int(idnum[0]))
        elif entry_num == 0:
                raise Exception, SampleID+' has not been entered into received_samples.csv!'

	'''
	elif entry_num == 3:
		gaf_curs.execute("SELECT Idnum FROM Sample Where SampleID=%s and SeqType='RNAseq'", (SampleID))
		RNAseq_idnum = gaf_curs.fetchone()
		return str(int(RNAseq_idnum[0]))
	'''
