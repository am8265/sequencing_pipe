#release_GAFdb_v2.01_beta.py
#Updates status for SequenceDB with a sample list

import getopt
import glob
import logging
import os
import sys
import traceback
from commands import getoutput
from datetime import datetime
from MySQLdb import Connect
from CHGV_mysql import getSequenceDB
from CHGV_mysql import getUserID
from checkRelease import checkUpdateDB

#Checks sample list
def Samples():
    path = sys.argv[1]
    if os.path.isfile(path) == False:
        print '/nfs/genotyping/bridgers/'+path
        raise Exception, "sample list file does not exist"
    return path

def checkS2R(sequenceDB,SampleID,prepID):
    '''check that all lanes are complete in S2R'''
    sequenceDB.execute("SELECT complete FROM samplesTOrun s2r JOIN poolMembers pm ON s2r.dbID=pm.poolid WHERE pm.prepID=%s AND complete=0 AND fail=0", prepID)
    S2Rsamp = sequenceDB.fetchall()
    #print S2Rsamp,SampleID,prepID
    #print "SELECT complete FROM samplesTOrun s2r JOIN poolMembers pm ON s2r.dbID=pm.poolid WHERE pm.prepID=%s AND complete=0", prepID
    if S2Rsamp != () :
        raise Exception, "%s still in S2R with uncompleted lanes" % (SampleID)

def check_Lane(sequenceDB,prepID,CHGVID):
	#check to see if LnYield is populated for another check
	sequenceDB.execute("SELECT prepID,LnYield FROM Lane WHERE prepID=%s and FailR1 is NULL and FailR2 is NULL and LnYield is NULL", prepID)
	sequencing_lanes = sequenceDB.fetchall()
	#print sequencing_lanes
	if sequencing_lanes != ():
		raise Exception, "Sample %s is still sequencing or LnYield is NULL!" % CHGVID

#Determines what Platforms the sample was on
def getPlatformChemVer(sequenceDB,prepID):
    sequenceDB.execute("SELECT DISTINCT f.Machine,f.FCillumID FROM Lane l JOIN Flowcell f ON l.FCID=f.FCID WHERE l.prepID=%s", prepID)
    Machines = sequenceDB.fetchall()
    ChemVer = 'v4'
    query="SELECT DISTINCT chemver FROM Lane l JOIN Flowcell f ON l.FCID=f.FCID WHERE l.prepID=%s" % prepID

    if verbose == True:
        print query
    sequenceDB.execute(query)
    ChemVer = sequenceDB.fetchall()[0]

    if type(ChemVer) == tuple:
        ChemVer = ChemVer[0]

    #print Machines,prepID
    H2500 = 0
    H2000 = 0
    HX = 0 #placeholder for HiSeq X
    rapid = 0

    for HiSeq in Machines:
        #print HiSeq,HiSeq[0][0]
        if HiSeq[1][0] == 'H':
            rapid = 1

    if HiSeq[0] is None:
        H2500 = 1
    elif HiSeq[0] == 'H8A' or HiSeq[0] == 'H8B' or len(HiSeq[0]) > 3:
        H2500 = 1
    else:
        H2000 = 1

    #print rapid,H2500,H2000
    if H2500 == 1 and H2000 == 1:
        platform = 'HiSeq2000 and HiSeq2500'
    elif H2000 == 1 and H2500 == 0:
        platform = 'HiSeq2000'
    elif H2000 == 0 and H2500 == 1:
        if rapid == 1:
            platform = 'HiSeq2500'
        else:
            platform = 'HiSeq2500'
    else:
        raise Exception, "Check Machine for prepID %s!" % (prepID)
    return platform,ChemVer

#Processes Fastq Location information from the Flowcell and Lane tables
def getFastqLoc(sequenceDB,SampleID,prepID,seqtype):
    sequenceDB.execute("SELECT GROUP_CONCAT(DISTINCT(CONCAT('/nfs/',f.SeqSataLoc,'/',%s,'/',%s,'/',fcillumid)) SEPARATOR ' ') FROM Flowcell f JOIN Lane l ON f.FCID=l.FCID WHERE SeqSataLoc IS NOT NULL and SeqSataLoc != '' and prepID=%s", (seqtype,SampleID,prepID))
    FastqLoc = sequenceDB.fetchone()[0]
    #print FastqLoc
    return FastqLoc

def getPoolID(sequenceDB,prepID):
    sql=("SELECT "
            "(SELECT CHGVID FROM prepT "
            "WHERE dbid="
                "(SELECT poolID FROM Lane l "
                "JOIN Flowcell f ON l.FCID=f.FCID "
                "WHERE l.prepID={0} ORDER BY FCtime DESC LIMIT 1)) "
        "poolID,exomeKit,libKit FROM prepT "
        "WHERE prepID={0};"
        ).format(prepID)

    #print sql
    #raw_input('test')
    sequenceDB.execute(sql)
    info = sequenceDB.fetchone()
    return info

def sequenceDB_exist_check(sequenceDB,prepID):
	logger = logging.getLogger('sequenceDB_exist_check')
	logger.info("SELECT idnum FROM seqdbClone WHERE prepID=%s" % prepID)
	sequenceDB.execute("SELECT idnum FROM seqdbClone WHERE prepID=%s", prepID)
	samp = sequenceDB.fetchone()
	logger.debug(samp,prepID)

	if samp == None:
		return 0
	else:
		return 1

def checkExomeKit(sequenceDB,ExomeSamPrepKit):
    #print ExomeSamPrepKit
    sequenceDB.execute("SELECT * FROM captureKit where name=%s and chr='all'"\
            ,(ExomeSamPrepKit))
    isCaptureKit = sequenceDB.fetchone()

    #print isCaptureKit
    if isCaptureKit is None:
        raise Exception, "%s: ExomeSamPrepKit information does not exist" % sequencing_info[1]


def process_exomekit(sequencing_info,seqtype,sequenceDB):

    if seqtype == 'exome':
        if sequencing_info[1] == () or sequencing_info[1] == '':
            raise Exception, '%s does not have Exome Prep Kit information' % SampleID
        elif 'EZ V3' in sequencing_info[1]:
            ExomeSamPrepKit = 'Roche'
        elif sequencing_info[1].isdigit() == True:
           ExomeSamPrepKit = sequencing_info[1]+'MB'
        else:
           ExomeSamPrepKit =  sequencing_info[1]
        checkExomeKit(sequenceDB,ExomeSamPrepKit)

    elif seqtype == 'Custom Capture' or seqtype == 'Custom_Capture':
        if sequencing_info[1] == () or sequencing_info[1] == '':
            raise Exception, '%s does not have Exome Prep Kit information' % SampleID
        elif 'SczEpi' in sequencing_info[1] or 'SchizEpi' in sequencing_info[1]:
            ExomeSamPrepKit = 'SchizoEpi'
        elif 'Roche SeqCap Custom - EpiMIR' in sequencing_info[1]:
            ExomeSamPrepKit = 'EpiMIR'
        elif 'Agilent SureSelect Custom - DILIN' in sequencing_info[1]:
            ExomeSamPrepKit = 'DILIN'
        elif 'MattHaloplex2015' in sequencing_info[1]:
            ExomeSamPrepKit = 'MattHaloplex2015'
        elif 'MTOR' in sequencing_info[1]:
            ExomeSamPrepKit = 'MCDMTOR'
        else:
            ExomeSamPrepKit =  sequencing_info[1]

        checkExomeKit(sequenceDB,ExomeSamPrepKit)
    else:
        ExomeSamPrepKit = 'N/A'
    return ExomeSamPrepKit

def getSampleTinfo(sequenceDB,prepID):
    #Note that CHGVID comes from the prepT.  This is due to redmine issue #2038
    sql = ("SELECT pt.CHGVID,Phenotype,FamilyID,DNASrc,GwasID,"
        "FinalRepLoc,Topstrandfile,AKA,ProjName,Protocol,AvaiContUsed,"
        "SelfDeclEthnic,SelfDeclEthnicDetail,SelfDeclGender,GwasGender,"
        "GwasEthnic,GenoChip,GAFbin,FundCode,RepConFamMem,"
        "FamilyRelationProband,priority,OrigID,CurrProjLeader,AlignedBuild36,"
        "priority,BroadPhenotype,DetailedPhenotype "
        "from SampleT s "
        "JOIN prepT pt on s.DBID=pt.DBID "
        "WHERE pt.prepID={0}"
        ).format(prepID)
    sequenceDB.execute(sql)
    SampleT_info = sequenceDB.fetchone()
    #print SampleT_info
    return SampleT_info

def constructMySql(qType,tableName,columns):
    if qType != 'UPDATE' and qType != 'INSERT':
        raise Exception, "Unhandled qType: %s" % qType
    else:
        if qType == 'INSERT':
            qType += ' INTO'
        query = "%s %s SET " % (qType,tableName)
        for item in columns:
            if len(item) != 3:
                raise Exception, "Improper column item: %s" % item
            if item[2] == 0:
                query += "%s='%s'," % (item[0],item[1])
            elif item[2] == 1:
                query += "%s=%s," % (item[0],item[1])
            else:
                raise Exception, "Improper value for columns list: %s" % item[2]
    return query

def getplateInfo(sequenceDB,DBID,prepID):
    sequenceDB.execute("SELECT plateName,Well,Concentration,VolWell,PicogreenValue FROM plateLocation pl JOIN plates p ON p.plateID=pl.plateID WHERE p.plateID=(SELECT parentPlateID FROM plateLocation pl JOIN plates p on pl.plateID=p.plateID WHERE PREPID=%s AND (plateType='RNA QC' or plateType='Fragmentation') ORDER BY locationDate DESC LIMIT 1) AND DBID=%s", (prepID,DBID))
    plateInfo = sequenceDB.fetchone()
    #print plateInfo

    if plateInfo:
        return plateInfo
    else:
        return ['','','','','']

def updateDB(sequenceDB,prepID,SampleID,seqtype,DBID):
    logger = logging.getLogger('updateDB')
    #In SequenceDB seqtype "Custom Capture" is represented as "Custom_Capture".  Note the underscore
    if seqtype == 'Custom Capture':
            seqtype = 'Custom_Capture'

    if old == False:
        sequencing_info = getPoolID(sequenceDB,prepID)
        PoolID = sequencing_info[0]
        Platform,ChemVer = getPlatformChemVer(sequenceDB,prepID)
        FastqLoc = getFastqLoc(sequenceDB,SampleID,prepID,seqtype.upper())
        ExomeSamPrepKit = process_exomekit(sequencing_info,seqtype,sequenceDB)
        PrepKit = sequencing_info[2]
        SamMultiplex = getSamMultiplex(sequenceDB,prepID)
        CoreSeqSoftware = 'CASAVA 1.8'
        sInfo = getSampleTinfo(sequenceDB,prepID)
        plateInfo = getplateInfo(sequenceDB,DBID,prepID)

        #Check whether there are already seqdbClone entries
        sequenceDB_exist = sequenceDB_exist_check(sequenceDB,prepID)

        #list of tuples containing keys and values for a MySQL UPDATE or INSERT statement.  The third item is whether to wrap the value in single quotes for the MySQL statement
        columns = []
        columns.append(('Status','Released to Bioinformatics Team',0))
        columns.append(('FastqLoc',FastqLoc,0))
        columns.append(('SeqCoreRelDate','curdate()',1))
        columns.append(('PoolID',PoolID,0))
        columns.append(('Platform',Platform,0))
        columns.append(('ChemVer',ChemVer,0))
        columns.append(('ExomeSamPrepKit',ExomeSamPrepKit,0))
        columns.append(('PrepKit',PrepKit,0))
        columns.append(('SamMultiplex',SamMultiplex,0))
        columns.append(('StatusChangeDate','curdate()',1))
        columns.append(('Phenotype',sInfo[1],0))
        columns.append(('BroadPhenotype',sInfo[26],0))
        try:
            columns.append(('DetailedPhenotype',sInfo[27].replace("'",''),0))
        except:
            columns.append(('DetailedPhenotype',sInfo[27],0))
        columns.append(('FamilyID',sInfo[2],0))
        columns.append(('DNASrc',sInfo[3],0))
        columns.append(('GwasID',sInfo[4],0))
        columns.append(('FinalRepLoc',sInfo[5],0))
        columns.append(('Topstrandfile',sInfo[6],0))
        columns.append(('AKA',sInfo[7],0))
        columns.append(('ProjName',sInfo[8],0))
        columns.append(('Protocol',sInfo[9],0))
        columns.append(('AvaiContUsed',sInfo[10],0))
        columns.append(('SelfDeclEthnic',sInfo[11],0))
        columns.append(('SelfDeclEthnicDetail',sInfo[12],0))
        columns.append(('SelfDeclGender',sInfo[13],0))
        columns.append(('GwasGender',sInfo[14],0))
        columns.append(('GwasEthnic',sInfo[15],0))
        columns.append(('GenoChip',sInfo[16],0))
        columns.append(('GAFbin',sInfo[17],0))
        columns.append(('FundCode',sInfo[18],0))
        columns.append(('RepConFamMem',sInfo[19],0))
        columns.append(('FamilyRelationProband',sInfo[20],0))
        columns.append(('OrigID',sInfo[22],0))
        columns.append(('prepID',prepID,0))
        columns.append(('FastqLoc',FastqLoc,0))
        columns.append(('SeqCoreRelDate','curdate()',1))
        columns.append(('PoolID',PoolID,0))
        columns.append(('Platform',Platform,0))
        columns.append(('ChemVer',ChemVer,0))
        columns.append(('ExomeSamPrepKit',ExomeSamPrepKit,0))
        columns.append(('PrepKit',PrepKit,0))
        columns.append(('SamMultiplex',SamMultiplex,0))
        columns.append(('CoreSeqSoftware',CoreSeqSoftware,0))
        columns.append(('CurrProjLeader',sInfo[23],0))
        columns.append(('AlignedBuild36',sInfo[24],0))
        columns.append(('BioinfoPriority',sInfo[25],0))
        columns.append(('PlateName',plateInfo[0],0))
        columns.append(('Well',plateInfo[1],0))
        columns.append(('Concentration',plateInfo[2],0))
        columns.append(('VolWell',plateInfo[3],0))
        columns.append(('PicogreenValue',plateInfo[4],0))

        if sequenceDB_exist == 1:
            UpdateQuery = constructMySql('UPDATE','seqdbClone',columns)
            UpdateQuery = UpdateQuery[:-1] + " WHERE prepID='%s'" % prepID
            if verbose == True:
                print UpdateQuery

            logger.info(UpdateQuery)
            sequenceDB.execute(UpdateQuery)

        if sequenceDB_exist != 1:
            columns = []
            columns.append(('Status','Released to Bioinformatics Team',0))
            columns.append(('CHGVID',sInfo[0],0))
            columns.append(('Seqtype',seqtype,0))
            columns.append(('Phenotype',sInfo[1],0))
            columns.append(('BroadPhenotype',sInfo[26],0))
            try:
                columns.append(('DetailedPhenotype',sInfo[27].replace("'",''),0))
            except:
                columns.append(('DetailedPhenotype',sInfo[27],0))
            columns.append(('FamilyID',sInfo[2],0))
            columns.append(('DNASrc',sInfo[3],0))
            columns.append(('GwasID',sInfo[4],0))
            columns.append(('FinalRepLoc',sInfo[5],0))
            columns.append(('Topstrandfile',sInfo[6],0))
            columns.append(('AKA',sInfo[7],0))
            columns.append(('ProjName',sInfo[8],0))
            columns.append(('Protocol',sInfo[9],0))
            columns.append(('AvaiContUsed',sInfo[10],0))
            columns.append(('SelfDeclEthnic',sInfo[11],0))
            columns.append(('SelfDeclEthnicDetail',sInfo[12],0))
            columns.append(('SelfDeclGender',sInfo[13],0))
            columns.append(('GwasGender',sInfo[14],0))
            columns.append(('GwasEthnic',sInfo[15],0))
            columns.append(('GenoChip',sInfo[16],0))
            columns.append(('GAFbin',sInfo[17],0))
            columns.append(('FundCode',sInfo[18],0))
            columns.append(('RepConFamMem',sInfo[19],0))
            columns.append(('FamilyRelationProband',sInfo[20],0))
            columns.append(('OrigID',sInfo[22],0))
            columns.append(('idnum',DBID,0))
            columns.append(('prepID',prepID,0))
            columns.append(('FastqLoc',FastqLoc,0))
            columns.append(('SeqCoreRelDate','curdate()',1))
            columns.append(('PoolID',PoolID,0))
            columns.append(('Platform',Platform,0))
            columns.append(('ChemVer',ChemVer,0))
            columns.append(('ExomeSamPrepKit',ExomeSamPrepKit,0))
            columns.append(('PrepKit',PrepKit,0))
            columns.append(('SamMultiplex',SamMultiplex,0))
            columns.append(('CoreSeqSoftware',CoreSeqSoftware,0))
            columns.append(('CurrProjLeader',sInfo[23],0))
            columns.append(('AlignedBuild36',sInfo[24],0))
            columns.append(('BioinfoPriority',sInfo[25],0))
            columns.append(('StatusChangeDate','curdate()',1))
            columns.append(('PlateName',plateInfo[0],0))
            columns.append(('Well',plateInfo[1],0))
            columns.append(('Concentration',plateInfo[2],0))
            columns.append(('VolWell',plateInfo[3],0))
            columns.append(('PicogreenValue',plateInfo[4],0))

            InsertQuery = constructMySql('INSERT','seqdbClone',columns)
            InsertQuery = InsertQuery[:-1]
            if verbose == True:
                print InsertQuery
            logger.info(InsertQuery)
            sequenceDB.execute(InsertQuery)


def getSamMultiplex(sequenceDB,prepID):
	multiplexSamples = []
	sequenceDB.execute("SELECT LaneNum,FCID FROM Lane WHERE prepID=%s", prepID)
	info = sequenceDB.fetchall()
	for lanes in info:
		sequenceDB.execute("SELECT pt.CHGVID FROM prepT pt JOIN Lane l ON pt.prepID=l.prepID WHERE l.LaneNum=%s and l.FCID=%s;",(lanes[0],lanes[1]))
		allSamp = sequenceDB.fetchall()
		for samp in allSamp:
			multiplexSamples.append(samp[0])
	return ' '.join([str(x) for x in set(multiplexSamples)])

def getPL(sequenceDB,SampleID,release_list,prepID):
	'''Grab project lead's (PL) email address and add it to the release_list '''
	sequenceDB.execute("SELECT u.email FROM SampleT s JOIN SeqType st ON s.DBID=st.DBID JOIN users u on (u.name like CONCAT('%',s.CurrProjLeader,'%')) WHERE st.prepid='"+prepID+"';")
	PL = sequenceDB.fetchone()

	#print PL
	if PL == ():
		raise Exception, SampleID+' does not exist in SequenceDB!'
	else:
		release_list.append((SampleID,PL[0]))

def usage():
    print
    print '='*80
    print 'USAGE: python release_GAFdb.py [parameters] [sample list file]'
    print
    print '-c, --custom\t\tRelease custom capture samples'
    print '-e, --email\t\tEmail is sent out to project leads'
    print '-g, --genome\t\tRelease genome samples'
    print '-k, --capturekit\tSpecify Exomesamprepkit of samples'
    print '-h, --help\t\tShows this help message and exit'
    print '-o, --old\t\tUse for samples without sequencing information in sequenceDB'
    print '-n, --noemail\t\tNo emails are sent out to the Project Leads'
    print '-r, --rna\t\tRelease RNASeq samples'
    print '-t, --test\t\tUses jb3816\'s test MySQL database for queries and updates'
    print '-v, --verbose\t\tPrints out MySQL injection commands'
    print '-x, --exome\t\tRelease exome samples'
    print
    print '--ignoreMultiPrepIDs\tIgnores multiple prepID check for SeqType and chooses'
    print '\t\t\thighest prepID for release. Use with the -k, --capturekit option.'
    print '\t\t\tUse with caution.'
    print '--nostatus\t\tDoes not perform status updates on statusT'
    print "--skipStorageCheck\tSkips check that seqdbClone's FastqLoc is correct"
    print '='*80
    print
    sys.exit(2)

def convert_release_list(sequenceDB,release_list):
    email = {}
#get all PLs from seqdb
    sequenceDB.execute("SELECT email from users")
    PLs = sequenceDB.fetchall()
    for PL in PLs:
        email[PL[0]] = []

    for samp in list(set(release_list)):
        if samp[1] in email.keys():
            email[samp[1]].append(samp[0])
        else:
            raise Exception, samp[0] + ' does not have a PL!'

    for name in email.keys():
        if email[name] != [] and sendemail == True:
            email_PL(sequenceDB,email[name],name)

def email_PL(sequenceDB,samples,name):
    sequenceDB.execute("SELECT Email FROM users WHERE email='"+name+"';")
    address = sequenceDB.fetchone()[0]
    email = open('email.tmp','w')
    email.write('The following sample(s) are being released to the Bioinformatics core:\n\n')
    for samp in sorted(samples):
            email.write(samp+'\n')
    email.write('\nThis is an AUTOMATED message.  If you have any questions please direct them to jb3816@cumc.columbia.edu\n')
    email.close()
    print "Sending email to "+address
    #os.system("cat email.tmp")
    os.system("/nfs/goldstein/software/mutt-1.5.23 -s 'Sample Release' "+address+" < email.tmp")

def check_sequenceDB(sequenceDB,SampleID,prepID,SeqType):
    sequenceDB.execute("SELECT status from statusT WHERE (prepID=%s and Status='Sequencing Complete') or (prepID=%s and Status='External Data Submitted')", (prepID,prepID))
    complete = sequenceDB.fetchone()
    if complete is None:
        raise Exception, "Flowcell has not been completed SequenceDB for %s.  Old sample?" % SampleID
    sequenceDB.execute("SELECT SUM(LnYield) from Lane where prepID=%s and failr1 is NULL and failr2 is NULL", prepID)
    info = sequenceDB.fetchone()

    sequenceDB.execute("SELECT status from statusT WHERE prepID=%s ORDER BY status_time DESC LIMIT 1", (prepID))
    status = sequenceDB.fetchone()
    #sequenceDB.execute("SELECT status from statusT WHERE prepID=%s AND status='External Data Submitted'", prepID)
    #externalCheck = sequenceDB.fetchone()

    failYield = ''
    poolID = ''

    #-->code a better check for External Data Submitted <--#
    if status[0] == 'External Data Submitted':
        pass
    elif int(info[0]) < 100000 and SeqType.lower() == 'genome':
        fail_switch = raw_input("%s has yield < 100000 (%s).  Is this ok to release (Y)es or (N)o? " % (SampleID,int(info[0])))
        if fail_switch.lower() != 'y' and status[0] != 'External Data Submitted':
            failYield = SampleID
            #raise Exception, "%s has yield of %s!" % (SampleID,info[0])

    elif int(info[0]) < 7500 and SeqType.lower() == 'exome':
        fail_switch = raw_input("%s has yield < 7500 (%s).  Is this ok to release (Y)es or (N)o? " % (SampleID,int(info[0])))
        if fail_switch.lower() != 'y' and status[0] != 'External Data Submitted':
            failYield = SampleID
            poolID = getPoolID(sequenceDB,prepID)
            #print poolID[0]

            #raise Exception, "%s has yield of %s!" % (SampleID,info[0])

    elif int(info[0]) < 150 and SeqType == 'Custom Capture':
        fail_switch = raw_input("%s has yield < 150.  Is this ok to release (Y)es or (N)o? " % (SampleID))
        if fail_switch.lower() != 'y' and status[0] != 'External Data Submitted':
            failYield = SampleID
            poolID = getPoolID(sequenceDB,prepID)

            #raise Exception, "%s has yield of %s!" % (SampleID,info[0])

    #print SampleID,prepID
    #print status[0]
    return status[0],complete[0],failYield,poolID,info[0]

def fixReleaseStatusT(SequenceDB,SampleID,SeqType,prepID):
    logger = logging.getLogger('fixReleaseStatusT')

    print "Deleted Release status from %s %s %s on statusT" % (SampleID,SeqType,prepID)

    logger.info("DELETE FROM statusT WHERE prepID=%s and status='Released to Bioinformatics Team" % prepID)
    #SequenceDB.execute("DELETE FROM statusT WHERE prepID=%s and status='Released to Bioinformatics Team", prepID)


def check_Storage(sequenceDB,prepID,SeqType):
#check that all fastq files exist
    query = ("SELECT p.CHGVID,f.FCIllumID,l.LaneNum "
        "FROM prepT p "
        "JOIN Lane l on p.DBID=l.DBID "
        "JOIN Flowcell f on f.FCID=l.FCID "
        "WHERE FailR1 IS NULL AND "
        "FailR2 IS NULL "
        "AND prepID='{}'"
        ).format(prepID)

    sequenceDB.execute(query)
    Lanes = sequenceDB.fetchall()
    for lane in Lanes:
        SampleID = lane[0]
        FCID = lane[1]
        LaneNum = lane[2]
        #print "/nfs/seqsata*/seqfinal*/whole_genome/%s/%s/%s*L00%s*R1*fastq.gz" % (SampleID,FCID,SampleID,LaneNum)

        #checks how many files are there for that FC, LaneNum and Sample.  If 0 warning sent out

        SeqType=SeqType.upper().replace(' ','_')
        r1files = len(glob.glob("/nfs/fastq18/%s/%s/%s/%s*L00%s*R1*fastq.gz" % (SeqType,SampleID,FCID,SampleID,LaneNum)))
        r2files = len(glob.glob("/nfs/fastq18/%s/%s/%s/%s*L00%s*R2*fastq.gz" % (SeqType,SampleID,FCID,SampleID,LaneNum)))

        if r1files == 0 or r2files == 0:
            r1files = len(glob.glob("/nfs/seqsata*/seqfinal*/whole_genome/%s/%s/%s*L00%s*R1*fastq.gz" % (SampleID,FCID,SampleID,LaneNum)))
            r2files = len(glob.glob("/nfs/seqsata*/seqfinal*/whole_genome/%s/%s/%s*L00%s*R2*fastq.gz" % (SampleID,FCID,SampleID,LaneNum)))

        if r1files != r2files:
            #checks to see if Lane was SE
            sequenceDB.execute("SELECT LaneNum FROM Lane where MmR2='0' and PerQ30R2='0' and MnQscR2='0' and MedInsSize='0' and AbvStDev='0' and BelStDev='0' and FCID=%s and SampleID=%s and LaneNum=%s", (FCID,SampleID,LaneNum))
            SE_Lane = sequenceDB.fetchall()
            if SE_Lane == ():
                raise Exception, "%s for FCID %s does not have equal number of files for Read1 and Read2" % (SampleID,FCID)
        elif r1files == 0 or r2files == 0:
                raise Exception, "%s for FCID %s does not have any fastq.gz files for Lane %s" % (SampleID,FCID,LaneNum)

def subject_header(SeqType,today,today2,release_loc,Version):
    if SeqType == 'exome':
        Subject = 'Exome Release '+today2+' ('+release_loc+')'
        run_sum = today+'_exRun_Summary_'+release_loc+'_v'+Version+'.txt'
    elif SeqType == 'genome':
        Subject = 'Whole Genome Release '+today2+' ('+release_loc+')'
        run_sum = today+'_wgRun_Summary_'+release_loc+'_v'+Version+'.txt'
    elif SeqType == 'Custom Capture':
        Subject = 'Custom Capture Release '+today2+' ('+release_loc+')'
        run_sum = today+'_ccRun_Summary_'+release_loc+'_v'+Version+'.txt'
    elif SeqType == 'RNASeq':
        Subject = 'RNASeq Release '+today2+' ('+release_loc+')'
        run_sum = today+'_rsRun_Summary_'+release_loc+'_v'+Version+'.txt'
    return Subject,run_sum

def updateStatus(sequenceDB,prepID,SampleID,DBID):
    logger = logging.getLogger('updateStatus')

    global userID
    userID = getUserID()

    logger.info("INSERT INTO statusT(CHGVID,status_time,status,DBID,prepID,userID) VALUES(%s,unix_timestamp(),'Released to Bioinformatics Team',%s,%s,%s)" % (SampleID,DBID,prepID,userID))
    if verbose == True:
        print "INSERT INTO statusT(CHGVID,status_time,status,DBID,prepID,userID) VALUES('%s',unix_timestamp(),'Released to Bioinformatics Team','%s','%s','%s')" % (SampleID,DBID,prepID,userID)
    sequenceDB.execute("INSERT INTO statusT(CHGVID,status_time,status,DBID,prepID,userID) VALUES(%s,unix_timestamp(),'Released to Bioinformatics Team',%s,%s,%s)", (SampleID,DBID,prepID,userID))

    logger.info("UPDATE seqdbClone SET Status='Released to Bioinformatics Team', SeqCoreRelDate=curdate() WHERE prepID='%s'" % prepID)
    if statusOn == True:
        sequenceDB.execute("UPDATE seqdbClone SET Status='Released to Bioinformatics Team',SeqCoreRelDate=curdate() WHERE prepID=%s", prepID)
    if verbose == True:
        print "UPDATE seqdbClone SET Status='Released to Bioinformatics Team',SeqCoreRelDate=curdate() WHERE prepID='%s'" % prepID




def get_release_var():

    #release_loc = 'IGMC'
	release_loc = 'LSRC'

	Summary_Path = '/nfs/sva01/Summaries/'

	tVersion = raw_input('Version 1.7 or 1.8? ')
	if tVersion == '1.7':
		Version = '1.7'
		CoreSeqSoftware = 'CASAVA 1.7'
	elif tVersion == '1.8':
		Version = '1.8'
		CoreSeqSoftware = 'CASAVA 1.8'

	return release_loc,Summary_Path,Version,CoreSeqSoftware

def release_email(sequenceDB,SeqType,IDs,failedSamples):
    logger = logging.getLogger('release_email')
    release_loc,Summary_Path,Version,CoreSeqSoftware = get_release_var()

    today = datetime.today().strftime("%y.%m.%d")
    today2 = datetime.today().strftime("%m.%d.%y")
    info = subject_header(SeqType,today,today2,release_loc,Version)
    Subject = info[0]
    run_sum = info[1]

    sampleNumber = len(IDs)

    release_email = open('release_email.txt','w')
    release_email.write('The following %s %s sample(s) are ready for alignment (BUILD37 v%s):\n' % (sampleNumber,SeqType,Version))
    release_email.write('\n')

    RunSumPath = Summary_Path+run_sum
    #print RunSumPath,SeqType
    """
    if os.path.isfile(RunSumPath) == True:
        pass
    elif release_loc == 'DSCR':
        pass

    else:
        print RunSumPath
        raise Exception, 'Run Summary file does not exist!'
    """
    for samp in IDs.keys():
        prepID = IDs[samp][1]
        DBID = IDs[samp][0]
        release_email.write("%s\t%s\n" % (prepID,samp))
        SampleID = str(samp.strip())

        if statusOn == True:
            updateStatus(sequenceDB,prepID,SampleID,DBID)

        #print samp,prepID
        if verbose == True:
            print "UPDATE seqdbClone SET CoreSeqSoftware='%s', ReadSummaryFileLoc='%s' WHERE prepID='%s'" % (CoreSeqSoftware,RunSumPath,prepID)
        sequenceDB.execute("UPDATE seqdbClone SET CoreSeqSoftware=%s, ReadSummaryFileLoc=%s WHERE prepID=%s",(CoreSeqSoftware,RunSumPath,prepID))
        logger.info("UPDATE seqdbClone SET CoreSeqSoftware='%s', ReadSummaryFileLoc='%s' WHERE prepID='%s'" % (CoreSeqSoftware,RunSumPath,prepID))

    sequenceDB.execute("SELECT name FROM users WHERE userid=%s", userID)
    name = sequenceDB.fetchone()[0]

    print; print '==========Release Email=========='
    print Subject;print

    release_email.write('\n'+RunSumPath+'\n')
    release_email.write('\nThanks,\n')
    release_email.write('\n')
    release_email.write('%s\n' % name)

    release_email.close()


    os.system('cat release_email.txt')
    print
    print "Samples not released:"
    print
    print "IGMID\t\tPoolID\t\tSeqYield"
    print "=" * 80
    for failedSample in failedSamples:
        print "%s\t\t%s\t\t%s" % (failedSample[0],failedSample[1][0],failedSample[2])

    os.system('rm release_email.txt')

def opts(argv):
    global capturekit ; capturekit = ''
    global ignore ; ignore = False
    global old ; old = False
    global sendemail ; sendemail = False
    global SeqType ; SeqType = 'NA'
    global skipStorageCheck ; skipStorageCheck = False
    global statusOn ; statusOn = True
    global test ; test = False
    global verbose ; verbose = False

    try:
        opts,args = getopt.getopt(argv, "henotvxgcrk:", ['help','email','noemail','old','verbose','exome','genome','custom','rna','nostatus','capturekit','test','ignoreMultiPrepIDs','skipStorageCheck'])
    except getopt.GetoptError, err:
        #print str(err)
        usage()

    #options
    for o,a in opts:
        if o in ('-e','--email'): #Send email to PL
            sendemail = True
        elif o in ('-h','--help'): #Invoke help
            usage()
        elif o in ('--ignoreMultiPrepIDs'):
            ignore = True
        elif o in ('-k','--capturekit'): #Specify capture kit
            capturekit = a
        elif o in ('-n','--noemail'): #Don't send email to PL (Default)
            sendemail = False
        elif o in ('--nostatus'): #Do not perform status update
            statusOn = False
        elif o in ('-o','--old'): #Sample is from old database
            old = True
        elif o in ('--skipStorageCheck'):
            skipStorageCheck = True
        elif o in ('-t','--test'): #Update test database
            test = True
        elif o in ('-v','--verbose'): #Verbose output
            verbose = True

        #SeqTypes
        elif o in ('-r','--rna'): #RNA
            if SeqType!='NA':
                raise Exception, "Conflicting SeqTypes were specified"
            SeqType = 'RNASeq'
        elif o in ('-x','--exome'): #Exome
            if SeqType!='NA':
                raise Exception, "Conflicting SeqTypes were specified"
            SeqType = 'exome'
        elif o in ('-g','--genome'): #Genome
            if SeqType!='NA':
                raise Exception, "Conflicting SeqTypes were specified"
            SeqType = 'genome'
        elif o in ('-c','--custom'): #Custom
            if SeqType!='NA':
                raise Exception, "Conflicting SeqTypes were specified"
            SeqType = 'Custom Capture'


        else:
            assert False, "Unhandled argument present"

def getIDs(SampleID,SeqType,sequenceDB,capturekit):
    '''Checks SequenceDB then GAFDB to determine which database if any the sample resides in.'''

    #DBID check
    sql = ("SELECT DISTINCT p.DBID "
        "FROM prepT p "
        "JOIN SeqType st on p.DBID=st.DBID "
        "WHERE CHGVID='{0}' and st.SeqType='{1}'"
        ).format(SampleID,SeqType)
    sequenceDB.execute(sql)
    DBID = sequenceDB.fetchall()
    if len(DBID) != 1:
        print SampleID,SeqType,DBID

        raise Exception, "Incorrect number of DBID's found for Sample %s" % SampleID

    if capturekit == '':
        sql = ("SELECT p.prepID "
                    "FROM SeqType st "
                    "join prepT p on st.prepid=p.prepid "
                    "where p.DBID={} and "
                    "SeqType='{}' and "
                    "p.CHGVID='{}' AND "
                    "failedPrep!=1"
                    ).format(DBID[0][0],SeqType,SampleID)
        #Grabs most recent prepID
        #sequenceDB.execute("SELECT prepID FROM SeqType where DBID=%s and SeqType=%s order by prepID desc limit 1", (DBID[0][0],SeqType))

    else:
        if ignore == False:
            sql = ("SELECT p.prepID "
                        "FROM prepT p "
                        "JOIN SeqType s ON p.prepID=s.prepID "
                        "where p.DBID={} and "
                        "s.SeqType='{}' and "
                        "p.exomeKit='{}' and "
                        "p.CHGVID='{}' AND "
                        "failedPrep!=1").format(DBID[0][0],SeqType,capturekit,SampleID)
        else:
            #Grabs most recent prepID
            sql = ("SELECT max(p.prepID) "
                        "FROM prepT p "
                        "JOIN SeqType s ON p.prepID=s.prepID "
                        "where p.DBID={} and "
                        "s.SeqType='{}' and "
                        "p.CHGVID='{}' AND "
                        "p.exomeKit='{}'").format(DBID[0][0],SeqType,capturekit,SampleID)

    sequenceDB.execute(sql)
    prepID = sequenceDB.fetchall()
    if len(prepID) != 1:
        print "Incorrect number of prepID's found for Sample %s." % SampleID
        print "prepIDs SeqType, CaptureKit"
        print "="*30
        print SampleID,prepID,SeqType,capturekit
        print "="*30

        multiplePrepOK = raw_input("Is this ok? (y/n)").lower()

        if multiplePrepOK == 'n':
            raise Exception, "Incorrect number of prepID's found for Sample %s.  Try using -k or --ignoreMultiPrepIDs option." % SampleID
        elif multiplePrepOK == 'y':
            return str(max(prepID)[0])
        else:
            raise Exception, "Incorrect input"

    else:
        #print sequenceDB_sample
        #raw_input('')
        return str(prepID[0][0])

def getSeqType():
    tSeqType = raw_input('E(x)ome, (g)enome, (r)NASeq or (c)ustom Capture sample? ')
    if tSeqType[0].lower() == 'g':
        SeqType = 'genome'
    elif tSeqType[0].lower() == 'x':
        SeqType = 'exome'
    elif tSeqType[0].lower() == 'r':
        SeqType = 'RNASeq'
    elif tSeqType[0].lower() == 'c':
        SeqType = 'Custom Capture'
    else:
        raise Exception, "Sample type not specified!"

    return SeqType

def checkPoolingRelease(IDs,failedSamples):
    poolDict= {}
    for sample in failedSamples:
        poolDict[sample[1][0]] = []
    for sample in failedSamples:
        poolDict[sample[1][0]].append(sample[0])

    #print len(IDs)

    for pool in poolDict.keys():
        if len(poolDict[pool]) > 1:
            for sample in IDs.keys():
                #print IDs[sample][2][0],pool,sample
                if IDs[sample][2][0] == pool:
                    print "Sample {0} was not released due to pool {1} having two samples with low sequencing yields".format(sample,pool)
                    IDs.pop(sample)
    #print len(IDs)
    return IDs

def getDBID(sequenceDB,prepID):
    #print prepID
    sequenceDB.execute("select DBID from prepT where prepID=%s", prepID)
    DBID = sequenceDB.fetchone()
    #print DBID,prepID
    DBID = DBID[0]
    return DBID

def setup_logging():
	logging.basicConfig(level=logging.INFO,format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',filename='/nfs/sva01/Summaries/RELEASE_LOG.log')
	logger = logging.getLogger(__name__)

def getTestSequenceDB():
    sequencedb=Connect(host='192.168.1.59',db='jb371',user='jb371',passwd="jb371_password")
    sequence_curs=sequencedb.cursor()
    sequence_curs.execute('USE jb371')
    sequence_curs.execute('BEGIN;')
    return sequence_curs


def main():
    opts(sys.argv[2:])
    setup_logging()
    logger = logging.getLogger('main')
    if test == False:
        sequenceDB = getSequenceDB()
    else:
        sequenceDB = getTestSequenceDB()

    path = Samples()
    global SeqType

    if SeqType == 'NA':
        SeqType = getSeqType()

    sampleList = open(path)
    global release_list
    release_list = []
    logger.info('Starting Release Script')
    logger.debug('Initializing Parameters: SeqType:%s, Parameters:%s',(SeqType,sys.argv[1:]))

    try:
        IDs = {}
        failedSamples = []
        for SampleID in sampleList.readlines():

            SampleID = SampleID.strip()
            print 'getting IDs for {}'.format(SampleID)
            prepID = getIDs(SampleID,SeqType,sequenceDB,capturekit)
            logger.debug('Releasing sample %s %s...' % (SampleID,prepID))
            DBID = getDBID(sequenceDB,prepID)
            poolID = getPoolID(sequenceDB,prepID)

            #print prepID,SampleID
            IDs[SampleID] = (DBID,long(prepID),poolID)
            logger.debug("Sample %s's prepID and DBID: %s, %s" % (SampleID,prepID,DBID))
            #print SampleID,prepID,DBID
            if old == False:
                status,completeStatus,failYieldSample,poolID,seqYield = check_sequenceDB(sequenceDB,SampleID,prepID,SeqType)
                if failYieldSample:
                    failedSamples.append((failYieldSample,poolID,seqYield))
                #print completeStatus
                if completeStatus != 'External Data Submitted' and skipStorageCheck == False:
                    pass
                    #check_Storage(sequenceDB,prepID,SeqType)
                check_Lane(sequenceDB,prepID,SampleID)
                checkS2R(sequenceDB,SampleID,prepID)
            #fixReleaseStatusT(SequenceDB,SampleID,SeqType,prepID)

            if sendemail == True:
                getPL(sequenceDB,SampleID,release_list,prepID)

        #print SeqType

        #remove failed Samples from ID list
        for failedSample in failedSamples:
            IDs.pop(failedSample[0],None)

        for passedSamples in IDs.keys():
            prepID = IDs[passedSamples][1]
            DBID = IDs[passedSamples][0]

            updateDB(sequenceDB,prepID,passedSamples,SeqType,DBID)

        #check it multiple samples from the same pool need a re-prep
        IDs = checkPoolingRelease(IDs,failedSamples)

        release_email(sequenceDB,SeqType,IDs,failedSamples)
        if sendemail == True:
            raw_input('Please hit Enter to continue')
            convert_release_list(sequenceDB,release_list)
        sampleList.close()
        sequenceDB.execute('COMMIT')


        #Release checking
        for samps in IDs:
            sequenceDB.execute("SELECT status from statusT WHERE prepID=%s AND status='External Data Submitted'", IDs[samps][1])
            externalCheck = sequenceDB.fetchone()
            if externalCheck:
                pass
            else:
                checkUpdateDB(sequenceDB,IDs[samps][1],samps,SeqType,IDs[samps][0])

        sequenceDB.close()
        logger.info('Release updates completed')

    except Exception,e:
        traceback.print_exc()
        sampleList.close()
        sequenceDB.execute('ROLLBACK')
        sequenceDB.close()
        logger.info('Release updates canceled/failed')


main()
