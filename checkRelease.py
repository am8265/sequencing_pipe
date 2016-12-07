#checkRelease.py

import sys
from MySQLdb import Connect
from CHGV_mysql import getSequenceDB

def checkUpdateDB(sequenceDB,prepID,SampleID,SeqType,DBID):
    if SeqType == 'Custom_Capture':
        SeqType == 'Custom Capture'
    sequenceDB.execute("SELECT CHGVID,prepID,idnum,ExomeSamPrepKit,ReadSummaryFileLoc,SeqType FROM seqdbClone WHERE prepID='%s'", int(prepID))
    samInfo = sequenceDB.fetchone()
    #print samInfo,prepID,SampleID,SeqType,DBID

    if samInfo:
        """
        print samInfo[0] != SampleID
        print str(int(samInfo[1])) != prepID
        print type(prepID)
        print samInfo[2] != DBID
        print samInfo[4] == ''
        print samInfo[5] == ''
        """
        if samInfo[0] != SampleID or samInfo[1] != prepID or samInfo[2] != DBID or samInfo[4] == '' or samInfo[5] == '':
            raise Exception, 'Sample information does not match or is missing!'
    else:
        raise Exception, 'Sample information was not found for prepID: %s' % prepID

    if SeqType == 'Custom Capture' or SeqType == 'Custom_Capture':
        if samInfo[4] == '' or samInfo[4] == 'na' or samInfo[4] == 'N/A' or samInfo[4] == 'n/a':
            raise Exception, 'ExomeSamPrepKit information is missing for sample: %s' % CHGVID

    #check lane table
    sequenceDB.execute("SELECT mnQscR1 from Lane WHERE prepID=%s and FailR1 is NULL and mnQscR1 is NULL", prepID)
    mnQscR1Missing = sequenceDB.fetchone()
    #print prepID,SampleID,SeqType,DBID,mnQscR1Missing
    if mnQscR1Missing == None:
        pass
    elif mnQscR1Missing[0] == None:
        pass
    elif mnQscR1Missing:
        raise Exception, '%s is missing mnQscR1 information for prepID: %s' % (SampleID,prepID)

def getIDs(sequenceDB,SampleID,SeqType):
    sequenceDB.execute("SELECT DISTINCT p.DBID FROM prepT p JOIN SeqType st on p.DBID=st.DBID WHERE CHGVID=%s and st.SeqType=%s", (SampleID,SeqType))
    DBID = sequenceDB.fetchall()
    if len(DBID) != 1:
        raise Exception, "Incorrect number of DBID's found for Sample %s" % SampleID
    #Grabs most recent prepID
    sequenceDB.execute("SELECT PrepID FROM SeqType where DBID=%s and SeqType=%s ORDER BY PREPID DESC LIMIT 1", (DBID[0][0],SeqType))
    prepID = sequenceDB.fetchone()

    return DBID[0][0],prepID[0]

def main():
    sequenceDB = getSequenceDB()
    sampleFile = open(sys.argv[1])
    SeqType = sys.argv[2]
    for s in sampleFile.readlines():
        s = s.strip()
        print "Checking %s..." % s 
        DBID,prepID = getIDs(sequenceDB,s,SeqType)
        checkUpdateDB(sequenceDB,prepID,s,SeqType,DBID)
    print 'Release Check complete'

if __name__ == '__main__':
    main()
