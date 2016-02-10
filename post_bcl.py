#!/usr/bin/python
# bcl_mysql.py
# Joshua Bridgers
# 04/26/2012
# jb3816@cumc.columbia.edu
#
# Updates gafdb after bcl completion

import collections
import getopt
import logging
import os
import re
import sys
import traceback
from commands import getoutput
from datetime import datetime
from CHGV_mysql import getSequenceDB
from CHGV_mysql import setup_logging
from CHGV_mysql import getUserID

def UpdateFC(sequenceDB,FCID,Unaligned):
    '''Updates Date Align, Actual Lane Fraction for FCID in the GAFdb'''

    logger = logging.getLogger('UpdateFC')
    CasavaVer = getoutput(('grep bcl2fastq -i {0}/Basecall_Stats_{1}/Demultiplex_Stats.htm  | ' 
        'cut -d- -f2 | cut -d\< -f1').format(Unaligned,FCID))
    sequenceDB.execute("SELECT Sum(LnYield) FROM Lane l join Flowcell f on l.fcid=f.fcid where FCillumID='%s'" % FCID)
    fcYield = sequenceDB.fetchone()
    sql = ("UPDATE Flowcell f "
        "JOIN Lane l ON f.FCID=l.FCID "
        "SET fcYield={0},CasavaVer='{1}',DateStor=now() "
        "WHERE f.FCillumID='{2}'"
        ).format(fcYield[0],CasavaVer,FCID)

    if verbose == True:
        print sql
    logger.info(sql)
    sequenceDB.execute(sql)

def UpdateSampleLane(sequenceDB,Unaligned,FCID):
    '''Gets Actual Lane Fraction FROM Demultiplex_Stats.htm, Updates sample status'''

    logger = logging.getLogger('UpdateSampleLane')
    os.system('cat '+Unaligned+'/Basecall_Stats_'+FCID+'/Demultiplex_Stats.htm | /usr/bin/w3m -dump -T text/html | egrep "^[12345678] " > %s/Demultiplex_Stats.txt' % Unaligned)
    os.system('chmod 775 %s/Demultiplex_Stats.txt' % Unaligned)
    Demulti = open('%s/Demultiplex_Stats.txt' % Unaligned,'r')
    for D in Demulti.readlines():
        info = D.split()
        LaneNum = info[0]
        SampleID = info[1]
        if len(info) == 15:
            LnFractionAct = info[10]
            LnYield = re.sub('\D','',info[7])
        else:
            LnFractionAct = '0'
            LnYield = '0'

        if info[3] == 'Undetermined' and float(LnFractionAct) > 7:
            print 'Percent of Undetermined Indexes is '+LnFractionAct+' for Lane '+LaneNum+'!'
            logger.info('Percent of Undetermined Indexes is '+LnFractionAct+' for Lane '+LaneNum+'!')
        elif info[1][0:4] == 'lane':
            pass
        else:
            sql = ("UPDATE Lane l "
                "JOIN prepT pt ON l.prepID=pt.prepID "
                "JOIN Flowcell f ON l.FCID=f.FCID "
                "SET l.LnYield={0}, l.LnFractionAct={1} "
                "WHERE pt.chgvid='{2}' AND f.FCillumID='{3}' AND l.LaneNum={4}"
                ).format(LnYield,LnFractionAct,SampleID,FCID,LaneNum)

            if verbose == True:
                print sql
            logger.info(sql)
            sequenceDB.execute(sql)
        Demulti.close()

def checkStatus(sequenceDB,FCID):
    """Check Status of all samples on the flowcell before the status update"""
    logger = logging.getLogger('checkStatus')

    sequenceDB.execute('SELECT DISTINCT prepID FROM Lane l join Flowcell f on l.FCID=f.FCID WHERE f.FCillumID=%s', FCID)
    samples = sequenceDB.fetchall()
    except_swt = 0
    for prepID in samples:
        sequenceDB.execute('SELECT CHGVID,status FROM statusT WHERE prepID=%s ORDER BY status_time DESC LIMIT 1', (prepID))
        status = sequenceDB.fetchone()
        sequenceDB.execute('SELECT SeqType from SeqType where prepID=%s', prepID)
        SeqType = sequenceDB.fetchone()
        if status[1] != 'BCL' and SeqType[0] != 'Genome':
            except_swt = 1
            logger.warn("%s does not have the correct status!  Current status: '%s'" % (status[0],status[1]))
            print "%s does not have the correct status!  Current status: '%s'" % (status[0],status[1])

    if except_swt == 1:
        logger.warn("Samples in the %s does not have the correct status!" % (FCID))
        raise Exception, 'Samples in the %s does not have the correct status!' % (FCID)

def UpdateSample(sequenceDB,FCID):
    """gafdb Sample Update"""
    logger = logging.getLogger('UpdateSample')
    if noStatusCheck == False:
        checkStatus(sequenceDB,FCID)
    userID = getUserID()

    #Status update for entire flowcell
    if noStatusUpdate == False:

        if verbose == True:
            print "INSERT INTO statusT (CHGVID,status_time,status,DBID,prepID,userID) SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'Storage',pt.DBID,pt.prepID,'%s' FROM Flowcell f join Lane l on l.FCID=f.FCID join prepT pt on pt.prepID=l.prepID where FCillumid='%s'" % (userID,FCID)
        logger.info("INSERT INTO statusT (CHGVID,status_time,status,DBID,prepID,userID) SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'Storage',pt.DBID,pt.prepID,'%s' FROM Flowcell f join Lane l on l.FCID=f.FCID join prepT pt on pt.prepID=l.prepID where FCillumid='%s'" % (userID,FCID))
        sequenceDB.execute("INSERT INTO statusT (CHGVID,status_time,status,DBID,prepID,userID) SELECT DISTINCT(pt.CHGVID),unix_timestamp(),'Storage',pt.DBID,pt.prepID,%s FROM Flowcell f join Lane l on l.FCID=f.FCID join prepT pt on pt.prepID=l.prepID where FCillumid=%s", (userID,FCID))

def usage():
    print '-i, --input\t\tParameter is the location of the Unaligned folder'
    print '-h, --help\t\tShows this help message and exit'
    print '-v, --verbose\t\tPrints out MYSQL injection commands'
    print '--noStatusCheck\t\tDoes not perform a status check of all samples in the flowcell'
    print '--noStatusUpdate\tDoes not update the status of all samples in the flowcell'
    sys.exit(2)

def getTotalLanes(FCID):
	if FCID[0] == 'H':
		return 2
	else:
		return 8

def checkLaneFractions(sequenceDB,FCID,Machine,Unaligned):
    logger = logging.getLogger('checkLaneFractions')
    Demulti = open('%s/Demultiplex_Stats.txt' % Unaligned,'r')
    email = open('%s/LnFractionEmail.txt' % Unaligned,'wb')
    lane = 1
    email_switch = 0
    demulti_d = collections.defaultdict(list)

    email.write("<style>\n")
    email.write("th, td {padding: 5px;}\n")
    email.write("th, td {text-align: center;}\n")
    email.write("th ding: 5px;}\n")

    email.write("</style>\n")



    email.write('<table border="1" style="border:1px solid black;border-collapse:collapse;width:95%">\n')

    email.write("<tr>\n")
    email.write("<th>SampleID</th>\n")
    email.write("<th>SeqType</th>\n")
    email.write("<th>LnFrac_pM</th>\n")
    email.write("<th>LnFracAct</th>\n")
    email.write("<th>ClustDen</th>\n")
    email.write("<th>Pool</th>\n")
    email.write("<th>Kapa</th>\n")
    #email.write("<th>Machine</th>\n")
    email.write("<th>cbot </th>\n")
    email.write("</tr>\n")

    #iterate over Demultiplex.txt
    for d in Demulti.readlines():
        d_info = d.split()
        SampleID = d_info[1]
        Lane = d_info[0]
        if d_info[7] != '0':
            lnFracAct = float(d_info[10])
        else:
            lnFracAct = 0
        if d_info[4] == 'unmatched':
            lnFrac = 0
            SeqType = 'N/A'
        else:
            lnFrac = float(d_info[4].split('_')[0])
            sequenceDB.execute("SELECT SeqType FROM SeqType s where prepID=(SELECT Distinct l.prepID FROM Lane l join Flowcell f on l.FCID=f.FCID join SampleT s on l.DBID=s.DBID where CHGVID=%s and FCillumID=%s and LaneNum=%s)", (SampleID,FCID,Lane))
            SeqType = sequenceDB.fetchone()
            #print SampleID,FCID,Lane,SeqType
            #print SampleID,FCID,SeqType
            SeqType = SeqType[0]
        if lnFracAct != 0:
            lnFracDiff = lnFrac/lnFracAct
        else:
            lnFracDiff = lnFrac
        #create hash of Demultiplex.txt
        demulti_d[Lane].append((SampleID,lnFracAct,lnFrac,lnFracDiff,SeqType,Lane,d_info[4]))
        Demulti.close()
    TotalLanes = getTotalLanes(FCID)

    for num in range(1,TotalLanes+1):
        lanes = demulti_d[str(num)]
        lane_switch = 0
        
        highlightRowNumber = []
        for samp in lanes:
            #print samp
            ClusterDensity = getClusterDensity(sequenceDB,FCID,lanes[0][5])
            #print ClusterDensity,ClusterDensity == None,ClusterDensity == ''
            """
            if ClusterDensity == '' or ClusterDensity == None:
                email_switch = 1
                lane_switch = 1
                highlightRowNumber = 999
            """
            if samp[0][0:4] == 'lane':
                if float(samp[1]) > 5:
                    email_switch = 1
                    lane_switch = 1
                    highlightRowNumber.append(10)

            elif (float(samp[3])*100 > 1.25 or float(samp[3])*100 < 0.80):
                email_switch = 1
                lane_switch = 1
                highlightRowNumber.append(3)

            elif float(ClusterDensity) > 1375 or float(ClusterDensity) < 1000:
                email_switch = 1
                lane_switch = 1
                highlightRowNumber.append(4)

            #print samp,lane_switch,samp[0][0:4] == 'lane' and samp[1] > 5 
        if lane_switch == 1:
            EmailLane(sequenceDB,FCID,lanes,email,highlightRowNumber)

    email.close()
    send_email(email_switch,FCID,Machine,Unaligned)

def send_email(email_switch,FCID,Machine,Unaligned):
    logger = logging.getLogger('send_email')

    if email_switch ==1 and os.path.isfile('%s/EmailSent.txt' % Unaligned) == False:
        address = "igm-hts@columbia.edu"
        #address = 'jb3816@cumc.columbia.edu'
        emailProgramLocation = '/nfs/goldstein/software/mutt-1.5.23/bin/mutt '
        emailCmd = emailProgramLocation + '-e "set content_type=text/html" '
        emailCmd += '-s \"Problem with Lane Fractions for flowcell %s %s\" ' % (FCID,Machine)
        emailCmd += address
        emailCmd += " < %s/LnFractionEmail.txt" % Unaligned

        logger.info(emailCmd)
        print emailCmd
        os.system(emailCmd)
        os.system("touch %s/EmailSent.txt" % Unaligned)

def getClusterDensity(sequenceDB,FCID,LaneNum):
	sequenceDB.execute("SELECT DISTINCT ClustDen FROM Lane l JOIN Flowcell f on l.FCID=f.FCID WHERE FCillumID=%s and LaneNum=%s", (FCID,LaneNum))
	ClustDen = sequenceDB.fetchone()
	return ClustDen[0]

def getCbot(sequenceDB,FCID):
    sequenceDB.execute("SELECT cbot from Flowcell where FCillumID=%s", FCID)
    cbot = sequenceDB.fetchone()
    return cbot[0]

def getKapaPicoDBID(sequenceDB,FCID,CHGVID):
    sequenceDB.execute("SELECT kapa_conc,picomoles,s2r.dbid FROM Lane l join Flowcell f\
            ON f.FCID=l.FCID join prepT pt ON l.prepID=pt.prepID join\
            samplesTOrun s2r ON s2r.seqID=l.seqID join SampleT s ON\
            s.DBID=pt.DBID LEFT JOIN poolMembers pm ON (pm.DBID=pt.DBID AND\
                pm.poolID=l.poolID)  where FCIllumID=%s and\
                pt.CHGVID=%s", (FCID,CHGVID))
    kapaPicoDBID = sequenceDB.fetchone()
    #print kapaPicoDBID
    return kapaPicoDBID

def getPoolName(sequenceDB,DBID):
    sequenceDB.execute("select CHGVID from SampleT where dbid=%s", DBID)
    poolName = sequenceDB.fetchone()
    return poolName[0]

def EmailLane(sequenceDB,FCID,lanes,email,highlightRowNumber):
    logger = logging.getLogger('EmailLane')
    logger.info('Lane fraction problem with lane: %s' % lanes[0][5])
    email.write('<tr><td colspan="8">Lane %s</td></tr>\n' %lanes[0][5]+'\n')
    ClusterDensity = getClusterDensity(sequenceDB,FCID,lanes[0][5])
    cbot = getCbot(sequenceDB,FCID)

    for samp in lanes:
        if samp[0][0:4] == 'lane':
            pass
        else:
            kapaPicoDBID = getKapaPicoDBID(sequenceDB,FCID,samp[0])
           
            poolName = getPoolName(sequenceDB,kapaPicoDBID[2])
            logging.info('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %
                    (samp[0],samp[4],samp[6],samp[1],ClusterDensity,poolName,kapaPicoDBID[0],cbot))


            email.write("<tr>\n")
            highlightRowCount = 0
            for column in (samp[0],samp[4],samp[6],samp[1],ClusterDensity,poolName,kapaPicoDBID[0],cbot):
                #print highlightRowNumber,highlightRowCount
                if highlightRowCount not in  highlightRowNumber:
                    email.write("<td>%s</td>\n" % column)
                else:
                    email.write('<td bgcolor="#FFFF00">%s</td>\n' % column)
                highlightRowCount += 1

                 
            email.write("</tr>\n")


    for samp in lanes:
        if samp[0][0:4] == 'lane':
            logging.info("Lane %s's unmatched reads percent: %s" % (lanes[0][5],samp[1]))
            if 10 in highlightRowNumber:
                email.write('<tr><td colspan="8" align="center" bgcolor="#FFFF00">Lane %s\'s unmatched reads percent: %s</td></tr>\n' %
                    (lanes[0][5],samp[1]))

            else:
                email.write('<tr><td colspan="8" align="center">Lane %s\'s unmatched reads percent: %s</td></tr>\n' %
                    (lanes[0][5],samp[1]))
            email.write('<tr><td colspan="8" align="center">&nbsp</td></tr>\n')
          
def opts(argv):
    global verbose
    verbose = False
    global runFolder
    runFolder = os.getcwd()
    global seqsata_drive
    sata_loc = ''
    global noStatusUpdate
    noStatusUpdate = False
    global noStatusCheck
    noStatusCheck = False
    global BCLDrive
    BCLDrive = ''


    try:
        opts,args = getopt.getopt(argv, "ab:hvi:s:", ['bcl=','input=','help','seqsata=','verbose','noStatusUpdate','noStatusCheck'])
    except getopt.GetoptError, err:
        print str(err)
        usage()

    for o,a in opts:
        if o in ('-v','--verbose'):
            verbose = True
        elif o in ('-h','--help'):
            usage()
        elif o in ('-i','--input'):
            runFolder =  a
        elif o in ('--noStatusUpdate'):
            noStatusUpdate = True
        elif o in ('-b','--bcl'):
            BCLDrive = a
        elif o in ('--noStatusCheck'):
            noStatusCheck = True
        elif o in ('-s','--seqsata'):
            seqsata_drive = a
        else:
            assert False, "Unhandled argument present"


def main():
    #accessing mysql gaf database
    sequenceDB = getSequenceDB()
    opts(sys.argv[1:])
    info = runFolder.split('_')
    Date = info[0]
    FCID = info[3]
    Machine = info[1]

    Unaligned = '{0}/{1}_{2}_{3}_Unaligned'.format(BCLDrive,Date,Machine,FCID)
    setup_logging(Machine,FCID,seqsata_drive)
    logger = logging.getLogger('main')
    logger.info('Starting post BCL script')


    try:
        UpdateSampleLane(sequenceDB,Unaligned,FCID)
        UpdateFC(sequenceDB,FCID,Unaligned)
        checkLaneFractions(sequenceDB,FCID,Machine,Unaligned)
        UpdateSample(sequenceDB,FCID)
       
        sequenceDB.execute('COMMIT;')
        sequenceDB.close()
        logger.info('Done')
        print 'Done'
    except:
        traceback.print_exc()
        sequenceDB.execute('ROLLBACK;')
        sequenceDB.close()
        logger.info('Post BCL Failure')
        print 'Post BCL Failure'
        sys.exit(1)

main()
