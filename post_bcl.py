#!/usr/bin/python
# post_bcl.py
# Joshua Bridgers
# 04/26/2012
# jb3816@cumc.columbia.edu
#
# Updates SequenceDB after bcl completion

import collections
import getopt
import logging
import os
import re
import subprocess
import sys
import traceback
from datetime import datetime
from CHGV_mysql_3_6 import getSequenceDB
from CHGV_mysql_3_6 import getTestSequenceDB
from CHGV_mysql_3_6 import setup_logging
from CHGV_mysql_3_6 import getUserID
from xml.etree.ElementTree import fromstring

def UpdateFC(sequenceDB,FCID,Unaligned,sampleInfo):
    '''Updates Flowcell total Yield, Casava version for FCIllumID'''

    logger = logging.getLogger('UpdateFC')

    fcYield = 0
    for i in sampleInfo.keys():
        LaneNum,SampleID,adapter = i.split('_')
        LnFractionAct = sampleInfo[i][0]
        LnYield = int(sampleInfo[i][1].replace(',',''))
        fcYield += LnYield

    casavaCmd = ['grep','bcl2fastq','-m1','{}/nohup.sge'.format(Unaligned)]
    proc = subprocess.Popen(casavaCmd, stdout=subprocess.PIPE)
    CasavaVer = proc.stdout.read().decode().split()[1].strip()
    sql = ("UPDATE Flowcell f "
        "JOIN Lane l ON f.FCID=l.FCID "
        "SET fcYield={0},CasavaVer='{1}' "
        "WHERE f.FCillumID='{2}'"
        ).format(fcYield,CasavaVer,FCID)

    if verbose == True:
        print(sql)
    logger.info(sql)
    sequenceDB.execute(sql)

def getHTML(Unaligned,FCID):
    '''Gets the HTML code specifically for the second table'''
    logger = logging.getLogger('getHTML')
    htmlLoc = ('{}/Reports/html/{}/all/all/all/laneBarcode.html'
            ).format(Unaligned,FCID)

    with open(htmlLoc) as f:
        laneSummaryFlag=0
        laneBarcodeHTML = ''
        for line in f:
            line = line.replace('<br>','')
            if line == '<h2>Top Unknown Barcodes</h2>\n':
                laneSummaryFlag = 0
            if laneSummaryFlag:
                laneBarcodeHTML += line
            if line == '<h2>Lane Summary</h2>\n':
                laneSummaryFlag = 1
    return laneBarcodeHTML

def UpdateSampleLane(sequenceDB,Unaligned,FCID):
    '''Gets Actual Lane Fraction FROM laneBarcode.htm, Updates sample status'''

    logger = logging.getLogger('UpdateSampleLane')
    laneBarcodeHTML = getHTML(Unaligned,FCID)

    #use xml.etree for HTML table parsing
    tree = fromstring(laneBarcodeHTML)
    rows = tree.findall("tr")
    headrow = rows[0]
    datarows = rows[1:]

    #convert xml elements into dict with the actual lane fraction and lane yield
    sampleInfo = {}
    for sampleNum,i in enumerate(datarows):
        key = ''
        for num, h in enumerate(headrow):
            if num == 0:
                key = datarows[sampleNum][num].text
            elif num == 2:
                key += '_' + datarows[sampleNum][num].text
            elif num == 3:
                key += '_' + datarows[sampleNum][num].text
                sampleInfo[key] = []
            elif num == 5 or num == 8:
                sampleInfo[key].append(datarows[sampleNum][num].text)

    for i in sampleInfo.keys():
        #print(sampleInfo[i],i)
        LaneNum,SampleID,adapter = i.split('_')
        LnFractionAct = sampleInfo[i][0]
        LnYield = sampleInfo[i][1].replace(',','')

        if SampleID == 'Undetermined' and float(LnFractionAct) > 5:
            print('Percent of Undetermined Indexes is '+LnFractionAct+' for Lane '+LaneNum+'!')
            logger.info('Percent of Undetermined Indexes is '+LnFractionAct+' for Lane '+LaneNum+'!')
        else:
            sql = ("UPDATE Lane l "
                "JOIN prepT pt ON l.prepID=pt.prepID "
                "JOIN Flowcell f ON l.FCID=f.FCID "
                "SET l.LnYield={0}, l.LnFractionAct={1} "
                "WHERE pt.chgvid='{2}' AND f.FCillumID='{3}' AND l.LaneNum={4}"
                ).format(LnYield,LnFractionAct,SampleID,FCID,LaneNum)

            if verbose == True:
                print(sql)
            logger.info(sql)
            sequenceDB.execute(sql)
    return sampleInfo

def checkStatus(sequenceDB,FCID):
    """Check Status of all samples on the flowcell before the status update"""
    logger = logging.getLogger('checkStatus')

    sequenceDB.execute('SELECT DISTINCT prepID FROM Lane l JOIN Flowcell f on l.FCID=f.FCID WHERE f.FCillumID=%s', FCID)
    samples = sequenceDB.fetchall()
    except_swt = 0
    for prepID in samples:
        sequenceDB.execute('SELECT CHGVID,status FROM statusT WHERE prepID=%s ORDER BY status_time DESC LIMIT 1',(prepID['prepID']))
        status = sequenceDB.fetchone()
        sequenceDB.execute('SELECT SeqType from SeqType WHERE prepID=%s', prepID['prepID'])
        SeqType = sequenceDB.fetchone()
        # Sometimes multiple flowcells of Genomes finish bcl2fastq near the same time so its status =='storage'
        if status['status'] != 'BCL' and SeqType['SeqType'] != 'Genome':
            except_swt = 1
            logger.warn("{} does not have the correct status!  Current status: '{}'".format(status[0],status[1]))
            print("{} does not have the correct status!  Current status: '{}'".format(status[0],status[1]))
        else:
            if status['status'] != 'BCL':
                except_swt = 1
                logger.warn("%s does not have the correct status!  Current status: '%s'" % (status[0],status[1]))
                print("%s does not have the correct status!  Current status: '%s'" % (status[0],status[1]))

    if except_swt == 1:
        logger.warn("Samples in the {} does not have the correct status!".format(FCID))
        raise Exception('Samples in the {} does not have the correct status!'.format(FCID))


def usage():
    print('-i, --input\t\tParameter is the location of the Unaligned folder')
    print('-h, --help\t\tShows this help message and exit')
    print('-v, --verbose\t\tPrints out MYSQL injection commands')
    print('--noStatusCheck\t\tDoes not perform a status check of all samples in the flowcell')
    print('--noStatusUpdate\tDoes not update the status of all samples in the flowcell')
    sys.exit(2)

def getTotalLanes(sequenceDB,FCID):
    query = ("SELECT max(laneNum) "
                "from Lane l "
                "join Flowcell f on l.fcid=f.fcid "
                "where FCILLUMID='{}'").format(FCID)
    sequenceDB.execute(query)
    totalLanes = sequenceDB.fetchone()['max(laneNum)']
    return totalLanes

def checkLaneFractions(sequenceDB,FCID,Machine,Unaligned,sampleInfo):
    logger = logging.getLogger('checkLaneFractions')
    email = open('%s/LnFractionEmail.txt' % Unaligned,'w')
    lane = 1
    emailSwitch = 0
    demulti_d = collections.defaultdict(list)

    #email header
    email.write("<style>\n")
    email.write("th, td {padding: 5px;}\n")
    email.write("th, td {text-align: center;}\n")
    email.write("th ding: 5px;}\n")
    email.write("</style>\n")
    email.write('<table border="1" style="border:1px solid black;border-collapse:collapse;width:95%">\n')
    email.write("<tr>\n")
    email.write("<th>SampleID</th>\n")
    email.write("<th>SeqType</th>\n")
    email.write("<th>LnFrac</th>\n")
    email.write("<th>LnFracAct</th>\n")
    email.write("<th>ClustDen</th>\n")
    email.write("<th>Pool</th>\n")
    email.write("<th>Kapa</th>\n")
    email.write("</tr>\n")

    for sample in sampleInfo.keys():
        LaneNum,SampleID,adapter = sample.split('_')
        LnFractionAct = float(sampleInfo[sample][0])
        LnYield = sampleInfo[sample][1].replace(',','')

        if SampleID == 'Undetermined':
            LnFrac = 0
            SeqType = 'N/A'
        else:
            LnFracQuery = ("SELECT LnFraction "
                            "FROM Lane l "
                            "JOIN Flowcell f on l.FCID=f.FCID "
                            "JOIN prepT p on p.prepid=l.prepid "
                            "WHERE FCillumID='{}' "
                            "AND CHGVID='{}' "
                            "AND LaneNum={}"
                            ).format(FCID,SampleID,LaneNum)
            sequenceDB.execute(LnFracQuery)
            LnFrac = sequenceDB.fetchone()['LnFraction']
            query = ("SELECT SeqType "
                    "FROM SeqType s "
                    "WHERE prepID=(SELECT Distinct l.prepID "
                        "FROM Lane l "
                        "JOIN Flowcell f on l.FCID=f.FCID "
                        "JOIN prepT p on l.DBID=p.DBID "
                        "WHERE CHGVID='{}' "
                        "AND FCillumID='{}' "
                        "AND LaneNum={})"
                    ).format(SampleID,FCID,LaneNum)
            sequenceDB.execute(query)
            SeqType = sequenceDB.fetchone()
            #print(SampleID,FCID,Lane,SeqType)
            #print(SampleID,FCID,SeqType)
            SeqType = SeqType['SeqType']
        if LnFractionAct != 0:
            LnFracDiff = LnFrac/LnFractionAct
        else:
            LnFracDiff = 0

        demulti_d[LaneNum].append((SampleID,LnFractionAct,LnFrac,LnFracDiff,SeqType,LaneNum))

    TotalLanes = getTotalLanes(sequenceDB,FCID)
    for num in range(1,TotalLanes+1):
        lanes = demulti_d[str(num)]
        lane_switch = 0
        highlightRowNumber = []

        for samp in lanes:
            #print(samp)
            ClusterDensity = getClusterDensity(sequenceDB,FCID,lanes[0][5])
            #print(ClusterDensity,ClusterDensity == None,ClusterDensity == '')

            if samp[0] == 'Undetermined':
                if float(samp[1]) > 5:
                    emailSwitch = 1
                    lane_switch = 1
                    highlightRowNumber.append(10)

            elif (float(samp[3])*100 > 1.25 or float(samp[3])*100 < 0.80):
                emailSwitch = 1
                lane_switch = 1
                highlightRowNumber.append(3)

            elif float(ClusterDensity) > 1375 or float(ClusterDensity) < 1000:
                emailSwitch = 1
                lane_switch = 1
                highlightRowNumber.append(4)

            #print(samp,lane_switch,samp[0][0:4] == 'lane' and samp[1] > 5) 
        if lane_switch == 1:
            EmailLane(sequenceDB,FCID,lanes,email,highlightRowNumber)

    email.close()
    send_email(emailSwitch,FCID,Machine,Unaligned)

def send_email(emailSwitch,FCID,Machine,Unaligned):
    logger = logging.getLogger('send_email')
    if emailSwitch == 1 and os.path.isfile('%s/EmailSent.txt' % Unaligned) == False:
        address = "igm-hts@columbia.edu"
        #address = 'jb3816@cumc.columbia.edu'
        emailProgramLocation = '/nfs/goldstein/software/mutt-1.5.23/bin/mutt '
        emailCmd = emailProgramLocation + '-e "set content_type=text/html" '
        emailCmd += '-s \"Problem with Lane Fractions for flowcell %s %s\" ' % (FCID,Machine)
        emailCmd += address
        emailCmd += " < %s/LnFractionEmail.txt" % Unaligned

        logger.info(emailCmd)
        #print(emailCmd)
        os.system(emailCmd)
        os.system("touch %s/EmailSent.txt" % Unaligned)

def getClusterDensity(sequenceDB,FCID,LaneNum):
	sequenceDB.execute("SELECT DISTINCT ClustDen FROM Lane l JOIN Flowcell f on l.FCID=f.FCID WHERE FCillumID=%s AND LaneNum=%s", (FCID,LaneNum))
	ClustDen = sequenceDB.fetchone()
	return ClustDen['ClustDen']


def getKapaPicoDBID(sequenceDB,FCID,CHGVID):
    sequenceDB.execute("SELECT kapa_conc,picomoles,s2r.dbid FROM Lane l JOIN Flowcell f\
            ON f.FCID=l.FCID JOIN prepT pt ON l.prepID=pt.prepID JOIN\
            samplesTOrun s2r ON s2r.seqID=l.seqID JOIN SampleT s ON\
            s.DBID=pt.DBID LEFT JOIN poolMembers pm ON (pm.DBID=pt.DBID AND\
                pm.poolID=l.poolID)  WHERE FCIllumID=%s and\
                pt.CHGVID=%s", (FCID,CHGVID))
    kapaPicoDBID = sequenceDB.fetchone()
    #print(kapaPicoDBID)
    return kapaPicoDBID

def getPoolName(sequenceDB,DBID):
    sequenceDB.execute("select CHGVID from SampleT WHERE dbid=%s", DBID)
    poolName = sequenceDB.fetchone()
    return poolName['CHGVID']

def EmailLane(sequenceDB,FCID,lanes,email,highlightRowNumber):
    logger = logging.getLogger('EmailLane')
    logger.info('Lane fraction problem with lane: %s' % lanes[0][5])
    email.write('<tr><td colspan="7">Lane %s</td></tr>\n' %lanes[0][5]+'\n')
    ClusterDensity = getClusterDensity(sequenceDB,FCID,lanes[0][5])

    for samp in lanes:
        if samp[0] == 'Undetermined':
            pass
        else:
            kapaPicoDBID = getKapaPicoDBID(sequenceDB,FCID,samp[0])
            poolName = getPoolName(sequenceDB,kapaPicoDBID['dbid'])
            logging.info('%s\t%s\t%s\t%s\t%s\t%s\t%s' %
                    (samp[0],samp[4],samp[2],samp[1],ClusterDensity,poolName,kapaPicoDBID['kapa_conc']))

            email.write("<tr>\n")
            highlightRowCount = 0
            for column in (samp[0],samp[4],samp[2],samp[1],ClusterDensity,poolName,kapaPicoDBID['kapa_conc']):
                #print(highlightRowNumber,highlightRowCount)
                if highlightRowCount not in  highlightRowNumber:
                    email.write("<td>%s</td>\n" % column)
                else:
                    email.write('<td bgcolor="#FFFF00">%s</td>\n' % column)
                highlightRowCount += 1
            email.write("</tr>\n")

    #demulti_d[LaneNum].append((SampleID,LnFractionAct,LnFrac,LnFracDiff,SeqType,LaneNum))
    for samp in lanes:
        if samp[0] == 'Undetermined':
            logging.info("Lane %s's unmatched reads percent: %s" % (lanes[0][5],samp[1]))
            if 10 in highlightRowNumber:
                email.write('<tr><td colspan="7" align="center" bgcolor="#FFFF00">Lane %s\'s unmatched reads percent: %s</td></tr>\n' %
                    (lanes[0][5],samp[1]))
            else:
                email.write('<tr><td colspan="7" align="center">Lane %s\'s unmatched reads percent: %s</td></tr>\n' %
                    (lanes[0][5],samp[1]))
            email.write('<tr><td colspan="7" align="center">&nbsp</td></tr>\n')

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
    global test
    test = False

    try:
        opts,args = getopt.getopt(argv, "ab:hvi:s:",
            ['bcl=','input=','help','seqsata=','verbose','noStatusUpdate','noStatusCheck','test'])
    except getopt.GetoptError:
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
        elif o in ('--test'):
            test = True
        else:
            assert False, "Unhandled argument present"

def main():
    opts(sys.argv[1:])
    #accessing mysql gaf database
    if test:
        sequenceDB = getTestSequenceDB()
    else:
        sequenceDB = getSequenceDB()
    #often the experiment field has extra underscores
    Date,Machine,RunNum,FCID,*rest = runFolder.split('_')

    Unaligned = '{0}/{1}_{2}_{3}_Unaligned'.format(BCLDrive,Date,Machine,FCID)
    setup_logging(Machine,FCID,seqsata_drive)
    logger = logging.getLogger('main')
    logger.info('Starting post BCL script')


    try:
        sampleInfo = UpdateSampleLane(sequenceDB,Unaligned,FCID)
        if noStatusCheck == False:
            checkStatus(sequenceDB,FCID)
        UpdateFC(sequenceDB,FCID,Unaligned,sampleInfo)
        checkLaneFractions(sequenceDB,FCID,Machine,Unaligned,sampleInfo)

        sequenceDB.execute('COMMIT;')
        sequenceDB.close()
        logger.info('Done')
        print('Done')
    except:
        traceback.print_exc()
        sequenceDB.execute('ROLLBACK;')
        sequenceDB.close()
        logger.info('Post BCL Failure')
        print('Post BCL Failure')
        sys.exit(1)

main()
