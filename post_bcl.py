#!/usr/bin/python
# post_bcl.py
# Joshua Bridgers
# 04/26/2012
# jb3816@cumc.columbia.edu
#
# Updates SequenceDB after bcl completion

import argparse
import collections
import getopt
import logging
import os
import re
import subprocess
import sys
import traceback
from datetime import datetime
from utilities import *
from xml.etree.ElementTree import fromstring

def UpdateFC(fcillumid,unaligned_dir,sample_info,verbose,database):
    '''Updates Flowcell total Yield, Casava version for FCIllumID'''

    logger = logging.getLogger(__name__)

    fcYield = 0
    for i in sample_info.keys():
        info = i.split('_')
        LaneNum = info[0]
        SampleID = '_'.join(info[1:-1])
        adapter = info[-1]
        LnFractionAct = sample_info[i][0]
        LnYield = int(sample_info[i][1].replace(',',''))
        fcYield += LnYield

    casavaCmd = ['grep','bcl2fastq','-m1','{}/nohup.sge'.format(unaligned_dir)]
    proc = subprocess.Popen(casavaCmd, stdout=subprocess.PIPE)
    CasavaVer = proc.stdout.read().decode().split()[1].strip()
    sql = ("UPDATE Flowcell f "
        "JOIN Lane l ON f.FCID=l.FCID "
        "SET fcYield={0},CasavaVer='{1}' "
        "WHERE f.FCillumID='{2}'"
        ).format(fcYield,CasavaVer,fcillumid)

    if verbose == True:
        print(sql)
    logger.info(sql)
    run_query(sql,database)

def getHTML(unaligned_dir,fcillumid,verbose):
    '''Gets the HTML code specifically for the second table'''
    logger = logging.getLogger(__name__)
    htmlLoc = ('{}/Reports/html/{}/all/all/all/laneBarcode.html'
            ).format(unaligned_dir,fcillumid)
    if verbose:
        print(htmlLoc)
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

def update_lane(unaligned_dir,fcillumid,verbose,database):
    '''Gets Actual Lane Fraction FROM laneBarcode.htm, Updates sample status'''

    logger = logging.getLogger(__name__)
    laneBarcodeHTML = getHTML(unaligned_dir,fcillumid,verbose)

    #use xml.etree for HTML table parsing
    tree = fromstring(laneBarcodeHTML)
    rows = tree.findall("tr")
    headrow = rows[0]
    datarows = rows[1:]

    #convert xml elements into dict with the actual lane fraction and lane yield
    sample_info = {}
    for sampleNum,i in enumerate(datarows):
        key = ''
        for num, h in enumerate(headrow):
            if num == 0:
                key = datarows[sampleNum][num].text
            elif num == 2:
                key += '_' + datarows[sampleNum][num].text
            elif num == 3:
                key += '_' + datarows[sampleNum][num].text
                sample_info[key] = []
            elif num == 5 or num == 8:
                sample_info[key].append(datarows[sampleNum][num].text)
    for i in sample_info.keys():
        info = i.split('_')
        LaneNum = info[0]
        SampleID = '_'.join(info[1:-1])
        adapter = info[-1]
        LnFractionAct = sample_info[i][0]
        LnYield = sample_info[i][1].replace(',','')

        if SampleID == 'Undetermined' and float(LnFractionAct) > 5:
            print('Percent of Undetermined Indexes is '+LnFractionAct+' for Lane '+LaneNum+'!')
            logger.info('Percent of Undetermined Indexes is '+LnFractionAct+' for Lane '+LaneNum+'!')
        else:
            sql = ("UPDATE Lane l "
                "JOIN prepT pt ON l.prepID=pt.prepID "
                "JOIN Flowcell f ON l.FCID=f.FCID "
                "SET l.LnYield={0}, l.LnFractionAct={1} "
                "WHERE pt.chgvid='{2}' AND f.FCillumID='{3}' AND l.LaneNum={4}"
                ).format(LnYield,LnFractionAct,SampleID,fcillumid,LaneNum)

            if verbose == True:
                print(sql)
            logger.info(sql)
            run_query(sql,database)
    return sample_info

def checkStatus(fcillumid,database):
    """Check Status of all samples on the flowcell before the status update"""
    logger = logging.getLogger(__name__)

    samples = run_query("SELECT DISTINCT prepID "
                        "FROM Lane l "
                        "JOIN Flowcell f on l.FCID=f.FCID "
                        "WHERE f.FCillumID='{}'".format(fcillumid),database)
    raise_exception_switch = 0
    for prepID in samples:
        status = run_query("SELECT CHGVID,status "
                           "FROM statusT "
                           "WHERE prepID={} "
                           "ORDER BY status_time "
                           "DESC LIMIT 1".format(prepID['prepID']),database)
        SeqType = run_query(GET_SEQTYPE_FROM_PREPID.format(prepid=prepID['prepID']),database)
        # Sometimes multiple flowcells of Genomes finish bcl2fastq near the same time so its status =='storage'
        if status[0]['status'] != 'BCL' and SeqType[0]['SAMPLE_TYPE'] != 'Genome':
            raise_exception_switch = 1
            logger.warn("{} does not have the correct status!  Current status: '{}'".format(status[0]['CHGVID'],status[0]['status']))
            print("{} does not have the correct status!  Current status: '{}'".format(status[0]['CHGVID'],status[0]['status']))
        else:
            if status[0]['status'] != 'BCL':
                raise_exception_switch = 1
                logger.warn("{} does not have the correct status!  Current status: '{}'".format(status[0]['CHGVID'],status[0]['status']))
                print("{} does not have the correct status!  Current status: '{}'".format(status[0]['CHGVID'],status[0]['status']))

    if raise_exception_switch == 1:
        logger.warn("Samples in the {} does not have the correct status!".format(fcillumid))
        raise Exception('Samples in the {} does not have the correct status!'.format(fcillumid))


def getTotalLanes(fcillumid,database):
    query = ("SELECT max(laneNum) as LANECOUNT "
                "FROM Lane l "
                "JOIN Flowcell f on l.fcid=f.fcid "
                "WHERE FCILLUMID='{}'").format(fcillumid)
    totalLanes = run_query(query,database)
    return totalLanes[0]['LANECOUNT']

def addHeader(email):
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

    return email

def checkLaneFractions(fcillumid,Machine,unaligned_dir,sample_info,database):
    logger = logging.getLogger(__name__)
    email = open('%s/LnFractionEmail.txt' % unaligned_dir,'w')
    lane = 1
    emailSwitch = 0
    demulti_d = collections.defaultdict(list)
    email = addHeader(email)

    chemVer = run_query(GET_FLOWCELL_CHEMVER.format(fcillumid=fcillumid),database)

    for sample in sample_info.keys():
        #print(sample)
        info = sample.split('_')
        LaneNum = info[0]
        SampleID = '_'.join(info[1:-1])
        adapter = info[-1]

        LnFractionAct = float(sample_info[sample][0])
        LnYield = sample_info[sample][1].replace(',','')

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
                            ).format(fcillumid,SampleID,LaneNum)
            LnFrac = run_query(LnFracQuery,database)
            LnFrac = LnFrac[0]['LnFraction']
            query = ("SELECT SAMPLE_TYPE "
                    "FROM prepT p "
                    "WHERE prepID=(SELECT Distinct l.prepID "
                        "FROM Lane l "
                        "JOIN Flowcell f on l.FCID=f.FCID "
                        "JOIN prepT p on l.DBID=p.DBID "
                        "WHERE CHGVID='{}' "
                        "AND FCillumID='{}' "
                        "AND LaneNum={})"
                    ).format(SampleID,fcillumid,LaneNum)
            SeqType = run_query(query,database)
            #print(SampleID,fcillumid,Lane,SeqType)
            #print(SampleID,fcillumid,SeqType)
            SeqType = SeqType[0]['SAMPLE_TYPE']
        if LnFractionAct != 0:
            LnFracDiff = LnFrac/LnFractionAct
        else:
            LnFracDiff = 0

        demulti_d[LaneNum].append((SampleID,LnFractionAct,LnFrac,LnFracDiff,SeqType,LaneNum))

    TotalLanes = getTotalLanes(fcillumid,database)
    for num in range(1,TotalLanes+1):
        lanes = demulti_d[str(num)]
        lane_switch = 0
        highlightRowNumber = []

        for samp in lanes:
            #print(samp)
            cluster_density = run_query(GET_CLUSTER_DENSITY_FOR_LANE.format(fcillumid=fcillumid,lanenum=lanes[0][5]),database)
            cluster_density = cluster_density[0]['CLUSTDEN']
            #HiSeqs
            if chemVer[0]['CHEMVER'][0] == 'H' or chemVer[0]['CHEMVER'] == 'v4':
                if samp[0] == 'undetermined':
                    if float(samp[1]) > 3:
                        emailSwitch = 1
                        lane_switch = 1
                        highlightRowNumber.append(10)
                elif (float(samp[3])*100 > 1.25 or float(samp[3])*100 < 0.80):
                    emailSwitch = 1
                    lane_switch = 1
                    highlightRowNumber.append(3)
                elif float(cluster_density) > 1375 or float(cluster_density) < 1000: #HiSeq ClustDen Check
                    emailSwitch = 1
                    lane_switch = 1
                    highlightRowNumber.append(4)
            #NovaSeqs
            elif chemVer[0]['CHEMVER'][0] == 'N':
                if samp[0] == 'undetermined':
                    if float(samp[1]) > 3:
                        emailSwitch = 1
                        lane_switch = 1
                        highlightRowNumber.append(10)
                elif (float(samp[3])*100 > 8 or float(samp[3])*100 < 0.5):
                    emailSwitch = 1
                    lane_switch = 1
                    highlightRowNumber.append(3)
                elif float(cluster_density) != 2961.28: #NovaSeq N2 ClustDen Check
                    emailSwitch = 1
                    lane_switch = 1
                    highlightRowNumber.append(4)

            else:
                raise Exception("Unhandled chemVer type: {}!".format(chemVer[0]['CHEMVER']))

            #print(samp,lane_switch,samp[0][0:4] == 'lane' and samp[1] > 5) 
        if lane_switch == 1:
            EmailLane(fcillumid,lanes,email,highlightRowNumber,database)

    email.close()
    send_email(emailSwitch,fcillumid,Machine,unaligned_dir)

def send_email(emailSwitch,fcillumid,Machine,unaligned_dir):
    logger = logging.getLogger(__name__)
    if emailSwitch == 1 and os.path.isfile('%s/EmailSent.txt' % unaligned_dir) == False:
        address = "igm-hts@columbia.edu"
        emailProgramLocation = '/nfs/goldstein/software/mutt-1.5.23/bin/mutt '
        emailCmd = emailProgramLocation + '-e "set content_type=text/html" '
        emailCmd += '-s \"Problem with Lane Fractions for flowcell %s %s\" ' % (fcillumid,Machine)
        emailCmd += address
        emailCmd += " < %s/LnFractionEmail.txt" % unaligned_dir

        logger.info(emailCmd)
        #print(emailCmd)
        os.system(emailCmd)
        os.system("touch %s/EmailSent.txt" % unaligned_dir)


def getKapaPicoDBID(fcillumid,CHGVID,database):
    kapaPicoDBID = run_query("SELECT KAPA_CONC,s2r.dbid as DBID \
            FROM Lane l JOIN Flowcell f ON f.FCID=l.FCID \
            JOIN prepT pt ON l.prepID=pt.prepID \
            JOIN samplesTOrun s2r ON s2r.seqID=l.seqID \
            JOIN SampleT s ON s.DBID=pt.DBID \
            LEFT JOIN poolMembers pm ON \
                (pm.DBID=pt.DBID AND pm.poolID=l.poolID) \
            WHERE FCIllumID='{}' and\
                pt.CHGVID='{}'".format(fcillumid,CHGVID),database)
    return kapaPicoDBID

def EmailLane(fcillumid,lanes,email,highlightRowNumber,database):
    logger = logging.getLogger(__name__)
    logger.info('Lane fraction problem with lane: %s' % lanes[0][5])
    email.write('<tr><td colspan="7">Lane %s</td></tr>\n' %lanes[0][5]+'\n')
    cluster_density = run_query(GET_CLUSTER_DENSITY_FOR_LANE.format(fcillumid=fcillumid,lanenum=lanes[0][5]),database)

    for samp in lanes:
        if samp[0] == 'Undetermined':
            pass
        else:
            kapaPicoDBID = getKapaPicoDBID(fcillumid,samp[0],database)
            poolName = run_query(GET_POOLID_FROM_DBID.format(DBID=kapaPicoDBID[0]['DBID']),database)
            poolName= poolName[0]['CHGVID']
            logging.info('%s\t%s\t%s\t%s\t%s\t%s\t%s' %
                    (samp[0],samp[4],samp[2],samp[1],cluster_density[0]['CLUSTDEN'],
                     poolName,kapaPicoDBID[0]['KAPA_CONC']))

            email.write("<tr>\n")
            highlightRowCount = 0
            for column in (samp[0],samp[4],samp[2],samp[1],cluster_density[0]['CLUSTDEN'],
                           poolName,kapaPicoDBID[0]['KAPA_CONC']):
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


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", dest='fcillumid', required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-b','--bcl_drive', default='seqscratch_ssd', dest='bcl_drive',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('-a','--archive_dir', default='igmdata01', dest='archive_dir',
                        help="Specify scratch dir for bcl2fastq")
    parser.add_argument('--noStatusCheck', action='store_true', default=False,
                        help="Do not check status")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    return args

def main():
    config = get_config()
    args = parse_arguments()

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    run_info_dict = parse_run_parameters_xml(args.fcillumid,database)
    setup_logging(run_info_dict['machine'],args.fcillumid,config.get('locs','logs_dir'))
    logger = logging.getLogger(__name__)
    logger.info('Starting post BCL script')

    fcillumid=run_info_dict['FCIllumID']
    machine=run_info_dict['machine']
    unaligned_dir = '{}/{}_{}_{}_Unaligned'.format(config.get('locs','bcl2fastq_scratch_dir'),
                                             run_info_dict['runDate'],
                                             machine,fcillumid)

    try:
        sample_info = update_lane(unaligned_dir,fcillumid,args.verbose,database)
        if args.noStatusCheck == False:
            checkStatus(fcillumid,database)
        UpdateFC(fcillumid,unaligned_dir,sample_info,args.verbose,database)
        checkLaneFractions(fcillumid,machine,unaligned_dir,sample_info,database)

        logger.info('Done')
        print('Done')
    except:
        traceback.print_exc()
        logger.info('Post BCL Failure')
        print('Post BCL Failure')
        sys.exit(1)

main()
