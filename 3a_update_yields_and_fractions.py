#!/usr/bin/python
# post_bcl.py
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

def update_flowcell_yield_and_casava_verision(fcillumid,fc_yield,unaligned_dir,run_info_dict,sample_info,verbose,database):
    '''Updates Flowcell total Yield, Casava version for FCIllumID'''
    logger = logging.getLogger(__name__)
    casavaCmd = ['grep','bcl2fastq','-m1','{}/nohup.sge'.format(unaligned_dir)]
    if verbose:
        print(' '.join(casavaCmd))
    proc = subprocess.Popen(casavaCmd, stdout=subprocess.PIPE)
    CasavaVer = proc.stdout.read().decode().split()[1].strip()
    sql = ("UPDATE Flowcell f "
        "JOIN Lane l ON f.FCID=l.FCID "
        "SET fcYield={0},CasavaVer='{1}' "
        "WHERE f.FCillumID='{2}'"
        ).format(fc_yield,CasavaVer,fcillumid)

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

def update_per_rg_aka_lane_table_yield_and_fractions_by_flowcell(config,machine,unaligned_dir,fcillumid,verbose,test,database):
    '''Gets Actual Lane Fraction FROM laneBarcode.htm, Updates sample status'''
    logger = logging.getLogger(__name__)
    laneBarcodeHTML = getHTML(unaligned_dir,fcillumid,verbose)

    #use xml.etree for HTML table parsing
    tree = fromstring(laneBarcodeHTML)
    rows = tree.findall("tr")
    headrow = rows[0]
    datarows = rows[1:]
    sample_info = {}
    max_lane = getTotalLanes(fcillumid,database)
    for i in range(1,max_lane+1):
        sample_info[str(i)] = []
    fc_yield = 0
    for sampleNum,i in enumerate(datarows):
        for num, h in enumerate(headrow):
            if num == 0:
                lanenum = datarows[sampleNum][num].text
            elif num == 2:
                chgvid = datarows[sampleNum][num].text
            elif num == 3:
                adapter = datarows[sampleNum][num].text
            elif num == 5:
                LnFracAct = datarows[sampleNum][num].text
            elif num == 8:
                LnYield = datarows[sampleNum][num].text
                LnYield = LnYield.replace(',','')
                fc_yield += int(LnYield)
                if chgvid != 'Undetermined':
                    LnFracQuery = ("SELECT SAMPLE_TYPE,LNFRACTION "
                                    "FROM Lane l "
                                    "JOIN Flowcell f ON l.FCID=f.FCID "
                                    "JOIN prepT p ON p.prepid=l.prepid "
                                    "WHERE FCillumID='{}' "
                                    "AND CHGVID='{}' "
                                    "AND LANENUM={}"
                                    ).format(fcillumid,chgvid,lanenum)
                    results = run_query(LnFracQuery,database)
                    LnFrac = results[0]['LNFRACTION']
                    SeqType = results[0]['SAMPLE_TYPE']
                    sql = ("UPDATE Lane l "
                        "JOIN prepT pt ON l.prepID=pt.prepID "
                        "JOIN Flowcell f ON l.FCID=f.FCID "
                        "SET l.LnYield={0}, l.LnFractionAct={1} "
                        "WHERE pt.chgvid='{2}' AND f.FCillumID='{3}' AND l.LaneNum={4}"
                        ).format(LnYield,LnFracAct,chgvid,fcillumid,lanenum)
                    if verbose == True:
                        print(sql)
                    logger.info(sql)
                    run_query(sql,database)
                else:
                    LnFrac = 0
                    SeqType = ''
                email_highlight_pos = check_clustden_lnfrac(config,chgvid,
                    LnFrac,LnFracAct,fcillumid,lanenum,database)
                sample_info[lanenum].append([chgvid,SeqType,LnFrac,LnFracAct,email_highlight_pos])
    makeHTMLEmail(config,unaligned_dir,machine,fcillumid,sample_info,test,database)
    return fc_yield,sample_info

def check_clustden_lnfrac(config,chgvid,LnFrac,LnFractionAct,fcillumid,LaneNum,database):
    chemver = run_query(GET_FLOWCELL_CHEMVER.format(fcillumid=fcillumid),database)[0]['CHEMVER']
    query = GET_CLUSTER_DENSITY_FOR_LANE.format(fcillumid=fcillumid,lanenum=LaneNum)
    cluster_density = run_query(query,database)[0]['CLUSTDEN']
    email_highlight_pos = []
    if chgvid == 'Undetermined':
        if float(LnFractionAct) > 3:
            email_highlight_pos.append(10)
    else:
        LnFracDiff = float(LnFractionAct)/(float(LnFrac) * 100)
        #print(LnFracDiff,LnFractionAct,LnFrac)
        if float(LnFracDiff) > 1.15 or float(LnFracDiff) < 0.85:
            email_highlight_pos.append(3)
        if float(cluster_density) > float(config.get(chemver,'high_cluster_den_threshold')):
            email_highlight_pos.append(4)
        elif float(cluster_density) < float(config.get(chemver,'low_cluster_den_threshold')):
            email_highlight_pos.append(4)

    return email_highlight_pos

def check_db_status_or_exit(fcillumid,database):
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
        status = status[0]['status']
        SeqType = run_query(GET_SEQTYPE_FROM_PREPID.format(prepid=prepID['prepID']),database)
        SeqType = SeqType[0]['SAMPLE_TYPE']
        # Sometimes multiple flowcells of Genomes finish bcl2fastq near the same time so its status =='storage'
        msg = ("{} does not have the correct status!  Current status: '{}'"
              ).format(status,status)
        if status != 'BCL' and status != 'BCL Complete' and status != 'Archiving' and SeqType != 'Genome':
            raise_exception_switch = 1
            logger.warn(msg)
            print(msg)
        elif (status != 'BCL' and status != 'Storage' and status != 'BCL Complete' and status != 'Archiving') and SeqType == 'Genome':
            raise_exception_switch = 1
            logger.warn(msg)
            print(msg)

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

def send_email(config,fcillumid,Machine,unaligned_dir,test):
    logger = logging.getLogger(__name__)
    if os.path.isfile('%s/EmailSent.txt' % unaligned_dir) == False:
        subject = '\"Lane Fraction Report for Flowcell {} {}'.format(fcillumid,Machine)
        if test == True:
            subject += ' TEST EMAIL\"'
        else:
            subject += '\"'
        address = config.get('emails','post_bcl_email')
        emailProgramLocation = '/nfs/goldstein/software/mutt-1.5.23/bin/mutt '
        emailCmd = emailProgramLocation + '-e "set content_type=text/html" '
        emailCmd += '-s {} '.format(subject)
        emailCmd += address
        emailCmd += " < %s/LnFractionEmail.txt" % unaligned_dir

        logger.info(emailCmd)
        print(emailCmd)
        op = os.system(emailCmd)
        if op != 0:
            raise Exception("{} gave incorrect return value: {}".format(emailCmd,op))
        op = os.system("touch %s/EmailSent.txt" % unaligned_dir)
        if op != 0:
            raise Exception("touch {}/EmailSent.txt gave incorrect return value: {}".format(unaligned_dir,op))
    else:
        raise Exception("Lane fraction email has already been sent")

def getKapaPicoPoolName(fcillumid,CHGVID,database):
    logger = logging.getLogger(__name__)
    cmd = "SELECT DISTINCT KAPA_CONC, pool.name as PoolName \
            FROM Flowcell f JOIN Lane l ON f.fcid=l.fcid JOIN prepT pt ON pt.prepid = l.prepid \
            JOIN samplesTOrun s2r ON s2r.seqid=l.seqid JOIN pool ON pool.id = pt.poolID \
            WHERE pt.chgvid='{}' and f.fcillumid='{}'".format(CHGVID,fcillumid)
    logger.info(cmd)
    kapaPicoPoolName = run_query(cmd,database)
    if len(kapaPicoPoolName) != 1:
        raise Exception("get_kapa_conc query gave incorrect number of results: {}".format(str(len(kapaPicoPoolName))))
    return kapaPicoPoolName

def makeHTMLEmail(config,unaligned_dir,machine,fcillumid,sample_info,test,database):
    email = open('%s/LnFractionEmail.txt' % unaligned_dir,'w')
    email = addHeader(email)
    logger = logging.getLogger(__name__)
    for lane in sample_info.keys():
        query = GET_CLUSTER_DENSITY_FOR_LANE.format(fcillumid=fcillumid,lanenum=lane)
        cluster_density = run_query(query,database)
        for samp in sample_info[lane]:
            if samp[0] == 'Undetermined':
                logging.info("Lane %s's unmatched reads percent: %s" % (lane,samp[1]))
                if 10 in samp[4]:
                    email.write('<tr><td colspan="7" align="center" bgcolor="#FFFF00">Lane %s\'s unmatched reads percent: %s</td></tr>\n' %
                        (lane,samp[3]))
                else:
                    email.write('<tr><td colspan="7" align="center">Lane %s\'s unmatched reads percent: %s</td></tr>\n' %
                        (lane,samp[3]))
                email.write('<tr><td colspan="7" align="center">&nbsp</td></tr>\n')
            else:
                kapaPicoPoolName = getKapaPicoPoolName(fcillumid,samp[0],database)
                #pool_query = GET_POOLID_FROM_DBID.format(DBID=kapaPicoDBID[0]['DBID'])
                #poolName = run_query(pool_query,database)
                #poolName= poolName[0]['CHGVID']
                logging.info('%s\t%s\t%s\t%s\t%s\t%s\t%s' %
                        (samp[0],samp[1],samp[2],samp[3],cluster_density[0]['CLUSTDEN'],
                         kapaPicoPoolName[0]['PoolName'],kapaPicoPoolName[0]['KAPA_CONC']))

                email.write("<tr>\n")
                for i,column in enumerate((samp[0],samp[1],samp[2],samp[3],
                                        cluster_density[0]['CLUSTDEN'],kapaPicoPoolName[0]['PoolName'],
                                        kapaPicoPoolName[0]['KAPA_CONC'])):
                    if i in samp[4]:
                        if float(samp[3]) < 0.05: #super low lane fract
                            email.write('<td bgcolor="##FF4500">%s</td>\n' % column)
                        else:
                            email.write('<td bgcolor="#FFFF00">%s</td>\n' % column)
                    else:
                        email.write("<td>%s</td>\n" % column)
                email.write("</tr>\n")
    email.close()
    send_email(config,fcillumid,machine,unaligned_dir,test)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", dest='fcillumid', required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('--noStatusCheck', action='store_true', default=False,
                        help="Do not check status")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    return args

def set_bcl_complete(fcillumid,database):

    samples = run_query("SELECT DISTINCT prepID FROM Lane l "
                        "JOIN Flowcell f on l.FCID=f.FCID "
                        "WHERE f.FCillumID='{}'".format(fcillumid),database)
    userID = get_user_id(database)
    for prepID in samples:
        status = run_query("SELECT status FROM prepT WHERE prepID={}".format(prepID['prepID']),database)
        status = status[0]['status']
        if status == 'BCL Started':
            #update
            run_query("UPDATE prepT p SET STATUS='BCL Complete',status_time=unix_timestamp() WHERE prepID={} AND STATUS='BCL Started'".format(prepID['prepID']),database) 
    #insert into statusT
    sample_status_insert_query = ("INSERT INTO statusT "
                                  "(CHGVID,STATUS_TIME,STATUS,DBID,PREPID,USERID,POOLID,SEQID,PLATENAME) "
                                  "SELECT DISTINCT(pt.CHGVID),unix_timestamp(),"
                                  "'BCL Complete',pt.DBID,pt.prepID,{},0,0,' ' "
                                  "FROM Flowcell f "
                                  "JOIN Lane l ON l.FCID=f.FCID "
                                  "JOIN prepT pt ON pt.prepID=l.prepID "
                                  "WHERE FCillumid='{}'").format(userID,fcillumid)
    run_query(sample_status_insert_query,database)



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

    ###### check for post-bcl2fastq proper exit checkpoint touch file
    check_bcl_complete(unaligned_dir)

    set_bcl_complete(fcillumid,database)

    #cmd="chmod -R 770 {0}".format(unaligned_dir)
    #op = os.system(cmd)
    #if op != 0:
    #    raise Exception("{} cmd incorrectly returned {}".format(cmd,op))
    
    try:

        fc_yield,sample_info = update_per_rg_aka_lane_table_yield_and_fractions_by_flowcell(
          config, machine, unaligned_dir, fcillumid,
          args.verbose, args.test, database
        )

        if args.noStatusCheck == False:
            check_db_status_or_exit(fcillumid,database)

        update_flowcell_yield_and_casava_verision(fcillumid,fc_yield,unaligned_dir,run_info_dict,sample_info,args.verbose,database)

        logger.info('Done')
        print('Done')

    except:

        traceback.print_exc()
        logger.info('Post BCL Failure')
        print('Post BCL Failure')
        sys.exit(1)

main()
