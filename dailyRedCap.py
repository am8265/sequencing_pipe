# bcl.py
# Joshua Bridgers
# 10/02/2015
# jb3816@cumc.columbia.edu
#
# Script designed to run daily to search for samples that have been updated
# with a specified status and then update the redcap database and email
# collaborators.

from requests import post
from json import dumps
import copy
import datetime
import itertools
import os
import pymysql.cursors
import sys
import traceback

def main():
    
    try:
        sequenceDB = getSequenceDB()
        flagableStatus = ["Sequencing Queue","External Data Submitted","QC Review Needed","In Annotation DB"]
        #flagableStatus = ["Passed Bioinfo QC"]
        samplePrefixes = ['EGI']
        emailProgramLocation = '/nfs/goldstein/software/mutt-1.5.23/bin/mutt '

        emails = getEmails()

        allSamples = getSamples(sequenceDB,flagableStatus,samplePrefixes)
        allSamples = fixSamples(allSamples)
        #print(allSamples)
        #allSamples = ('test','Sequencing Queue')
        sampleSites = updateRedcap(allSamples)
        emailCollaborators(sampleSites,emails,allSamples,emailProgramLocation)
    except:
        failureText = traceback.print_exc()
        emailFailure(emailProgramLocation,failureText)

def emailFailure(emailProgramLocation,failureText):
    today_date = datetime.datetime.now().strftime("%Y-%m-%d")
    print(failureText)
    emailCmd = emailProgramLocation
    emailCmd += '-s \"%s %s\" ' % ('dailyRedcap failure',today_date)
    emailCmd += 'jb3816@cumc.columbia.edu '
    emailCmd += "< <(echo test)" 
    print(emailCmd)
    os.system(emailCmd)


def getEmails():

    email = {}
    email['Melbourne'] = ['aschneider@unimelb.edu.au','06']
    email['Boston'] = ['beth.sheidley@childrens.harvard.edu','04']
    email['NYU'] = ['patricia.tolete@nyumc.org','05']
    email['UCSF'] = ['joseph.sullivan@ucsf.edu','02']
    email['CHOP'] = ['dubbsh@email.chop.edu','03']
    email['Lurie'] = ['DMiazga@luriechildrens.org','07']
    email['Columbia'] = ['lb2993@cumc.columbia.edu','01']
    
    return email
def fixSamples(allSamples):
    '''Fixes double status update for externally submitted samples.  Samples
    that are externally submitted will have both the "External Data Submitted"
    status and also the "Sequencing Queue" status.  This function will remove
    the "Sequencing Queue" status from the dictionary'''

    allSamplesCopy = copy.deepcopy(allSamples)
    for items in allSamples:
        if items['status'] == 'External Data Submitted':
            #print(items)
            allSamplesCopy.remove({'CHGVID': items['CHGVID'], 'status':
                'Sequencing Queue'})
    
    #print(allSamplesCopy)
    return allSamplesCopy

def getSamples(sequenceDB,flagableStatus,samplePrefixes):
    allSamples = []
    timeDelta = 24
    #timeDelta = 168
    yesterday_date = datetime.datetime.now() - datetime.timedelta(hours = timeDelta)
    yesterday_date = yesterday_date.strftime("%Y-%m-%d")
    #print yesterday_date
    for status in flagableStatus:
        for prefix in samplePrefixes:

            sql = ("SELECT distinct CHGVID,status "
                "FROM statusT "
                "WHERE from_unixtime(status_time) > '{0}' "
                "AND CHGVID LIKE '{1}%' "
                "AND status = '{2}' "
                ).format(yesterday_date,prefix,status)
            sequenceDB.execute(sql)
            #print(sequenceDB.fetchall())
            allSamples.append(sequenceDB.fetchall())

    mergedSamples = list(itertools.chain(*allSamples))
    #print mergedSamples
    return mergedSamples

def updateRedcap(allSamples):
    #TOKEN = sys.argv[1]
    TOKEN = '3B7C39D9720D72EB8B03CD468B060C1A'
    URL = "https://wchredcap.cumc.columbia.edu/redcap/api/"

    #allSamples = [{'status': 'Josh was here', 'IGMID': 'EGI99.999TEST1'}]
    payload = {'token': TOKEN, 'format': 'json', 'content': 'metadata'}
    #response = post(URL, data=payload)
    #print(response.status_code)
    payload['content'] = 'record'
    payload['type'] = 'flat'
    response = post(URL, data=payload)
    data = response.json()

    sampleSites = []
    for sample in allSamples:
        IGMID = sample['CHGVID']
        status = sample['status']
        for record in data:
            #print(record)

            if record['igm_seq'] == IGMID:
                print("Updating sample {0} with status '{1}'".format(IGMID,status))
                #print(record['chgv_seq_status'])
                record['chgv_seq_status'] = status
                to_import_json = dumps([record], separators=(',',':'))
                payload['data'] = to_import_json
                response = post(URL, data=payload)

                #print(response.json()['count'])
                #print(record)

                sampleSites.append((record['igm_seq'],'Columbia'))
    
    return sampleSites
def emailCollaborators(sampleSites,emails,allSamples,EmailProgramLocation):

    #emails = []
    #emails = ['jb3816@cumc.columbia.edu','joshbridgers@gmail.com']
    today = datetime.datetime.now().strftime("%Y-%m-%d")
    subject = 'Samples in SequenceDB pipeline %s' % today

  
    for email in emails:
    
        sampleList = []
        for sample in allSamples:
            #The two digits after 'EGI' denote the site
            #print(sample['IGMID'][3:5])
            if sample['CHGVID'][3:5] == emails[email][1]:
                sampleList.append(sample)

        sampleNumber = len(sampleList)

        release_email = open('/home/jb3816/email.tmp','w')
        release_email.write('The following %s sample(s) are being processed through SequenceDB\n' % (sampleNumber))
        release_email.write('\n')
        release_email.write('IGM SEQ\t\tStatus\n')
        release_email.write('='*80+'\n')

        #print(sampleList,emails[email])
        if len(sampleList) > 0:

            #addresses = [emails[email][0],'jb3816@cumc.columbia.edu','lb2993@cumc.columbia.edu']
            addresses = ['jb3816@cumc.columbia.edu','lb2993@cumc.columbia.edu']
            #addresses = ['joshbridgers@gmail.com','jb3816@cumc.columbia.edu']

            for sample in sampleList:
                release_email.write("{0}\t{1}\n".format(sample['CHGVID'],sample['status']))
            
            release_email.write('\nFor any questions regarding the status of your'
                    ' samples or this message please email Louise Bier at lb2993@cumc.columbia.edu\n')

            release_email.write('\nEmail intended for %s\n' % emails[email][0])
            release_email.write('\nThanks,\n')
            release_email.write('\n')
            release_email.write('%s\n' % 'Joshua Bridgers')

            release_email.close()
            os.system('cat /home/jb3816/email.tmp')

            for address in addresses:

                #print("Emailing {0}".format(address))
                emailCmd = emailProgramLocation
                emailCmd += '-s \"%s\" ' % (subject)
                emailCmd += address
                emailCmd += " < /home/jb3816/email.tmp"
                print(emailCmd)
                os.system(emailCmd)

        release_email.close()
        os.system('rm /home/jb3816/email.tmp')


def getSequenceDB():
    connection = pymysql.connect(host='10.73.50.38',
        user='sequence_connect',
        passwd='g3n3t1c5213',
        db='sequenceDB',
        charset='utf8mb4',
        cursorclass=pymysql.cursors.DictCursor)
    return connection.cursor()

if __name__ == '__main__':
    main()

