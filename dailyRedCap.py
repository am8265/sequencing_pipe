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
import datetime
import itertools
import os
import sys
import pymysql.cursors

def main():
    
    sequenceDB = getSequenceDB()
    flagableStatus = ["Sequencing Queue","External Data Submitted"]
    samplePrefixes = ['EGI','IGM']
    emails = ['jb3816@cumc.columbia.edu']

    allSamples = getSamples(sequenceDB,flagableStatus,samplePrefixes)
    #allSamples = ('test','Sequencing Queue')
    updateRedcap(allSamples)
    #emailCollaborators(emails,allSamples)

def getSamples(sequenceDB,flagableStatus,samplePrefixes):
    allSamples = []
    timeDelta = 24
    timeDelta = 144
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
    TOKEN = sys.argv[1]
    URL = "https://wchredcap.cumc.columbia.edu/redcap/api/"

    allSamples = [{'status': 'Josh was here', 'CHGVID': 'EGI99.999TEST1'}]
    payload = {'token': TOKEN, 'format': 'json', 'content': 'metadata'}
    #response = post(URL, data=payload)
    #print(response.status_code)
    payload['content'] = 'record'
    payload['type'] = 'flat'
    response = post(URL, data=payload)
    data = response.json()
    for sample in allSamples:
        CHGVID = sample['CHGVID']
        status = sample['status']
        for record in data:
            #print(record)

            if record['igm_seq'] == CHGVID:
                print("Updating sample {0} with status '{1}'".format(CHGVID,status))
                print(record['chgv_seq_status'])
                record['chgv_seq_status'] = status
                to_import_json = dumps([record], separators=(',',':'))
                payload['data'] = to_import_json
                response = post(URL, data=payload)

                #print(response.json()['count'])


def emailCollaborators(emails,allSamples):
    today = datetime.datetime.now().strftime("%Y-%m-%d")
    subject = 'Samples in SequenceDB pipeline %s' % today
    emailProgramLocation = '/nfs/goldstein/software/mutt-1.5.23/bin/mutt '
    sampleNumber = len(allSamples)

    release_email = open('email.tmp','w')
    release_email.write('The following %s sample(s) are being processed through SequenceDB\n' % (sampleNumber))
    release_email.write('\n')

    for sample in allSamples:
        release_email.write("%s\n" % (sample))

    release_email.write('\nThanks,\n')
    release_email.write('\n')
    release_email.write('%s\n' % 'Joshua Bridgers')

    release_email.close()

    #os.system('cat email.tmp')

    for address in emails:
        emailCmd = emailProgramLocation
        emailCmd += '-s \"%s\" ' % (subject)
        emailCmd += address
        emailCmd += " < email.tmp"

    os.system(emailCmd)
    os.system('rm email.tmp')


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


