#!/usr/bin/python
# weekly_reports.py
# Joshua Bridgers
# jb3816@cumc.columbia.edu
#
# Creates weekly report for directors overviewing informatics operations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import re
from collections import Counter
from datetime import datetime
from CHGV_mysql import getSequenceDB
from dateutil.relativedelta import relativedelta, MO

def main():

    sequenceDB = getSequenceDB()
    # Date of the last Monday @ 1:00pm.  Note that relativedelta = -2. 
    # Designed to run on the Monday the report is due
    lastMonday = datetime.today() + relativedelta(weekday=MO(-2))
    lastMonday = lastMonday.replace(hour=13,minute=0,second=0,microsecond=0)
    fig = plt.figure(figsize={10,44})
    ReleasedSamples(sequenceDB,lastMonday,fig)
    DragenAlignedSamples(sequenceDB,lastMonday,fig)
    PassedSamples(sequenceDB,lastMonday,fig)
    #OverallActivity(sequenceDB,lastMonday,fig)
    SamplesInQueue(sequenceDB)
    #QCReviewNeeded(sequenceDB)
    print 'Saving...'
    fig.savefig('test.png')
    #Notes

def SamplesInQueue(sequenceDB):
    queueQuery = ("SELECT COUNT(*) FROM tempStatus ts "
                  "WHERE (status='QC review needed' OR "
                  "status like 'Pipeline launch%' OR "
                  "status = 'Processing error OR' ")

def ReleasedSamples(sequenceDB,lastMonday,fig):

    internalQuery = ("SELECT CHGVID FROM statusT st "
              "JOIN Lane l ON st.prepID=l.prepID "
              "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
              "FROM_UNIXTIME(status_time) > '{}' AND "
              "FCILLUMID NOT LIKE 'X%' AND "
              "status='Released to Bioinformatics Team'"
              ).format(lastMonday)
    sequenceDB.execute(internalQuery)
    internalReleasedSamples = sequenceDB.fetchall()
    title = 'Internal Released Samples (n={})'.format(len(internalReleasedSamples))
    getPieChart(internalReleasedSamples,title,fig,'421')

    externalQuery =  ("SELECT CHGVID FROM statusT st "
              "JOIN Lane l ON st.prepID=l.prepID "
              "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
              "FROM_UNIXTIME(status_time) > '{}' AND "
              "FCILLUMID LIKE 'X%' AND "
              "status='Released to Bioinformatics Team'"
              ).format(lastMonday)
    sequenceDB.execute(externalQuery)
    externalReleasedSamples = sequenceDB.fetchall()
    title = 'External Released Samples (n={})'.format(len(externalReleasedSamples))
    getPieChart(externalReleasedSamples,title,fig,'422')

def DragenAlignedSamples(sequenceDB,lastMonday,fig):
    internalDragenQuery = ("SELECT CHGVID FROM Lane l "
                          "JOIN prepT p on p.prepid=l.prepID "
                          "JOIN pseudo_prepid pp ON p.prepid=pp.prepid "
                          "JOIN dragen_pipeline_step dps on dps.pseudo_prepid=pp.pseudo_prepid "
                          "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
                          "FCILLUMID NOT LIKE 'X%' AND "
                          "pipeline_step_id = 11 AND "
                          "finish_time > '{}' "
                          ).format(lastMonday)
    sequenceDB.execute(internalDragenQuery)
    internalDragenAlignedSamples = sequenceDB.fetchall()
    title = 'Internal Dragen Aligned Samples (n={})'.format(len(internalDragenAlignedSamples))
    getPieChart(internalDragenAlignedSamples,title,fig,'423')

    externalDragenQuery = ("SELECT CHGVID FROM Lane l "
                          "JOIN prepT p on p.prepid=l.prepID "
                          "JOIN pseudo_prepid pp ON p.prepid=pp.prepid "
                          "JOIN dragen_pipeline_step dps on dps.pseudo_prepid=pp.pseudo_prepid "
                          "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
                          "FCILLUMID LIKE 'X%' AND "
                          "pipeline_step_id = 11 AND "
                          "finish_time > '{}' "
                          ).format(lastMonday)
    sequenceDB.execute(externalDragenQuery)
    externalDragenAlignedSamples = sequenceDB.fetchall()
    title = 'External Dragen Align Samples (n={})'.format(len(externalDragenAlignedSamples))
    getPieChart(externalDragenAlignedSamples,title,fig,'424')


def PassedSamples(sequenceDB,lastMonday,fig):
    internalQuery =  ("SELECT CHGVID FROM statusT st "
              "JOIN Lane l ON st.prepID=l.prepID "
              "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
              "FROM_UNIXTIME(status_time) > '{}' AND "
              "FCILLUMID NOT LIKE 'X%' AND "
              "(status='Passed Bioinfo QC' OR "
              "status='QC review needed')"
              ).format(lastMonday)
    sequenceDB.execute(internalQuery)
    internalPassedSamples = sequenceDB.fetchall()
    title = 'Internal Prod. Passed Samples (n={})'.format(len(internalPassedSamples))
    getPieChart(internalPassedSamples,title,fig,'425')

    externalQuery =  ("SELECT CHGVID FROM statusT st "
              "JOIN Lane l ON st.prepID=l.prepID "
              "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
              "FROM_UNIXTIME(status_time) > '{}' AND "
              "FCILLUMID LIKE 'X%' AND "
              "(status='Passed Bioinfo QC' OR "
              "status='QC review needed')"
              ).format(lastMonday)
    sequenceDB.execute(externalQuery)
    externalPassedSamples = sequenceDB.fetchall()
    title = 'External Prod. Passed Samples (n={})'.format(len(externalPassedSamples))
    getPieChart(externalPassedSamples,title,fig,'426')

    internalDragenQuery = ("SELECT CHGVID FROM Lane l "
                          "JOIN prepT p on p.prepid=l.prepID "
                          "JOIN pseudo_prepid pp ON p.prepid=pp.prepid "
                          "JOIN dragen_pipeline_step dps on dps.pseudo_prepid=pp.pseudo_prepid "
                          "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
                          "FCILLUMID NOT LIKE 'X%' AND "
                          "pipeline_step_id = 31 AND "
                          "finish_time > '{}' "
                          ).format(lastMonday)
    sequenceDB.execute(internalDragenQuery)
    internalDragenPassedSamples = sequenceDB.fetchall()
    title = 'Internal Dragen Passed Samples (n={})'.format(len(internalDragenPassedSamples))
    getPieChart(internalDragenPassedSamples,title,fig,'427')

    externalDragenQuery = ("SELECT CHGVID FROM Lane l "
                          "JOIN prepT p on p.prepid=l.prepID "
                          "JOIN pseudo_prepid pp ON p.prepid=pp.prepid "
                          "JOIN dragen_pipeline_step dps on dps.pseudo_prepid=pp.pseudo_prepid "
                          "JOIN Flowcell f ON f.FCID=l.FCID WHERE "
                          "FCILLUMID LIKE 'X%' AND "
                          "pipeline_step_id = 31 AND "
                          "finish_time > '{}' "
                          ).format(lastMonday)
    sequenceDB.execute(externalDragenQuery)
    externalDragenPassedSamples = sequenceDB.fetchall()
    title = 'External Dragen Passed Samples (n={})'.format(len(externalDragenPassedSamples))
    getPieChart(externalDragenPassedSamples,title,fig,'428')

def OverallActivity(sequenceDB,lastMonday,fig):
    activityQuery = ("SELECT status")
    pass

def getPieChart(sampleList,plotTitle,fig,subplot):
    sampleAlphaPrefixList = []
    sampleNumPrefixList = []
    for sample in sampleList:
        tmpPrefix = ''
        #Get prefixes of samples starting with alphabetic chars.
        if sample[0][0].isalpha():
            for char in sample[0]:
                if char.isalpha():
                    tmpPrefix += char
                else:
                    sampleAlphaPrefixList.append(tmpPrefix)
                    break
        #Get prefixes of samples starting with numeric chars.
        else:
            for char in sample[0]:
                if char.isalpha() == False:
                    tmpPrefix += char
                else:
                    sampleNumPrefixList.append(tmpPrefix)
                    break

    prefixCount = Counter(sampleAlphaPrefixList + sampleNumPrefixList)

    k = []
    v = []
    miscCount = 0
    for key in prefixCount.keys():
        if prefixCount[key] < 10:
            miscCount += prefixCount[key]
        else:
            k.append(key)
            v.append(prefixCount[key])
    if miscCount > 0:
        k.append('MISC')
        v.append(miscCount)
    v = numpy.array(v)

    def absolute_value(val):
        a = numpy.round(val/100.*v.sum(), 0)
        return int(a)
    ax = fig.add_subplot(subplot)
    ax.set_title(plotTitle)
    plt.pie(v, labels=k,autopct=absolute_value,shadow=True,
            labeldistance=1.2,pctdistance=0.9)

main()


