from CHGV_mysql import getSequenceDB
import glob
import os
import sys

def main():
    sequenceDB = getSequenceDB()
    
    #grab unfailed, completed flowcells that haven't completed the sequenceDB pipeline
    query = ("select FCILLUMID "
            "from Flowcell f "
            "join Lane l on f.fcid=l.fcid "
            "join SampleT s on s.dbid=l.dbid "
            "where (pipelineComplete != 1 and complete = 1 and fail != 1) "
            "group by l.FCID "
            "order by min(priority) asc , from_unixtime(seqend) asc")
    
    #print query
    sequenceDB.execute(query)
    completeFlowcell = sequenceDB.fetchall()
    if len(completeFlowcell) != 0: 
        checkBCLFolderExist(completeFlowcell,sequenceDB)
    else:
        sys.exit(0)
   


def checkBCLFolderExist(completeFlowcell,sequenceDB):
    '''Checks for the exisitance of a bcl2fastq folder.  The folder has to be
       deleted before proper running of the sequence pipeline'''

    for FCIllumID in completeFlowcell:
        #check seqscratch or fastq folders for
        BCLFolder = glob.glob("/nfs/[sf][ea][qs][st][cq]*[0-9]/BCL/*%s_Unaligned" % FCIllumID)
        if BCLFolder == []:
            cmd = "/nfs/goldstein/software/python2.7/bin/python2.7 /nfs/goldstein/software/sequencing_pipe/production/GAF_PIPELINE.py --FCID %s -r" % FCIllumID
            cmd = "/nfs/goldstein/software/python2.7/bin/python2.7 /home/jb3816/github/sequencing_pipe/GAF_PIPELINE.py --FCID %s -r" % FCIllumID

            print cmd
            query = ("UPDATE Flowcell "
                "SET pipelineComplete = pipelineComplete + 1 "
                "WHERE FCIllumID = '{0}'"
                ).format(FCIllumID[0])
            print query

            os.system(cmd)
            #sequenceDB.execute(query)
        else:
            print "BCL folder already exists for %s" % FCIllumID
            print BCLFolder
        sequenceDB.execute('COMMIT')



if __name__ == '__main__':
    main()
