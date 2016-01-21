from CHGV_mysql import getSequenceDB
import glob
import os

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
    
    print query
    sequenceDB.execute(query)
    completeFlowcell = sequenceDB.fetchall()

    checkBCLFolderExist(completeFlowcell)
   


def checkBCLFolderExist(completeFlowcell):
    '''Checks for the exisitance of a bcl2fastq folder.  The folder has to be
       deleted before proper running of the sequence pipeline'''

    for FCIllumID in completeFlowcell:

        BCLFolder = glob.glob("/nfs/[sf][ea][qs][st][cq]*[0-9]/BCL/*%s_Unaligned" % FCIllumID)
        if BCLFolder == []:
            cmd = "python2.7 /nfs/goldstein/software/sequencing_pipe/production/GAF_PIPELINE.py --FCID %s -r" % FCIllumID
            print cmd
            #os.system(cmd)

if __name__ == '__main__':
    main()
