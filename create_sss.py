from bcl import create_sss
from CHGV_mysql import getSequenceDB,MachineCheck
import logging
import os
import sys

def main():
    pwd = os.getcwd()
    sequenceDB = getSequenceDB()
    
    if len(sys.argv) == 1:
        
        info = pwd.split('/')[4].split('_')
        Date = info[0]
        FCID = info[3]
        Machine = MachineCheck(sequenceDB,info[1],FCID)
    else:
        Date = sys.argv[1]
        Machine = sys.argv[2]
        FCID = sys.argv[3]

    create_sss(FCID,Machine,Date,sequenceDB)


main()
