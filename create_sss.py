from bcl import create_sss
from CHGV_mysql import getSequenceDB,MachineCheck
import logging
import os

def main():
    pwd = os.getcwd()
    sequenceDB = getSequenceDB()

    info = pwd.split('/')[4].split('_')
    Date = info[0]
    FCID = info[3]
    Machine = MachineCheck(sequenceDB,info[1],FCID)


    create_sss(FCID,Machine,Date,sequenceDB)


main()
