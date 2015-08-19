from bcl import create_sss
from CHGV_mysql import getSequenceDB
import logging
import os

def main():
	pwd = os.getcwd()
	info = pwd.split('/')[4].split('_')
	Date = info[0]
	FCID = info[3]
	Machine = info[1]

	sequenceDB = getSequenceDB()
	create_sss(FCID,Machine,Date,sequenceDB)


main()
