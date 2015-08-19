#!/usr/bin/python
#Joshua Bridgers
#jb371@dm.duke.edu
#11/10/2013
#Creates Lane and FCID entries into the GAF database from Sequencing Sample Sheets.

import glob
import os
import sys
import traceback
from commands import getoutput

def check_dup_index(FCID,lanenum):
	all_indices = getoutput("grep ,%s, *%s.csv | cut -d, -f5 | wc -l" % (lanenum,FCID))
	sort_indices = getoutput("grep ,%s, *%s.csv | cut -d, -f5 | sort -u | wc -l" % (lanenum,FCID))

	if all_indices != sort_indices:
		raise Exception, "Duplicate entry for index exists for lane %s" % lanenum

#performs a checks on sss
def sss_qc(FCID):
	#opening Sequencing Sample Sheet
	sss = glob.glob('/nfs/genotyping/Sequencing_SampleSheets/*%s*.csv' % FCID)
	input1 = open(sss[0],'r')
	samples = input1.readlines()

	if samples[0][0:4] == 'FCID':
		del samples[0] #removes header
	else:
		raise Exception, 'Warning: Header is not formatted correctly'
	start_swt = 0
	#need to remove blank lines
	sss_lanes = getoutput("cut -d, -f2 /nfs/genotyping/Sequencing_SampleSheets/*%s* | sort -u | grep -v Lane" % FCID).split('\n')
	adict = {}
	for l in sss_lanes:
		adict[l] = 0
	for s in samples:
		#check if there are duplicate index entries
		if ' ' in s:
			#print s
			raise Exception, 'Space(s) are present in SSS!'
		elif '/' in s or '\\' in s or '"' in s or "'" in s:
			raise Exception, 'Special characters are present in SSS!'
		else:
			splitted = s.strip(',').split(',')
			if len(splitted) <10:
				raise Exception, 'Warning: Line is not formatted correctly!'
			adict[splitted[1]] += float(splitted[5].split('_')[0])

			sampleID = splitted[2]
			index = splitted[4]
			if index == '':
				index = None

			if start_swt == 0:
				old_sssFCID = splitted[0]
				old_reads = splitted[7]
				old_TechFC = splitted[8]
				start_swt = 1
			else:
				sssFCID = splitted[0]
				if len(sssFCID) != 9:
					raise Exception, 'Warning: FCID is not correct!'
				LnFraction = splitted[5].split('_')[0]
				reads = splitted[7]
				TechFC = splitted[8]
				if sssFCID == old_sssFCID and sssFCID == FCID:
					old_sssFCID = sssFCID[:]
				else:
					raise Exception, 'Warning: FCID does not match!'
				if TechFC == old_TechFC:
					old_TechFC = TechFC[:]
				else:
					raise Exception, 'Warning: Operator does not match throughout!'
				if reads == old_reads:
					old_reads = reads[:]
				else:
					raise Exception, 'Warning: Recipe does not match throughout'

	for key in adict:
		check_dup_index(FCID,key)
		TotalLnFrac = adict[key]
		if key != '':
			#print TotalLnFrac
			if TotalLnFrac < 0.97:
				raise Exception, 'Warning: Lane Fraction do not add approximately to 1 for %s for lane %s!' % (FCID,key)
			elif TotalLnFrac > 1.002:
				raise Exception, 'Warning: Lane Fraction adds up to greater than 1!'
	return len(adict)

def main():
	email = 'jb3816@cumc.columbia.edu'
	#parsing values from filename
	info = sys.argv[1].split('/')[-1].split('_')
	print info
	Machine	= info[0]
	date = info[1]
	FCID = info[2].split('.')[0]

	total_lane_num = sss_qc(FCID)

if __name__ == '__main__':
	main()
