#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Define a class describing samples going through the dragen pipeline
"""

import argparse
import pymysql
import os
import pdb
import subprocess
import sys
import time
from db_statements import *
from glob import glob
from utilities import check_number_query_results,get_connection,is_external_or_legacy_sample

class FastqError(Exception):
    """Fastq not found"""
    pass

def run_query(query,database):
    connection = get_connection(database)
    try:
        with connection.cursor() as cursor:
            cursor.execute(query)
            results = cursor.fetchall()
            connection.commit()
            return results
    finally:
        connection.close()

def get_prepid(database,sample):
    """Retrieve qualifying prepids"""
    query = ("SELECT prepid FROM prepT WHERE p_prepid={0}"
            ).format(sample["pseudo_prepid"])
    prepids = run_query(query,database)
    prepids = [x['prepid'] for x in prepids]
    return prepids

def get_dbid(database,sample):
    query = """SELECT DBID FROM prepT WHERE p_prepid={}
            """.format(sample['pseudo_prepid'])
    results = run_query(query,database)
    check_number_query_results(results,1)
    return results[0]['DBID']

def get_priority(database,sample):
    query = ("SELECT PRIORITY "
            "FROM SampleT s "
            "JOIN prepT p ON s.DBID=p.DBID "
            "WHERE P_PREPID ={pseudo_prepid} "
            "ORDER BY PRIORITY DESC "
            "LIMIT 1"
            ).format(**sample)
    priority = run_query(query,database)[0]['PRIORITY']
    return priority

def get_pseudo_prepid(database,sample):
    query = ("SELECT DISTINCT p_prepid "
            "FROM prepT p "
            "WHERE CHGVID='{sample_name}' "
            "AND SAMPLE_TYPE='{sample_type}' "
            "AND FAILEDPREP=0 "
            ).format(**sample)
    if sample_type.lower() != 'genome':
        query += ("AND EXOMEKIT='{capture_kit}' "
                 ).format(**sample)
    query += "ORDER BY PRIORITY DESC "

    pseudo_prepid = run_query(query)
    if len(pseudo_prepid) == 1:
        return pseudo_prepid[0]['p_prepid']
    else:
        print(pseudo_prepid)
        raise ValueError('Too many pseudo_prepids found!')

def get_bed_file_loc(database,capture_kit):
    """Retrieve the bed file location for a given capture_kit name"""
    query = (("SELECT region_file_lsrc FROM captureKit WHERE name='{0}' "
        "and chr = 'all'"
        ).format(capture_kit))
    #print(query)
    bed_file_loc = run_query(query,database)[0]['region_file_lsrc']
    return bed_file_loc

def get_fastq_loc(database, sample):
    found_locs = []
    #correcting the downstream consequences of "custom capture" as the sample_type
    corrected_sample_type = sample['sample_type'].upper().replace(" ","_")
    sample_name=sample['sample_name']
    for prepid in sample["prepid"]:
        external_or_legacy_flag = is_external_or_legacy_sample(prepid,database)
        """For cases where there is not flowell information the sequeuncing
        will have to found via brute force.  There will be three types of samples
        that under this catergory: Old Casava 1.8 sample (pre-seqDB),
        Casava 1.7 samples sequenced by the Illumina GAII and external samples.

        For externally submitted sample they have a fake Flowcell entry
        in the database.  The Flowcell.FCIllumID begins wt'X' always.
        Typically when processing the external sample each read group is
        enumerated 1,2,3...etc.
        Secondly when external samples are archived sometimes the FCIllumID
        is preserved otherwise its enumerated."""

        potential_locs = ['/nfs/seqscratch_ssd/{}/'.format(corrected_sample_type),
                          '/nfs/fastq_temp2/{}/'.format(corrected_sample_type),
                          '/nfs/seqscratch10/SRR/',
                          '/nfs/sequencing/tx_2390/CGND*/Project*/',
                          '/nfs/fastq_temp/tx_temp/tx_2390/CGND*/Project*/',
                          '/nfs/seqscratch*/tx_temp/tx_2390/CGND*/Project*/',
                          '/nfs/stornext/seqfinal/casava1.8/whole_{}/'.format(corrected_sample_type),
                          '/nfs/fastq1[568]/{}/'.format(corrected_sample_type),
                          '/nfs/igmdata01/{}/'.format(corrected_sample_type),
                          '/nfs/seqsata*/seqfinal/whole_genome/']

        if external_or_legacy_flag == True:
            for potential_loc in potential_locs:
                potential_path = '{}/{}/*[0-9xXyY]'.format(potential_loc,sample_name)
                #print('Checking {} for fastqs....'.format(potential_loc))
                fastq_loc = glob(potential_path)
                if fastq_loc != []:
                    print('FOUND: {}'.format(potential_path))
                    for folder in fastq_loc:
                        found_locs.append(folder)
                    break
            if fastq_loc == []:
                raise FastqError('{} Fastq files not found'.format(sample_name))
        else:
            query = ("SELECT DISTINCT SEQSATALOC,FCILLUMID FROM Flowcell f "
                "JOIN Lane l ON f.FCID=l.FCID "
                "JOIN prepT p ON l.prepID=p.prepID "
                "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
                "AND l.prepID={0} AND f.fail=0 "
                "GROUP BY f.fcid"
                ).format(prepid)
            seqsatalocs = run_query(query,database)
            print(seqsatalocs)
            for flowcell in seqsatalocs:
                if flowcell['SEQSATALOC'] == '':
                    msg = ('Sample {} on flowcell {} is missing seqsataloc!  Still sequencing?'
                          ).format(sample_name,flowcell['FCILLUMID'])

                    raise ValueError(msg)
                for potential_loc in potential_locs:
                    potential_path = '{}/{}/{}'.format(potential_loc,sample_name,flowcell['FCILLUMID'])
                    fastq_loc = glob(potential_path)
                    if fastq_loc != []:
                        print('FOUND: {}'.format(potential_path))
                        for folder in fastq_loc:
                            found_locs.append(folder)
                        break
                if fastq_loc == []:
                    raise FastqError('{} Fastq files not found'.format(sample_name))


    """For samples in the database we can exclude any samples that only have
    R1 data however for sampels that predate the database we have to manually
    check for R2 existance"""
    final_locs = check_fastq_locs(found_locs)
    return final_locs

def check_fastq_locs(locs):
    """Determine if loc has read 2 fastq files"""

    valid_locs = []
    for loc in locs:
        read2 = glob("{loc}/*R2_[0-9]*fastq*".format(loc=loc))
        if read2 != []:
            valid_locs.append(loc)
        else:
            print('Did not find fastq mate pair for: {}!'.format(loc))
    return valid_locs


def get_lanes(database,sample):
    """retrieve all qualifying lanes,flowcells for the prepids associated
    with the sample while creating PE and SE conf files for the sample"""
    lanes = []
    """For cases where there is not flowell information the sequeuncing
    will have to be manually.  There will be two types of samples that
    under this catergory: Old Casava 1.8 sample (pre-seqDB) and Casava 1.7
    samples sequenced by the Illumina GAII."""

    for prepid in sample['prepid']:
        """ Queries Lane table for lanes that do not have a failing lane for read1
        or read two since the dragen cannot use a mix of single and paired end reads
        """
        query = ("SELECT DISTINCT lanenum,FCIllumID,p.prepID from Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            "AND l.prepID={0} AND p.failedPrep=0 "
            ).format(prepid)
        lane = run_query(query,database)
        leg_or_ext = is_external_or_legacy_sample(prepid,database)

        if lane and leg_or_ext == False:
            for entry in lane:
                lanes.append((entry['lanenum'],entry['FCIllumID'],entry['prepID']))
        else: #In case of no database entry or external sample
            for flowcell in sample['fastq_loc']:
                fastqs = glob(flowcell + '/*fastq.gz')
                lane = []
                for fastq in fastqs:
                    lane_num = fastq[fastq.find('L00')+3]
                    lane.append(lane_num)
                lane_nums = set(sorted(lane))
                for lane_num in lane_nums:
                    lanes.append((lane_num,flowcell.split('/')[-1],prepid))
        lanes = (lanes,)
    return lanes

def isInternalSample(laneFCIDTuples):
    internalFlag = 1
    for laneFCID in laneFCIDTuples:
        lane,FCID,prepid = laneFCID
        if FCID[0] != 'X':
            internalFlag = 0

    if internalFlag == 0:
        return True
    else:
        return False

class dragen_sample:
    # store all relevant information about a sample in a dictionary
    def __init__(self, sample_name, sample_type, pseudo_prepid, capture_kit, database):
        self.metadata = {}
        self.metadata['sample_name'] = sample_name
        self.metadata['sample_type'] = sample_type.lower()
        self.metadata['pseudo_prepid'] = pseudo_prepid
        if sample_type.lower() != 'genome':
            self.metadata['capture_kit'] = capture_kit
            self.metadata['bed_file_loc'] = get_bed_file_loc(database,self.metadata['capture_kit'])
        else:
            self.metadata['capture_kit'] = ''
            self.metadata['bed_file_loc'] = '/nfs/goldsteindata/refDB/captured_regions/Build37/65MB_build37/SeqCap_EZ_Exome_v3_capture.bed'
        self.metadata['dbid'] = get_dbid(database,self.metadata)
        self.metadata['prepid'] = get_prepid(database,self.metadata)
        self.metadata['priority'] = get_priority(database,self.metadata)
        self.metadata['fastq_loc'] = get_fastq_loc(database,self.metadata)
        self.metadata['lane'] = get_lanes(database,self.metadata)

    def get_attribute(self, attribute):
        """return the value requested if present, otherwise raise a TypeError"""
        if attribute in self.metadata:
            return self.metadata[attribute]
        else:
            raise TypeError("{attribute} is not defined for {CHGVID}".format(
                attribute=attribute, CHGVID=self.CHGVID))

    def get_dict(self):
        """return a dict of the values of this sample"""
        return self.metadata

    def set(self, attribute, value):
        """set the specified attribute to the given value"""
        self.metadata[attribute] = value
