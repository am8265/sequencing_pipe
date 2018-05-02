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

def get_sample_id(database,sample):
    query = """SELECT sample_id FROM prepT JOIN Experiment ON prepT.experiment_id=Experiment.id WHERE p_prepid={}
            """.format(sample['pseudo_prepid'])
    results = run_query(query,database)
    check_number_query_results(results,1)
    return results[0]['sample_id']

def get_priority(database,sample):
    query = ("SELECT PRIORITY "
            "FROM prepT "
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

        from pprint import pprint as pp
        pp(sample)

        """For cases where there is not flowell information the sequeuncing
        will have to be found via brute force.  There will be three types of samples
        that under this catergory: Old Casava 1.8 sample (pre-seqDB),
        Casava 1.7 samples sequenced by the Illumina GAII and external samples.

        Typically when processing the external sample previously each read group is
        enumerated 1,2,3...etc.
        Now external samples are archived with the fcillumid."""

        ####### this is all an abomination!?!
        potential_locs = []

        if corrected_sample_type != 'GENOME':
            potential_locs = ['/nfs/seqscratch_ssd/{}/'.format(corrected_sample_type)]

        if sample['pseudo_prepid']<20000:
            potential_locs.insert(0,'/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type))

        potential_locs.extend([
            '/nfs/fastq_temp/{}/'.format(corrected_sample_type),
            '/nfs/fastq_temp5/tx_temp/tx_3036/',
            '/nfs/fastq_temp2/{}_ssd/'.format(corrected_sample_type),       ##### truelly evil!?!
            '/nfs/fastq_temp2/{}/'.format(corrected_sample_type),
            '/nfs/seqscratch*/tx_temp/tx_*/',
            '/nfs/sequencing/tx_2390/CGND*/Project*/',
            '/nfs/fastq_temp/tx_temp/tx_2390/CGND*/Project*/',
            '/nfs/seqscratch*/tx_temp/tx_2390/CGND*/Project*/',
            '/nfs/igmdata01/{}/'.format(corrected_sample_type),
            '/nfs/stornext/seqfinal/casava1.8/whole_{}/'.format(corrected_sample_type),
            '/nfs/fastq1[568]/{}/'.format(corrected_sample_type),
            '/nfs/fastq_temp5/tx_2995/Project_CGND_12737_B01_GRM_WGS.fastq.2017-08-29/'
            ###### this is causing so many issues that better to remove and add back later!?!
            ###### i.e. just check ALL fastq error!?!
            # '/nfs/seqsata*/seqfinal/whole_genome/'
        ])

        if external_or_legacy_flag == True:

            if sample["prepid"][0]<50000:
            # if sample["prepid"][0]<22000:
                # from pprint import pprint as pp
                # pp(sample)
                # print("avoid back contamination of older sequencing by newer")
                # potential_locs = ['/nfs/stornext/seqfinal/casava1.8/whole_exome/'.format(corrected_sample_type)]
                ##### fs isn't even mounted anymore so doing a copy
                # potential_locs = ['/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type)]
                ##### needs to be compatible with other old samples - why aren't they run?!?!
                potential_locs.insert(0,'/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type))

            print ("using globs for location")
            for potential_loc in potential_locs:

                ###### use the above with sample name alone - i.e. need to make sure legacy and external samples are first AND visible...
                ###### or we get the newer sapmle if sequence more than once - it doesn't use sample_type/capture_kit or anything!?!
                potential_path = '{}/{}/*[0-9xXyY]'.format(potential_loc,sample_name)
                print("globbing for legacy {}, {}".format(potential_loc,sample_name))
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

            print ("using query for location")
            query = ("SELECT DISTINCT SEQSATALOC,FCILLUMID FROM Flowcell f "
                "JOIN Lane l ON f.FCID=l.FCID "
                "JOIN prepT p ON l.prepID=p.prepID "
                "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
                "AND l.prepID={0} AND f.fail=0 "
                "GROUP BY f.fcid"
                ).format(prepid)

            seqsatalocs = run_query(query,database)
            if len(seqsatalocs)==0:
                raise ValueError("no db registered fastq locations")

            print('we have the following locations from db to check : {}'.format(seqsatalocs))

            print(seqsatalocs)
            for flowcell in seqsatalocs:

                if flowcell['SEQSATALOC'] == '':
                    msg = ('Sample {} on flowcell {} is missing seqsataloc! Still sequencing?' ).format(sample_name,flowcell['FCILLUMID'])
                    raise ValueError(msg)

                for potential_loc in potential_locs:

                    ##### use the above along with sample name and flowcell id
                    potential_path = '{}/{}/{}'.format(potential_loc,sample_name,flowcell['FCILLUMID'])
                    print("globbing for non-legacy {}".format(potential_path))
                    fastq_loc = glob(potential_path)
                    if fastq_loc != []:
                        print('FOUND: {}'.format(potential_path))
                        ##### do we only get the first one and not check for ambiguity?!?
                        for folder in fastq_loc:
                            found_locs.append(folder)
                        break

                print('sample={} and loc'.format(sample_name,fastq_loc))
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
            raise ValueError('Did not find fastq mate pair for: {}!'.format(loc))
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
            if capture_kit == '':
                raise ValueError("must supply valid capture kit")
            self.metadata['capture_kit'] = capture_kit
            self.metadata['bed_file_loc'] = get_bed_file_loc(database,self.metadata['capture_kit'])
        else:
            self.metadata['capture_kit'] = ''
            self.metadata['bed_file_loc'] = '/nfs/goldsteindata/refDB/captured_regions/Build37/65MB_build37/SeqCap_EZ_Exome_v3_capture.bed'
        self.metadata['sample_id'] = get_sample_id(database,self.metadata)
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
