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
from glob import glob
from utilities import get_connection
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
    query = ("SELECT prepid FROM pseudo_prepid WHERE pseudo_prepid={0}"
            ).format(sample["pseudo_prepid"])
    #print(query)
    prepids = run_query(query,database)
    prepids = [x['prepid'] for x in prepids]
    return prepids

def get_priority(database,sample):
    query = ("SELECT priority "
            "FROM SampleT s "
            "JOIN prepT p ON s.CHGVID=p.CHGVID "
            "JOIN SeqType st ON st.prepid=p.prepid "
            "WHERE p.CHGVID='{sample_name}' "
            "AND st.seqtype='{sample_type}' "
            "AND p.exomekit='{capture_kit}' "
            "ORDER BY priority ASC "
            "LIMIT 1"
            ).format(sample_name=sample['sample_name'],
                    sample_type=sample['sample_type'],
                    capture_kit=sample['capture_kit'])
    #print(query)
    priority = run_query(query,database)[0]['priority']
    return priority

def get_pseudo_prepid(database,sample):
    query = ("SELECT DISTINCT pseudo_prepid "
            "FROM pseudo_prepid pp "
            "JOIN prepT p ON pp.prepid=p.prepid "
            "JOIN SeqType s ON s.prepid=p.prepid "
            "WHERE CHGVID='{sample_name}' "
            "AND seqtype='{sample_type}' "
            "AND p.exomekit='{capture_kit}' "
            "AND failedprep=0"
            ).format(sample_name=sample['sample_name'],
                    sample_type=sample['sample_type'],
                    capture_kit=sample['capture_kit'])
    pseudo_prepid = run_query(query)[0]['pseudo_prepid']
    return pseudo_prepid

def get_bed_file_loc(database,capture_kit):
    """Retrieve the bed file location for a given capture_kit name"""
    query = (("SELECT region_file_lsrc FROM captureKit WHERE name='{0}' "
        "and chr = 'all'"
        ).format(capture_kit))
    bed_file_loc = run_query(query,database)[0]['region_file_lsrc']
    return bed_file_loc

def get_fastq_loc(database, sample):
    locs = []
    #correcting the downstream consequences of "custom capture" as the sample_type
    corrected_sample_type = sample['sample_type'].upper().replace(" ","_")
    sample_name=sample['sample_name']
    #pdb.set_trace()
    for prepid in sample["prepid"]:
        query = ("SELECT seqsataloc,FCIllumID FROM Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            "AND l.prepID={0} AND f.fail=0 "
            "GROUP BY f.fcid"
            ).format(prepid)
        seqsatalocs = run_query(query,database)
        """For cases where there is not flowell information the sequeuncing
        will have to found be manually.  There will be two types of samples
        that under this catergory: Old Casava 1.8 sample (pre-seqDB) and
        Casava 1.7 samples sequenced by the Illumina GAII."""
        if seqsatalocs:
            """For externally submitted sample they have a fake Flowcell entry
            in the database.  The Flowcell.FCIllumID begins wt'X' always.
            Typically when processing the external sample each read group is
            enumerated 1,2,3...etc.

            Secondly when external samples are archived sometimes the FCIllumID
            is preserved otherwise its enumerated."""
            #print(sample)
            if seqsatalocs[0]['FCIllumID'][0] == 'X':
                print("Externally submitted samples")
                if 'SRR' in sample_name: #specifically for SRR samples
                    loc = ('/nfs/seqscratch10/SRA/{}/*X[XY]').format(sample_name)
                    fastq_loc = glob(loc)
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif 'CGNDHDA' in sample_name:
                    fastq_loc = glob(('/nfs/seqscratch09/tx_temp/tx_2390/CGND_11418-fastq/Project_CGND_11418_B01_GRM_WGS.2016-03-30/{}/1').format(sample_name))
                    if fastq_loc == []:
                        fastq_loc = glob(('/nfs/seqscratch09/tx_temp/tx_2390/CGND_11645-fastq/Project_CGND_11645_B01_GRM_WGS.2016-08-17/{}/*[XY]').format(sample_name))
                    if fastq_loc ==[]:
                        fastq_loc = glob(('/nfs/seqscratch12/tx_temp/tx_2390/CGND_11911-fastq/Project_CGND_11911_B01_GRM_WGS.2016-09-29/{}/*[XY]').format(sample_name))
                    if fastq_loc ==[]:
                        fastq_loc = glob(('/nfs/seqscratch09/tx_temp/tx_2390/CGND_11772-fastq/Project_CGND_11772_B01_GRM_WGS.2016-06-09/{}/*[XY]').format(sample_name))
                    if fastq_loc ==[]:
                        fastq_loc = glob(('/nfs/seqscratch09/tx_temp/tx_2390/CGND_11772-fastq/Project_CGND_11772_B02_GRM_WGS.fastq_bam.2017-04-29{}/*[XY]').format(sample_name))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/igmdata0[0-9]/{}/{}/[0-9]'
                    ).format(corrected_sample_type,sample_name)) != []: #external igmdata samples wt enumerated folders
                    fastq_loc = glob(('/nfs/igmdata0[0-9]/{}/{}/[0-9]'
                                ).format(corrected_sample_type,sample_name))
                    if fastq_loc:
                        for flowcell in fastq_loc:
                            locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/igmdata0[0-9]/{}/{}/*X[XY]'
                    ).format(corrected_sample_type,sample_name)) != []: #external igmdata samples wt actual flowcell name
                    fastq_loc = glob(('/nfs/igmdata0[0-9]/{}/{}/*X[XY]'
                            ).format(corrected_sample_type,sample_name))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/fastq1[0-9]/{}/{}/[0-9]'
                    ).format(corrected_sample_type,sample_name)) != []: #external fastq16 samples wt enumerated folders
                    fastq_loc = glob(('/nfs/fastq1[0-9]/{}/{}/[0-9]'
                                ).format(corrected_sample_type,sample_name))
                    if fastq_loc:
                        for flowcell in fastq_loc:
                            locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/fastq1[0-9]/{}/{}/*X[XY]'
                    ).format(corrected_sample_type,sample_name)) != []: #external fastq16 samples wt actual flowcell name
                    fastq_loc = glob(('/nfs/fastq1[0-9]/{}/{}/*X[XY]'
                            ).format(corrected_sample_type,sample_name))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/[0-9]'
                    ).format(sample_name)) != []: #external seqsata samples wt enumerated folders
                    fastq_loc = glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/[0-9]'
                        ).format(sample_name))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/*X[XY]'
                    ).format(sample_name)) != []: #external seqsata samples wt actual flowcell name
                    fastq_loc = glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/*X[XY]'
                        ).format(sample_name))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                else:
                    raise FastqError('{} Fastq files not found'.format(sample_name))
            else:
                print("Internally sequenced sample in sequenceDB")
                for flowcell in seqsatalocs:
                    print(flowcell)
                    if 'igmdata' in flowcell['seqsataloc'] or 'fastq' in flowcell['seqsataloc']: #igmdata## or fastq_temp##
                        fastq_loc = ('/nfs/{0}/{1}/{2}/{3}'
                                ).format('seqscratch_ssd',corrected_sample_type,
                                        sample_name,flowcell['FCIllumID'])
                        if glob(fastq_loc) == []:
                            fastq_loc = ('/nfs/{0}/{1}/{2}/{3}'
                                    ).format(flowcell['seqsataloc'],corrected_sample_type,
                                            sample_name,flowcell['FCIllumID'])
                    elif 'seqsata' in flowcell['seqsataloc']: #seqsata## drives
                        fastq_loc = ('/nfs/{0}/seqfinal/whole_genome/{1}/{2}'
                                ).format(flowcell['seqsataloc'],sample_name,flowcell['FCIllumID'])
                    else:
                        raise FastqError('{} Fastq files not found'.format(sample_name))
                    read = glob(os.path.realpath(fastq_loc))
                    if read != []:
                        locs.append(os.path.realpath(fastq_loc))
                    elif glob(('/nfs/igmdata01/{}/{}/{}'
                              ).format(corrected_sample_type,sample_name,flowcell['FCIllumID'])) != []:
                        """ this is a hack for samples suppose to be on a
                            RO avere volume downloaded from AWS to igmdata01
                            instead"""
                        fastq_loc = '/nfs/igmdata01/{}/{}/{}'.format(corrected_sample_type,sample_name,flowcell['FCIllumID'])
                        locs.append(os.path.realpath(fastq_loc))

                    else:
                        # For samples in the database but stored on the quantum and 
                        # have not had their location properly restored
                        fastq_loc = glob('/nfs/stornext/seqfinal/casava1.8/whole_{0}/{1}/*X[XY]'.format(
                            corrected_sample_type.lower(),sample_name))
                        if fastq_loc:
                            for path in fastq_loc:
                                locs.append(os.path.realpath(path))

        else:
            print("Internally sequenced sample  but not in SequenceDB")
            if glob('/nfs/fastq_temp2/{}/{}/*[0-9XY]'.format(corrected_sample_type,sample_name)):
                fastq_loc = glob('/nfs/fastq_temp2/{}/{}/*[0-9XY]'.format(corrected_sample_type,sample_name))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/igmdata01/{}/{}/*X[XY]'.format(corrected_sample_type,sample_name)):
                fastq_loc = glob('/nfs/fastq16/{}/{}/*X[XY]'.format(corrected_sample_type,sample_name))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/igmdata01/{}/{}/[0-9]'.format(corrected_sample_type,sample_name)) != []:
                fastq_loc = glob('/nfs/fastq16/{}/{}/[0-9]'.format(corrected_sample_type,sample_name))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/fastq16/{}/{}/*X[XY]'.format(corrected_sample_type,sample_name)):
                fastq_loc = glob('/nfs/fastq16/{}/{}/*X[XY]'.format(corrected_sample_type,sample_name))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/fastq16/{}/{}/[0-9]'.format(corrected_sample_type,sample_name)) != []:
                fastq_loc = glob('/nfs/fastq16/{}/{}/[0-9]'.format(corrected_sample_type,sample_name))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/stornext/seqfinal/casava1.8/whole_{0}/{1}/*X[XY]'.format(
                corrected_sample_type.lower(),sample_name)):
                fastq_loc = glob('/nfs/stornext/seqfinal/casava1.8/whole_{0}/{1}/*X[XY]'.format(corrected_sample_type.lower(),sample_name))
                for flowcell in fastq_loc:
                    locs.append(flowcell)
            elif glob('/nfs/seqsata*/seqfinal/whole_genome/{}/*X[XY]'.format(sample_name)) != []:
                fastq_loc = glob('/nfs/seqsata*/seqfinal/whole_genome/{}/*X[XY]'.format(sample_name))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            else:
                raise FastqError('{} Fastq files not found'.format(sample_name))

    """For samples in the database we can exclude any samples that only have
    R1 data however for sampels that predate the database we have to manually
    check for R2 existance"""
    locs = check_fastq_locs(list(set(locs)))
    return locs

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
        if lane and isInternalSample(lane) == True:
            lanes.append((lane[0]['lanenum'],lane[0]['FCIllumID'],lane[0]['prepID']))
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
            #Genome samples are set using the most current capture kit for any case which requires a target region.
            self.metadata['bed_file_loc'] = '/nfs/goldsteindata/refDB/captured_regions/Build37/65MB_build37/SeqCap_EZ_Exome_v3_capture.bed'
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
