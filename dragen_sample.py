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
import json

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
    # query = ("SELECT prepid FROM prepT WHERE p_prepid={0} and failedprep=0"
    query = ("SELECT prepid FROM prepT WHERE experiment_id={} and (failedprep=0 or failedprep>=100)").format(sample["experiment_id"])
    prepids = run_query(query,database)
    prepids = [x['prepid'] for x in prepids]
    return prepids

def get_sample_id(database,sample):
    query = """SELECT distinct(sample_id) FROM prepT JOIN Experiment ON prepT.experiment_id=Experiment.id WHERE experiment_id={} and (failedprep=0 or failedprep>=100)
            """.format(sample['experiment_id'])
    results = run_query(query,database)
    check_number_query_results(results,1)
    return results[0]['sample_id']

def get_priority(database,sample):
    query = ("SELECT PRIORITY "
            "FROM prepT "
            "WHERE P_PREPID ={pseudo_prepid} and (failedprep=0 or failedprep>=100) "
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
            "AND (failedprep=0 or failedprep>=100)"
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

def get_fastq_loc(database, sample, rarp ):

    total = 0
    
    ##### people have been breaking the dir/file name combos wrt., the Lane.data json field so simpler to downgrade and patch up...
    NEW_GLOBS_FOR_DOWNGRADING_PATCH = []

    for prepid in sample["prepid"]:

        query = ("SELECT rg_status,lanenum,fcillumid, data FROM Flowcell f JOIN Lane l ON f.FCID=l.FCID JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            "AND l.prepID={0} AND f.fail=0 and (failedprep=0 or failedprep>=100) "
            # "GROUP BY f.fcid"
            ).format(prepid)
        info = run_query(query,database)
        for d in info:
            total = total+1
            if d['data'] != None:
                j=(json.loads(d['data']))
                if 'bam' in j:
                    bam=j['bam']
                    rarp.append( [ d['lanenum'], d['fcillumid'], '{}/{}'.format(bam['path']['scratch'],bam['basename']) ] )
                else:
                    raise ValueError("since this only happens post alignment invocation this should be possible")
                NEW_GLOBS_FOR_DOWNGRADING_PATCH.append( 
                        '/'.join( j['fastq']['path']['archive'].split('/')[:-3] )
                )
    ############# since we are retreiving multiple prepids in init we can actually do things this way!?!

    if sample['sample_name'][0:6]=="igm160" and ( int(sample['experiment_id'])>=180922 and int(sample['experiment_id'])<=182983 ):
        print("> SUMMARY : this should have been were we re-staged db bams for merging with topup SEs but instead is a BS hack for Gilead")
        rarp.clear() # del rarp[:]
    elif len(rarp)>0:
        if len(rarp)!=total:
            rarp.clear() # del rarp[:]

#            if len(sample["prepid"])>1:
#                query = ("update dragen_sample_metadata set is_merged = 80121 where pseudo_prepid = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                query = ("update Experiment set is_released = 'legacy_multiprep_error' where id = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                exit(1)
#            else:
#                query = ("update dragen_sample_metadata set is_merged = 80122 where pseudo_prepid = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                query = ("update Experiment set is_released = 'legacy_hybrid_error' where id = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                run_query(query,database)
#                exit(1)
            ########## implement hybrid procedure directly with Experiment table along with metrics - i.e. at each bam finishing check if all bams for

        else:
            return []

    else:
        print("> SUMMARY : LEGACY")

    found_locs = []
    #correcting the downstream consequences of "custom capture" as the sample_type
    corrected_sample_type = sample['sample_type'].upper().replace(" ","_")
    sample_name=sample['sample_name']

    p_c=0
    for prepid in sample["prepid"]:
        p_c+=1
        print ("> (get_fastq_loc) CHECKING PREPID = '{}', '{}' of '{}'".format(prepid,p_c,len(sample['prepid'])))
        external_or_legacy_flag = is_external_or_legacy_sample(prepid,database)
        ####### this is all an abomination!?!
        potential_locs = []

        if corrected_sample_type != 'GENOME':
            potential_locs = ['/nfs/seqscratch_ssd/{}/'.format(corrected_sample_type)]

        if sample['pseudo_prepid']<20000:
            potential_locs.insert(0,'/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type))

        potential_locs.extend([
            '/nfs/archive/p2018/FASTQ/EXOME','/nfs/archive/p2018/FASTQ/EXOME',
            '/nfs/fastq_temp/{}/'.format(corrected_sample_type),
            '/nfs/seqscratch_ssd/transfer_finished_symlinks/{}/'.format(corrected_sample_type),
            '/nfs/fastq_temp2/{}/'.format(corrected_sample_type),
            '/nfs/seqscratch*/tx_temp/tx_*/',
            '/nfs/igmdata01/{}/'.format(corrected_sample_type),
            '/nfs/fastq1[568]/{}/'.format(corrected_sample_type),
        ])

        if external_or_legacy_flag == True:

            print(" > pure, evil name-coliding nonsense - in general will mix up multiply sequenced samples and truncate multiprep...")
            if sample["prepid"][0]<50000:
                # potential_locs = ['/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type)]
                ##### needs to be compatible with other old samples - why aren't they run?!?!
                potential_locs.insert(0,'/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type))

            for potential_loc in potential_locs:

                potential_path = '{}/{}/*[0-9xXyY]'.format(potential_loc,sample_name)
                fastq_loc = glob(potential_path)
                if fastq_loc != []:
                    for folder in fastq_loc:
                        found_locs.append(folder)
                    break

            if fastq_loc == []:
                raise FastqError('{} Fastq files not found'.format(sample_name))
        else:

            query = ("SELECT DISTINCT SEQSATALOC,FCILLUMID,SAMPLE_TYPE FROM Flowcell f "
                "JOIN Lane l ON f.FCID=l.FCID "
                "JOIN prepT p ON l.prepID=p.prepID "
                "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
                "AND l.prepID={0} AND f.fail=0 and (failedprep=0 or failedprep>=100) "
                "GROUP BY f.fcid"
            ).format(prepid)

            flowcell_seqsatalocs = run_query(query,database)

            if len(flowcell_seqsatalocs)==0:

                #### do single joined update?!?
                query = ("update dragen_sample_metadata set is_merged = -2 where experiment_id = {}".format(sample['experiment_id']))
                run_query(query,database)
                query = ("update prepT set is_released = 0, status = 'Failed/Low-Qual Sample; Has no passing sequence events - will require deprecation', status_time = UNIX_TIMESTAMP(CURRENT_TIMESTAMP()) where experiment_id = {}".format(sample['experiment_id']))
                run_query(query,database)
                query = ("update Experiment set is_released = 'release_rejected' where id = {}".format(sample['experiment_id']))
                run_query(query,database)
                print("whatever, just don't let exception handling catch this one")
                os._exit(1)

                raise ValueError("prepid {} of experiment {} does not appear to have been sequenced - sequence or fail it to proceed".format(
                  prepid,sample['experiment_id']))

            print('we have the following locations from db to check : {}'.format(flowcell_seqsatalocs))

            from pprint import pprint as pp

            f_c=0
            for flowcell in flowcell_seqsatalocs:

                f_c+=1

                if flowcell['SEQSATALOC'] == '':

                    msg = ('Sample {} on flowcell {} is missing seqsataloc! Still sequencing?' ).format(sample_name,flowcell['FCILLUMID'])
                    raise ValueError(msg)

                if sample['sample_name'][0:6]=="ALSNEU":
                    potential_locs.insert(0,'/nfs/{}/{}/'.format(flowcell['SEQSATALOC'],flowcell['SAMPLE_TYPE'].upper()))
                    pp(potential_locs)

                if len(NEW_GLOBS_FOR_DOWNGRADING_PATCH):

                    ##### duh, need to put list into list!?!?
                    # potential_locs.append(NEW_GLOBS_FOR_DOWNGRADING_PATCH)
                    potential_locs.extend(NEW_GLOBS_FOR_DOWNGRADING_PATCH)

                for potential_loc in potential_locs:

                    # potential_path = '{}/{}/{}'.format(potential_loc,sample_name,flowcell['FCILLUMID'])
                    potential_path = '{}/{}/{}/*.gz'.format(potential_loc,sample_name,flowcell['FCILLUMID'])
                    if sample['sample_name'][0:6]=="ALSNEU":
                        print("This is offensitve!?!? : using even nastier glob!?!")
                        potential_path = '{}/{}*/{}/*.gz'.format(potential_loc,sample_name,flowcell['FCILLUMID'])

                    fastq_loc = glob(potential_path)
                    if fastq_loc != []:

                        # potential_path.replace('/*.gz','')
                        potential_path=potential_path[0:-4]

                        ##### do we only get the first one and not check for ambiguity?!?
                        for folder in fastq_loc:
                            folder=os.path.dirname(folder)
                            found_locs.append(folder)
                        break

                print('sample={} and loc'.format(sample_name,fastq_loc))

                if fastq_loc == []:

                    # potential_path = '/nfs/{}/{}/{}/{}'.format(flowcell['SEQSATALOC'],flowcell['SAMPLE_TYPE'].upper(),sample_name,flowcell['FCILLUMID'])
                    potential_path = '/nfs/{}/{}/{}/{}/*gz'.format(flowcell['SEQSATALOC'],flowcell['SAMPLE_TYPE'].upper(),sample_name,flowcell['FCILLUMID'])

                    fastq_loc = glob(potential_path)
                    if fastq_loc != []:
                        potential_path=potential_path[0:-4]
                        # potential_path.replace('/*.gz','')
                        for folder in fastq_loc:
                            folder=os.path.dirname(folder)
                            found_locs.append(folder)
                    else:
                        raise FastqError('{} Fastq files not found'.format(sample_name))


        #### evil!?!
        found_locs=list(set(found_locs))
        print("prepid = '{}' has '{}'".format(prepid,found_locs))
    """For samples in the database we can exclude any samples that only have
    R1 data however for sampels that predate the database we have to manually
    check for R2 existance"""

    ###### oh, just horrible
    found_locs=list(set(found_locs))
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
            raise ValueError('Did not find fastq mate pair for: {} i.e. {}/*R2_[0-9]*fastq* !!!'.format(loc,loc))
    return valid_locs


def get_lanes(database,sample):

    lanes = []

    dist_fc=[]
    ffs=False
    for prepid in sample['prepid']:

        print("> (get_lanes) RUNNING PREPID = '{}'".format(prepid))

        query = ("SELECT DISTINCT lanenum,FCIllumID,p.prepID from Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            ###### don't require pipelinecomplete/complete as want to see them
            "and f.fail=0 "
            "AND l.prepID={0} AND (p.failedPrep=0 or failedprep>=100) "
        ).format(prepid)

        lane = run_query(query,database)

        leg_or_ext = is_external_or_legacy_sample(prepid,database)

        if lane and leg_or_ext == False:
            for entry in lane:
                lanes.append((entry['lanenum'],entry['FCIllumID'],entry['prepID']))
                dist_fc.append(entry['FCIllumID'])
        else: #In case of no database entry or external sample
            ffs=True
            for flowcell in sample['fastq_loc']:
                print("this has got to be a joke : {}".format(flowcell))
                fastqs = glob(flowcell + '/*fastq.gz')
                lane = []
                for fastq in fastqs:
                    lane_num = fastq[fastq.find('L00')+3]
                    lane.append(lane_num)
                lane_nums = set(sorted(lane))
                for lane_num in lane_nums:
                    lanes.append((lane_num,flowcell.split('/')[-1],prepid))
                    dist_fc.append(flowcell.split('/')[-1])
        ###### cretin - this makes it a tuple!?! what's it even for?!?
        # lanes = (lanes,)
    if ffs and len(sample['prepid'])>len(set(dist_fc)):
        from pprint import pprint as pp
        pp(sample['prepid'])
        pp(set(dist_fc))
        # query = ("update dragen_sample_metadata set is_merged = 80140 where experiment_id = {}".format(sample['experiment_id']))
        # run_query(query,database)
        # raise ValueError("legacy globbing procedure is unable to locate distinct flowcells")
    ### evil?!?
    lanes=list(set(lanes))
    # pp(lanes) 
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

        #### will need to allow fro multiple preps here too!?!
        from pprint import pprint as pp
        pp(self) ; pp(sample_name) ; pp(sample_type) ; pp(pseudo_prepid) ; pp(capture_kit)

        self.metadata = {}
        self.metadata['sample_name'] = sample_name
        self.metadata['sample_type'] = sample_type.lower()
        self.metadata['pseudo_prepid'] = pseudo_prepid
        self.metadata['experiment_id'] = pseudo_prepid
        ############
        xxx=run_query("select experiment_id from dragen_sample_metadata where pseudo_prepid = {}".format(pseudo_prepid),database);
        if len(xxx)!=1 or xxx[0]['experiment_id']!=pseudo_prepid:
            raise ValueError("experiment_id legacy pseudo_prepid mismatch")
        ############
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
        tracked_files=[]; # {}
        self.metadata['fastq_loc'] = get_fastq_loc(database,self.metadata,tracked_files)
        from pprint import pprint as pp
        pp(self.metadata['fastq_loc'])
        if len(self.metadata['fastq_loc'])==0 and len(tracked_files)==0:
            raise ValueError("this is very, very, very wrong")
        elif len(tracked_files)>0:
            if len(self.metadata['fastq_loc'])>0:
                raise ValueError("these are mutually exclusive")
            self.metadata['bams'] = tracked_files
        self.metadata['lane'] = get_lanes(database,self.metadata)
        from pprint import pprint as pp
        pp(self.metadata['lane'])
        print("> handing back silly object thing")

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
