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
    # rarp['arse']='fdsafds';
    # https://stackoverflow.com/questions/986006/how-do-i-pass-a-variable-by-reference
    # If you pass a mutable object into a method, the method gets a reference to that same object and you can mutate it to your heart's delight, but if you rebind the reference in the method, the outer scope will know nothing about it,
    # hence 
    #   rarp.append('fdfd')
    # not 
    #   rarp='fff';
    # print ("select count(1) count, data using query for location")
    # new_fq = 0
    total = 0

    
    ##### people have been breaking the dir/file name combos wrt., the Lane.data json field so simpler to downgrade and patch up...
    NEW_GLOBS_FOR_DOWNGRADING_PATCH = []
    

    print("we have '{}' preps for this experiment".format(len(sample['prepid'])))

    for prepid in sample["prepid"]:
        ###### grab absolutely ALL data for sample - BUT - needs to be changed to join via experiment_id and not prepid!?!
        ###### grab absolutely ALL data for sample - BUT - needs to be changed to join via experiment_id and not prepid!?!
        ###### grab absolutely ALL data for sample - BUT - needs to be changed to join via experiment_id and not prepid!?!
        query = ("SELECT rg_status,lanenum,fcillumid, data FROM Flowcell f JOIN Lane l ON f.FCID=l.FCID JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            # let's do it with or without bam - i.e. just use the fastq from scratch location directly
            # "WHERE data->'$.bam' is not null and (FailR1 IS NULL and FailR2 IS NULL) "
######################### strictly should be checking is_released here!?!
######################### strictly should be checking is_released here!?!
######################### strictly should be checking is_released here!?!
            "AND l.prepID={0} AND f.fail=0 and (failedprep=0 or failedprep>=100) "
            # "GROUP BY f.fcid"
            ).format(prepid)
        info = run_query(query,database)
        # from pprint import pprint as pp
        for d in info:
            total = total+1
            # pp(d['data'])
            ##### should we allow for fastq direct - i.e. align directly when in hurry?!?
            if d['data'] != None:
                j=(json.loads(d['data']))
                if 'bam' in j:
                    print(d['rg_status'])
                    bam=j['bam']
                    print("new flavour with {} @ {}".format(d['rg_status'],j['bam']['path']['scratch']))
                    # pp(bam)
                    rarp.append( [ d['lanenum'], d['fcillumid'], '{}/{}'.format(bam['path']['scratch'],bam['basename']) ] )
                else:
                    raise ValueError("since this only happens post alignment invocation this should be possible")
                NEW_GLOBS_FOR_DOWNGRADING_PATCH.append( 
                        '/'.join( j['fastq']['path']['archive'].split('/')[:-3] )
                )
        # pp(info)
    ############# since we are retreiving multiple prepids in init we can actually do things this way!?!
    ############# we can also hack it to allow for hybrids to work via ond procedure
    # if ( len(sample['prepid'])>1 and len(rarp)==0 ) or len(rarp)>0:
    if sample['sample_name'][0:6]=="igm160" and ( int(sample['experiment_id'])>=181922 and int(sample['experiment_id'])<=181983 ):
        print("> SUMMARY : F0CKW1T, CUN7 HANDLING. this should have been were we re-staged db bams for merging with topup SEs but instead is a BS hack for Gilead")
        from pprint import pprint as pp
        pp(sample)
        print ("cunts")
    elif len(rarp)>0:
        if len(rarp)!=total:
            print("> SUMMARY : HYBRID : JUST HACK THIS BY ADDING IN KNOWN PATH TO NASTY GLOB STUFF BELOW!?!?")
            print("annoying {}".format(len(rarp)))
            # this actually causes a purely local change :
            # rarp=[] # force legacy procedure
            rarp.clear() # del rarp[:]
            print("annoying {}".format(len(rarp)))
#            from pprint import pprint as pp
#            pp(sample)
#            if len(sample["prepid"])>1:
#                query = ("update dragen_sample_metadata set is_merged = 80121 where pseudo_prepid = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                query = ("update Experiment set is_released = 'legacy_multiprep_error' where id = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                ### or fastq/bam ...
#                print("We aren't allowing legacy_multiprep_error mixes - please patch all SEs")
#                exit(1)
#                raise ValueError("We aren't allowing hybrids or bam/fastq mixes at this time (should we make release require bam - we shouldn't get here without a bam present anyway)")
#            else:
#                query = ("update dragen_sample_metadata set is_merged = 80122 where pseudo_prepid = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                query = ("update Experiment set is_released = 'legacy_hybrid_error' where id = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                query = ("update prepT set status = 'Not allowing legacy-hybrid mixes - please patch all SEs' where experiment_id = {}".format(sample['pseudo_prepid']))
#                run_query(query,database)
#                print("We aren't allowing legacy_hybrid_error mixes - please patch all SEs")
#                exit(1)
#                raise ValueError("We aren't allowing hybrids or bam/fastq mixes at this time (should we make release require bam - we shouldn't get here without a bam present anyway)")
            ########## implement hybrid procedure directly with Experiment table along with metrics - i.e. at each bam finishing check if all bams for
            ########## expt_id are mapped - if they're old trigger them for dragen_align_se else directly merge them and do stats...!?!?
        else:
            print("> SUMMARY : NEW PRE-MAPPED PROCEDURE")
            return []
    else:
        print("> SUMMARY : LEGACY")
        print("\n\t>legacy procedure : total={}, registered_fastq={}, preps={}\n\n".format(total,len(rarp),len(sample['prepid'])))

    found_locs = []
    #correcting the downstream consequences of "custom capture" as the sample_type
    corrected_sample_type = sample['sample_type'].upper().replace(" ","_")
    sample_name=sample['sample_name']

    p_c=0
    for prepid in sample["prepid"]:
        p_c+=1
        print ("> (get_fastq_loc) CHECKING PREPID = '{}', '{}' of '{}'".format(prepid,p_c,len(sample['prepid'])))
        external_or_legacy_flag = is_external_or_legacy_sample(prepid,database)
        # pp(sample)
        ####### this is all an abomination!?!
        potential_locs = []

        if corrected_sample_type != 'GENOME':
            potential_locs = ['/nfs/seqscratch_ssd/{}/'.format(corrected_sample_type)]

        if sample['pseudo_prepid']<20000:
            potential_locs.insert(0,'/nfs/seqscratch_ssd/dsth/APPALING_ALS_ISSUE/'.format(corrected_sample_type))

        potential_locs.extend([
            '/nfs/fastq_temp/{}/'.format(corrected_sample_type),
            '/nfs/seqscratch_ssd/transfer_finished_symlinks/{}/'.format(corrected_sample_type),
            # '/nfs/fastq_temp5/tx_temp/tx_3036/',
            # '/nfs/fastq_temp2/{}_ssd/'.format(corrected_sample_type),       ##### truelly evil!?!
            '/nfs/fastq_temp2/{}/'.format(corrected_sample_type),
            '/nfs/seqscratch*/tx_temp/tx_*/',
            '/nfs/igmdata01/{}/'.format(corrected_sample_type),
            '/nfs/fastq1[568]/{}/'.format(corrected_sample_type),
            # '/nfs/stornext/seqfinal/casava1.8/whole_{}/'.format(corrected_sample_type),
            # '/nfs/fastq_temp5/tx_2995/Project_CGND_12737_B01_GRM_WGS.fastq.2017-08-29/'
            ###### this is causing so many issues that better to remove and add back later!?!
            ###### i.e. just check ALL fastq error!?!
            # '/nfs/seqsata*/seqfinal/whole_genome/'
        ])

        if external_or_legacy_flag == True:

            print(" > pure, evil name-coliding nonsense - in general will mix up multiply sequenced samples and truncate multiprep...")
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
                    print('LEGACY FOUND: {}'.format(potential_path))
                    for folder in fastq_loc:
                        found_locs.append(folder)
                    break

            if fastq_loc == []:
                raise FastqError('{} Fastq files not found'.format(sample_name))
        else:

            print ("> using, but largely ignoring, query for location")

            query = ("SELECT DISTINCT SEQSATALOC,FCILLUMID,SAMPLE_TYPE FROM Flowcell f "
                "JOIN Lane l ON f.FCID=l.FCID "
                "JOIN prepT p ON l.prepID=p.prepID "
                "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
                "AND l.prepID={0} AND f.fail=0 and (failedprep=0 or failedprep>=100) "
                "GROUP BY f.fcid"
            ).format(prepid)

            print("using\n{}".format(query))
            flowcell_seqsatalocs = run_query(query,database)
            if len(flowcell_seqsatalocs)==0:
                query = ("update dragen_sample_metadata set is_merged = 80130 where experiment_id = {}".format(sample['experiment_id']))
                run_query(query,database)
                raise ValueError("prepid {} of experiment {} does not appear to have been sequenced - sequence or fail it to proceed".format(
                  prepid,sample['experiment_id']))

            print('we have the following locations from db to check : {}'.format(flowcell_seqsatalocs))

            from pprint import pprint as pp

            # print(seqsatalocs)
            f_c=0
            for flowcell in flowcell_seqsatalocs:
                print(flowcell)
                f_c+=1
                print (" > CHECKING FLOWCELL = '{}', '{}' of '{}'".format(flowcell['FCILLUMID'],f_c,len(flowcell_seqsatalocs)))

                ######### FUCK ME!?! this only actually checks the seqsata location when NONE of the hard-coded locations work?!?
                if flowcell['SEQSATALOC'] == '':

                    msg = ('Sample {} on flowcell {} is missing seqsataloc! Still sequencing?' ).format(sample_name,flowcell['FCILLUMID'])
                    raise ValueError(msg)

                # pp(sample)
                if sample['sample_name'][0:6]=="ALSNEU":
                    print("This is offensitve!?!? : prepending 'actual' archive location!?!")
                    potential_locs.insert(0,'/nfs/{}/{}/'.format(flowcell['SEQSATALOC'],flowcell['SAMPLE_TYPE'].upper()))
                    pp(potential_locs)

                if len(NEW_GLOBS_FOR_DOWNGRADING_PATCH):
                    print("we patch in yet more paths for hybrids?!? : {}".format(NEW_GLOBS_FOR_DOWNGRADING_PATCH))
                    ##### duh, need to put list into list!?!?
                    # potential_locs.append(NEW_GLOBS_FOR_DOWNGRADING_PATCH)
                    potential_locs.extend(NEW_GLOBS_FOR_DOWNGRADING_PATCH)

                for potential_loc in potential_locs:

                    ##### use the above along with sample name and flowcell id

                    ####### this BS always gets screwed by empty dirs!?!
                    # potential_path = '{}/{}/{}'.format(potential_loc,sample_name,flowcell['FCILLUMID'])
                    potential_path = '{}/{}/{}/*.gz'.format(potential_loc,sample_name,flowcell['FCILLUMID'])
                    if sample['sample_name'][0:6]=="ALSNEU":
                        print("This is offensitve!?!? : using even nastier glob!?!")
                        ###### need to account for 'rep' and '_LIMS_CHANGE_OVER_DEPRECATION' suffixes?!?
                        potential_path = '{}/{}*/{}/*.gz'.format(potential_loc,sample_name,flowcell['FCILLUMID'])

                    print("globbing for non-legacy {}".format(potential_path))
                    fastq_loc = glob(potential_path)
                    if fastq_loc != []:

                        ####### dipshit!?!?
                        # potential_path.replace('/*.gz','')
                        potential_path=potential_path[0:-4]
                        print('SEQSATA FOUND: {}'.format(potential_path))

                        ########## dipshit!?!
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
                        print('FOUND: {}'.format(potential_path))
                        for folder in fastq_loc:
                            ### fuckwit
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
    ######### the 'legacy' ones may well find the non-legacy ones and end up with junk - i.e. like the exome/capture messup?!?
    found_locs=list(set(found_locs))
    final_locs = check_fastq_locs(found_locs)
    print("final list of unique locations for all preps = '{}'".format(final_locs))
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
    # print("wtf='{}'".format(type(lanes)))

    dist_fc=[]
    ffs=False
    for prepid in sample['prepid']:

        print("> (get_lanes) RUNNING PREPID = '{}'".format(prepid))

        query = ("SELECT DISTINCT lanenum,FCIllumID,p.prepID from Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            "AND l.prepID={0} AND (p.failedPrep=0 or failedprep>=100) "
        ).format(prepid)

        lane = run_query(query,database)

        from pprint import pprint as pp
        # print('raw=')
        # pp(lane) 

        leg_or_ext = is_external_or_legacy_sample(prepid,database)

        # print("wtf='{}'".format(type(lanes)))
        if lane and leg_or_ext == False:
            for entry in lane:
                # print("wtf='{}'".format(type(lanes)))
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
    print("we have {} distinct prepids and {} distinct fc".format(len(sample['prepid']),len(set(dist_fc))))
    if ffs and len(sample['prepid'])>len(set(dist_fc)):
        query = ("update dragen_sample_metadata set is_merged = 80140 where experiment_id = {}".format(sample['experiment_id']))
        run_query(query,database)
        raise ValueError("legacy globbin procedure is unable to locate distinct flowcells")
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
        # print("before={}".format(tracked_files))
        self.metadata['fastq_loc'] = get_fastq_loc(database,self.metadata,tracked_files)
        from pprint import pprint as pp
        pp(self.metadata['fastq_loc'])
        # print("after={}".format(tracked_files))
        if len(self.metadata['fastq_loc'])==0 and len(tracked_files)==0:
            raise ValueError("this is very, very, very wrong")
        elif len(tracked_files)>0:
            print("> WE HAVE {} TRACKED FILES AND {} LOCATIONS".format(len(tracked_files),len(self.metadata['fastq_loc'])))
            if len(self.metadata['fastq_loc'])>0:
                raise ValueError("these are mutually exclusive")
            self.metadata['bams'] = tracked_files
        print("> WE HAVE {} TRACKED FILES AND {} LOCATIONS".format(len(tracked_files),len(self.metadata['fastq_loc'])))
        # from pprint import pprint as pp
        # pp(vars(self))
        # exit(1)
        self.metadata['lane'] = get_lanes(database,self.metadata)
        from pprint import pprint as pp
        pp(self.metadata['lane'])
        print("> handing back silly object thing")
        # exit(1)

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
