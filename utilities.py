import configparser
import logging
import os
import pymysql
import sys
import traceback
from db_statements import *
from glob import glob
from xml.etree import ElementTree

def parse_run_parameters_xml(fcillumid,database):
    config = get_config()
    run_folder = '{}/*{}*'.format(config.get('locs','bcl_dir'),fcillumid)
    xml_type = 'HiSeq'
    run_parameters_xml = glob('{}/runParameters.xml'.format(run_folder))
    if run_parameters_xml == []:
        xml_type = 'NovaSeq'
        run_parameters_xml_loc = '{}/RunParameters.xml'.format(run_folder)
        run_parameters_xml = glob(run_parameters_xml_loc)
    if run_parameters_xml == []:
        update_pipelinecomplete(fcillumid,'-404',database)
        raise Exception(("Could not find run parameters xml files: {}!"
                        ).format(run_parameters_xml_loc))
    if len(run_parameters_xml) != 1:
        update_pipelinecomplete(fcillumid,'-6',database)
        raise Exception("Too many run folders found!")
    else:
        run_parameters_xml = run_parameters_xml[0]

    with open(run_parameters_xml,'r') as xml:
        run_info_dict = {}
        run_info_dict['type'] = xml_type
        tree = ElementTree.parse(xml)
        if xml_type == 'NovaSeq':
            run_info_dict['RtaVersion'] = tree.findall('RtaVersion')[0].text
            run_info_dict['Side'] = tree.findall('Side')[0].text
            if run_info_dict['Side'] == 'None':
                run_info_dict['Side'] = tree.findall('PreRunFolder')[0].text[-1]
            run_info_dict['experimentName'] = tree.findall('ExperimentName')[0].text
            run_info_dict['runFolder'] = tree.findall('RunId')[0].text
            run_info_dict['runDate'] = run_info_dict['runFolder'].split('_')[0]
            run_info_dict['machineName'] = tree.findall('InstrumentName')[0].text
            run_info_dict['ControlSoftwareVer'] = tree.findall('ApplicationVersion')[0].text
            for FCInfo in tree.findall('RfidsInfo'):
                run_info_dict['FCIllumID'] = FCInfo.find('FlowCellSerialBarcode').text
            run_info_dict['machine'] = run_info_dict['machineName'] + run_info_dict['Side'][0]
        elif xml_type == 'HiSeq':
            run_info_dict['RtaVersion'] = tree.findall('Setup')[0].find('RTAVersion').text
            run_info_dict['experimentName'] = tree.findall('Setup')[0].find('ExperimentName').text
            run_info_dict['runFolder'] = tree.findall('Setup')[0].find('RunID').text
            run_info_dict['runDate'] = tree.findall('Setup')[0].find('RunStartDate').text
            ComputerName = tree.findall('Setup')[0].find('ComputerName').text
            if '-' in ComputerName:
                run_info_dict['machineName'] = ComputerName.split('-')[-1]
            else:
                run_info_dict['machineName'] = 'External'
            run_info_dict['ControlSoftwareVer'] = tree.findall('Setup')[0].find('ApplicationVersion').text
            run_info_dict['FCIllumID'] = tree.findall('Setup')[0].find('Barcode').text
            run_info_dict['Side'] = tree.findall('Setup')[0].find('FCPosition').text
            run_info_dict['machine'] = run_info_dict['machineName'] + run_info_dict['Side'][0]
        else:
            update_pipelinecomplete(fcillumid,'-100',database)
            raise Exception("Undefined xml_type:{}!".format(xml_type))
        run_info_dict['out_dir'] = '{}/{}_{}_{}_Unaligned'.format(config.get('locs','bcl2fastq_scratch_dir'),
        run_info_dict['runDate'],
        run_info_dict['machine'],
        run_info_dict['FCIllumID'])

    return run_info_dict


def get_connection(database):

    try:
        reader = configparser.RawConfigParser()
        my_cnf='{}/bcl_user.cnf'.format(os.path.dirname(os.path.realpath(__file__)))
        if not os.path.exists(my_cnf):
            raise Exception("\n\nSeriously?!? where's your '{}' for '{}' access".format(my_cnf,database))
        # else:
        #    print("using db config file '{}'".format(my_cnf))
        # hilarious!?!
        # my_cnf = glob(os.path.expanduser(conf))
        # reader.read(my_cnf[0])
        reader.read(my_cnf)
        db = 'client' + database; db_host = reader.get(db, 'host'); 
        db_user = reader.get(db, 'user'); db_pass = reader.get(db,'password')

        connection = pymysql.connect( host=db_host, user=db_user, password=db_pass,
          db='sequenceDB', cursorclass=pymysql.cursors.DictCursor )

        return connection

    except pymysql.err.OperationalError:
        traceback.print_exc()
        sys.exit("Wrong username/database or password found, please try again")

def is_external_or_legacy_sample(prepid,database):
    """In cases of legecy samples (prepid < 22000) or external samples. There are no cases of legecy samples with
    multiple preps"""
    # if prepid < 22000:
        # print("Sample has a prepID < 22000")
    if prepid < 50000:
        print("Sample has a prepID < 50000")
        return True

    else:
        IS_SAMPLE_EXTERNAL_FROM_PREPID = """
            SELECT CASE
                WHEN EXTERNALDATA = 1
                    THEN 1
                    ElSE 0
                END AS IS_EXTERNAL
            FROM prepT
            WHERE PREPID = {prepid}
            """
        query = IS_SAMPLE_EXTERNAL_FROM_PREPID.format(prepid=prepid)
        is_external_sample = int(run_query(query,database)[0]['IS_EXTERNAL'])
        if is_external_sample == 1:
            #print("Sample is an external sample!")
            return True
        elif is_external_sample == 0:
            #print("Sample is not an external sample")
            return False
        else: #Returns None
            raise Exception("No value found.  Does prepID exist?")


def check_number_query_results(results,expected):
    if len(results) == 0:
        raise Exception("No results found!".format(results,expected))
    elif len(results) > 1 and expected == 'many':
        pass
    elif len(results) <= 1 and expected == 'many':
        raise Exception("Number of results: {} < expected: {}".format(results,expected))
    elif len(results) == expected:
        pass
    elif len(results) < expected:
        raise Exception("Number of results: {} < expected: {}".format(results,expected))
    elif len(results) > expected:
        raise Exception("Number of results: {} > expected: {}".format(results,expected))
    else:
        raise Exception("Unhandled {} condition. Results: {}  Expected: {}").format(__name__,results,expected)

def check_machine(machine,fcillumid,database):
    """Confirms the machine in the database matches machine value in the
       runParameters.xml. If they differ we update the Flowcell.machine
       entry.  Differeneces can exist in cases of re-hybridization of a
       flowcell on a different machine.

       Note that NS1 is the HTS name for the novaseq and the database
       but on the NovaSeq itself its called A00116 from the runParameters.xml
    """
    machine_from_db = run_query(
        GET_MACHINE_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    check_number_query_results(machine_from_db,1)
    machine_from_db = machine_from_db[0]['MACHINE']

    if machine_from_db == 'NS1':
        machine_from_db = 'A00116'
    elif machine_from_db == 'NS2':
        machine_from_db = 'A00123'

    if machine_from_db != machine:
        update_machine_query = """UPDATE Flowcell
            SET Machine='{0}'
            WHERE FCIllumID = '{1}'
            """.format(machine,fcillumid)
        run_query(update_machine_query,database)

    machine_complete = run_query(
        GET_MACHINE_COMPLETE_STATUS_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    check_number_query_results(machine_complete,1)
    if machine_complete[0]['COMPLETE'] != 1:
        update_pipelinecomplete(fcillumid,'-3',database)
        raise Exception("Flowcell {} has not yet been completed".format(fcillumid))

    machine_failed = run_query(
        GET_MACHINE_FAILED_STATUS_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    if machine_failed[0]['FAIL'] == 1:
        update_pipelinecomplete(fcillumid,'-99',database)
        raise Exception("Flowcell {} has been failed".format(fcillumid))

def getSSSLaneFraction(prepid,fcillumid,chem_version,lane_num,config,database):
    #get seqtype
    query = GET_SEQTYPE_FROM_PREPID.format(prepid=prepid)
    seqtype = run_query(query,database)[0]['SAMPLE_TYPE'].lower()
    #Values currently hard coded at Colin's request.  See Redmine#2320
    if seqtype == 'exome':
        return 1.0/int(config.get(chem_version,'exome_per_lane'))
    elif seqtype == 'genome':
        return 1.0/int(config.get(chem_version,'genome_per_lane'))
    elif seqtype == 'custom capture': #note not custom_capture
        return 1.0/int(config.get(chem_version,'custom_per_lane'))
    elif seqtype == 'rnaseq':
        return 1.0/int(config.get(chem_version,'rna_per_lane'))
    else:
        raise Exception("Unhandled seqtype: {}!".format(seqtype))

def update_pipelinecomplete(fcillumid,code,database):
    query = UPDATE_PIPELINE_COMPLETE.format(fcillumid=fcillumid,code=code)
    return run_query(query,database)

def create_sss_from_database(fcillumid,machine,run_info_dict,config,database):
    """Newer version of bcl2fastq v2 require a new sample sheet format
       usually generated by the Illumina Experiment Manager"""
    sss_loc = '{}/{}/{}_{}_{}.csv'.format(config.get('locs','bcl_dir'),
                                          run_info_dict['runFolder'],
                                          machine,
                                          run_info_dict['runDate'],
                                          fcillumid
                                         )
    print("Creating SSS: {}".format(sss_loc))
    recipe = run_query(GET_FLOWCELL_RECIPE.format(fcillumid=fcillumid),database)
    chem_version = run_query(GET_FLOWCELL_CHEMVER.format(fcillumid=fcillumid),database)
    all_projects = run_query(GET_FLOWCELL_PROJECTS.format(fcillumid=fcillumid),database)
    all_projects = ','.join(list(map((lambda x:x['GAFBIN']),all_projects)))
    total_num_lanes = run_query(GET_TOTAL_NUM_LANES_FROM_FLOWCELL.format(fcillumid=fcillumid),database)
    flowcell_creator = run_query(GET_FLOWCELL_CREATOR.format(fcillumid=fcillumid),database)
    outfile=open(sss_loc,'w')
    #Sample Sheet Header
    outfile.write("[Header]\n")
    outfile.write("IEMFileVersion,4\n")
    outfile.write("Investigator Name,{}\n".format(flowcell_creator[0]['FULLNAME']))
    outfile.write("Experiment Name,{}\n".format(fcillumid))
    outfile.write("Date,{}\n".format(run_info_dict['runDate']))
    outfile.write("Workflow,GenerateFASTQ\n")
    outfile.write("Application,HiSeq FASTQ Only\n")
    outfile.write("Assay,TruSeq HT\n")
    outfile.write("Description,{}\n".format(all_projects))
    outfile.write("Chemistry,Default\n")
    outfile.write("\n")
    outfile.write("[Reads]\n")
    outfile.write("{}\n".format(recipe[0]['LenR1']))
    outfile.write("{}\n".format(recipe[0]['LenR2']))
    outfile.write("\n")
    outfile.write("[Data]\n")
    outfile.write("\n")
    outfile.write("Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n")
    for lane_num in range(1,total_num_lanes[0]['TOTAL_NUM_LANES']+1):
        prepid_on_lane = run_query(GET_PREPID_ON_LANE.format(
            fcillumid=fcillumid,
            lane_num=lane_num),database)
        check_number_query_results(prepid_on_lane,'many')
        for prepid in prepid_on_lane:
            prepid = prepid['PREPID']
            laneFraction = getSSSLaneFraction(prepid,fcillumid,chem_version[0]['CHEMVER'],lane_num,config,database)

            # empty strings are for 'Sample_Plate,Sample_Well,I7_Index_ID'
            # Description contains expected lane fraction, picomolar amount, flowcell user
            sql = "SELECT l.LaneNum Lane,\
                pt.CHGVID SampleID,\
                pt.CHGVID SampleID,\
                '',\
                '',\
                '',\
                (case \
                    WHEN f.recipe=2 THEN '' \
                    ELSE pt.adapterLet \
                END) 'Index',\
                replace(s.GAFbin,' ','') Project, \
                CONCAT(round('{0}',4),'_',round(s2r.picomoles,1),'pM') Description \
                FROM Lane l \
                    JOIN Flowcell f ON f.FCID=l.FCID \
                    JOIN prepT pt ON l.prepID=pt.prepID \
                    JOIN samplesTOrun s2r ON s2r.seqID=l.seqID \
                    JOIN SampleT s ON s.sample_id=pt.sample_id \
                WHERE \
                    l.prepid='{1}' AND \
                    f.FCillumID='{2}' AND \
                    LaneNum='{3}'".format(laneFraction,prepid,fcillumid,lane_num)
            #print sql
            sss_line = run_query(sql,database)
            check_number_query_results(sss_line,1)
            str_desc = sss_line[0]['Description'] #why did I use decode?
            #str_desc = sss_line[0]['Description']
            outfile.write("{Lane},{SampleID},{SampleID},,,,{Index},{Project},{str_desc}\n".format(str_desc=str_desc,**sss_line[0]))
    outfile.close()
    #copies sequencing sample sheet to genotyping location
    print("Finish creating SSS: {}".format(sss_loc))
    cmd="chmod -R 770 {0}".format(sss_loc)
    if os.system(cmd) != 0:
        raise Exception("error occured in {}".format(cmd))
    return sss_loc

def setup_logging(machine,fcillumid,logs_dir):
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        filename=('{}/{}_{}.log'
                 ).format(logs_dir,machine,
                          fcillumid))
    logger = logging.getLogger(__name__)

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

def get_config():
    config = configparser.ConfigParser()
    config.read('{}/config.ini'.format(os.path.dirname(os.path.realpath(__file__))))
    return config

def check_exist_bcl_dir(fcillumid,BCLDrive,database):
    #### is this a joke? : os.path.isdir(path) : https://docs.python.org/2/library/os.path.html
    dir_path = glob(BCLDrive)

    if dir_path != []:
        update_pipelinecomplete(fcillumid,'-7',database)
        raise Exception('\n\nBCL directory already exists: {}. So, me thinks me should check to see if \
there is a checkpoint file bcl_complete and if so re-run next step or complain and ask if \
i should delete and re-run the whole thing. if there is not, then perhaps i should check.\
that funky magic job name and see if it exists and if not we could safely assume it failed...\n\n\
However, i really should not be using a very busy cluster when i have fancy dragen machines \
who are very bored heating the marvelous state of NJ. Come to think of it, if we then mapped the \
fastq RG immediately we could even only ever store bam and make the release step simply trigger \
merging and not mapping cos we would already have that lovely, juicy bam!'.format(BCLDrive))

def check_bcl_complete(bcl2fastq_dir):
    bcl_complete_flag_loc = glob(bcl2fastq_dir + '/bcl_complete')
    if bcl_complete_flag_loc == []:
        print(bcl2fastq_dir + '/bcl_complete')
        raise Exception('bcl_complete flag file not found! Check if BCL2Fastq completed successfully')

def check_flowcell_complete(fcillumid,bcl_dir,run_folder_path,machine_type,database):
    rta_complete_loc = bcl_dir + run_folder_path + '/RTAComplete.txt'
    if os.path.isfile(rta_complete_loc) == False:
        print(rta_complete_loc)
        update_pipelinecomplete(fcillumid,'-4',database)
        raise Exception("RTA has not completed!")
    else:
        print("RTAComplete.txt check: OK!")
    if machine_type == 'NovaSeq':
        storage_complete_loc = bcl_dir + run_folder_path + '/CopyComplete.txt'
        if os.path.isfile(storage_complete_loc) == False:
            update_pipelinecomplete(fcillumid,'-5',database)
            raise Exception("Copy has not completed!")
        else:
            print("CopyComplete.txt check: OK!")
        #run_complete_loc = bcl_dir + run_folder_path + '/RunComplete.txt'
        #if os.path.isfile(run_complete_loc) == False:
        #    update_pipelinecomplete(fcillumid,'-8',database)
        #    raise Exception("Run has not completed!")
        #else:
        #    print("RunComplete.txt check: OK!")

def get_user():
    p = os.popen('echo $USER')
    userName = p.readline().strip()
    p.close()
    return userName

def get_user_id(database):
    uni = get_user()
    get_userID_query = GET_USERID_FROM_UNI.format(uni=uni)
    userID = run_query(get_userID_query,database)
    check_number_query_results(userID,1)
    userID = userID[0]['USERID']
    return userID
