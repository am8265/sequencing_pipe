import configparser
import logging
import os
import pymysql
import sys
from db_statements import *
from glob import glob
from xml.etree import ElementTree

def parse_run_parameters_xml(fcillumid):
    run_folder = '/nfs/seqscratch1/Runs/*{}*'.format(fcillumid)
    xml_type = 'HiSeq'
    run_parameters_xml = glob('{}/runParameters.xml'.format(run_folder))
    if run_parameters_xml == []:
        xml_type = 'NovaSeq'
        run_parameters_xml = glob('{}/RunParameters.xml'.format(run_folder))
    if run_parameters_xml == []:
        raise ValueError("Could not find run parameters xml files!")
    if len(run_parameters_xml) != 1:
        raise ValueError("Too many run folders found!")
    else:
        run_parameters_xml = run_parameters_xml[0]

    with open(run_parameters_xml,'r') as xml:
        run_info_dict = {}
        tree = ElementTree.parse(xml)
        if xml_type == 'NovaSeq':
            run_info_dict['RTAVersion'] = tree.findall('RTAVersion')[0].text
            run_info_dict['experimentName'] = tree.findall('ExperimentName')[0].text
            run_info_dict['runFolder'] = tree.findall('RunId')[0].text
            run_info_dict['runDate'] = runFolder.split('_')[0]
            run_info_dict['machineName'] = tree.findall('InstrumentName')[0].text
            run_info_dict['ControlSoftwareVer'] = tree.findall('ApplicationVersion')[0].text
            for FCInfo in tree.findall('RfidsInfo'):
                run_info_dict['FCIllumID'] = FCInfo.find('FlowCellSerialBarcode')
                run_info_dict['machineSide'] = FCInfo.find('FlowCellPartNumber')
        elif xml_type == 'HiSeq':
            run_info_dict['RTAVersion'] = tree.findall('Setup')[0].find('RTAVersion').text
            run_info_dict['experimentName'] = tree.findall('Setup')[0].find('ExperimentName').text
            run_info_dict['runFolder'] = tree.findall('Setup')[0].find('RunID').text
            run_info_dict['runDate'] = tree.findall('Setup')[0].find('RunStartDate').text
            run_info_dict['machineName'] = tree.findall('Setup')[0].find('ComputerName').text.split('-')[1]
            run_info_dict['ControlSoftwareVer'] = tree.findall('Setup')[0].find('ApplicationVersion').text
            run_info_dict['FCIllumID'] = tree.findall('Setup')[0].find('Barcode').text
            run_info_dict['machineSide'] = tree.findall('Setup')[0].find('FCPosition').text

        else:
            raise ValueError("Undefined xml_type:{}!".format(xml_type))
    return run_info_dict


def get_connection(database):
    try:
        reader = configparser.RawConfigParser()
        reader.read('/home/jb3816/.my.cnf')
        db = 'client' + database
        db_host = reader.get(db, 'host')
        db_user = reader.get(db, 'user')
        db_pass = reader.get(db,'password')
        connection = pymysql.connect(host=db_host,
                                     user=db_user,
                                     password=db_pass,
                                     db='sequenceDB',
                                     cursorclass=pymysql.cursors.DictCursor)
        return connection
    except pymysql.err.OperationalError:
            sys.exit("Wrong username/database or password found, please try again")

def check_number_query_results(results,expected):
    if len(results) == 0:
        raise ValueError("No results found!".format(results,expected))
    elif len(results) > 1 and expected == 'many':
        pass
    elif len(results) <= 1 and expected == 'many':
        raise ValueError("Number of results: {} < expected: {}").format(results,expected)
    elif len(results) == expected:
        pass
    elif len(results) < expected:
        raise ValueError("Number of results: {} < expected: {}").format(results,expected)
    elif len(results) > expected:
        raise ValueError("Number of results: {} > expected: {}").format(results,expected)
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
            """.format(Machine,fcillumid)
        run_query(update_machine_query,database)

    machine_complete = run_query(
        GET_MACHINE_COMPLETE_STATUS_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    check_number_query_results(machine_complete,1)
    if machine_complete[0]['COMPLETE'] != 1:
        raise ValueError("Flowcell {} has not yet been completed".format(fcillumid))

    machine_failed = run_query(
        GET_MACHINE_FAILED_STATUS_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    if machine_failed[0]['FAIL'] == 1:
        raise ValueError("Flowcell {} has been failed".format(fcillumid))

def getSSSLaneFraction(prepid,fcillumid,chem_version,lane_num,config,database):
    #get seqtype
    seqtype = run_query(GET_SEQTYPE_FROM_PREPID.format(prepid=prepid),database)[0]['SEQTYPE'].lower()
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
        raise ValueError("Unhandled seqtype: {}!".format(seqtype))

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
                    WHEN st.SeqType='Exome' THEN pm.adapterlet \
                    WHEN st.Seqtype='RNAseq' THEN s2r.adapterlet \
                    WHEN st.Seqtype='Genome' THEN s2r.adapterlet \
                    WHEN st.SeqType='Custom Capture' THEN pm.adapterlet \
                END) 'Index',\
                replace(s.GAFbin,' ','') Project, \
                CONCAT(round('{0}',4),'_',round(s2r.picomoles,1),'pM') Description \
                FROM Lane l \
                    JOIN Flowcell f ON f.FCID=l.FCID \
                    JOIN prepT pt ON l.prepID=pt.prepID \
                    JOIN samplesTOrun s2r ON s2r.seqID=l.seqID \
                    JOIN SeqType st ON l.prepID=st.prepID \
                    JOIN SampleT s ON s.DBID=pt.DBID \
                    LEFT JOIN poolMembers pm ON \
                        (pm.DBID=pt.DBID AND pm.poolID=l.poolID) \
                WHERE \
                    l.prepid='{1}' AND \
                    f.FCillumID='{2}' AND \
                    LaneNum='{3}'".format(laneFraction,prepid,fcillumid,lane_num)
            #print sql
            sss_line = run_query(sql,database)
            check_number_query_results(sss_line,1)
            outfile.write("{Lane},{SampleID},{SampleID},,,,{Index},{Project},{str_desc}\n".format(**sss_line[0],str_desc=sss_line[0]['Description'].decode('UTF-8')))
    outfile.close()
    #copies sequencing sample sheet to genotyping location
    cp_cmd= ('cp {} /nfs/igmdata01/Sequencing_SampleSheets/'
             ).format(sss_loc)
    os.system(cp_cmd)
    return sss_loc

def setup_logging(machine,fcillumid,fastq_archive_drive,bcl_dir):
    bcl_drive = bcl_dir.split('/')[2]
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        filename=('/nfs/{}/summary/GAF_PIPELINE_LOGS/{}_{}_{}.log'
                 ).format(fastq_archive_drive,machine,
                          fcillumid,bcl_drive))
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
    config.read('config.ini')
    return config

def check_exist_bcl_dir(BCLDrive):
    dir_path = glob(BCLDrive)
    if dir_path != []:
        raise Exception('BCL directory already exists: {}!'.format(BCLDrive))

def check_rta_complete(bcl_dir,run_folder_path):
    rta_complete_loc = bcl_dir + run_folder_path + '/RTAComplete.txt'
    if os.path.isfile(rta_complete_loc) == False:
        raise Exception("RTA has not completed!")
    else:
        print("RTA has already completed")

def get_user_id(database):
    p = os.popen('echo $USER')
    userName = p.readline().strip()
    p.close()
    get_userID_query = "SELECT userID FROM users WHERE netid='{}'".format(userName)
    userID = run_query(get_userID_query,database)

    if len(userID) == 1:
        return userID[0]['userID']
    elif len(userID) > 1:
        raise ValueError("Too many userIDs were found with netid {}!".format(userName))
    else: #len == 0
        raise ValueError("No userID was found with netid ()!".format(userName))
