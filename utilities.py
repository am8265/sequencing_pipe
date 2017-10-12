import configparser
import pymysql
import os
import sys
from db_statements import *

def get_connection(database):
    try:
        reader = configparser.RawConfigParser()
        reader.read('/home/jb3816/.my.cnf')
        db = 'client' + database
        db_host = reader.get(db, 'host')
        db_user = reader.get(db, 'user')
        db_pass = reader.get(db,'password')
        print(db_host,db_user,db_pass,database)
        connection = pymysql.connect(host=db_host,
                                     user=db_user,
                                     password=db_pass,
                                     db=database,
                                     cursorclass=pymysql.cursors.DictCursor)
        return connection
    except pymysql.err.OperationalError:
            sys.exit("Wrong username/database or password found, please try again")

def check_number_query_results(results,expected):
    if len(results) == 0:
        raise ValueError("No results found!").format(results,expected)
    elif len(results) > 1 and expected == 'many':
        pass
    elif len(results) =< 1 and expected == 'many':
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

    if machine_from_db != machine
        update_machine_query = """UPDATE Flowcell
            SET Machine='{0}'
            WHERE FCIllumID = '{1}'
            """.format(Machine,FCID)
        run_query(update_machine_query,database)

def check_machine_status(machine,fcillumid,database):
    machine_complete = run_query(
        GET_MACHINE_COMPLETE_STATUS_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    check_number_query_results(machine_complete,1)
    if machine_complete[0]['COMPLETE'] != '1':
        raise ValueError("Flowcell {} has not yet been completed".format(fcillumid))

    machine_failed = run_query(
        GET_MACHINE_FAILED_STATUS_FROM_FCILLUMID_QUERY.format(
            fcillumid=fcillumid),database)
    if machine_failed[0]['FAIL'] != '1':
        raise ValueError("Flowcell {} has been failed".format(fcillumid))

def setup_logging(machine,FCID,fastq_archive_drive,bcl2fastq_scratch_dir):
    bcl2fastq_drive = bcl2fastq_scratch_dir.split('/')[2]
	logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        filename=('/nfs/{}/summary/GAF_PIPELINE_LOGS/{}_{}_{}.log'
                                 ).format(fastq_archive_drive,machine,
                                          fcillumid,bcl2fastq_scratch_dir)
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
        logger.warn('BCL directory already exists!')
        raise Exception('BCL directory already exists!')

def check_rta_complete(run_folder_path):
    if os.path.isfile('{}/RTAComplete.txt'.format(run_folder_path)) == False:
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
