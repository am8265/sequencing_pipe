import configparser
import pymysql
import os
import sys

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

def check_machine(machine,fcillumid,database):
    get_machine_query = ("SELECT MACHINE "
                          "FROM Flowcell "
                          "WHERE FCIllumID = '{0}'"
                         ).format(fcillumid)
    machine_from_db = run_query(get_machine_query)[0]['machine']
    if len(machine_from_db) != 1:
        raise Exception("Incorrect number of FCIDs found!")
    else:
        if machine_from_db != Machine:

            #Updates database entry if incorrect.  Need to add logging
            sql = """UPDATE Flowcell
                    SET Machine='{0}' 
                    WHERE FCIllumID = '{1}'
                    """.format(Machine,FCID)

            sequenceDB.execute(sql)
            sequenceDB.execute('COMMIT;')

    return Machine


def setup_logging(machine,FCID,seqsata_drive):
	logging.basicConfig(level=logging.INFO,format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',filename='/nfs/%s/summary/GAF_PIPELINE_LOGS/%s_%s_%s.log' % (seqsata_drive,machine,FCID,seqsata_drive))
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
