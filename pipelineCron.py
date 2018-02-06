#!/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python
import configparser
import os
from glob import glob
from subprocess import Popen, PIPE
from utilities import *

def main(config):
    database='sequenceDB'

    """grab unfailed flowcells that haven't completed sequenceDB
       order by priority of samples on flowcell"""

    query = """SELECT FCILLUMID
            FROM Flowcell f
            JOIN Lane l ON f.fcid=l.fcid
            JOIN prepT p ON p.dbid=l.dbid
            WHERE (PIPELINECOMPLETE = 0 AND COMPLETE = 1 AND FAIL != 1)
            GROUP BY l.FCID
            ORDER BY MIN(PRIORITY) ASC , FROM_UNIXTIME(SEQEND) ASC
            """
    #print(query)
    # print("checking for flowcells")
    complete_flowcells = run_query(query,database)
    for fcillumid in complete_flowcells:
        fcillumid = fcillumid['FCILLUMID']
        check_and_run_flowcell(fcillumid,config,database)
    if complete_flowcells == ():
        #print('No flowcells were found')
        pass

def check_and_run_flowcell(fcillumid,config,database):
    print("running {}".format(fcillumid))
    run_info_dict = parse_run_parameters_xml(fcillumid,database)
    bcl_dir = '{}/'.format(config.get('locs','bcl_dir'))
    check_flowcell_complete(fcillumid,bcl_dir,run_info_dict['runFolder'],
                            run_info_dict['type'],database)
    out_dir = run_info_dict['out_dir']

    if glob(out_dir) == []:
        src_dir=os.path.dirname(os.path.realpath(__file__))
        print("we'll use '{}'".format(src_dir))
        cmd = ("{}/bcl2fastq_pipeline_wrapper.py --fcillumid {} -r").format( src_dir, fcillumid)
        print(cmd)
        print()
        print('BCL2FASTQ STARTING')
        print()
        if os.system(cmd) != 0:
            raise Exception("\n\ncouldn't run '{}'!".format(cmd))
    else:
        print('BCL DIR EXISTS!')
        print("***************************")
        print("* FLOWCELL NOT STARTED!!! *")
        print("***************************")

if __name__ == '__main__':
    config = get_config()
    try:
        main(config)
    except:
        import traceback
        print(traceback.format_exc())
