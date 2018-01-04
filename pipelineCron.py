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
            JOIN SampleT s ON s.dbid=l.dbid
            WHERE (PIPELINECOMPLETE = 0 AND COMPLETE = 1 AND FAIL != 1)
            GROUP BY l.FCID
            ORDER BY MIN(PRIORITY) ASC , FROM_UNIXTIME(SEQEND) ASC
            """
    #print query
    complete_flowcells = run_query(query,database)
    for fcillumid in complete_flowcells:
        fcillumid = fcillumid['FCILLUMID']
        check_and_run_flowcell(fcillumid,config,database)
    if complete_flowcells == ():
        pass
        #print('No flowcells were found')
def check_and_run_flowcell(fcillumid,config,database):
    run_info_dict = parse_run_parameters_xml(fcillumid)
    bcl_dir = '{}/'.format(config.get('locs','bcl_dir'))
    check_flowcell_complete(fcillumid,bcl_dir,run_info_dict['runFolder'],
                            run_info_dict['type'],database)
    out_dir = ('{}/{}_{}_{}_Unaligned'
              ).format(config.get('locs','bcl2fastq_scratch_dir'),
                                  run_info_dict['runDate'],
                                  run_info_dict['machine'],
                                  run_info_dict['FCIllumID'])
    if glob(out_dir) == []:
        cmd = ("{} {}/bcl2fastq_pipeline_wrapper.py --fcillumid {} -r"
              ).format(config.get('programs','python36_program'),
                       config.get('locs','scr_dir'),
                       fcillumid)
        print(cmd)
        print()
        print('BCL2FASTQ STARTED')
        print()
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
