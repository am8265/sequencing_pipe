import configparser
import os
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
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
            WHERE (PIPELINECOMPLETE != 1 AND COMPLETE = 1 AND FAIL != 1)
            GROUP BY l.FCID
            ORDER BY MIN(PRIORITY) ASC , FROM_UNIXTIME(SEQEND) ASC
            """
    #print query
    complete_flowcells = run_query(query,database)
    for fcillumid in complete_flowcells:
        fcillumid = fcillumid['FCILLUMID']
        check_and_run_flowcell(fcillumid,config,database)
        input(fcillumid)

def check_and_run_flowcell(fcillumid,config,database):
    run_info_dict = parse_run_parameters_xml(fcillumid)
    try:
        bcl_dir = '{}/'.format(config.get('locs','bcl_dir'))
        check_flowcell_complete(bcl_dir,run_info_dict['runFolder'],
                                run_info_dict['type'])
        out_dir = ('{}/{}_{}_{}_Unaligned'
                  ).format(config.get('locs','bcl2fastq_scratch_dir'),
                                      run_info_dict['runDate'],
                                      run_info_dict['machine'],
                                      run_info_dict['FCIllumID'])
    except Exception:
        status = 'RTA or copy not complete'
        print(status)
        set_fc_status(status,fcillumid,database)
    if glob(out_dir):
        status = 'BCL dir exists'
        print(status)
        set_fc_status(status,fcillumid,database)

    else:
        print('ok')
        cmd = ("{} {}/bcl2fastq_pipeline_wrapper.py --fcillumid {} -r"
              ).format(config.get('programs','python36_program'),
                       config.get('locs','scr_dir'),
                       fcillumid)
        print(cmd)
        p = Popen(cmd, shell=True, stdout=PIPE)
        output,error = p.communicate()
        print(output.replace('\\n','\n'))
        if p.returncode != 0:
            print(error)
            failure_email(config,fcillumid,error)
        else:
            status = 'bcl2fastq started'
            print(status)
            set_fc_status(status,fcillumid,database)

def set_fc_status(status,fcillumid,database):
    if status == 'bcl2fastq started':
        code = 1
    elif status == 'BCL dir exists':
        code = -2
    elif status == 'RTA or copy not complete':
        code = -3
    else:
        raise Exception('Unhandled status: {}!'.format(status))

    query = """UPDATE Flowcell
               SET PIPELINECOMPLETE = {}
               WHERE FCILLUMID = '{}'""".format(code,fcillumid)
    run_query(query,database)
    if code != 1:
       failure_email(config,fcillumid,status)

def failure_email(config,fcillumid,error):
    uni = get_user()
    sender = '{}@cumc.columbia.edu'.format(uni)
    recipient = [config.get('emails','failures_email')]
    recipient = ['joshbridgers@gmail.com']
    msg = MIMEMultipart()
    msg['Subject'] = "bcl2fastq wrapper failure: {}".format(fcillumid)
    msg['From'] = sender
    msg['To'] = recipient[0]
    if error is not None:
        msg.attach(MIMEText(error, 'plain'))
    server = smtplib.SMTP('localhost')
    server.sendmail(sender, recipient, msg.as_string())
    server.quit()

if __name__ == '__main__':
    config = get_config()
    try:
        main(config)
    except:
        import traceback
        print(traceback.format_exc())
        #failure_email(config,'CRON_JOB',traceback.format_exc())
