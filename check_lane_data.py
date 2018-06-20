#!/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python

from utils import *
import os
import warnings
import json
import datetime
import pwd
import grp
import stat

def main():
    do_md5 = 0
    if len(sys.argv) == 2 and sys.argv[1] == 'md5check':
        do_md5 = 1
    setup_logging('/nfs/seqscratch10/mml2204/check_lane/','check_lane_table')
    logger = logging.getLogger(__name__)
    conn_sdb = get_connection_mml('sequenceDB','sequenceDB')
    cmd="select distinct prepid,data from Lane where data is not NULL and rg_status='fastq_archived'"
    cursor = conn_sdb.cursor()
    cursor.execute(cmd)
    row_info=cursor.fetchone()
    if row_info is None:
        raise Exception("No data entries in Lane?")
    f = open('/nfs/seqscratch10/mml2204/check_lane/error_file.txt','w')
    fout = open('/nfs/seqscratch10/mml2204/check_lane/output_file.txt','w')
    #femail = open('/nfs/seqscratch10/mml2204/check_lane/email_file.txt','w')
    fout.write("Filename\tExists\tSz_archive\tSz_json\tSz_mtime\tSz_mtime_json\tPerm\tOwner\tGroup\n")
    email_list = []
    while row_info is not None:
        logger.info("prepid {} is being checked".format(row_info['prepid']))
        #print(json.loads(row_info['data']))
        email_list = check_prepid(json.loads(row_info['data']),conn_sdb,f,fout,email_list)
        row_info=cursor.fetchone()

    f.close()
    fout.close()
    if len(email_list) == 0:
        send_email("No problems in names and permissions with archived read groups")
    else:
        send_email("MISMATCH between database and archived info in filenames/permissions",''.join(list(set(email_list))))

def check_prepid(row_info,conn_sdb,f,fout,email_list):
    archive_dir = row_info['fastq']['path']['archive']
    if not os.path.exists(archive_dir):
        f.write("{} not found\n".format(archive_dir))
        fout.write("{}\tNo\n".format(archive_dir))
        email_list.append("{} not found\n".format(archive_dir))
        return email_list
    
    fout.write("{}\tYes".format(archive_dir))
    file_stats = os.stat(archive_dir)
    perm =  stat.filemode(file_stats.st_mode)
    if perm != 'dr-xr-s---':
            f.write("Directory {} Problem with permissions: Expected:'dr-xr-s---', actual:'{}'\n".format(archive_dir,perm))
            email_list.append("Directory {} Problem with permissions: Expected:'dr-xr-s---', actual:'{}'\n".format(archive_dir,perm))
    owner = pwd.getpwuid(file_stats.st_uid)[0]
    if owner != 'dragen':
            f.write("Directory {} Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,owner))
            email_list.append("Directory {} Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,owner))
    group = grp.getgrgid(file_stats.st_gid)[0]
    fout.write("\t{}\t\t{}\t\t{}\t{}\t{}\n".format(file_stats.st_size,file_stats.st_mtime,perm,owner,group)) 
    
    num_fastqs = 0
    for entry in row_info['fastq']:
        if entry in ['r1','r2']: #'basename' in row_info['fastq'][entry]
            #check file exists
            if not os.path.exists('{}/{}'.format(archive_dir,row_info['fastq'][entry]['basename'])):
                 f.write('File {}/{} does not exist\n'.format(archive_dir,row_info['fastq'][entry]['basename']))
                 fout.write("{}/{}\tNo\n".format(archive_dir,row_info['fastq'][entry]['basename']))
                 email_list.append('File {}/{} does not exist\n'.format(archive_dir,row_info['fastq'][entry]['basename']))
                 continue
            else:
                fout.write("{}/{}\tYes".format(archive_dir,row_info['fastq'][entry]['basename']))
            #check sz, modification,permissions
            file_stats = os.stat('{}/{}'.format(archive_dir,row_info['fastq'][entry]['basename']))
            if file_stats.st_size != row_info['fastq'][entry]['size']:
                f.write("File: {}/{} Size mismatch: actual-{}, in_json-{}\n".format(archive_dir,row_info['fastq'][entry]['basename'],file_stats.st_size,row_info['fastq'][entry]['size']))
                email_list.append('File {}/{} does not exist\n'.format(archive_dir,row_info['fastq'][entry]['basename']))

            if file_stats.st_mtime != row_info['fastq'][entry]['modification']:
              f.write("File: {}/{} modification time mismatch: actual-{}, in_json-{}\n".format(archive_dir,row_info['fastq'][entry]['basename'],file_stats.st_mtime,row_info['fastq'][entry]['modification']))
              email_list.append('File: {}/{} modification time mismatch: actual-{}, in_json-{}\n'.format(archive_dir,row_info['fastq'][entry]['basename'],file_stats.st_mtime,row_info['fastq'][entry]['modification']))
            perm =  stat.filemode(file_stats.st_mode)
            if perm != '-r--r-----':
                f.write("File {}/{} Problem with permissions: expected:'-r--r-----', actual:'{}'\n".format(archive_dir,row_info['fastq'][entry]['basename'],perm))
                email_list.append("File {}/{} Problem with permissions: expected:'-r--r-----', actual:'{}'\n".format(archive_dir,row_info['fastq'][entry]['basename'],perm))

            owner = pwd.getpwuid(file_stats.st_uid)[0]
            if owner != 'dragen':
                f.write("File {}/{} Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,row_info['fastq'][entry]['basename'],owner))
                email_list.append("File {}/{} Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,row_info['fastq'][entry]['basename'],owner))
            
            group = grp.getgrgid(file_stats.st_gid)[0]
            fout.write("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file_stats.st_size,row_info['fastq'][entry]['size'],file_stats.st_mtime,row_info['fastq'][entry]['modification'],perm,owner,group))
            
            #md5sum file
            if not os.path.exists('{}/{}.md5sum'.format(archive_dir,row_info['fastq'][entry]['basename'])):
                f.write('File {}/{}.md5sum does not exist\n'.format(archive_dir,row_info['fastq'][entry]['basename']))
                fout.write("{}/{}.md5sum\tNo\n".format(archive_dir,row_info['fastq'][entry]['basename']))
                email_list.append("{}/{}.md5sum\tDoesnt exist\n".format(archive_dir,row_info['fastq'][entry]['basename']))
            else:
                fout.write("{}/{}.md5sum\tYes".format(archive_dir,row_info['fastq'][entry]['basename']))
                file_stats = os.stat('{}/{}.md5sum'.format(archive_dir,row_info['fastq'][entry]['basename']))
                perm =  stat.filemode(file_stats.st_mode)
                if perm != '-r--r-----':
                    f.write("File {}/{}.md5sum Problem with permissions: expected:'-r--r-----', actual:'{}'\n".format(archive_dir,row_info['fastq'][entry]['basename'],perm))
                    email_list.append("File {}/{}.md5sum Problem with permissions: expected:'-r--r-----', actual:'{}'\n".format(archive_dir,row_info['fastq'][entry]['basename'],perm))
                owner = pwd.getpwuid(file_stats.st_uid)[0]
                if owner != 'dragen':
                    f.write("File {}/{}.md5sum Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,row_info['fastq'][entry]['basename'],owner))
                    email_list.append("File {}/{}.md5sum Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,row_info['fastq'][entry]['basename'],owner))
                group = grp.getgrgid(file_stats.st_gid)[0]
                fout.write("\t{}\t\t{}\t\t{}\t{}\t{}\n".format(file_stats.st_size,file_stats.st_mtime,perm,owner,group))
                    
    #check laneBarcode and md5s.
    if not os.path.exists('{}/laneBarcode.html'.format(archive_dir)):
        f.write("{}/laneBarcode.html not found\n".format(archive_dir))
        fout.write("{}/laneBarcode.html\tNo\n".format(archive_dir))
        email_list.append("{}/laneBarcode.html not found\n".format(archive_dir))
    else:
        fout.write("{}/laneBarcode.html\tYes".format(archive_dir))
        file_stats = os.stat("{}/laneBarcode.html".format(archive_dir))
        perm =  stat.filemode(file_stats.st_mode)
        if perm != '-r--r-----':
            f.write("File {}/laneBarcode.html Problem with permissions: expected:'-r--r-----', actual:'{}'\n".format(archive_dir,perm))
            email_list.append("File {}/laneBarcode.html Problem with permissions: expected:'-r--r-----', actual:'{}'\n".format(archive_dir,perm))
        owner = pwd.getpwuid(file_stats.st_uid)[0]
        if owner != 'dragen':
            f.write("File {}/laneBarcode.html Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,owner))
            email_list.append("File {}/laneBarcode.html Incorrect owner: Reqd-dragen, actual-{}\n".format(archive_dir,owner))
        group = grp.getgrgid(file_stats.st_gid)[0]
        fout.write("\t{}\t\t{}\t\t{}\t{}\t{}\n".format(file_stats.st_size,file_stats.st_mtime,perm,owner,group))
    
    return email_list

main()

