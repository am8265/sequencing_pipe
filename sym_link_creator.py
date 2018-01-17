""" EX: [CHGVID]/raw/cg1096416_FASTQ13.fastq.gz --> [CHGVID]/[fcillumid]/[IGMID]_AAAAAA_L001_R1_001.fastq.gz"""
import argparse
import gzip
import io
import os
import re
import subprocess
import sys
from glob import glob


def main():
    args = parse_arguments()
    sample_path = args.fastq_folder
    verbose = args.verbose

    while  sample_path[-1] == '/':
        sample_path = sample_path[0:-1]
    if os.path.isabs(sample_path) == False:
        raise Exception("Location is not absolute path!")
    sample = sample_path.split('/')[-1]
    if os.path.isdir('{}/raw'.format(sample_path))== False:
        raise Exception("Raw fastq folder not found!")
    check_for_fastqs(sample_path)
    check_for_symlinks(sample_path,verbose)

    fastqs = glob('{}/raw/*fastq.gz'.format(sample_path))
    if fastqs == []:
        fastqs = glob('{}/raw/*/*fastq.gz'.format(sample_path))
        if fastqs == []:
            fastqs = glob('{}/raw/*txt.gz'.format(sample_path))
            if fastqs == []:
                raise Exception("No fastq.gz were found!")
    total_fastqs = len(fastqs)
    for fastq in fastqs:
        fcillumid,lane,read,adapter = get_rg_info(fastq,verbose)
        fastq_counter = fastq.split('.')[0].split('_')[-1]
        if len(fastq_counter) != 3 or fastq_counter.isdigit() == False:
            fastq_counter = '001'
        sym_link_fastq = ('{}/{}/{}_{}_L00{}_R{}_{}.fastq.gz'
                         ).format(sample_path,fcillumid,sample,adapter,
                                  lane,read,fastq_counter)
        fastq_dir = '{}/{}'.format(sample_path,fcillumid)
        make_fcillumid_dir(fastq_dir,verbose)

        ln_cmd = ['ln','-s',fastq,sym_link_fastq]

        if os.path.exists(sym_link_fastq):

            raise Exception('Symlink already exists!: {}'.format(ln_cmd))
        else:
            if verbose:
                print(' '.join(ln_cmd))
            subprocess.check_call(ln_cmd)
    symlinks = glob('{}/*/*fastq.gz'.format(sample_path))
    total_symlinks = []
    for symlink in symlinks:
        if 'raw' not in symlink:
            total_symlinks.append(symlink)
    total_symlinks_count = len(total_symlinks)
    if total_symlinks_count != total_fastqs:
        print(total_symlinks_count,total_fastqs)
        raise Exception('Number of raw fastqs and symlinks do not match!')
    else:
        print("Sample {}'s symlink creation done".format(sample))

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-f','--fastq-folder', required=True,
                        help="Specify scratch dir for bcl2fastq")
    args=parser.parse_args()
    return args

def check_for_fastqs(sample_path):
    fastqs = glob('{}/raw/*fastq')
    if fastqs != []:
        raise Exception("Fastqs not gzip!")

def check_for_symlinks(sample_path,verbose):
    possible_folders = glob('{}/*'.format(sample_path))
    for folder in possible_folders:
        if os.path.isdir(folder) == True and 'raw' not in folder:
            clean_up_symlink_folder(folder,verbose)

def clean_up_symlink_folder(folder,verbose):
    files = glob('{}/*'.format(folder))
    for file in files:
        if os.path.islink(file) == True:
            if verbose:
                print('removing symlink: {}'.format(file))
            os.remove(file)
    files = glob('{}/*'.format(folder))
    if files != []:
        raise Exception("Files still exist in folder")
    else:
        if verbose:
            print('removing folder: {}'.format(folder))
        os.rmdir(folder)


def make_fcillumid_dir(fastq_dir,verbose):
    if not os.path.exists(fastq_dir):
        try:
            if verbose:
                print('mkdir {}'.format(fastq_dir))
            os.makedirs(fastq_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                    raise

def check_fcillumid(rg_info):
    tmp_fcillumid = rg_info[0]
    if len(rg_info) == 10:
        fcillumid = rg_info[2]
    else:
        search_obj = re.search('[CHD].*[XY]',tmp_fcillumid)
        fcillumid = search_obj.group()
    if len(fcillumid) != 9 and fcillumid[-1] in 'xyXY':
        raise Exception('Check fcillumid: {}'.format(fcillumid))
    else:
        return fcillumid

def get_rg_info(fastq,verbose):
    gz = gzip.open(fastq,'rt')
    fastq_read_header = gz.readline()
    rg_info = fastq_read_header.strip().split(':')
    fcillumid = check_fcillumid(rg_info)
    gz.close()
    if verbose:
        print(fastq_read_header.strip())
    if len(rg_info) == 10:
        """ IGM Illumina read head example
        @K00347:44:HL7K2BBXX:5:1101:13352:1384 1:N:0:NCTATCCT+NGGATAGG
        @Instrument:RunID:FlowCellID:Lane:Tile:X:Y:UMI ReadNum:FilterFlag:0:IndexSequence or SampleNumber"""
        machine,run_number,fcillumid2,lane,tile,X,read,control,something,adapter = rg_info

        read = read.split(' ')[1]
    elif len(fastq_read_header.split(':')) == 7:
        """@HISEQ:549:C6PE0ANXX:2:2305:17233:17109/1
        @Instrument:RunID:FlowCellID:Lane:Tile:X:Y IndexSequence"""
        machine,run_number,fcillumid2,lane,tile,X,Y_and_read = rg_info
        read = Y_and_read.split('/')[1]
        adapter = 'AAAAAA'
    elif len(fastq_read_header.split(':')) == 5:
        if re.search('[ACGNT]',rg_info[4]):
            """@FCHNYHLBCXX:1:1207:6982:25608#TGACCAAN/1"""
            fcillumid2,lane,tile,X,Y_and_read = rg_info
            fcillumid2 = fcillumid2[3:]
            tmp = Y_and_read.split('#')[1]
            adapter,read = tmp.split('/')
        else:
            fcillumid2,lane,tile,X,Y_and_read = rg_info
            fcillumid2 = fcillumid2[3:]
            read = Y_and_read.split('/')[1]
            adapter = 'AAAAAA'
    else:

        print(fastq_read_header)
        raise Exception('Incorrect read header format!')

    if read not in '123':
        raise Exception('read is not formatted correctly: {}!'.format(read))

    adapter = adapter.replace('+','-')
    if verbose:
        print(fcillumid,lane,read,adapter)
    #rg_info_sanity_check(fcillumid,lane,read,adapter)
    return fcillumid,lane,read,adapter

if __name__ == '__main__':
    main()
