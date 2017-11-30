import argparse
import configparser
import logging
from dragen_sample import dragen_sample
from utilities import *


def main():
    config = get_config()
    args = parse_arguments()

    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'

    if args.fcillumid:
        auto_release_flag = True
        get_sample_info_from_flowcell = ("SELECT DISTINCT CHVID,p.PREPID,SEQTYPE,p.EXOMEKIT"
                                         "FROM Lane l "
                                         "JOIN Flowcell f ON f.fcid=l.fcid "
                                         "JOIN SeqType st on st.prepid=l.prepid "
                                         "JOIN prepT p ON p.prepid=l.prepid "
                                         "WHERE FCILLUMID='{}'"
                                        ).format(args.fcillumid)

        for sample_info in run_query(get_sample_info_from_flowcell)
            ppid = check_ppid(sample_info)
            sample = dragen_sample(sample_info['CHVID'],sample_info['SeqType'],
                                   ppid,sample_info['EXOMEKIT'],get_curs())

    print("Starting sample release for: {}".format(args.fcillumid))



def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fcillumid", dest='fcillumid',
                        help="Specify Illumina's Flowcell ID")
    group.add_argument('-s', '--sample_name',dest='sample_name',
                        help="Specify sample name")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-t','--sample_type', dest='sample_type',
                        help="Specify sample type (exome,genome,etc)")
    parser.add_argument('-k','--capture_kit', dest='capture_kit',
                        help="Specify capture kit ")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', action='version',
                        version='%(prog)s v3.0')
    args=parser.parse_args()
    if args.sample_name and (bool(args.sample_type) == False or  bool(args.capture_kit) == False) :
        parser.error('--sample_name and --sample_type and --capture_kit must be given together')
    return args

if __name__ == '__main__':
    main()
