"""Create a Illumina bcl2fastq v2.X.X compatibile sequencing
   sample sheet via queries from sequenceDB"""
   
import argparse
from utilities import create_sss_from_database,get_config,parse_run_parameters_xml

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fcillumid", dest='fcillumid', required=True,
                        help="Specify Illumina's Flowcell ID")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    args=parser.parse_args()
    return args

def main():
    config = get_config()
    args = parse_arguments()
    run_info_dict = parse_run_parameters_xml(args.fcillumid)
    if args.test:
        database = 'testDB'
    else:
        database = 'sequenceDB'
    create_sss_from_database(args.fcillumid,run_info_dict['machine'],run_info_dict,config,database)

if __name__ == '__main__':
    main()
