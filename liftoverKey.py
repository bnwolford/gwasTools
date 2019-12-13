#!/usr/bin/python3

#Script by Brooke Wolford
## December 13, 2019

import gzip
import sys
import re
import argparse
import linecache
from collections import OrderedDict

### Written for Python 3

argparser = argparse.ArgumentParser('This script creates a key for old reference genome coords to new reference genome coordinates after liftover\n')
argparser.add_argument("--old_bed",type=str,help="Old coordinate bed file that was submitted to liftOver")
argparser.add_argument("--new_bed",type=str,help="New coordinate bed file from liftOver")
argparser.add_argument("--unmapped",type=str,help="Unmapped file from liftOver")
argparser.add_argument("--prefix",type=str,help="Prefix to add to .txt for output")

def open_zip(f):
    """open files"""
    if ".gz" in f:
        command=gzip.open(f,"rt")
        print("Opening gzipped file %s\n" % f,file=sys.stderr)
    elif f == "-":
        sys.exit("Cannot read file from stdin\n")
    else:
        command=open(f,"rt")
        print("Opening file %s\n" % f, file=sys.stderr)
        return command
    
def read_bed(f):
    command=open_zip(f)
    my_dict=OrderedDict()
    with command as bed_file:
        for line in bed_file:
            ls=line.rstrip()
            if not ls.startswith("#"): #ignore comment lines 
                line_list=ls.split("\t")
                my_dict[":".join([line_list[0],line_list[1],line_list[2]])]=line_list
    bed_file.close()
    return(my_dict)


def make_key(new,old,unmap,prefix):
    command=open_zip(old)
    ucount=0 #unampped count
    lcount=0 #line count
    results_filename=".".join([prefix,"txt"])
    with command as new_file:
        for line in new_file:
            ls=line.rstrip()
            line_list=ls.split("\t")
            old_coord=":".join([line_list[0],line_list[1],line_list[2]])
            with open(results_filename, "w") as out:
                out.write("\t".join("old_chrom","old_posStart","old_posEnd","old_coord","new_chrom","new_posStart","new_posEnd","new_coord"))
                if old_coord in unmap: #old coordiante didn't get mapped so NA for key
                    ucount+=1
                    lcount+=1
                    out.write("\t".join(["\t".join(line_list),old_coord,"\t".join(["NA","NA","NA","NA"])]))
                else:
                    lcount+=1
                    nl=linecache.getline(new,lcount-ucount) #new line
                    nl_list=nl.rstrip().split("\t")
                    out.write("\t".join(["\t".join(line_list),old_coord,"\t".join(nl_list),":".join([nl_list[0],nl_list[1],nl_list[2]])]))

    new_file.close()
    out.close()

###############################
######## MAIN ###############
############################
            
if __name__ == '__main__':
    args = argparser.parse_args()
    print(args, file=sys.stderr)
    
    unmapped=read_bed(args.unmapped)
    make_key(args.new_bed,args.old_bed,unmapped,args.prefix)

