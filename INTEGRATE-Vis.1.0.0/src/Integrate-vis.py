#!/usr/bin/python
import sys
import os
import getopt
from subprocess import Popen, PIPE, STDOUT
import vis_a_structure as pa_wrapper
import vis_b_domain as pb_wrapper
import vis_c_exon_exp as pc_wrapper
import vis_d_gene_exp as pd_wrapper

def usage(argv):
    print """
    Usages:

    Integrate-vis structure <parameters>
    Integrate-vis domain    <parameters>
    Integrate-vis exon-exp  <parameters>
    Integrate-vis gene-exp  <parameters>

    Version 1.0.0
          """
    print "    (A) Parameters for structure:"
    pa_wrapper.main([argv[0]]+["structure"]+argv[1:])
    print "    (B) Parameters for domain:"
    pb_wrapper.main([argv[0]]+["domain"]+argv[1:])
    print "    (C) Parameters for exon-exp:"
    pc_wrapper.main([argv[0]]+["exon-exp"]+argv[1:])
    print "    (D) Parameters for gene-exp:"
    pd_wrapper.main([argv[0]]+["gene-exp"]+argv[1:])

def main(argv):
    if len(argv)==1:
        usage(argv)
    elif len(argv)==2 and (argv[1]=="-h" or argv[1]=="--help"):
        usage(argv[0:1])
    else:
        if argv[1]=="structure":
            pa_wrapper.main(argv)
        elif argv[1]=="domain":
            pb_wrapper.main(argv)
        elif argv[1]=="exon-exp":
            pc_wrapper.main(argv)
        elif argv[1]=="gene-exp":
            pd_wrapper.main(argv)
        else:
            usage(argv[0:1])


if __name__ == '__main__':
    sys.exit(main(sys.argv))
