import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
import time

def usage():
    print """
    pd_expression_converter -s <comma separated sample names> -f <comma separated expression file list>
                            -g <column number for gene Id> -e <column number for expression value>

    Requested parameters:

        -s/--sample-list                [string:    Comma separated sample names.                               ]
        
        -f/--file-list                  [string:    Comma separated expression files (ordered the same as -s).  ]

        -g/--gene-id-column             [int:       The column of gene Id in a expression file.                 ]

        -e/--expression-value-column    [int:       The column of experssion value in a expression file.        ]
    
    Optional parameters:

        -o/--output-dir                 [string:    output directory.         Default: current directory        ]

        -c/--cohort-name                [string:    default: cohort                                             ]

    Use the following instead of -s and -f as an alternative:

        -a/--all-in-one-file            [string:    a tsv file, first column sample names, second column path   ]
                                        [           to expression files                                         ]

    Version:                1.0.0
          """

#parameters
sample_list = ''
file_list = ''
output_dir = ''
cur = ''
gene_id_column=0
expression_column=0
all_file=''
cohort_name=''

def setDefault():
    global output_dir
    output_dir = './'
    global cohort_name
    cohort_name='cohort'

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir);

def initialSetupFile():
    global cur
    cur=os.path.dirname(os.path.abspath(__file__))

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hs:f:g:e:o:c:a:",["help",
                                                           "sample-list=",
                                                           "file-list=",
                                                           "gene-id-column=",
                                                           "expression-value-column=",
                                                           "output-dir=",
                                                           "cohort-name=",
                                                           "all-in-one-file="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(0)
        elif opt in ("-s", "--sample-list"):
            global sample_list
            sample_list = arg
        elif opt in ("-f", "--file-list"):
            global file_list
            file_list = arg
        elif opt in ("-g", "--gene-id-column"):
            global gene_id_column
            gene_id_column = int(arg)
        elif opt in ("-e", "--expression-column"):
            global expression_column
            expression_column = int(arg)
        elif opt in ("-a", "--all-in-one-file"):
            global all_file
            all_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-c", "--cohort-name"):
            global cohort_name
            cohort_name = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0775 )

s_list = []
f_list = []

def get_lists_from_file():
    f=open(all_file, "r")
    while True:
        line=f.readline()
        if line=="":
            break
        if line[0]=="#":
            continue
        else:
            tmp=line.split("\t")
            s_list.append(tmp[0])
            f_list.append(tmp[1][0:len(tmp[1])-1])
    f.close()

def get_lists_from_str():
    tmp1=sample_list.split(",")
    tmp2=file_list.split(",") 
    for x in range(len(tmp1)):
        s_list.append(tmp1[x])
        f_list.append(tmp2[x])

def generate_expression_file():
    #get all the gene ids from the first file
    id_list = []
    f=open(f_list[0],"r")
    while True:
        line=f.readline()
        if line=="":
            break
        if line[0]=="#":
            continue
        else:
            tmp=line.split("\t")
            if gene_id_column-1==len(tmp):
                tmp[gene_id_column-1]=tmp[gene_id_column-1][0:len(tmp[gene_id_column-1]-1)]
            id_list.append(tmp[gene_id_column-1])
    f.close()
    #get each column of values, append in matrix
    mat = []
    for x in range(len(f_list)):
        f=open(f_list[x],"r")
        h = {}
        while True:
            line=f.readline()
            if line=="":
                break
            if line[0]=="#":
                continue
            else:
                tmp=line.split("\t")
                if expression_column > len(tmp) or gene_id_column > len(tmp):
                    continue
                tmp[len(tmp)-1]=tmp[len(tmp)-1][0:len(tmp[len(tmp)-1])-1]
                h[tmp[gene_id_column-1]] = tmp[expression_column-1]
        f.close()
        #get 
        v=[]
        for t in range(len(id_list)):#here assume all files were from the runs with the same parameters, so all id exist
            id1=id_list[t]
            v.append(h[id1])
        mat.append(v)
    #write to file
    f=open(output_dir+"/"+cohort_name+".gene_expression.tsv","w")
    f.write("#%s\n" % ('\t'.join(map(str, s_list)))) 
    for r in range(len(id_list)):#row
        f.write("%s" % (id_list[r]))
        for c in range(len(mat)):#column
            f.write("\t%s" % (mat[c][r]))
        f.write("\n")
    f.close()

def main(argv):
    t=time.time()
    setDefault()
    initialSetupFile()
    getParameters(argv[1:])
    #print sample_list,file_list,all_file,gene_id_column,expression_column
    if not ((sample_list!='' and file_list!='') or all_file!='') or gene_id_column==0 or expression_column==0:
        usage()
        exit(1);

    print "Running expression_converter..."
    use_real_path()
    make_dir(output_dir)
    if all_file!='':
        get_lists_from_file()
    else:
        get_lists_from_str()
    generate_expression_file()
    print "Time elapsed:", ("%0.2f" % (time.time()-t)), "seconds"

if __name__ == '__main__':
    sys.exit(main(sys.argv))
