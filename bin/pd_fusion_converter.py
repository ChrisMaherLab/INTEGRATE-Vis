import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
import time

def usage():
    print """
    pd_fusion_converter -s <comma separated sample names> -f <comma separated fusion bedpe files>
                        -r <reference> -g <genes.gtf>
    Requested parameters:

        -s/--sample-list      [string:    Comma separated sample names.                                 ]
        
        -f/--file-list        [string:    Comma separated fusion bedpe files (ordered the same as -s).  ]

        -r/--reference        [string:    path to refrence genome in FASTA                              ]

        -g/--gene-model       [string:    path to GTF gene file                                         ]

    Optional parameters:

        -o/--output-dir       [string:    output directory.         Default: current directory          ]

        -c/--cohort-name      [string:    default: cohort                                               ]

        -k/--keep-tmp         [      :    keep tmp directory.       Default: not keeping tmp            ]

    Use the following instead of -s and -f as an alternative:

        -a/--all-in-one-file  [string:    a tsv file, first column sample names, second column path     ]
                              [           to fusion bedpe files                                         ]

    Version:                  1.0.0
          """

#parameters
sample_list = ''
file_list = ''
output_dir = ''
cur = ''
all_file=''
cohort_name=''
reference=''
gene_model = ''
is_rm_tmp=True

def setDefault():
    global output_dir
    output_dir = './'

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir);
    global gene_model
    global reference
    output_dir=os.path.realpath(output_dir);
    gene_model=os.path.realpath(gene_model)
    reference=os.path.realpath(reference)

def initialSetupFile():
    global cur
    cur=os.path.dirname(os.path.abspath(__file__))
    global cohort_name
    cohort_name='cohort'

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hks:f:o:c:r:g:a:",["help",
                                                            "keep-tmp",
                                                            "sample-list=",
                                                            "file-list=",
                                                            "output-dir=",
                                                            "cohort-name=",
                                                            "reference=",
                                                            "gene-model=",
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
        elif opt in ("-a", "--all-in-one-file"):
            global all_file
            all_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-c", "--cohort-name"):
            global cohort_name
            cohort_name = arg
        elif opt in ("-k","--keep-tmp"):
            global is_rm_tmp
            is_rm_tmp = False
        elif opt in ("-r", "--reference"):
            global reference
            reference = arg
        elif opt in ("-g", "--gene-model"):
            global gene_model
            gene_model = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0775 )

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def make_annot_one(sample_name, fusion_bedpe):
    path, filename = os.path.split(gene_model)
    global genePred_file
    genePred_file=output_dir +'/tmp/'+filename+'.genePred'

    if os.path.exists(genePred_file)==False:
        #GTF to GenePred
        cmd = 'gtfToGenePred -genePredExt -geneNameAsName2 '+gene_model+' '+genePred_file # require gtfToGenePred
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        #print output

    path, filename = os.path.split(fusion_bedpe)
    bedpeAnnot=output_dir+"/tmp/"+filename+"."+"annot"
    
    if os.path.exists(bedpeAnnot):
        return

    print "    Annotating fusion bedpe file for ",sample_name,"..."
    di_file=cur+'/difile.txt'
    cmd = cur+"/"+"fusionBedpeAnnotator"+" --reference-file "+reference+" --gene-annotation-file "+genePred_file+" --di-file "+di_file \
    +" --input-file "+fusion_bedpe+" --output-file "+bedpeAnnot
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    #print output

    if not os.path.exists(bedpeAnnot):
        cmd = 'touch '+ bedpeAnnot
        Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

    return

def rm_quoat(aa):
    if aa[0]=='\"' and aa[len(aa)-1]==";":
        return aa[1:len(aa)-2]
    else:
        return aa

trans_gene_dict = {}
gene_name_dict = {}

def set_trans_gene_dict():
    print "  --Setting gene name dictionaries...",
    t2=time.time()
    dict_file=output_dir +'/tmp/dict.txt'
    if not os.path.exists(dict_file):
        f=open(gene_model, "r")
        line=f.readline()
        while True:
            line=f.readline()
            if line=="":
                break
            if line[0]=="#":
                continue
            else:
                tmp=line.split("\t")
                tmp2=tmp[8].split(" ")
                isT=False
                isG=False
                isN=False
                tt=''
                gg=''
                nn=''
                for x in range(len(tmp2)):
                    if tmp2[x]=='gene_id':
                        isG=True
                        continue
                    if tmp2[x]=='transcript_id':
                        isT=True
                        continue
                    if tmp2[x]=='gene_name':
                        isN=True
                        continue
                    if isG:
                        gg=rm_quoat(tmp2[x])
                        isG=False
                    if isT:
                        tt=rm_quoat(tmp2[x])
                        isT=False
                    if isN:
                        nn=rm_quoat(tmp2[x])
                        isN=False
                if tt!='' and gg!='':
                     trans_gene_dict[tt]=gg
                if gg!='' and nn!='':
                     gene_name_dict[gg]=nn
        f.close()
        f=open(dict_file, "w")
        for x in trans_gene_dict:
            f.write("%s\t%s\t%s\n" % (gene_name_dict[trans_gene_dict[x]],x,trans_gene_dict[x]))
        f.close()
    else:
        f=open(dict_file, "r")
        while True:
            line=f.readline()
            if line=="":
                break
            else:
                tmp=line.split("\t")
                trans_gene_dict[tmp[1]]=tmp[2][0:len(tmp[2])-1]
                gene_name_dict[tmp[2][0:len(tmp[2])-1]]=tmp[0]
        f.close()
    print "\t", ("%0.2f" % (time.time()-t2)), "seconds"

def is_on_boundary(fusionAnnot):
    if fusionAnnot[11]=='1':
        return True
    else:
        return False

def gen_col17_2(fusionAnnot):
    s1=fusionAnnot[8]
    s2=fusionAnnot[9]
    chr1=fusionAnnot[0]
    chr2=fusionAnnot[3]
    p1,p2=get_pos(fusionAnnot)
    tmp_file=gen_name(output_dir+"/tmp/",sample_name,"tmp",fusionAnnot,"txt")
    cmd = 'awk \''+p1+'>$4 && '+p1+'<$5 && $2=="'+chr1+'" && $3=="'+s1+'"{print $1}\' '+genePred_file+' > '+tmp_file
    #print cmd
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    #print output
    new_col17=get_trans_str(tmp_file,"(0)")
    cmd = 'awk \''+p2+'>$4 && '+p2+'<$5 && $2=="'+chr2+'" && $3=="'+s2+'"{print $1}\' '+genePred_file+' > '+tmp_file
    #print cmd
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    #print output
    new_col17=new_col17+";"+get_trans_str(tmp_file,"")+";"
    cmd="rm -f "+tmp_file
    #print cmd
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    #print output
    return new_col17

def get_gene_ids(fusionAnnot):
    tmp=fusionAnnot[16].split(";")
    tmp2=tmp[0].split("|")[0]
    t5=tmp2.split("(")[0]
    t3=''
    if tmp[1]!='':
        t3=tmp[1].split("|")[0]
    else:
        t3=tmp[2].split("|")[0]
    if t5 in trans_gene_dict and t3 in trans_gene_dict:
        return trans_gene_dict[t5],trans_gene_dict[t3]
    else:
        return '',''    

def get_gene_ids_2(fusionAnnot):
    fusionAnnot[16]=gen_col17_2(fusionAnnot)
    return get_gene_ids(fusionAnnot)

def get_pos(fusionAnnot):
    s1=fusionAnnot[8]
    s2=fusionAnnot[9]
    c1=''
    if s1=="+":
       c1=2
    else:
       c1=1
    c2=''
    if s2=="+":
       c2=4
    else:
       c2=5
    return fusionAnnot[c1],fusionAnnot[c2]

def get_real_pos(fusionAnnot):
    s1=fusionAnnot[8]
    s2=fusionAnnot[9]
    c1=''
    if s1=="+":
       c1=2
    else:
       c1=1
    c2=''
    if s2=="+":
       c2=4
    else:
       c2=5
    p5=int(fusionAnnot[c1])
    p3=int(fusionAnnot[c2])
    if c1==1:
        p5=p5+1
    if c2==4:
        p3=p3+1
    return str(p5),str(p3)

def get_trans_str(tmp_file, suffix):
    tran_str=""
    f=open(tmp_file, "r")
    while True:
        line=f.readline()
        if line=="":
           break
        else:
           line=line[0:len(line)-1]
           tran_str=tran_str+line+suffix+"|"
    f.close()
    tran_str=tran_str[0:len(tran_str)-1]
    return tran_str

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

def gen_name(output_dir, sample_name, fusion_name, fusionAnnot, suffix):
    s1=fusionAnnot[8]
    s2=fusionAnnot[9]
    chr1=fusionAnnot[0]
    chr2=fusionAnnot[3]
    p1,p2=get_real_pos(fusionAnnot)
    return output_dir+"/"+sample_name+"."+fusion_name+"."+chr1+s1+p1+"."+chr2+s2+p2+"."+suffix

sample_name = ''
id_sample_has_fusion_dict = {}
id_pairs = {}

def make_annot_all():
    print "  --Making annotation files..."
    for i in range(len(s_list)):
        make_annot_one(s_list[i],f_list[i])

def get_fusion_matrix():
    print "  --Processing all gene fusion bedpe files..."
    for i in range(len(s_list)):
        t2=time.time()
        all_fusions = []
        path, filename = os.path.split(f_list[i])
        bedpeAnnot=output_dir+"/tmp/"+filename+"."+"annot"
        f=open(bedpeAnnot, "r")
        while True:
            line=f.readline()
            if line=="":
               break
            else:
               tmp=line.split("\t")
               all_fusions.append(tmp)
        f.close()
        sample_name=s_list[i]#global to use older code with it as global
        for x in range(len(all_fusions)):
           id1,id2='',''
           if all_fusions[x][16]!="NA":
               id1,id2=get_gene_ids(all_fusions[x])
           else:
               id1,id2=get_gene_ids_2(all_fusions[x])
           if id1!='' and id2!='':
               id_sample_has_fusion_dict[id1+"__"+id2+"__"+s_list[i]] = 1
               id_pairs[id1+"__"+id2] = 1
        print "    Sample:",s_list[i],"took\t", ("%0.2f" % (time.time()-t2)), "seconds"
    t2=time.time()
    print "  --Creating fusion matrix...",
    #write to file
    f=open(output_dir+"/"+cohort_name+".fusions.tsv","w")
    f.write("#%s\t%s\t" % ("5p_id","3p_id"))
    f.write("%s\n" % ('\t'.join(map(str, s_list))))
    for r in id_pairs:#row
        tmp=r.split("__")
        id1=tmp[0]
        id2=tmp[1]
        f.write("%s\t%s\t%s" % (gene_name_dict[id1]+"--"+gene_name_dict[id2],id1,id2))
        for c in range(len(s_list)):#column
            value="0"
            if id1+"__"+id2+"__"+s_list[c] in id_sample_has_fusion_dict:
                value="1"
            f.write("\t%s" % (value))
        f.write("\n")
    f.close()       
    print "\t", ("%0.2f" % (time.time()-t2)), "seconds"   
 
def main(argv):
    t=time.time()
    setDefault()
    initialSetupFile()
    getParameters(argv[1:])
    if not ((sample_list!='' and file_list!='') or all_file!='') or gene_model=='' or reference=='':
        usage()
        exit(1);
    print "Running fusion_converter..."
    use_real_path()
    make_dir(output_dir)
    make_dir(output_dir+'/tmp')
    t3=time.time()
    
    print "  --Getting lists of samples and gene fusion bedpes...",
    use_real_path()
    make_dir(output_dir)
    if all_file!='':
        get_lists_from_file()
    else:
        get_lists_from_str()
    print "\t", ("%0.2f" % (time.time()-t3)), "seconds"
    
    set_trans_gene_dict()
    make_annot_all()
    get_fusion_matrix()    
    remove_tmp()
    print "Time elapsed:", ("%0.2f" % (time.time()-t)), "seconds"

if __name__ == '__main__':
    sys.exit(main(sys.argv))
