#!/usr/bin/python

import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
import time

prefix_usage = "vis_b_domain"

def usage():
    print
    print "    "+prefix_usage+" -b <fusions.bedpe> -s <sample-name> -d <ideogram.txt> -r <reference.fasta> -g <genes.gtf>"
    print """
    Requested parameters:

        -b/--fusion-bedpe   [string:    path to gene fusion BEDPE file (SMC-RNA 11 columns)    ]

        -s/--sample-name    [string:    sample name for the fusion bedpe file                  ]

        -d/--domain-table   [string:    path to domain information table file                  ]
      
        -r/--reference      [string:    path to refrence genome in FASTA                       ] 
       
        -g/--gene-model     [string:    path to GTF gene file                                  ]
       
    Optional parameters:

        -o/--output-dir     [string:    output directory.         Default: current directory   ]

        -k/--keep-tmp       [      :    keep tmp directory.       Default: not keeping tmp     ]

    Use the following to only plot one gene fusion:    

        -5/--5-prime        [string:    5' gene partner Id                                     ]
    
        -3/--3-prime        [string:    3' gene partner Id                                     ]

    Version:                1.0.0
          """

#parameters
fusion_bedpe = ''
domain_table = ''
gene_model = ''
output_dir = ''
is_rm_tmp=True
id_5p=''
id_3p=''
cur = ''
reference=''
sample_name = ''
bedpeAnnot = ''
genePred_file = ''

def setDefault():
    global output_dir
    output_dir = './'

def use_real_path():
    global output_dir
    global fusion_bedpe
    global domain_table
    global gene_model
    global reference
    output_dir=os.path.realpath(output_dir);
    fusion_bedpe=os.path.realpath(fusion_bedpe)
    domain_table=os.path.realpath(domain_table)
    gene_model=os.path.realpath(gene_model)
    reference=os.path.realpath(reference)


def initialSetupFile():
    global cur
    cur=os.path.dirname(os.path.abspath(__file__))

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hkb:d:g:o:s:r:5:3:",["help",
                                                              "keep-tmp",
                                                              "fusion-bedpe=",
                                                              "domain-table=",
                                                              "gene-model=",
                                                              "output-dir=",
                                                              "sample-name=",
                                                              "reference=",
                                                              "5-prime=",
                                                              "3-prime="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        #if opt in ("-h","--help"):
        #    usage()
        if opt in ("-k","--keep-tmp"):
            global is_rm_tmp
            is_rm_tmp = False
        elif opt in ("-5", "--5-prime"):
            global id_5p
            id_5p = arg
        elif opt in ("-3", "--3-prime"):
            global id_3p
            id_3p = arg
        elif opt in ("-b", "--fusion-bedpe"):
            global fusion_bedpe
            fusion_bedpe = arg
        elif opt in ("-d", "--domain-table"):
            global domain_table
            domain_table = arg
        elif opt in ("-r", "--reference"):
            global reference
            reference = arg
        elif opt in ("-g", "--gene-model"):
            global gene_model
            gene_model = arg
        elif opt in ("-s", "--sample-name"):
            global sample_name
            sample_name = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0775 ) 

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def make_annot():
    path, filename = os.path.split(gene_model)
    global genePred_file
    genePred_file=output_dir +'/tmp/'+filename+'.genePred'

    if os.path.exists(genePred_file)==False:
        #GTF to GenePred
        cmd = 'gtfToGenePred -genePredExt -geneNameAsName2 '+gene_model+' '+genePred_file # require gtfToGenePred 
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        print output

    path, filename = os.path.split(fusion_bedpe)
    global bedpeAnnot
    bedpeAnnot=output_dir+"/tmp/"+filename+"."+"annot"
    
    if os.path.exists(bedpeAnnot):
        return

    print "[Domain] Annotating fusion bedpe file..."
    di_file=cur+'/difile.txt'
    cmd = cur+"/"+"fusionBedpeAnnotator"+" --reference-file "+reference+" --gene-annotation-file "+genePred_file+" --di-file "+di_file \
    +" --input-file "+fusion_bedpe+" --output-file "+bedpeAnnot
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

def plot_one(g5,g3,p5,p3,is_in_frame,frame_string,domain_file,gene_model,output_dir,file_name):
    cmd = 'python '+cur+'/pb.py -5 '+g5+' -3 '+g3+' -f '+p5+' -t '+p3+' -s '+is_in_frame+' -p '+frame_string
    cmd = cmd+' -d '+domain_file+' -m '+gene_model+' -o '+file_name
    print "Command:",cmd
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

def rm_quoat(aa):
    if aa[0]=='\"' and aa[len(aa)-1]==";":
        return aa[1:len(aa)-2]
    else:
        return aa

trans_gene_dict = {}
gene_name_dict = {}

def set_trans_gene_dict():
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

def is_on_boundary(fusionAnnot):
    if fusionAnnot[11]=='1':
        return True
    else:
        return False

def get_gene_ids(fusionAnnot):
    #print fusionAnnot[16]
    tmp=fusionAnnot[16].split(";")
    tmp2=tmp[0].split("|")[0]
    t5=tmp2.split("(")[0]
    t3=''
    if tmp[1]!='':
        t3=tmp[1].split("|")[0]
    else:
        t3=tmp[2].split("|")[0]
    return trans_gene_dict[t5],trans_gene_dict[t3]

def get_gene_ids_2(fusionAnnot):
    fusionAnnot[16]=gen_col17(fusionAnnot)
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

def get_in_frame(fusionAnnot):
    return fusionAnnot[12]

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

def gen_col17(fusionAnnot):
    s1=fusionAnnot[8]
    s2=fusionAnnot[9]
    chr1=fusionAnnot[0]
    chr2=fusionAnnot[3]
    p1,p2=get_pos(fusionAnnot)
    tmp_file=gen_name(output_dir+"/tmp/",sample_name,"tmp",fusionAnnot,"txt")
    cmd = 'awk \'($9~"'+p1+'" || $10~"'+p1+'") && $2=="'+chr1+'" && $3=="'+s1+'"{print $1}\' '+genePred_file+' > '+tmp_file
    #print cmd
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    #print output
    new_col17=get_trans_str(tmp_file,"(0)")
    cmd = 'awk \'($9~"'+p2+'" || $10~"'+p2+'") && $2=="'+chr2+'" && $3=="'+s2+'"{print $1}\' '+genePred_file+' > '+tmp_file
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

def get_frame_str(fusionAnnot):
    if fusionAnnot[13]!='NA':
        return "\""+fusionAnnot[13]+" "+fusionAnnot[14]+" "+fusionAnnot[15]+" "+fusionAnnot[16]+"\""  
    else:
        return "\" 0 0 "+gen_col17(fusionAnnot)+"\""
    

def gen_name(output_dir, sample_name, fusion_name, fusionAnnot, suffix):
    s1=fusionAnnot[8]
    s2=fusionAnnot[9]  
    chr1=fusionAnnot[0]
    chr2=fusionAnnot[3]
    p1,p2=get_real_pos(fusionAnnot)
    return output_dir+"/"+sample_name+"."+fusion_name+"."+chr1+s1+p1+"."+chr2+s2+p2+"."+suffix

def plot_all():
    all_fusions = []  
    f=open(bedpeAnnot, "r")
    while True:
        line=f.readline()
        if line=="":
           break
        else:
           tmp=line.split("\t")
           all_fusions.append(tmp)
    f.close()
    for x in range(len(all_fusions)):
       #get gene ids
       if is_on_boundary(all_fusions[x]):
           id1,id2='',''
           if all_fusions[x][16]!="NA":
               id1,id2=get_gene_ids(all_fusions[x])
           else:
               id1,id2=get_gene_ids_2(all_fusions[x])
           p5,p3=get_pos(all_fusions[x])
           is_in_frame=get_in_frame(all_fusions[x])
           frame_string=get_frame_str(all_fusions[x])
           fusion_name=gene_name_dict[id1]+"--"+gene_name_dict[id2]
           file_name=gen_name(output_dir, sample_name, fusion_name, all_fusions[x], "domain.pdf")
           if id_5p!='' and id_3p!='':
               if id1==id_5p and id2==id_3p:
                   plot_one(id1,id2,p5,p3,is_in_frame,frame_string,domain_table,gene_model,output_dir,file_name)     
           else:
                   plot_one(id1,id2,p5,p3,is_in_frame,frame_string,domain_table,gene_model,output_dir,file_name)

def main(argv):
    t=time.time()
    setDefault()
    initialSetupFile()
    if len(argv)>1 and argv[1]=="domain":
        global prefix_usage
        prefix_usage="Integrate-vis domain"
        if(len(argv)==2):
            usage()
            return    
        getParameters(argv[2:])
    else:    
        getParameters(argv[1:])
    #print fusion_bedpe,gene_model,domain_table,reference,sample_name
    if fusion_bedpe=='' or gene_model=='' or domain_table=='' or reference=='' or sample_name=='':
        usage()
        exit(1);
    print "[Domain] Running pb_wrapper on",
    print sample_name,"..."
    use_real_path()
    make_dir(output_dir)
    make_dir(output_dir+'/tmp')
    print "[Domain] Making annotation file..."
    make_annot()
    print "[Domain] Setting gene name dictionaries..."
    set_trans_gene_dict()
    print "[Domain] Plotting figures..."
    print
    plot_all()
    remove_tmp()
    print "[Domain] Total time elapsed:", ("%0.2f" % (time.time()-t)), "seconds."

if __name__ == '__main__':
    sys.exit(main(sys.argv))
