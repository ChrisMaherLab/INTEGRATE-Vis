#!/usr/bin/python

import sys
import os
import os.path
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
import time
import re
import geneclass
import shutil

prefix_usage = "vis_e_fcirc"

def usage():
    print
    print "    "+prefix_usage+" -f <fcircrna.txt> -s <sample-name> -d <ideogram.txt> -r <reference.fasta> -g <genes.gtf>"  
    print """
    Requested parameters:

        -f/--fcircrna       [string:    path to fcircRNA file                                  ]
        
        -s/--sample-name    [string:    sample name for the Fcirc file                         ]
        
        -r/--reference      [string:    reference fasta                                        ]
        
        -g/--gene-model     [string:    path to GTF gene file                                  ]
        
    Optional parameters:
        
        -o/--output-dir     [string:    output directory.      Default: current directory      ]
        
        -k/--delete-tmp     [      :    delete tmp directory (flag).    Default: keep tmp      ]
        
        -b/--bam-file       [string:    path to bam file for plotting read count information   ]
        
        -l/--legend         [      :    Read support legend (flag).   Default: no legend       ]
        
                             
    Version:                1.1.0 
        """


def setDefault():
    global output_dir
    output_dir = './'
    
def use_real_path():
    global output_dir
    global fcircrna_file
    global ideogram
    global gene_model
    global reference
    output_dir = os.path.realpath(output_dir)
    fcircrna_file = os.path.realpath(fcircrna_file)
    ideogram = os.path.realpath(ideogram)
    gene_model = os.path.realpath(gene_model)
    reference = os.path.realpath(reference)

    
def initialSetupFile():
    global cur
    cur=os.path.dirname(os.path.abspath(__file__))
    
def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hklf:d:g:o:s:r:5:3:b:",["help",
                                                              "delete-tmp",
                                                              "fcircrna=",
                                                              "gene-model=",
                                                              "output-dir=",
                                                              "sample-name=",
                                                              "reference=",
                                                              "5-prime=",
                                                              "3-prime=",
                                                              "bam-file="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
        
    global rm_tmp 
    rm_tmp = False
    #is_tmp_tmp = True
    #global show_legend 
    #show_legend = False
    
    for opt, arg in opts:
        if opt in ("-k","--delete-tmp"):
            rm_tmp = True
        elif opt in ("-l","--legend"):
            global show_legend
            show_legend = True
        elif opt in ("-5","--5-prime"):
            global id_5p
            id_5p = arg
        elif opt in ("-3","--3-prime"):
            global id_3p
            id_3p = arg
        elif opt in ("-f", "--fcircrna"):
            global fcircrna_file
            fcircrna_file = arg
        elif opt in ("-d", "--ideogram"):
            global ideogram
            ideogram = arg
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
        elif opt in ("-b", "--bam-file"):
            global bam_file
            bam_file = arg
    
def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path, 0775)
        
def parse_fcirc(pos5, pos3):
    chrm5 = re.findall("([0-9A-Z]*):", pos5)[0]
    chrm5strand = pos5[-1]
    chrm5posa = re.findall(":([0-9]*):", pos5)[0]
    chrm5posb = chrm5posa
    if chrm5strand == "+":
        chrm5posa = "-1"
    else:
        chrm5posb = "-1"
    chrm3 = re.findall("([0-9A-Z]*):", pos3)[0]
    chrm3strand = pos3[-1]
    chrm3posa = re.findall(":([0-9]*):", pos3)[0]
    chrm3posb = chrm3posa
    if chrm3strand == "+":
        chrm3posb = "-1"
    else:
        chrm3posa = "-1"
    return chrm5,chrm5posa,chrm5posb,chrm5strand,chrm3,chrm3posa,chrm3posb,chrm3strand
        
def fcircrna_to_bedpe():
    # Creates a standard bedpe file, with an id in the score column
    global fcircrna_file
    global tmp_bedpe
    tmp_bedpe = output_dir + '/tmp/' + sample_name + '.fcirc.bedpe'
    if(not os.path.isfile(tmp_bedpe)):
        bedpe_lines = []
        i = 1 # Used as an id for the fcircRNA
        with open(fcircrna_file, "r") as f:
            for line in f:
                if line[0] != "#":
                    line = line.rstrip()
                    fields = line.split("\t")
                    pos5 = fields[4]
                    pos3 = fields[5]
                    backstart = fields[2]
                    backend = fields[3]
                    fusion = fields[1]
                    reads = "0" # Not used for plot
                    genea = ""
                    geneb = ""
                    if "--" in fusion:
                        genea = re.findall("(.*)--", fusion)[0]
                        geneb = re.findall("--(.*)", fusion)[0]
                    elif "::" in fusion:
                        genea = re.findall("(.*)::", fusion)[0]
                        geneb = re.findall("::(.*)", fusion)[0]
                    else:
                        genea = re.findall("(.*)>>", fusion)[0]
                        geneb = re.findall(">>(.*)", fusion)[0]
                    chrm5,chrm5posa,chrm5posb,chrm5strand,chrm3,chrm3posa,chrm3posb,chrm3strand = parse_fcirc(pos5,pos3)
                    fusion_string = "\t".join([chrm5,chrm5posa,chrm5posb,
                                           chrm3,chrm3posa,chrm3posb,
                                           genea+">>"+geneb,str(i)+"-fusion",
                                           chrm5strand,chrm3strand,reads]) + "\n"
                    bedpe_lines.append(fusion_string)
                    chrm5,chrm5posa,chrm5posb,chrm5strand,chrm3,chrm3posa,chrm3posb,chrm3strand = parse_fcirc(backend,backstart)
                    backsplice_string = "\t".join([chrm5,chrm5posa,chrm5posb,
                                               chrm3,chrm3posa,chrm3posb,
                                               geneb+">>"+genea,str(i)+"-backsplice",
                                               chrm5strand,chrm3strand,reads]) + "\n"
                    bedpe_lines.append(backsplice_string)
                    i += 1
        f.close()
    
        with open(tmp_bedpe,"w") as f:
            for line in bedpe_lines:
                f.write(line)
                f.close()
    else:
        print "Modified bedpe already exists"
    
def create_gene(line):
    fields = line.split("\t")
    geneid = fields[0]
    genename = fields[11]
    exonstarts = fields[8]
    exonends = fields[9]
    strand = fields[2]
    gene = geneclass.gene(genename,strand,geneid)
    gene.addExons(exonstarts,exonends)
    return gene
    
def getBestExon(pos, gene_name, genes, is5p):
    # Get the exon that is 1) closest to breakpoint and
    # 2) shortest in length
    if gene_name not in genes.keys():
        return -1,-1
    transcripts = genes[gene_name]
    best = 100000000
    shortestExon = 100000000
    bestTranscript = None
    bestExon = None
    for transcript in transcripts:
        exons = transcript.exons
        strand = transcript.strand
        exonCount = 1
        for exon in exons:
            #if is5p and is +, then you want the pos near the right side of the exon
            #if is 3p and is -, then you also want the pos near the right side
            if (is5p and strand=="+") or (not is5p and strand=="-"):
                dist = int(exon[1]) - int(pos)
                if(dist<0):
                    dist = dist * -1
                if dist < best:
                    best = dist
                    exonLength = int(exon[1]) - int(exon[0])
                    if exonLength < shortestExon:
                        shortestExon = exonLength
                        bestTranscript = transcript.geneid
                        bestExon = exonCount
            #If is5p and is -, then you want the post to be near the left side of the exon
            elif (is5p and strand=="-") or (not is5p and strand=="+"):
                dist = int(exon[0]) - int(pos)
                if(dist<0):
                    dist = dist * -1
                if dist < best:
                    best = dist
                    exonLength = int(exon[1]) - int(exon[0])
                    if exonLength < shortestExon:
                        shortestExon = exonLength
                        bestTranscript = transcript.geneid
                        bestExon = exonCount
            exonCount += 1
    return bestExon, bestTranscript
        
                

def make_annot():
    ###
    #Create genePred_file
    path, filename = os.path.split(gene_model)
    global genePred_file
    genePred_file=output_dir + "/tmp/" + filename + '.genePred'
    
    if os.path.exists(genePred_file)==False:
        print "[Circ] Create genePred file..."
        cmd = 'gtfToGenePred -genePredExt -geneNameAsName2 '+gene_model+' '+genePred_file
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        print output
    
    # Get list of fused genes
    
    fused_genes = []
    with open(tmp_bedpe, "r") as f:
        for line in f:
            fields = line.split("\t")
            genea = re.findall("(.*)>>", fields[6])[0]
            geneb = re.findall(">>(.*)", fields[6])[0]
            fused_genes.append(genea)
            fused_genes.append(geneb)
    f.close()
    
    # Gather exon info on fused genes
    genes = {}
    with open(genePred_file, "r") as f:
        for line in f:
            fields = line.split("\t")
            genename = fields[11]
            if genename in fused_genes:
                gene = create_gene(line)
                if genename in genes.keys():
                    transcripts = genes[genename]
                    transcripts.append(gene)
                    genes[genename] = transcripts
                else:
                    genes[genename] = [gene]
    f.close()
    
    # Create custom annotation file
    print "[Circ] Creating custom annotation file..."
    outstrings = []
    with open(tmp_bedpe, "r") as f:
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            posstarta = fields[1]
            posenda = fields[2]
            posa = posstarta
            if posenda > posstarta:
                posa = posenda
            posstartb = fields[4]
            posendb = fields[5]
            posb = posstartb
            if posendb > posstartb:
                posb = posendb
            exonA,transcriptA = getBestExon(posa,re.findall("(.*)>>",fields[6])[0],genes,True)
            exonB,transcriptB = getBestExon(posb,re.findall(">>(.*)",fields[6])[0],genes,False)
            line += "\t" + transcriptA + ":" + str(exonA) + "\t" + transcriptB + ":" + str(exonB) + "\n"
            outstrings.append(line)
    f.close()
    
    global tmp_annot
    tmp_annot = output_dir + "/tmp/"  + sample_name + ".fcirc.annot"
    if(not os.path.isfile(tmp_annot)):
        with open(tmp_annot, "w") as f:
            for s in outstrings:
                f.write(s)
                f.close()
    else:
        print "Custom annotation file already exists"


def plot_all():
    all_events = []
    with open(tmp_annot, "r") as f:
        while True:
            fusion = f.readline()
            if fusion=="":
                break
            else:
                backsplice = f.readline()
                all_events.append([fusion.rstrip(),backsplice.rstrip()])
    f.close()
    counter = 1
    for event in all_events:
        fusion = event[0]
        backsplice = event[1]
        fields = fusion.split("\t")
        fusionname = fields[6]
        genea = re.findall("(.*)>>", fusionname)[0]
        geneb = re.findall(">>(.*)", fusionname)[0]
        chrma = fields[0]
        chrmb = fields[3]
        posa = fields[1]
        if int(posa) <= 0:
            posa = fields[2]
        posb = fields[4]
        if int(posb) <= 0:
            posb = fields[5]
        stranda = fields[8]
        annota = fields[11]
        strandb = fields[9]
        annotb = fields[12]
        fields = backsplice.split("\t")
        bposa = fields[1]
        if int(bposa) <= 0:
            bposa = fields[2]
        bposb = fields[4]
        if int(bposb) <= 0:
            bposb = fields[5]
        bannota = fields[11]
        bannotb = fields[12]
        shortname = genea +  "-" + geneb  + ".fcircRNA." + str(counter) + ".pdf"
        filename = output_dir + "/" + shortname
        mybamfile = "none"
        if 'bam_file' in globals():
            mybamfile = bam_file
        showlegend = "0"
        if 'show_legend' in globals() and 'bam_file' in globals():
            showlegend = "1"
        print "[Circ] Creating fcircRNA plot number " + str(counter)
        plot_one(filename,genea,geneb,
                 stranda,strandb,
                 chrma,posa,chrmb,posb,
                 annota,annotb,
                 bposa,bposb,
                 bannota,bannotb,
                 mybamfile,showlegend)
        counter += 1
        

def plot_one(filename, geneA, geneB, strandA, strandB,
             chrmA, fusionPosA, chrmB, fusionPosB,
             annotA, annotB, backposA, backposB,
             backannotA, backannotB, mybamfile, showlegend):
    
    cmd = 'python '+cur+"/pe.py "+filename+" " + geneA + " " + geneB
    cmd += " " + strandA + " " + strandB
    cmd += " " + chrmA + " " + chrmB
    cmd += " " + fusionPosA + " " + fusionPosB
    cmd += " " + annotA + " " + annotB
    cmd += " " + backposA + " " + backposB
    cmd += " " + backannotA + " " + backannotB + " " + ideogram + " " + genePred_file
    cmd += " " + mybamfile + " " + showlegend
    
    #print cmd
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output
        

def main(argv):
    t=time.time()
    setDefault()
    initialSetupFile()
    if len(argv)>1 and argv[1]=="fcircrna":
        global prefix_usage
        prefix_usage="Integrate-vis fcircrna"
        if(len(argv)==2):
            usage()
            return
        getParameters(argv[2:])
    else:
        getParameters(argv[1:])
    use_real_path()
    if fcircrna_file=='' or gene_model=='' or ideogram=='' or reference=='' or sample_name=='':
        usage()
        exit(0)
    print "[Circ] Running Integrate-vis fcirc on",
    print sample_name,"..."
    
    make_dir(output_dir)
    make_dir(output_dir+'/tmp')
    print "[Circ] Converting fcircRNA file to modified bedpe..."
    fcircrna_to_bedpe()
    print "[Circ] Making annotation file..."
    make_annot()
    print "[Circ] Plotting figures..."
    plot_all()
    if rm_tmp:
        print "[Circ] Deleting tmp directory..."
        tmp_dir = output_dir + "/tmp"
        shutil.rmtree(tmp_dir)


    print
    print "[Circ] Total time elapsed:", ("%0.2f" % (time.time()-t)), "seconds."



if __name__ == '__main__':
    sys.exit(main(sys.argv))
