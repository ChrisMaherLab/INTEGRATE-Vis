import subprocess
import sys
import os
from subprocess import Popen, PIPE, STDOUT

def num_of_inter(ll1,ll2):
    h={}
    inter = []
    for ll in ll1:
        h[ll] = 1
    for ll in ll2:
        if ll in ll1:
            inter.append(ll)
    return len(set(inter))-1 # because \n

def get_number(bam,chr1,strand1,start1,end1,junction1,chr2,strand2,start2,end2,junction2):
    #change start or end
    if strand1 == "+":
        end1=junction1
    else:
        start1=str(int(junction1)-1000)
    if strand2 == "+":
        start2=str(int(junction2)-1000)
    else:
        end2=junction2
    #encompassing
    remote_chr2 = "=" if chr1==chr2 else chr2
    en_from1=subprocess.check_output("samtools view %s %s:%s-%s | awk '$7==\"%s\" && $8>%s && $8<%s {print $1}'" % (bam, chr1, start1, end1, remote_chr2, start2, end2), shell=True)
    en1=en_from1.split("\n")
    remote_chr1 = "=" if chr2==chr1 else chr1
    en_from2=subprocess.check_output("samtools view %s %s:%s-%s | awk '$7==\"%s\" && $8>%s && $8<%s {print $1}'" % (bam, chr2, start2, end2, remote_chr1, start1, end1), shell=True)
    en2=en_from2.split("\n")
    n_en=num_of_inter(en1,en2)
    #spanning and normal
    self_from1=subprocess.check_output("samtools view %s %s:%s-%s | awk '$7==\"%s\" && $8>%s && $8<%s {print $1}'" % (bam, chr1, start1, end1, "=", start1, end1), shell=True)
    sf1=self_from1.split("\n")
    self_from2=subprocess.check_output("samtools view %s %s:%s-%s | awk '$7==\"%s\" && $8>%s && $8<%s {print $1}'" % (bam, chr2, start2, end2, "=", start2, end2), shell=True)
    sf2=self_from2.split("\n")
    #spanning looks like encompassing
    spen_from1=subprocess.check_output("samtools view %s %s:%s-%s | awk '$7==\"%s\" && $8>%s && $8<%s {print $1}'" % (bam, chr1, str(int(junction1)-10 if int(junction1)-10>=0 else 0),str(int(junction1)+10), remote_chr2, start2, end2), shell=True)
    spen1=spen_from1.split("\n")
    remote_chr1 = "=" if chr2==chr1 else chr1
    spen_from2=subprocess.check_output("samtools view %s %s:%s-%s | awk '$7==\"%s\" && $8>%s && $8<%s {print $1}'" % (bam, chr2, str(int(junction2)-10 if int(junction2)-10>=0 else 0),str(int(junction2)+10), remote_chr1, start1, end1), shell=True)
    spen2=spen_from2.split("\n")
    n_sp1=num_of_inter(spen1,sf2)
    n_sp2=num_of_inter(spen2,sf1)
    return n_en,n_sp1,n_sp2


#TMPRSS2 chr21:41464551-41508159
#ERG     chr21:38380027-38498483
#41508080        -1      21      -1      38445621
#/gscmnt/gc2601/maherlab/jzhang/PRAD2/19_BAMs/e530a9089495423cbe28143097cd1114.bam
#get_number("/gscmnt/gc2601/maherlab/jzhang/PRAD2/19_BAMs/e530a9089495423cbe28143097cd1114.bam","21","41464551","41508159","41508080","21","38380027","38498483","38445621")
#41498118        -1      21      -1      38445621

#print get_number("/gscmnt/gc6127/research/jzhang/INTEGRATE-Vis-code/running_aligners/STAR2/star.Chimeric.sort.bam","21","-","41464551","41508159","41498118","21","-","38380027","38498483","38445621")
