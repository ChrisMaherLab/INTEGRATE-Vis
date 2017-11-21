import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import subprocess
import re
import getopt
import sys
import os
import time
import pdb
import readsNum_getter as rng


def usage():
    print """
    panel_a -5 <five_p_gene_ID_or_transcript_ID> -3 <three_p_gene_ID_or_transcript_ID> -f <five_prime_gene_fusion_point> -t <three_prime_gene_fusion_point> -d <ideogram_dir> -p <transcript_frames> -m <gene_model_dir> -o <output_dir>
    Parameters:

        -5/--five_p_gene_ID      [Requested]
                                 [string:   5' gene ID.                                                             ]
        -3/--three_p_gene_ID     [Requested]
                                 [string:   3' gene ID.                                                             ]
        -f/--five_prime_gene_fusion_point   [Requested]
                                 [string:   5' gene fusion coordinate                                               ]
        -t/--three_prime_gene_fusion_point   [Requested]
                                 [string:   3' gene fusion coordinate                                               ]
        -d/--ideogram_dir        [Requested]
                                 [string:   ideogram file                                                           ]
        -m/--gene_model_dir      [Requested]
                                 [string:   provide your gene model                                                 ]
        -o/--output-dir          [Optional]
                                 [string:   output directory.                                                       ]
        -1/--five_p_trans_ID     [Optional]
                                 [string:   5' transcript ID. Only choose this transcrtipt.                         ]
        -2/--three_p_trans_ID    [Optional]
                                 [string:   3' transcript ID. Only choose this transcrtipt.                         ]
        -p/--transcript_frames   [Requested]
                                 [string:   stop info and shift frame status for different transcripts,
                                            a sub line of the annotated bedpe, seperated by space, comma, and '|'   ]
        -b/--bam-file            [Optional]
                                 [string:   path to bam file for plotting number of encompassing and spanning reads ]

    Version:                     1.0.0
          """

scale_global = 1.31 #make the chrom size shorter
scale_5 = 1.0
scale_3 = 1.0

class PANEL_A:
    # provided
    gene_ID_3 = ''
    gene_ID_5 = ''
    ideogram_dir = ''
    gene_model_dir =''
    gene_model_dir_2 = ''
    fusion_point_5 =''
    fusion_point_3 = ''
    output_dir = './'

    trans_ID_5 = 'None'
    trans_ID_3 = 'None'

    # fixed
    centrolmere_x = 1
    centrolmere_y = 0.2
    ch_regions = []
    pad_size = 0.025 #changed from 0.02 to 0.025
    centrolmere_radius = 0.025 # changed from 0.03 to 0.025
    arm_width = 0.01
    line_width = 1.5
    ch_graph_start_coord = 0.075 # 0.5 changed to 0.075: move to left
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    exon_strands_y = 0.75
    exon_strand_x_5 = 0.5
    exon_strand_x_3 = 1.5
    exon_strand_radius = 0.4
    exon_box_width = 0.05
    exon_box_length = 0.18
    exon_seperation = 0.25
    fusion_trans_x = 1
    fusion_trans_y = 1.1

    # reads

    reads_y_bottom = 1.25
    reads_x_middle = 0.995
    read_len = 0.15
    read_insert = 0.1
    read_width = 0.015
    increase = 0.04

    # computed
    TR_frames = ''
    in_frame_TRs_5 = []
    in_frame_TRs_3 = []
    out_frame_TRs_3 = []
    gene_name_3 = ''
    gene_name_5 = ''
    TR_ID_3 = ''
    TR_ID_5 = ''
    strand_3 = ''
    strand_5 = ''
    exons_3 = []
    exons_5 = []
    TR_start_5 = 0
    TR_start_3 = 0
    TR_end_3 = 0
    TR_end_5 = 0
    gene_band_code_3 = ''
    gene_band_code_5 = ''
    gene_region_3 = []
    gene_region_5 = []
    acen_length_5 = 0
    acen_length_3 = 0

    ch_length = 0
    arm_length_left = 0
    arm_length_right = 0
    ch_3 = 0
    ch_5 = 0
    fusion_exons_3 = []
    fusion_exons_5 = []
    fusion_center_exon_3 = 0
    fusion_center_exon_5 = 0

    bam_file = ''

    def grep_gene(self, gene_ID_5, gene_ID_3, output_dir, suffix):
        if self.trans_ID_5!='None':
            subprocess.check_output("grep %s %s | awk '$3==\"gene\"{print}' > %s" % (gene_ID_5, self.gene_model_dir, output_dir+"."+suffix), shell = True)
            subprocess.check_output("grep %s %s >> %s" % (self.trans_ID_5, self.gene_model_dir, output_dir+"."+suffix), shell = True)
        else:
            subprocess.check_output("grep %s %s > %s" % (gene_ID_5, self.gene_model_dir, output_dir+"."+suffix), shell = True)
        if self.trans_ID_3!='None':
            subprocess.check_output("grep %s %s | awk '$3==\"gene\"{print}' >> %s" % (gene_ID_3, self.gene_model_dir, output_dir+"."+suffix), shell = True)
            subprocess.check_output("grep %s %s >> %s" % (self.trans_ID_3, self.gene_model_dir, output_dir+"."+suffix), shell = True)
        else:
            subprocess.check_output("grep %s %s >> %s" % (gene_ID_3, self.gene_model_dir, output_dir+"."+suffix), shell = True)
        self.gene_model_dir_2=output_dir+"."+suffix

    def rm_tmp(self, output_dir, suffix):
        os.remove(output_dir+"."+suffix)

    def __init__(self,argv):
        try:
            opts, args = getopt.getopt(argv[1:],"h5:3:f:t:d:m:o:1:2:p:b:", ["help", "gene_ID_5=", "gene_ID_3=", "fusion_point_5=", "fusion_point_3=", "ideogram_dir=", "gene_model_dir=", "output_dir=", "five_p_trans_ID=", "three_p_trans_ID=", "transcript_frames=", "bam-file="])
            #print(opts, args)

        except getopt.GetoptError:
             usage()
             sys.exit(1)

        for opt, arg in opts:
            if opt in ("-h","--help"):
                usage()
                sys.exit(0)
            elif opt in ("-5", "--gene_ID_5"):
                self.gene_ID_5 = arg
            elif opt in ("-3", "--gene_ID_3"):
                self.gene_ID_3 = arg
            elif opt in ("-f", "--fusion_point_5"):
                self.fusion_point_5 = int(arg)
            elif opt in ("-t", "--fusion_point_3"):
                self.fusion_point_3 = int(arg)
            elif opt in ("-d", "--ideogram_dir"):
                self.ideogram_dir = arg
            elif opt in ("-m", "--gene_model_dir"):
                self.gene_model_dir = arg
            elif opt in ("-o", "--output_dir"):
                self.output_dir = arg
            elif opt in ("-1", "--five_p_trans_ID"):
                self.trans_ID_5 = arg
            elif opt in ("-2", "--five_p_trans_ID"):
                self.trans_ID_3 = arg
            elif opt in ("-p", "--transcript_frames"):
                self.TR_frames = arg
            elif opt in ("-b", "--bam-file"):
                self.bam_file = arg

        self.TR_frames = self.TR_frames.split(' ')


        transcript_lists = self.TR_frames[3].split(';')
        transcript_lists[:] = [TR_list.split('|') for TR_list in transcript_lists]

        self.in_frame_TRs_5 = [TR[:15] for TR in transcript_lists[0]]
        self.in_frame_TRs_3 = transcript_lists[1]
        self.out_frame_TRs_3 = transcript_lists[2]



        self.grep_gene(self.gene_ID_5, self.gene_ID_3, self.output_dir, "tmp.g.txt")

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_xlim([0, 2])
        self.ax.set_ylim([0, 2])
        plt.axis('on')
        self.get_gene_names(self.gene_ID_3, self.gene_ID_5, self.gene_model_dir_2)
        self.get_directionality(self.gene_ID_3, self.gene_ID_5, self.gene_model_dir_2)
        self.get_transcript(self.gene_ID_3, self.gene_ID_5)
        self.get_exons('5', self.TR_ID_5)
        self.get_exons('3', self.TR_ID_3)
        self.set_scales(self.ch_5, self.ch_3 ,self.ideogram_dir)
        self.get_chrom_band_code('5',self.ch_5, self.ideogram_dir, self.TR_start_5, self.TR_start_3)
        if self.ch_5!=self.ch_3:
            self.draw_chromosome('5')
            self.highlight_genes('5','5')
        self.get_chrom_band_code('3',self.ch_3, self.ideogram_dir, self.TR_start_5, self.TR_start_3)
        if self.ch_5!=self.ch_3:
            self.draw_chromosome('3')
            self.highlight_genes('3','3')
        if self.ch_5==self.ch_3:
            self.draw_chromosome('4')
            self.highlight_genes('4','5')
            self.highlight_genes('4','3')

        self.fusion_exons_3 = self.get_fusion_exons(self.gene_ID_3, self.strand_3, '3', self.fusion_point_3, self.exons_3)
        self.fusion_exons_5 = self.get_fusion_exons(self.gene_ID_5, self.strand_5, '5', self.fusion_point_5, self.exons_5)
        self.draw_fusion_exons('5', self.fusion_exons_5)
        self.draw_fusion_exons('3', self.fusion_exons_3)
        self.draw_fusion_trans('5', self.strand_5, self.fusion_center_exon_5, self.exons_5)
        self.draw_fusion_trans('3', self.strand_3, self.fusion_center_exon_3, self.exons_3)
        self.draw_combination_lines()
        self.plot_reads()
        self.draw(self.output_dir)
        self.rm_tmp(self.output_dir,"tmp.g.txt")

    def draw_one_sp1_read(self, percent, position): # anchor at 3'
        sp11_box = patches.Rectangle(
        (self.reads_x_middle-self.read_len*percent,self.reads_y_bottom+position*self.increase), self.read_len*percent, self.read_width, ec = "none", alpha = 1, fc = "blue"
        )
        self.ax.add_patch(sp11_box)
        sp12_box = patches.Rectangle(
        (self.reads_x_middle,self.reads_y_bottom+position*self.increase), self.read_len*(1.0-percent), self.read_width, ec = "none", alpha = 1, fc = "red"
        )
        self.ax.add_patch(sp12_box)
        sp1_an1_box = patches.Rectangle(
        (self.reads_x_middle+self.read_len*(1.0-percent),self.reads_y_bottom+self.read_width/2.0-self.read_width/5.0/2.0+position*self.increase), self.read_insert, self.read_width/5.0, ec = "none", alpha = 1, fc = "black"
        )
        self.ax.add_patch(sp1_an1_box)
        an1_box = patches.Rectangle(
        (self.reads_x_middle+self.read_len*(1.0-percent)+self.read_insert,self.reads_y_bottom+position*self.increase), self.read_len, self.read_width, ec = "none", alpha = 1, fc = "red"
        )
        self.ax.add_patch(an1_box)
    def draw_one_sp2_read(self, percent, position): # anchor at 5'
        sp22_box = patches.Rectangle(
        (self.reads_x_middle,self.reads_y_bottom+position*self.increase), self.read_len*percent, self.read_width, ec = "none", alpha = 1, fc = "red"
        )
        self.ax.add_patch(sp22_box)
        sp22_box = patches.Rectangle(
        (self.reads_x_middle-self.read_len*(1.0-percent),self.reads_y_bottom+position*self.increase), self.read_len*(1.0-percent), self.read_width, ec = "none", alpha = 1, fc = "blue"
        )
        self.ax.add_patch(sp22_box)
        sp2_an2_box = patches.Rectangle(
        (self.reads_x_middle-self.read_len*(1.0-percent)-self.read_insert,self.reads_y_bottom+self.read_width/2.0-self.read_width/5.0/2.0+position*self.increase), self.read_insert, self.read_width/5.0, ec = "none", alpha = 1, fc = "black"
        )
        self.ax.add_patch(sp2_an2_box)
        an2_box = patches.Rectangle(
        (self.reads_x_middle-self.read_len*(1.0-percent)-self.read_insert-self.read_len,self.reads_y_bottom+position*self.increase), self.read_len, self.read_width, ec = "none", alpha = 1, fc = "blue"
        )
        self.ax.add_patch(an2_box)

    def draw_one_en_read(self, percent, position): # encompassing
        mid_box = patches.Rectangle(
        (self.reads_x_middle-self.read_insert*percent,self.reads_y_bottom+self.read_width/2.0-self.read_width/5.0/2.0+position*self.increase), self.read_insert, self.read_width/5.0, ec = "none", alpha = 1, fc = "black"
        )
        self.ax.add_patch(mid_box)
        left_box = patches.Rectangle(
        (self.reads_x_middle-self.read_insert*percent-self.read_len,self.reads_y_bottom+position*self.increase), self.read_len, self.read_width, ec = "none", alpha = 1, fc = "blue"
        )
        self.ax.add_patch(left_box)
        right_box = patches.Rectangle(
        (self.reads_x_middle+self.read_insert*(1-percent),self.reads_y_bottom+position*self.increase), self.read_len, self.read_width, ec = "none", alpha = 1, fc = "red"
        )
        self.ax.add_patch(right_box)

    def plot_reads(self):
        if self.bam_file=='':
            return
        en,sp2,sp1 = rng.get_number(self.bam_file,self.ch_5,self.strand_5,self.TR_start_5,self.TR_end_5,self.fusion_point_5,self.ch_3,self.strand_3,self.TR_start_3,self.TR_end_3,self.fusion_point_3)

        #self.draw_one_sp1_read(0.25, 0)
        #self.draw_one_sp1_read(0.5,  1)
        #self.draw_one_sp1_read(0.75, 2)
        #self.draw_one_sp2_read(0.75, 3.5)
        #self.draw_one_sp2_read(0.5,  4.5)
        #self.draw_one_sp2_read(0.25, 5.5)
        #self.draw_one_en_read( 0.25, 7.5)
        #self.draw_one_en_read( 0.5,  8.5)
        #self.draw_one_en_read( 0.75, 9.5)

        pos_en = 0
        pos_sp2 = 0
        pos_sp1 = 0

        if en>0:
            self.draw_one_en_read( 0.25, 7.5)
            pos_en = 7.5
        if en>1:
            self.draw_one_en_read( 0.5,  8.5)
            pos_en = 8
        if en>=10:
            self.draw_one_en_read( 0.75, 9.5)
            pos_en = 8.5

        if sp2>0:
            self.draw_one_sp2_read(0.75, 3.5)
            pos_sp2 = 3.5
        if sp2>1:
            self.draw_one_sp2_read(0.5,  4.5)
            pos_sp2 = 4
        if sp2>=10:
            self.draw_one_sp2_read(0.25, 5.5)
            pos_sp2 = 4.5

        if sp1>0:
            self.draw_one_sp1_read(0.75, 2)
            pos_sp1 = 2
        if sp1>1:
            self.draw_one_sp1_read(0.5,  1)
            pos_sp1 = 1.5
        if sp1>=10:
            self.draw_one_sp1_read(0.25, 0)
            pos_sp1 = 1

        if en>0:
            plt.text(self.reads_x_middle+0.35, self.reads_y_bottom + pos_en*self.increase, "Encompassing reads "+str(en), fontsize = 8, ha = 'left')
        if sp2>0:
            plt.text(self.reads_x_middle+0.35, self.reads_y_bottom + pos_sp2*self.increase, "Spanning reads "+str(sp2)+" (anchor at 5')", fontsize = 8, ha = 'left')
        if sp1>0:
            plt.text(self.reads_x_middle+0.35, self.reads_y_bottom + pos_sp1*self.increase, "Spanning reads "+str(sp1)+" (anchor at 3')", fontsize = 8, ha = 'left')

    def get_gene_names(self, gene_ID_3, gene_ID_5, gene_model_dir):
        # description = subprocess.check_output("grep %s %s |  awk '$3 == \"gene\"' | cut -f9" % (gene_ID_3, gene_model_dir), shell = True).rstrip('\n').split('; ')
        # description[:] = [x.split(' ') for x in description]
        # self.gene_name_3 = next(x for x in description if x[0] == 'gene_name')[1][1:-1]
        # description = subprocess.check_output("grep %s %s |  awk '$3 == \"gene\"' | cut -f9" % (gene_ID_5, gene_model_dir), shell = True).rstrip('\n').split('; ')
        # description[:] = [x.split(' ') for x in description]
        # self.gene_name_5 = next(x for x in description if x[0] == 'gene_name')[1][1:-1]
        description = subprocess.check_output("grep %s %s |  awk '$3 == \"gene\"' | cut -f9" % (gene_ID_3, gene_model_dir), shell = True).rstrip('\n').rstrip(';').replace('"', '').split('; ')
        description = dict([x.split(' ') for x in description])
        self.gene_name_3 = description['gene_name']
        description = subprocess.check_output("grep %s %s |  awk '$3 == \"gene\"' | cut -f9" % (gene_ID_5, gene_model_dir), shell = True).rstrip('\n').rstrip(';').replace('"', '').split('; ')
        description = dict([x.split(' ') for x in description])
        self.gene_name_5 = description['gene_name']
        #print("got gene names: " + self.gene_name_3 + " " + self.gene_name_5)

    def get_directionality(self, gene_ID_3, gene_ID_5, gene_model_dir):
        self.strand_3 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f7" % (gene_ID_3, gene_model_dir), shell = True).rstrip('\n')[0]
        self.strand_5 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f7" % (gene_ID_5, gene_model_dir), shell = True).rstrip('\n')[0]
        #print("strands detected: " + self.strand_3 + " " +  self.strand_5)

    def get_transcript(self, gene_ID_3, gene_ID_5):
        if (gene_ID_5[:4] == 'ENSG'):
            self.TR_ID_5 = self.find_best_transcript(gene_ID_5, self.strand_5, '5')
        else:
            self.TR_ID_5 = gene_ID_5
        if (gene_ID_3[:4] == 'ENSG'):
            self.TR_ID_3 = self.find_best_transcript(gene_ID_3, self.strand_3, '3')
        else:
            self.TR_ID_3 = gene_ID_3
        #print("transcript ID: " + self.TR_ID_5)
        #print("transcript ID: " + self.TR_ID_3)

    def find_best_transcript(self, gene_ID, strand, side):
        description = subprocess.check_output("grep %s %s | awk '$3 == \"transcript\"' | cut -f4,5,9" % (gene_ID, self.gene_model_dir_2), shell = True).rstrip('\n').split('\n')
        description[:] = [transcript.split('\t') for transcript in description]
        description = [[transcript[0], transcript[1], dict([info.replace('"', '').split(' ') for info in transcript[2].rstrip(';').split('; ')])] for transcript in description]
        # description[:] = [re.split('; |\t', x) for x in description]
        # for transcript in description:
        #     transcript[:] = [x.split(' ') for x in transcript]

        TR_ID_list = [transcript[2]['transcript_id'] for transcript in description]

        min_dist_from_fusion = []
        length = []
        #to pass GRCh37 test only need change!!!!!1
        #biotype = [transcript[13][1][1:-1] for transcript in description]
        # biotype = ["protein_coding" for transcript in description]
        biotype = [transcript[2]['gene_biotype'] for transcript in description]
        # tsl = [int(transcript[-1][1][1:-2]) if isinstance(transcript[-1][1][1:-2], (int, long)) else 6  for transcript in description]
        tsl = []
        for transcript in description:
            if ('transcript_support_level' in transcript[2]):
                if isinstance(transcript[2]['transcript_support_level'], (int, long)):
                    tsl.append(int(transcript[2]['transcript_support_level']))
                else:
                    tsl.append(6)
            else:
                tsl.append(6)

        for ID in TR_ID_list:
            exon_junctions = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f4,5" % (ID, self.gene_model_dir_2), shell = True).rstrip('\n').split('\n')
            exon_junctions[:] = [x.split('\t') for x in exon_junctions]
            TR_length = sum([(int(junction_pair[1]) - int(junction_pair[0])) for junction_pair in exon_junctions])
            length.append(TR_length)
            if ((side == '5') and (strand == '+')):
                distances = [abs(int(junction_pair[1]) - int(self.fusion_point_5)) for junction_pair in exon_junctions]
            elif ((side == '3') and (strand == '-')):
                distances = [abs(int(junction_pair[1]) - int(self.fusion_point_3)) for junction_pair in exon_junctions]
            elif ((side == '5') and (strand == '-')):
                distances = [abs(int(junction_pair[0]) - int(self.fusion_point_5)) for junction_pair in exon_junctions]
            else:
                distances = [abs(int(junction_pair[0]) - int(self.fusion_point_3)) for junction_pair in exon_junctions]
            min_dist_from_fusion.append(min(distances))

        in_frame = [int(ID in self.in_frame_TRs_5) for ID in TR_ID_list] if side == '5' else [int(ID in self.in_frame_TRs_3) for ID in TR_ID_list]

        TR_list = zip(TR_ID_list, biotype, in_frame, min_dist_from_fusion, tsl, length)
        #TR_list.sort(key = lambda x: ((x[1] != 'trans_coding'), x[2], x[3], -x[4]))
        TR_list.sort(key = lambda x: ((x[1] != 'trans_coding'), -x[2], x[3], x[4], -x[5]))
        #print("Transcript metrics: ")
        #print(TR_list)
        return TR_list[0][0]

    def get_exons(self, side, TR_ID):
        gene_name = self.gene_name_5 if side == '5' else self.gene_name_3
        # exon_list = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f1,4,5,9" % (TR_ID, self.gene_model_dir_2), shell = True).rstrip("\n").split("\n")
        # exon_list[:] = [re.split('; |\t', exon) for exon in exon_list]
        # for exon in exon_list: exon[:] = [x.split(' ') for x in exon]
        # exon_list[:] = [(exon[0][0], exon[1][0], exon[2][0], exon[7][1][1:-1]) for exon in exon_list]

        # chr number, start, end, exon number
        exon_list = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f1,4,5,9" % (TR_ID, self.gene_model_dir_2), shell = True).rstrip("\n").split("\n")
        exon_list[:] = [exon.split('\t') for exon in exon_list]
        exon_list = [[exon[0], exon[1], exon[2], dict([info.replace('"', '').split(' ') for info in exon[3].rstrip(';').split('; ')])] for exon in exon_list]
        exon_list[:] = [(exon[0], int(exon[1]), int(exon[2]), exon[3]['exon_number']) for exon in exon_list]
        exon_list.sort(key = lambda x: x[1])

        if side == '5':
            self.exons_5 = exon_list
            self.ch_5 = exon_list[0][0] # not cast to int
            self.TR_start_5 = int(exon_list[0][1])
            self.TR_end_5 = int(exon_list[-1][2])
        else:
            self.exons_3 = exon_list
            self.ch_3 = exon_list[0][0] # not cast to int
            self.TR_start_3 = int(exon_list[0][1])
            self.TR_end_3 = int(exon_list[-1][2])

        #print("detected exons for " + gene_name + ":")
        #print(exon_list)

    def set_scales(self, chr5, chr3, ideogram_dir):
        if chr5==chr3:
            global scale_5
            scale_5 = 1.0*scale_global
            global scale_3
            scale_3 = 1.0*scale_global
            return
        ch_regions_5 = subprocess.check_output("grep -w chr%s %s" % (chr5, ideogram_dir), shell = True).rstrip('\n').split('\n')#add a chr to the name, becasue the code only uses GRCh38
        ch_regions_5[:] = [line.split('\t') for line in ch_regions_5]
        ch_length_5 = int(ch_regions_5[-1][2])

        ch_regions_3 = subprocess.check_output("grep -w chr%s %s" % (chr3, ideogram_dir), shell = True).rstrip('\n').split('\n')#add a chr to the name, becasue the code only uses GRCh38
        ch_regions_3[:] = [line.split('\t') for line in ch_regions_3]
        ch_length_3 = int(ch_regions_3[-1][2])

        if ch_length_5 > ch_length_3:
            scale_3 = 1.0*ch_length_5/ch_length_3*scale_global
            scale_5 = 1.0*scale_global
        else:
            scale_5 = 1.0*ch_length_3/ch_length_5*scale_global
            scale_3 = 1.0*scale_global

    def get_chrom_band_code(self, side, chromosome, ideogram_dir, TR_start_5, TR_start_3):
        self.ch_regions = subprocess.check_output("grep -w chr%s %s" % (chromosome, ideogram_dir), shell = True).rstrip('\n').split('\n')#add a chr to the name, becasue the code only uses GRCh38

        self.ch_regions[:] = [line.split('\t') for line in self.ch_regions]

        if side =='3':
            self.gene_region_3 = next(region for region in self.ch_regions if (int(TR_start_3) > int(region[1]) and int(TR_start_3) <  int(region[2])))
            self.gene_band_code_3 = str(self.ch_3)+str(self.gene_region_3[3])
        else:
            self.gene_region_5 = next(region for region in self.ch_regions if (int(TR_start_5) > int(region[1]) and int(TR_start_5) <  int(region[2])))
            self.gene_band_code_5 = str(self.ch_5)+str(self.gene_region_5[3])
        self.ch_length = int(self.ch_regions[-1][2])

    def draw_chromosome(self, side):
        tmp=self.ch_graph_start_coord
        if side=='3':
            self.ch_graph_start_coord=self.ch_graph_start_coord+1.0
        if side=='4':
            self.ch_graph_start_coord=self.ch_graph_start_coord+0.5

        #changed add /scale make it shorter
        scale = 1.0
        if side =='5':
            scale = scale_5
        else:
            scale = scale_3

        acen_start = int(next(region for region in self.ch_regions if region[4] == 'acen')[1])
        acen_end = int(next(region for region in list(reversed(self.ch_regions)) if region[4] == 'acen')[2])

        self.arm_length_left = float(acen_start)/self.ch_length/scale
        self.arm_length_right = (self.ch_length - float(acen_end))/self.ch_length/scale

        #changed
        if side =='5':
            self.acen_length_5 = acen_end - acen_start
        else:
            self.acen_length_3 = acen_end - acen_start

        self.centrolmere_x = self.ch_graph_start_coord + self.arm_length_left + self.centrolmere_radius + self.pad_size

        centrolmere = patches.Circle((self.centrolmere_x, self.centrolmere_y), self.centrolmere_radius, ec="none", color = "black")

        chrom_arm_1 = patches.FancyBboxPatch(
            (self.centrolmere_x + 0.05, self.centrolmere_y - self.arm_width/2), self.arm_length_right, self.arm_width,
            facecolor= "gray", ec = "black",
            boxstyle=patches.BoxStyle("Round", pad = self.pad_size))

        chrom_arm_2 = patches.FancyBboxPatch(
            (self.centrolmere_x - 0.05 - self.arm_length_left, self.centrolmere_y - self.arm_width/2), self.arm_length_left, self.arm_width,
            facecolor= "gray", ec = "black",
            boxstyle=patches.BoxStyle("Round", pad = self.pad_size))

        self.ax.add_patch(chrom_arm_1)
        self.ax.add_patch(chrom_arm_2)
        self.ax.add_patch(centrolmere)
        self.ch_graph_start_coord=tmp

    def highlight_genes(self, chrcode, side):
        tmp=self.ch_graph_start_coord
        if chrcode=='3':
            self.ch_graph_start_coord=self.ch_graph_start_coord+1.0
        if chrcode=='4':
            self.ch_graph_start_coord=self.ch_graph_start_coord+0.5

        global lines
        TR_start = self.TR_start_5 if side == '5' else self.TR_start_3
        gene_band_code = self.gene_band_code_5 if side == '5' else self.gene_band_code_3
        label_y = 0.5
        highlight_box_width = 0.02
        label_x = 0.5 if side == '5' else 1.5
        label_radius = 0.45
        color = "blue" if side == '5' else "red"
        line_width = 1.5
        ha_side = 'center'
        if chrcode=='4' and side =='5':
            if self.TR_start_5 < self.TR_start_3:
                ha_side = 'right'
            else:
                ha_side = 'left'
        if chrcode=='4' and side =='3':
            if self.TR_start_3 < self.TR_start_5:
                ha_side = 'right'
            else:
                ha_side = 'left'

        TR_ID = self.TR_ID_5 if side == '5' else self.TR_ID_3
        number_of_exons = len(self.exons_5) if side == '5' else len(self.exons_3)
        gene_name = self.gene_name_5 if side == '5' else self.gene_name_3

        #changed added /scale
        scale = 1.0
        if side =='5':
            scale = scale_5
        else:
            scale = scale_3

        highlight_box_offset = float(TR_start)/self.ch_length/scale
        #changed + self.pad_size * 2

        acen_length = self.acen_length_5 if side == '5' else self.acen_length_3

        #changed, need to subtract acen region bp length
        if highlight_box_offset > self.arm_length_left:
            highlight_box_offset = highlight_box_offset + self.centrolmere_radius * 2 + self.pad_size * 2 - float(acen_length)/self.ch_length/scale

        highlight_box = patches.Rectangle(
        (self.ch_graph_start_coord + highlight_box_offset, (self.centrolmere_y - self.arm_width/2) - self.pad_size * 2), highlight_box_width, self.arm_width + 4 * self.pad_size, ec = color, alpha = 1, fill = False, lw = line_width
        )

        line_1 = lines.Line2D(#lines connect chromosomes and transcript names
        [(self.ch_graph_start_coord + highlight_box_offset + highlight_box_width/2), label_x],
        [((self.centrolmere_y - self.arm_width/2) - self.pad_size * 2) + (self.arm_width + 4 * self.pad_size), label_y],
        c = color, alpha = 1, lw = line_width*0.5
        )

        line_2 = lines.Line2D(#horizontal line under transcript names
        [label_x - label_radius, label_x + label_radius],
        [label_y, label_y],
        c = color, alpha = 1, lw = line_width*0.5
        )

        line_3 = lines.Line2D(#left vertial for transcript names
        [label_x - label_radius, label_x - label_radius],
        [label_y, label_y + 0.02],
        c = color, alpha = 1, lw = line_width*0.5
        )

        line_4 = lines.Line2D(#right vertial for transcript names
        [label_x + label_radius, label_x + label_radius],
        [label_y, label_y + 0.02],
        c = color, alpha = 1, lw = line_width*0.5
        )

        add_space = 0.0
        if ha_side == 'left':
            add_space = 0.02
        plt.text(add_space+self.ch_graph_start_coord + highlight_box_offset, (self.centrolmere_y - self.arm_width/2) - self.pad_size * 2 - 0.05, gene_band_code, fontsize = 9, ha = ha_side) # Transcript name
        plt.text(label_x, label_y + 0.01, "(%s, %s exons)" % (TR_ID, number_of_exons), fontsize = 9, ha = 'center')
        plt.text(label_x, label_y + 0.07, gene_name, fontsize = 11, ha = 'center') #gene name

        highlight_lines = [line_1, line_2, line_3, line_4]
        for line in highlight_lines:
            self.ax.add_line(line)
        self.ax.add_patch(highlight_box)

        self.ch_graph_start_coord=tmp

    def get_fusion_exons(self, gene_ID, strand, side, fusion_point, exon_list):
        fusion_exon_list = []
        junction_side = 2 if ((side == '5') and (strand == '+')) or ((side == '3') and (strand == '-')) else 1
        distances = [abs(int(fusion_point) - int(exon_junctions[junction_side])) for exon_junctions in exon_list]
        i = distances.index(min(distances))
        if side == '5':
            self.fusion_center_exon_5 = i
        else:
            self.fusion_center_exon_3 = i
        ex_size=len(exon_list)
        if i-1>=0:#assume sorted
            fusion_exon_list.append(exon_list[i-1])
        fusion_exon_list.append(exon_list[i])
        if i+1<ex_size:
            fusion_exon_list.append(exon_list[i+1])
        return fusion_exon_list

    def draw_fusion_exons(self, side, fusion_exons):

        DNA_strand_5 = lines.Line2D(
        [self.exon_strand_x_5 - self.exon_strand_radius, self.exon_strand_x_5 + self.exon_strand_radius],
        [self.exon_strands_y, self.exon_strands_y],
        c = "blue", alpha = 1, lw = self.line_width
        )

        DNA_strand_3 = lines.Line2D(
        [self.exon_strand_x_3 - self.exon_strand_radius, self.exon_strand_x_3 + self.exon_strand_radius],
        [self.exon_strands_y, self.exon_strands_y],
        c = "red", alpha = 1, lw = self.line_width
        )
        self.ax.add_line(DNA_strand_5)
        self.ax.add_line(DNA_strand_3)

        #calculate for the position for the boxes in draw_exon_boxes
        if len(fusion_exons) == 3:
            show = [True, False, True, False, True]
            exon_numbers = [fusion_exons[0][3], 0, fusion_exons[1][3], 0, fusion_exons[2][3]]
        if len(fusion_exons) == 2:
            #if (side == '5' and self.strand_5=='-') or (side == '3' and self.strand_3=='-'):
            #    show = [True, False, True, False, False]
            #    exon_numbers = [fusion_exons[0][3], 0, fusion_exons[1][3], 0, 0]
            #if (side == '5' and self.strand_5=='+') or (side == '3' and self.strand_3=='+'):
            #    show = [False, False, True, False, True]
            #    exon_numbers = [0, 0, fusion_exons[0][3], 0, fusion_exons[1][3]]
            if (side == '5' and self.strand_5=='-') or (side == '3' and self.strand_3=='+'):
                show = [True, False, True, False, False]
                exon_numbers = [fusion_exons[0][3], 0, fusion_exons[1][3], 0, 0]
            if (side == '5' and self.strand_5=='+') or (side == '3' and self.strand_3=='-'):
                show = [False, False, True, False, True]
                exon_numbers = [0, 0, fusion_exons[0][3], 0, fusion_exons[1][3]]
               
        if len(fusion_exons) == 1:
            show = [False, False, True, False, False]
            exon_numbers = [0, 0, fusion_exons[0][3], 0, 0]

        if side == '5':
            box_color = "blue"
            center_x_coord = 0.5
            direction = self.strand_5
            fusion_coord = self.fusion_point_5
            strand = self.strand_5

        else:
            box_color = "red"
            center_x_coord = 1.5
            direction = self.strand_3
            fusion_coord = self.fusion_point_3
            strand = self.strand_3

        self.draw_exon_boxes(box_color, center_x_coord, show)
        self.label_exon_boxes(center_x_coord, exon_numbers, show)
        self.label_strand_direction(box_color, center_x_coord, direction, show)
        self.highlight_fusion_junction(side, strand, center_x_coord, direction, fusion_coord, fusion_exons)

    def draw_exon_boxes(self, box_color, center_x_coord, show):
        exon_boxes = []
        for i in range(-2, 3):
            x_offset = i * self.exon_seperation/2
            exon_box = patches.Rectangle(
            (center_x_coord - self.exon_box_length/2.0 + x_offset, self.exon_strands_y - self.exon_box_width/2.0), self.exon_box_length, self.exon_box_width, ec = "none", alpha = 1, fc = box_color
            )
            exon_boxes.append(exon_box)

        for box, show_status in zip(exon_boxes, show):
            if show_status:
                self.ax.add_patch(box)

    def label_exon_boxes(self, center_x_coord, exon_numbers, show):
        for number, show_status, i in zip(exon_numbers, show, range(-2, 3)):
            if show_status:
                x_offset = i * self.exon_seperation/2
                plt.text(center_x_coord + x_offset, self.exon_strands_y, number, fontsize = 9, ha = 'center', va = 'center', color = "white")

    def label_strand_direction(self, box_color, center_x_coord, direction, show):
        show2 = []
        if show[1]:
            show2.append(True)
        else:
            show2.append(False)
        for x in range(3):
            y=x+1
            if show[y-1] and show[y+1]:
                show2.append(True)
            else:
                show2.append(False)
        if show[4]:
            show2.append(True)
        else:
            show2.append(False)

        show = show2
        #show = [not i for i in show]
        h = 0.01
        w = 0.01 if direction == '-' else -0.01

        for show_status, i in zip(show, range(-2, 3)):
            if show_status:
                x_offset = i * self.exon_seperation/2
                arrow = lines.Line2D(
                [center_x_coord + w + x_offset, center_x_coord + x_offset, center_x_coord + w + x_offset],
                [self.exon_strands_y + h, self.exon_strands_y, self.exon_strands_y - h],
                c = box_color, alpha = 1, lw = self.line_width
                )
                self.ax.add_line(arrow)

    def highlight_fusion_junction(self, side, strand, center_x_coord, direction, fusion_coord, fusion_exons):
        color = "black"
        xx = 0
        if side == '5':
            adjust = self.exon_box_length/2 if self.strand_5=="+" else 0.0-self.exon_box_length/2
            xx=self.exon_strand_x_5 + adjust
        else:
            adjust = 0.0-self.exon_box_length/2 if self.strand_3=="+" else self.exon_box_length/2
            xx=self.exon_strand_x_3 + adjust
        yy=self.exon_strands_y - self.exon_box_width*1.5

        #junction_side = 2 if ((side == '5') and (strand == '+')) or ((side == '3') and (strand == '-')) else 1
        ch_number = self.ch_5 if side == '5' else self.ch_3
        junction_mark_width = 0.01
        junction_mark_length = self.exon_box_width

        #diff = int(fusion_coord) - int(fusion_exons[1][junction_side]) # need to accomodate 2 exon/1 exon case
        #if diff < 0:
        #    offset = self.exon_seperation * diff/(int(fusion_exons[1][junction_side]) - int(fusion_exons[0][junction_side]))
        #else:
        #    offset = self.exon_seperation * diff/(int(fusion_exons[2][junction_side]) - int(fusion_exons[1][junction_side]))

        #if junction_side == 2: # left
        #    fusion_x = (center_x_coord + self.exon_box_length/2 + offset)
        #else:
        #    fusion_x = (center_x_coord - self.exon_box_length/2 + offset)

        #if side == '5' and self.strand_5=="+":
        #    fusion_x=fusion_x + self.exon_seperation
        #if side == '3' and self.strand_3=="+": # this part need test
        #    fusion_x=fusion_x + self.exon_seperation

        #fusion_mark = patches.Rectangle(
        #(fusion_x, self.exon_strands_y - self.exon_box_width*1.5), junction_mark_width, junction_mark_length, ec = "none", alpha = 1, fc = 'black'
        #)

        #plt.text(fusion_x - 0.01, self.exon_strands_y - self.exon_box_width * 1.5, "chr%s:%s" % (ch_number, fusion_coord), fontsize = 6, ha = 'right', va = 'bottom', color = "black")

        fusion_mark = patches.Rectangle(
        (xx, yy), junction_mark_width, junction_mark_length, ec = "none", alpha = 1, fc = 'black'
        )
        plt.text(xx - 0.01, yy, "chr%s:%s" % (ch_number, fusion_coord), fontsize = 6, ha = 'right', va = 'bottom', color = "black")
        self.ax.add_patch(fusion_mark)

    def draw_fusion_trans(self, side, strand, fusion_center_exon, all_exons):
        fusion_trans_exon_numbers = []
        if ((side == '5') and (strand == '+')) or ((side == '3') and (strand == '-')):
            fusion_trans_exon_numbers[:] = [included_exon[3] for included_exon in all_exons[:(fusion_center_exon + 1)]]
        else:
            fusion_trans_exon_numbers[:] = [included_exon[3] for included_exon in all_exons[fusion_center_exon:]]

        if (strand == '-'):
            fusion_trans_exon_numbers.reverse()

        #print("exons included in the fusion trans: ")
        #print(fusion_trans_exon_numbers)

        if len(fusion_trans_exon_numbers) > 3:
            show = [1, 0, 1, 1, 1, 1, 1] if side == '5' else [1, 1, 0, 1, 1, 1, 1]
            numbers = [fusion_trans_exon_numbers[0], 0, fusion_trans_exon_numbers[-2], fusion_trans_exon_numbers[-1]] if side == '5' else [fusion_trans_exon_numbers[0], fusion_trans_exon_numbers[1], 0, fusion_trans_exon_numbers[-1]]
        elif len(fusion_trans_exon_numbers) == 3:
            #show = [0, 1, 1, 1, 0, 0, 0]
            show = [0, 1, 1, 1, 0, 0, 0] if side == '5' else [1, 1, 1, 0, 0, 0, 0]
            #numbers = [0, fusion_trans_exon_numbers[0], fusion_trans_exon_numbers[1], fusion_trans_exon_numbers[2]]
            numbers = [0, fusion_trans_exon_numbers[0], fusion_trans_exon_numbers[1], fusion_trans_exon_numbers[2]] if side =='5' else [fusion_trans_exon_numbers[0], fusion_trans_exon_numbers[1], fusion_trans_exon_numbers[2]]
        elif len(fusion_trans_exon_numbers) == 2:
            if side == '5':
                show = [0, 0, 1, 1, 0, 0, 0]
                numbers = [0, 0, fusion_trans_exon_numbers[0], fusion_trans_exon_numbers[1]]
            else:
                show = [1, 1, 0, 0, 0, 0, 0]
                numbers = [fusion_trans_exon_numbers[0], fusion_trans_exon_numbers[1], 0, 0, 0]
        else:
            if side  == '5':
                show = [0, 0, 0, 1, 0, 0, 0]
                numbers = [0, 0, 0, fusion_trans_exon_numbers[0]]
            else:
                show = [1, 0, 0, 0, 0, 0, 0]
                numbers = [fusion_trans_exon_numbers[0], 0, 0, 0, 0]
        self.draw_fusion_trans_boxes(side, show, numbers)

    def draw_fusion_trans_boxes(self, side, show, numbers):
        box_color = "blue" if side == '5' else "red"
        offset_multiplier = -1.05 if side == '5' else 1.05
        box_positions = [3.5, 3, 2, 1] if side == '5' else [0, 1, 2, 2.5]
        dot_positions = [2.125, 2.25, 2.375] if side == '5' else [2.125, 2.25, 2.375]
        fusion_trans_x = self.fusion_trans_x
        fusion_trans_y = self.fusion_trans_y
        graphics = []
        dot_radius = 0.005

        for i in box_positions:
            x_offset = i * offset_multiplier * self.exon_box_length
            exon_box = patches.Rectangle(
            (fusion_trans_x + x_offset, fusion_trans_y), self.exon_box_length, self.exon_box_width, ec = "none", alpha = 1, fc = box_color
            )
            graphics.append(exon_box)

        for i in dot_positions: #dots
            x_offset = i * offset_multiplier * self.exon_box_length
            dot = patches.Circle((fusion_trans_x + x_offset, fusion_trans_y), dot_radius, ec="none", color = "black")
            graphics.append(dot)

        for graphic, show_status in zip(graphics, show): #box on top: fusion transcript
            if show_status:
                self.ax.add_patch(graphic)

        for number, i, show_status in zip(numbers, box_positions, show[:4]): # exon name
            if show_status:
                x_offset = i * offset_multiplier * self.exon_box_length
                plt.text(fusion_trans_x + x_offset + self.exon_box_length/2, fusion_trans_y + self.exon_box_width/2, number, fontsize = 8, ha = 'center', va = 'center', color = "white")


    def draw_combination_lines(self): # thin black lines connects junctions
        color = "black"
        adjust_5 = self.exon_box_length/2 if self.strand_5=="+" else 0.0-self.exon_box_length/2
        line_5 = lines.Line2D(
        [self.exon_strand_x_5 + adjust_5, self.fusion_trans_x-0.01],
        [self.exon_strands_y + self.exon_box_width/2, self.fusion_trans_y],
        c = color, alpha = 1, lw = 0.5
        )
        adjust_3 = 0.0-self.exon_box_length/2 if self.strand_3=="+" else self.exon_box_length/2
        line_3 = lines.Line2D(
        [self.exon_strand_x_3 + adjust_3, self.fusion_trans_x],
        [self.exon_strands_y + self.exon_box_width/2, self.fusion_trans_y],
        c = color, alpha = 1, lw = 0.5
        )
        self.ax.add_line(line_5)
        self.ax.add_line(line_3)

    def draw(self, output_dir):
        self.fig.tight_layout()
        plt.axis('off')
        self.fig.savefig(self.output_dir,bbox_inches='tight')

def main(argv):
    t=time.time()
    panel = PANEL_A(argv)
    print "Time used: %0.2f" % (time.time()-t), "seconds."

if __name__ == '__main__':
    sys.exit(main(sys.argv))

