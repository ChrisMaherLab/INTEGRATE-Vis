import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import sys, subprocess
import getopt
from textwrap import wrap
import re
import time
import os
import pdb

# example: python pb.py -5 ENSG00000184012 -3 ENSG00000157554 -f 41498118 -t 38445621 -s 0 -p 'RHSNASQSLCEIVRLSRDQMLQIQNSTEPDPLLATLEKX 37 2 ENST00000332149(731)|ENST00000454499(731);ENST00000288319;ENST00000294304' -d domain_table.txt -m /gscmnt/gc2601/maherlab/jzhang/PRAD2/ensembl/Homo_sapiens.GRCh38.85.gtf -o plots/pb.pdf

def usage():
    print """
    panel_b -5 <five_p_gene_ID_or_transcript_ID> -3 <three_p_gene_ID_or_transcript_ID> -f <five_prime_gene_fusion_point> -t <three_prime_gene_fusion_point> -s <possible_in_frame> -p <transcript_frames> -d <domain_table_dir> -m <gene_model_dir> -o <output_dir>
    Parameters:

        -5/--five_p_gene_ID      [Requested]
                                 [string:   5' gene ID.                                                             ]
        -3/--three_p_gene_ID     [Requested]
                                 [string:   3' gene ID.                                                             ]
        -f/--five_prime_gene_fusion_point   [Requested]
                                 [string:   5' gene fusion coordinate                                               ]
        -t/--three_prime_gene_fusion_point   [Requested]
                                 [string:   3' gene fusion coordinate                                               ]
        -s/--possible_in_frame   [Requested]
                                 [string:   potential in frame shift or not                                         ]
        -p/--transcript_frames   [Requested]
                                 [string:   stop info and shift frame status for different transcripts,
                                            a sub line of the annotated bedpe, seperated by space, comma, and '|'   ]
        -d/--domain_table_dir    [Requested]
                                 [string:   domain information table director                                       ]
        -m/--gene_model_dir      [Requested]
                                 [string:   provide your own gene model or default is used                          ]
        -o/--output-dir          [Optional]
                                 [string:   output directory.                                                       ]
        -1/--five_p_trans_ID     [Optional]
                                 [string:   5' transcript ID. Only choose this transcrtipt.                         ]
        -2/--three_p_trans_ID    [Optional]
                                 [string:   3' transcript ID. Only choose this transcrtipt.                         ]

    Version:                     1.0.0
          """

def rm_quoat(aa):
    if aa[0]=='\"' and aa[len(aa)-1]==";":
        return aa[1:len(aa)-2]
    else:
        return aa

class PANEL_B:

    #requested
    gene_ID_5 = ''
    gene_ID_3 = ''
    fusion_point_5 = ''
    fusion_point_3 = ''
    domain_table_dir = ''
    gene_model_dir =''
    gene_model_dir_2 =''
    potential_in_frame = 1
    stop_info = ''

    #optional
    output_dir = './pb.pdf'
    trans_ID_5 = 'None'
    trans_ID_3 = 'None'

    #computed
    TR_frames = ''
    in_frame_TRs_5 = []
    in_frame_TRs_3 = []
    out_frame_TRs_3 = []
    gene_name_5 = ''
    gene_name_3 = ''
    display_ratio = 0 # in bp
    TR_start_5 = 0
    TR_ID_3 = ''
    TR_ID_5 = ''
    TR_length_3 = 0
    TR_length_5 = 0
    fusion_TR_length = 0
    fusion_protein_length = 0
    exons_3 = []
    exons_5 = []
    fusion_center_exon_5 = 0
    fusion_center_exon_3 = 0
    fusion_exon_list_5 = []
    fusion_exon_list_3 = []
    strand_3 = ''
    strand_5 = ''
    TR_fragment_length_5 = 0
    TR_fragment_length_3 = 0
    five_p_utr_length_5 = 0
    three_p_utr_length_5 = 0
    five_p_utr_length_3 = 0
    three_p_utr_length_3 = 0
    domain_list_3 = []
    domain_list_5 = []
    included_domain_list_3 = []
    included_domain_list_5 = []
    fusion_offset = 0

    #fixed
    font_size = 10
    graphic_protein_length = 1.7
    graphic_protein_width = 0.1
    graphic_protein_y = 0.3
    graphic_protein_x = 0.13
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']
    domain_colors = []
    protein_line_width = 0.25

    #def __init__(self, gene_ID_5, gene_ID_3, fusion_point_5, fusion_point_3, potential_in_frame, TR_frames, domain_table_dir, gene_model_dir, output_dir, trans_ID_5, trans_ID_3):
    def __init__(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:],"h5:3:f:t:s:p:d:m:o:1:2:", ["help", "gene_ID_5=", "gene_ID_3=", "fusion_point_5=", "fusion_point_3=", "domain_table_dir=", "potential_in_frame=", "transcript_frames=","gene_model_dir=", "output_dir=", "five_p_trans_ID=", "three_p_trans_ID="])
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
            elif opt in ("-s", "--possible_in_frame"):
                self.potential_in_frame = int(arg)
            elif opt in ("-p", "--transcript_frames"):
                self.TR_frames = arg
            elif opt in ("-d", "--domain_table_dir"):
                self.domain_table_dir = arg
            elif opt in ("-m", "--gene_model_dir"):
                self.gene_model_dir = arg
            elif opt in ("-o", "--output_dir"):
                self.output_dir = arg
            elif opt in ("-1", "--five_p_trans_ID"):
                self.trans_ID_5 = arg
            elif opt in ("-2", "--five_p_trans_ID"):
                self.trans_ID_3 = arg


    #    self.gene_ID_5 = gene_ID_5
    #    self.gene_ID_3 = gene_ID_3
    #    self.fusion_point_5 = fusion_point_5
    #    self.fusion_point_3 = fusion_point_3
    #    self.domain_table_dir = domain_table_dir
    #    self.gene_model_dir = gene_model_dir
    #    self.output_dir = output_dir
    #    self.potential_in_frame = potential_in_frame

    #    self.trans_ID_5=trans_ID_5
    #    self.trans_ID_3=trans_ID_3

        self.TR_frames = self.TR_frames.split(' ')
        self.stop_info = self.TR_frames[:3]

        transcript_lists = self.TR_frames[3].split(';')
        transcript_lists[:] = [TR_list.split('|') for TR_list in transcript_lists]

        self.in_frame_TRs_5 = [TR[:15] for TR in transcript_lists[0]]
        self.in_frame_TRs_3 = transcript_lists[1]
        self.out_frame_TRs_3 = transcript_lists[2]

        #print("in_frame_TRs_5")
        #print(self.in_frame_TRs_5)
        #print("in_frame_TRs_3")
        #print(self.in_frame_TRs_3)
        #print("out_frame_TRs_3")
        #print(self.out_frame_TRs_3)

    def grep_equal_t(self, tId, output_dir, suffix):
        name=output_dir + "." + suffix
        f.open(name,"r")
        f.write(tId_model[tId])
        f.close()

    def grep_equal_g(self, gId, output_dir, suffix):
        name=output_dir + "." + suffix
        f.open(name,"r")
        f.write(gId_model[tId])
        f.close()

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

    def get_gene_names(self, gene_ID_3, gene_ID_5, gene_model_dir):
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

        TR_ID_list = [transcript[2]['transcript_id'] for transcript in description]

        min_dist_from_fusion = []
        length = []

        biotype = [transcript[2]['gene_biotype'] for transcript in description]

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
        TR_list.sort(key = lambda x: ((x[1] != 'trans_coding'), -x[2], x[3], x[4], -x[5]))
        return TR_list[0][0]

    def get_exons(self, side, TR_ID):
        gene_name = self.gene_name_5 if side == '5' else self.gene_name_3
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
            self.TR_length_5 = sum([int(exon[2]) - int(exon[1]) for exon in exon_list])
        else:
            self.exons_3 = exon_list
            self.ch_3 = exon_list[0][0] # not cast to int
            self.TR_start_3 = int(exon_list[0][1])
            self.TR_end_3 = int(exon_list[-1][2])
            self.TR_length_3 = sum([int(exon[2]) - int(exon[1]) for exon in exon_list])

    def get_fusion_exons(self, strand, side, fusion_point, exon_list):
        fusion_exon_list = []
        junction_side = 2 if ((side == '5') and (strand == '+')) or ((side == '3') and (strand == '-')) else 1
        distances = [abs(int(fusion_point) - int(exon_junctions[junction_side])) for exon_junctions in exon_list]
        fusion_center_exon = distances.index(min(distances))

        if side == '5':
            self.fusion_center_exon = fusion_center_exon
        else:
            self.fusion_center_exon = fusion_center_exon

        if ((side == '5') and (strand == '+')) or ((side == '3') and (strand == '-')):
            fusion_exon_list[:] = [included_exon for included_exon in exon_list[:(fusion_center_exon + 1)]]
        else:
            fusion_exon_list[:] = [included_exon for included_exon in exon_list[fusion_center_exon:]]
        #print("exons included in the fusion protein: ")
        #print(fusion_exon_list)
        return fusion_exon_list

    def draw_fusion_protein(self):
        fusion_protein_box = patches.Rectangle(
        (self.graphic_protein_x, self.graphic_protein_y), self.graphic_protein_length, self.graphic_protein_width, ec = 'black', alpha = 1, fill = True, lw = self.protein_line_width, fc = '#b8b7bb'
        )
        self.ax.add_patch(fusion_protein_box)

    def compute_fusion_parameters(self):
        fusion_exon_lengths_5 = [int(exon[2]) - int(exon[1]) for exon in self.fusion_exon_list_5]
        fusion_exon_lengths_3 = [int(exon[2]) - int(exon[1]) for exon in self.fusion_exon_list_3]

        self.TR_fragment_length_5 = sum(fusion_exon_lengths_5)
        self.TR_fragment_length_3 = sum(fusion_exon_lengths_3)

        self.fusion_TR_length = self.TR_fragment_length_5 + self.TR_fragment_length_3
        self.fusion_protein_length = (self.TR_fragment_length_5 + self.TR_fragment_length_3)/3      # sometimes this isn't integer

        self.display_ratio = self.graphic_protein_length/self.fusion_TR_length

    def highlight_fusion_junction(self):
        self.fusion_offset = float(self.TR_fragment_length_5) * self.graphic_protein_length/self.fusion_TR_length

        #print("fusion_TR_length: " + str(self.fusion_TR_length))
        #print("self.fusion_offset: " + str(self.fusion_offset))

        fusion_mark_y = self.graphic_protein_y + self.graphic_protein_width + 0.01

        fusion_mark = lines.Line2D(
        [self.graphic_protein_x + self.fusion_offset - 0.05, self.graphic_protein_x + self.fusion_offset, self.graphic_protein_x + self.fusion_offset + 0.05],
        [fusion_mark_y + 0.05, fusion_mark_y, fusion_mark_y + 0.05],
        c = "black", alpha = 1, lw = 3
        )

        fusion_division = lines.Line2D(
        [self.graphic_protein_x + self.fusion_offset, self.graphic_protein_x + self.fusion_offset],
        [fusion_mark_y - self.graphic_protein_width - 0.01, fusion_mark_y],
        c = "black", alpha = 1, lw = self.protein_line_width * 4
        )
        self.ax.add_line(fusion_division)
        self.ax.add_line(fusion_mark)

    def get_utr_lengths(self, side, strand, TR_ID, exon_list):
        #CDS_list = subprocess.check_output("grep %s %s | awk '$3 == \"CDS\"' | cut -f4,5" % (TR_ID, self.gene_model_dir), shell = True).rstrip("\n").split("\n")
        CDS_list = subprocess.check_output("grep %s %s | awk '$3 == \"CDS\"' | cut -f4,5" % (TR_ID, self.gene_model_dir_2), shell = True).rstrip("\n").split("\n")
        CDS_list[:] = [CDS.split('\t') for CDS in CDS_list]
        if CDS_list[0][0] == '':
            if side == '5':
                # return self.TR_fragment_length_5, self.TR_fragment_length_5
                return self.TR_length_5, self.TR_length_5
            else:
                # return self.TR_fragment_length_3, self.TR_fragment_length_3
                return self.TR_length_3, self.TR_length_3
        CDS_list.sort(key = lambda x: x[0])
        #print(CDS_list)
        CDS_start = int(CDS_list[0][0])
        CDS_end = int(CDS_list[-1][1])

        #print("CDS for " + side + " prime gene")
        #print(CDS_list)

        five_p_utr_length = 0
        three_p_utr_length = 0

        for exon in exon_list:
            if exon[2] < CDS_start:
                five_p_utr_length += (exon[2] - exon[1])
            elif exon[1] < CDS_start and exon[2] > CDS_start:
                five_p_utr_length += (CDS_start - exon[1])
            elif exon[1] < CDS_end and exon[2] > CDS_end:
                three_p_utr_length += (exon[2] - CDS_end)
            elif exon[1] > CDS_end:
                three_p_utr_length += (exon[2] - exon[1])

        if strand == '-':
            five_p_utr_length, three_p_utr_length = three_p_utr_length, five_p_utr_length

        return five_p_utr_length, three_p_utr_length

    def draw_utrs(self):
        if self.five_p_utr_length_5 < self.TR_fragment_length_5:
            graphic_utr_length_5 = self.five_p_utr_length_5 * self.display_ratio
        else:
            graphic_utr_length_5 = self.TR_fragment_length_5 * self.display_ratio

        if self.three_p_utr_length_3 < self.TR_fragment_length_3:
            graphic_utr_length_3 = self.three_p_utr_length_3 * self.display_ratio
        else:
            graphic_utr_length_3 = self.TR_fragment_length_3 * self.display_ratio

        five_p_utr_box_5 = patches.Rectangle(
        (self.graphic_protein_x, self.graphic_protein_y), graphic_utr_length_5, self.graphic_protein_width, ec = 'black', fill = True, lw = self.protein_line_width, fc = '#808080', hatch = '////'
        )

        three_p_utr_box_3 = patches.Rectangle(
        (self.graphic_protein_x + self.graphic_protein_length - graphic_utr_length_3, self.graphic_protein_y), graphic_utr_length_3, self.graphic_protein_width, ec = 'black', fill = True, lw = self.protein_line_width, fc = '#808080', hatch = '////'
        )

        self.ax.add_patch(five_p_utr_box_5)
        self.ax.add_patch(three_p_utr_box_3)

        if self.five_p_utr_length_3 > (self.TR_length_3 - self.TR_fragment_length_3):
            five_p_utr_box_3 = patches.Rectangle(
            (self.graphic_protein_x + self.fusion_offset, self.graphic_protein_y), self.display_ratio * (self.five_p_utr_length_3 - (self.TR_length_3 - self.TR_fragment_length_3)), self.graphic_protein_width, ec = 'black', fill = True, lw = self.protein_line_width, fc = '#808080', hatch = '////')
            self.ax.add_patch(five_p_utr_box_3)

        if self.three_p_utr_length_5 > (self.TR_length_5 - self.TR_fragment_length_5):
            three_p_utr_length = self.display_ratio * (self.three_p_utr_length_5 - (self.TR_length_5 - self.TR_fragment_length_5))
            three_p_utr_box_5 = patches.Rectangle(
            (self.graphic_protein_x + self.fusion_offset - three_p_utr_length, self.graphic_protein_y), three_p_utr_length, self.graphic_protein_width, ec = 'black', fill = True, lw = self.protein_line_width, fc = '#808080', hatch = '////')
            self.ax.add_patch(three_p_utr_box_5)

    def get_domain_list(self, TR_ID):
        try:
            description = subprocess.check_output("grep %s %s" % (TR_ID, self.domain_table_dir), shell = True).rstrip('\n').split('\n')
        except:
            return []
        description[:] = [domain.split('\t') for domain in description]
        return description

    def gen_color(self):
        i = 0
        while i < len(self.colors):
            yield self.colors[i]
            i += 1

    def draw_domains(self, side, domain_list):
        included_domain_list = []

        fusion_cut = (self.TR_fragment_length_5 - self.five_p_utr_length_5) if side == '5' else (self.TR_length_3 - self.TR_fragment_length_3 - self.five_p_utr_length_3)    # in bp

        for domain in domain_list:
            domain_start = int(domain[2]) * 3  # in bp

            domain_end = int(domain[3]) * 3  # in bp
            is_included = (domain_end < fusion_cut) if side == '5' else (domain_start > fusion_cut)

            if is_included:
                included_domain_list.append(domain)
                graphic_offset = (self.five_p_utr_length_5 + domain_start) * self.display_ratio if side == '5' else (self.TR_length_3 - domain_start - self.five_p_utr_length_3) * self.display_ratio

                graphic_domain_length = (domain_end - domain_start) * self.display_ratio

                graphic_domain_x = (self.graphic_protein_x + graphic_offset) if side == '5' else (self.graphic_protein_x + self.graphic_protein_length - graphic_offset)

                domain_box = patches.Rectangle(
                (graphic_domain_x, self.graphic_protein_y), graphic_domain_length, self.graphic_protein_width, ec = 'black', alpha = 1, fill = True, lw = self.protein_line_width, fc = self.domain_colors.next()
                )

                self.ax.add_patch(domain_box)
        return included_domain_list

    def label_domains(self):
        utr_label_x = self.graphic_protein_x + 0.1
        utr_label_y = 1.75
        domain_lengend_figure_width = 0.05
        line_width = 0.07
        plt.text(utr_label_x, utr_label_y + domain_lengend_figure_width, "UTR", fontsize = self.font_size, ha = 'left', va = 'top')
        utr_label = patches.Rectangle(
        (utr_label_x - 0.1, utr_label_y), domain_lengend_figure_width, domain_lengend_figure_width, ec = 'black', alpha = 1, fill = True, lw = self.protein_line_width, fc = '#808080', hatch = '////'
        )
        self.ax.add_patch(utr_label)

        all_domain_list = self.included_domain_list_5 + self.included_domain_list_3 if self.potential_in_frame else self.included_domain_list_5
        domain_legend_x = utr_label_x
        domain_legend_y = utr_label_y - line_width

        for domain in all_domain_list:
            domain_name = domain[5]
            wrap_time = 1
            max_chars = 22
            if len(domain_name) > max_chars:
                wrap_time += (len(wrap(domain_name,max_chars)) - 1)
                domain_name = "\n".join(wrap(domain_name,max_chars))
            #print(domain_name)
            plt.text(domain_legend_x, domain_legend_y + domain_lengend_figure_width, domain_name, fontsize = self.font_size, ha = 'left', va = 'top')
            domain_lengend_figure = patches.Rectangle(
            (domain_legend_x - 0.1, domain_legend_y), domain_lengend_figure_width, domain_lengend_figure_width, ec = 'black', alpha = 1, fill = True, lw = self.protein_line_width, fc = self.domain_colors.next()
            )
            self.ax.add_patch(domain_lengend_figure)
            domain_legend_y -= line_width * wrap_time

    def check_out_of_frame(self):
        fusion_in_3p_utr = self.TR_fragment_length_3 < self.three_p_utr_length_3 # fusion falls in 3 prime utr region of 3 prime gene
        fusion_in_5p_utr = (self.TR_length_3 - self.TR_fragment_length_3) < self.five_p_utr_length_3 # fusion falls in 5 prime utr region of 3 prime gene
        three_p_gene_non_coding = (self.TR_fragment_length_3 == self.three_p_utr_length_3) # because I set three_p_utr_length_3 = TR_fragment_length_3 if no CDS is found
        return ((not self.potential_in_frame) and (not fusion_in_3p_utr) and (not three_p_gene_non_coding) and (not three_p_gene_non_coding))


    def draw_blank_box(self):
        # if (not self.potential_in_frame) and (self.TR_fragment_length_3 > self.three_p_utr_length_3) and (self.TR_fragment_length_3 != self.three_p_utr_length_3) and ((self.TR_length_3 - self.TR_fragment_length_3) > self.five_p_utr_length_3):

        # if (not self.potential_in_frame):
        if self.check_out_of_frame():
            white_box_x = self.graphic_protein_x + self.TR_fragment_length_5 * self.display_ratio
            white_box_length = self.TR_fragment_length_3 * self.display_ratio
            white_box = patches.Rectangle(
            (white_box_x, self.graphic_protein_y), white_box_length, self.graphic_protein_width, ec = 'black', alpha = 1, fill = True, lw = self.protein_line_width, fc = 'white'
            )
            self.ax.add_patch(white_box)

    def label_stop_codon(self, stop_info):
        if len(stop_info[0]) != 0:
            if (stop_info[0][-1] == "X") and self.check_out_of_frame():
                distance_from_fusion = (len(stop_info[0]) - 1 - int(stop_info[1])) * 3 - int(stop_info[2])

                mark_width = 0.03

                termination_mark_x = self.graphic_protein_x + (self.TR_fragment_length_5 + distance_from_fusion) * self.display_ratio - mark_width/2

                termination_mark_y = self.graphic_protein_y

                term_mark_1 = lines.Line2D(
                [termination_mark_x, termination_mark_x + mark_width],
                [termination_mark_y, termination_mark_y + mark_width],
                c = "red", alpha = 1, lw = 1.5
                )

                term_mark_2 = lines.Line2D(
                [termination_mark_x, termination_mark_x + mark_width],
                [termination_mark_y + mark_width, termination_mark_y],
                c = "red", alpha = 1, lw = 1.5
                )

                self.ax.add_line(term_mark_1)
                self.ax.add_line(term_mark_2)

    def draw_right_legend(self):

        right_lengend_x = 1.4
        right_legend_1_y = 1.75
        right_legend_2_y = right_legend_1_y - 0.2
        right_legend_3_y = right_legend_2_y - 0.2
        right_legend_4_y = 0.9
        right_legend_5_y = 0.7

        plt.text(right_lengend_x, right_legend_1_y, "Fusion Junction", fontsize = self.font_size, ha = 'left')
        plt.text(right_lengend_x, right_legend_2_y, "Termination", fontsize = self.font_size, ha = 'left')
        plt.text(right_lengend_x, right_legend_3_y, "Out of Frame", fontsize = self.font_size, ha = 'left')
        #plt.text(right_lengend_x, right_legend_4_y, "5' partner", fontsize = self.font_size, ha = 'center')
        #plt.text(right_lengend_x, right_legend_5_y, "3' partner", fontsize = self.font_size, ha = 'center')

        legend_figures_x = right_lengend_x - 0.1

        fusion_mark_x = legend_figures_x
        fusion_mark_lengend_y = right_legend_1_y

        fusion_mark = lines.Line2D(
        [fusion_mark_x - 0.05, fusion_mark_x, fusion_mark_x + 0.05],
        [fusion_mark_lengend_y + 0.05, fusion_mark_lengend_y, fusion_mark_lengend_y + 0.05],
        c = "black", alpha = 1, lw = 3
        )

        self.ax.add_line(fusion_mark)

        termination_mark_x = legend_figures_x - 0.025
        termination_mark_y = right_legend_2_y

        term_mark_1 = lines.Line2D(
        [termination_mark_x, termination_mark_x + 0.05],
        [termination_mark_y, termination_mark_y + 0.05],
        c = "red", alpha = 1, lw = 3
        )

        term_mark_2 = lines.Line2D(
        [termination_mark_x, termination_mark_x + 0.05],
        [termination_mark_y + 0.05, termination_mark_y],
        c = "red", alpha = 1, lw = 3
        )

        self.ax.add_line(term_mark_1)
        self.ax.add_line(term_mark_2)

        out_of_frame_mark_x = legend_figures_x - 0.025
        out_of_frame_mark_y = right_legend_3_y

        out_of_frame_mark = patches.Rectangle(
        (out_of_frame_mark_x, out_of_frame_mark_y), 0.05, 0.05, ec = 'black', alpha = 1, fill = False, lw = 1
        )

        self.ax.add_patch(out_of_frame_mark)

        #p5_mark_x = legend_figures_x - 0.025
        #p5_mark_y = right_legend_4_y

        #p5_mark = patches.Rectangle(
        #(p5_mark_x, p5_mark_y), 0.05, 0.05, ec = 'blue', alpha = 1, fill = False, lw = 1
        #)

        #self.ax.add_line(p5_mark)

        #p3_mark_x = legend_figures_x - 0.025
        #p3_mark_y = right_legend_5_y

        #p3_mark = patches.Rectangle(
        #(p3_mark_x, p3_mark_y), 0.05, 0.05, ec = 'red', alpha = 1, fill = False, lw = 1
        #)

        #self.ax.add_line(p3_mark)

    def add_tick_marks(self):
        # left side
        coding_start_left = min(self.TR_fragment_length_5, self.five_p_utr_length_5)
        coding_end_left = min(self.TR_fragment_length_5, self.TR_length_5 - self.three_p_utr_length_5)

        # right side
        five_p_utr_on_3p_fragment = self.five_p_utr_length_3 - (self.TR_length_3 - self.TR_fragment_length_3)
        coding_start_right = five_p_utr_on_3p_fragment + self.TR_fragment_length_5 if five_p_utr_on_3p_fragment > 0 else self.TR_fragment_length_5
        coding_end_right = self.TR_fragment_length_5 if self.check_out_of_frame() else self.fusion_TR_length - min(self.three_p_utr_length_3, self.TR_fragment_length_3)

        # conitinuous translated region
        if coding_end_left == coding_start_right:
            self.mark_helper(coding_start_left, coding_end_right, True)
        else:
            if self.check_out_of_frame():
                self.mark_helper(coding_start_left, coding_end_left, True)
            else:
                self.mark_helper(coding_start_left, coding_end_left, False)
                self.mark_helper(coding_start_right, coding_end_right, True)

    def mark_helper(self, coding_start, coding_end, aa_label):
        protein_length = int((coding_end - coding_start)/3)

        if protein_length > 0:
            number_of_digits = len(str(protein_length))
            first_digit = int(str(protein_length)[0])
            if first_digit >= 5:
                interval = int(1 * (10 ** (number_of_digits - 1))) * 2
            elif first_digit >= 2:
                interval = int(0.5 * (10 ** (number_of_digits - 1))) * 2
            else:
                interval = int(0.2 * (10 ** (number_of_digits - 1))) * 2
            if number_of_digits == 1:
                interval = 1
            if (protein_length * 3 * self.display_ratio < 0.3): # range too narrow, changed threashold from 0.25 to 0.3
                interval = protein_length

            number_of_ticks = int(protein_length/interval)

            for i in range(0, number_of_ticks + 1):
                ha = 'center' if i == 0 else 'right'
                offset = 3 * interval * self.display_ratio * i
                plt.text(self.graphic_protein_x + coding_start * self.display_ratio + offset, self.graphic_protein_y - 0.04, i * interval, fontsize = 6, ha = ha, va = 'bottom')

            # plt.text(self.graphic_protein_x + coding_start * self.display_ratio + aa_offset + 0.1, self.graphic_protein_y - 0.01, "aa", fontsize = 6, ha = 'center', va = 'top')
            if aa_label:
                aa_offset = 3 * interval * self.display_ratio * number_of_ticks
                plt.text(self.graphic_protein_x + coding_start * self.display_ratio + aa_offset + 0.05, self.graphic_protein_y - 0.04, "aa", fontsize = 6, ha = 'right', va = 'bottom')

    def label_genes(self):
        gene_label_x_5 = self.graphic_protein_x + self.TR_fragment_length_5 * self.display_ratio/2
        gene_label_x_3 = self.graphic_protein_x + self.graphic_protein_length - self.TR_fragment_length_3 * self.display_ratio/2
        gene_label_y = 0.2
        plt.text(gene_label_x_5, gene_label_y, self.gene_name_5, fontsize = 10, ha = 'center', va = 'top')
        plt.text(gene_label_x_3, gene_label_y, self.gene_name_3, fontsize = 10, ha = 'center', va = 'top')

    def label_terminals(self):
        offset = 0.08
        plt.text(self.graphic_protein_x - offset, self.graphic_protein_y + self.graphic_protein_width/2, "N", fontsize = 12, ha = 'right', va = 'center')
        plt.text(self.graphic_protein_x + self.graphic_protein_length + offset - 0.02, self.graphic_protein_y + self.graphic_protein_width/2, "C", fontsize = 12, ha = 'left', va = 'center')

    def highlight_53(self):
        highlight_box_5 = patches.Rectangle(
        (self.graphic_protein_x, self.graphic_protein_y), self.TR_fragment_length_5 * self.display_ratio, self.graphic_protein_width, ec = 'blue', alpha = 1, fill = False, lw = self.protein_line_width * 4)
        fusion_x = self.graphic_protein_x + self.TR_fragment_length_5 * self.display_ratio
        highlight_box_3 = patches.Rectangle(
        (fusion_x, self.graphic_protein_y), self.graphic_protein_x + self.graphic_protein_length - fusion_x, self.graphic_protein_width, ec = 'red', alpha = 1, fill = False, lw = self.protein_line_width * 4)
        self.ax.add_patch(highlight_box_3)
        self.ax.add_patch(highlight_box_5)

    def draw(self, output_dir):
        plt.axis('off')
        plt.savefig(output_dir, bbox_inches='tight')

def main(argv):

    t=time.time()

    #panel = PANEL_B(gene_ID_5, gene_ID_3, fusion_point_5, fusion_point_3, potential_in_frame, TR_frames, domain_table_dir, gene_model_dir, output_dir, trans_ID_5, trans_ID_3)
    panel = PANEL_B(argv)
    panel.grep_gene(panel.gene_ID_5, panel.gene_ID_3, panel.output_dir, "tmp.g.txt")
    panel.fig = plt.figure()
    panel.ax = panel.fig.add_subplot(111, aspect='equal')
    panel.ax.set_xlim([0, 2])
    panel.ax.set_ylim([0, 2])
    plt.axis('on')
    #print "TIME$$$$$$$$$$$$$$$$$$aaaa",time.time()-t
    #t=time.time()
    #panel.get_gene_names(panel.gene_ID_3, panel.gene_ID_5, panel.gene_model_dir)
    panel.get_gene_names(panel.gene_ID_3, panel.gene_ID_5, panel.gene_model_dir_2)
    #print "TIME$$$$$$$$$$$$$$$$$$bbbb",time.time()-t
    #t=time.time()
    #panel.get_directionality(panel.gene_ID_3, panel.gene_ID_5, panel.gene_model_dir)
    panel.get_directionality(panel.gene_ID_3, panel.gene_ID_5, panel.gene_model_dir_2)
    #print "TIME$$$$$$$$$$$$$$$$$$cccc",time.time()-t
    #t=time.time()
    panel.get_transcript(panel.gene_ID_3, panel.gene_ID_5)
    #print "TIME$$$$$$$$$$$$$$$$$$dddd",time.time()-t
    #t=time.time()
    panel.get_exons('5', panel.TR_ID_5)
    panel.get_exons('3', panel.TR_ID_3)
    #print "TIME$$$$$$$$$$$$$$$$$$eeee",time.time()-t
    #t=time.time()
    panel.fusion_exon_list_5 = panel.get_fusion_exons(panel.strand_5, '5', panel.fusion_point_5, panel.exons_5)
    panel.fusion_exon_list_3 = panel.get_fusion_exons(panel.strand_3, '3', panel.fusion_point_3, panel.exons_3)
    #print "TIME$$$$$$$$$$$$$$$$$$ffff",time.time()-t
    #t=time.time()
    panel.compute_fusion_parameters()
    panel.draw_fusion_protein()

    panel.highlight_fusion_junction()
    #print "TIME$$$$$$$$$$$$$$$$$gggg",time.time()-t
    #t=time.time()
    panel.five_p_utr_length_5, panel.three_p_utr_length_5 = panel.get_utr_lengths('5', panel.strand_5, panel.TR_ID_5, panel.exons_5)
    panel.five_p_utr_length_3, panel.three_p_utr_length_3 = panel.get_utr_lengths('3', panel.strand_3, panel.TR_ID_3, panel.exons_3)

    #print "TIME$$$$$$$$$$$$$$$$$$11111",time.time()-t
    #t=time.time()

    #print("5' UTR length for 5' partner: ")
    #print(panel.five_p_utr_length_5)

    #print("5' UTR length for 3' partner: ")
    #print(panel.five_p_utr_length_3)

    #print("3' UTR length for 3' partner: ")
    #print(panel.three_p_utr_length_3)

    panel.draw_utrs()

    panel.domain_list_5 = panel.get_domain_list(panel.TR_ID_5)
    panel.domain_list_3 = panel.get_domain_list(panel.TR_ID_3)

    #print "TIME$$$$$$$$$$$$$$$$$$2222222",time.time()-t
    #t=time.time()

    #print("all domains on the five prime gene:")
    #print(panel.domain_list_5)
    #print("all domains on the three prime gene:")
    #print(panel.domain_list_3)

    panel.domain_colors = panel.gen_color()
    panel.included_domain_list_5 = panel.draw_domains('5', panel.domain_list_5)
    panel.included_domain_list_3 = panel.draw_domains('3', panel.domain_list_3)

    #print "TIME$$$$$$$$$$$$$$$$$$33333",time.time()-t
    #t=time.time()

    panel.domain_colors = panel.gen_color()
    panel.label_domains()
    panel.label_stop_codon(panel.stop_info)
    panel.draw_blank_box()

    #print("domains on the five prime half:")
    #print(panel.included_domain_list_5)
    #print("domains on the three prime half:")
    #print(panel.included_domain_list_3)

    #print "TIME$$$$$$$$$$$$$$$$$$44444",time.time()-t
    #t=time.time()

    panel.draw_right_legend()
    panel.add_tick_marks()
    panel.label_genes()
    panel.label_terminals()
    panel.highlight_53()
    panel.draw(panel.output_dir)

    panel.rm_tmp(panel.output_dir,"tmp.g.txt")

    print "Time used: %0.2f" % (time.time()-t), "seconds."

if __name__ == '__main__':
    sys.exit(main(sys.argv))

