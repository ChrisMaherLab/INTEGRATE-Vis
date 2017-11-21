import sys, subprocess
import os
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import time
import pdb

def usage():
    print """
    panel_c -b <bam_file_directory> -5 <five_p_gene_ID_or_transcript_ID> -3 <three_p_gene_ID_or_transcript_ID> -m <gene_model_directory> -c <five_prime_gene_fusion_point> -d <three_prime_gene_fusion_point> -i <coverage_mode> -o <output-dir>
    Parameters:
        -b/--bam_file_directory      [Requested]
                                     [string:    bam file for a patient                                                  ]
        -5/--five_p_gene_ID          [Requested]
                                     [string:    5' gene ID.                                                             ]
        -3/--three_p_gene_ID         [Requested]
                                     [string:    3' gene ID.                                                             ]
        -c/--five_prime_gene_fusion_point   [Requested]
                                     [string:    5' gene fusion coordinate                                               ]
        -d/--three_prime_gene_fusion_point  [Requested]
                                     [string:    3' gene fusion coordinate                                               ]
        -m/--gene_model_directory    [Optional]
                                     [string:   provide your own gene model or default is used                          ]
        -i/--coverage_mode           [Optional]
                                     [string:   genome coverage (g) or transcript coverage (t)                          ]
        -o/--output-dir              [Optional]
                                     [string:   output directory.                                                       ]
        -1/--five_p_trans_ID         [Optional]
                                     [string:   5' transcript ID. Only choose this transcrtipt.                         ]
        -2/--three_p_trans_ID        [Optional]
                                     [string:   3' transcript ID. Only choose this transcrtipt.                         ]
        -x/--x-max                   [Optional]
                                     [value:    xlim for both gene partners                                             ]
        -y/--y-max                   [Optional]
                                     [value:    ylim for both gene partners                                             ]
        -q/--x-5p-max                [Optional]
                                     [value:    xlim for 5p gene, rule over -x/--x-max                                  ]
        -r/--y-5p-max                [Optional]
                                     [value:    ylim for 5p gene, rule over -y/--y-max                                  ]
        -s/--x-3p-max                [Optional]
                                     [value:    xlim for 3p gene, rule over -x/--x-max                                  ]
        -t/--y-3p-max                [Optional]
                                     [value:    ylim for 3p gene, rule over -y/--y-max                                  ]
        -u/--same-y                  [Optional]
                                     [          after the abovoe -max options, set y the same                           ]
        -n/--no-exon                 [Optional]
                                     [          no exon boundary lines                                                  ]
        -p/--transcript_frames       [Requested]
                                     [string:   stop info and shift frame status for different transcripts,
                                     [          a sub line of the annotated bedpe, seperated by space, comma, and '|'   ]

    Version:                         1.0.0
          """

#parameters
class PANEL_C:

    # requested:
    bam_file_directory=''
    five_prime_ID=''
    three_prime_ID=''
    five_prime_fusion=0
    three_prime_fusion=0

    # optional:
    gene_model_directory=''
    gene_model_dir_2 = ''
    output_dir='./pc.pdf'
    coverage_mode = 't'
    trans_ID_5 = 'None'
    trans_ID_3 = 'None'
    TR_frames = ''

    in_frame_TRs_5 = []
    in_frame_TRs_3 = []
    out_frame_TRs_3 = []


    m_xx=0.0
    m_yy=0.0
    m_xx_3=0.0
    m_yy_3=0.0
    m_xx_5=0.0
    m_yy_5=0.0
    is_same_y=False
    no_exon=False

    # computed:
    fig_size = [12,6]
    strand_3 = ''
    strand_5 = ''
    gene_name_5 = ''
    gene_name_3 = ''
    TR_3 = ''
    TR_5 = ''
    exons_3 = []
    exons_5 = []
    ch_3 = 0
    ch_5 = 0
    end_3 = 0
    end_5 = 0
    start_3 = 0
    start_5 = 0
    reads_3 = []
    reads_5 = []
    graph_junctions_3 = []
    graph_junctions_5 = []
    read_coords_3 = []
    read_coords_5 = []
    graph_coords_3 = []
    graph_coords_5 = []
    x_3 = []
    x_5 = []
    y_3 = []
    y_5 = []
    graph_fusion_3 = 0
    graph_fusion_5 = 0

    def __init__(self, argv):

        try:
            opts, args = getopt.getopt(argv,"hb:5:3:c:d:m:i:o:1:2:x:y:q:r:s:t:unp:", ["help", "bam_file_directory=", "five_prime_ID=", "three_prime_ID=", "five_prime_fusion=", "three_prime_fusion=", "gene_model_directory=", "coverage_mode=","output_dir=", "trans_ID_5=", "trans_ID_3=", "x-max=", "y-max=", "x-5p-max=", "y-5p-max=", "x-3p-max=", "y-3p-max=", "same-y=", "no-exon=", "transcript_frames="])
            #print(opts, args)

        except getopt.GetoptError:
            usage()
            sys.exit(1)

        for opt, arg in opts:
            if opt in ("-h","--help"):
                usage()
                sys.exit(0)
            elif opt in ("-b", "--bam_file_directory"):
                self.bam_file_directory = arg
            elif opt in ("-5", "--five_prime_ID"):
                self.five_prime_ID = arg
            elif opt in ("-3", "--three_prime_ID"):
                self.three_prime_ID = arg
            elif opt in ("-c", "--five_prime_fusion"):
                self.five_prime_fusion = int(arg)
            elif opt in ("-d", "--three_prime_fusion"):
                self.three_prime_fusion = int(arg)
            elif opt in ("-m", "--gene_model_directory"):
                self.gene_model_directory = arg
            elif opt in ("-i", "--coverage_mode"):
                if (arg != 't') and (arg != 'g'):
                    usage()
                    sys.exit(1)
                self.coverage_mode = arg
            elif opt in ("-o", "--output_dir"):
                self.output_dir = arg
            elif opt in ("-1", "--five_p_trans_ID"):
                self.trans_ID_5 = arg
            elif opt in ("-2", "--five_p_trans_ID"):
                self.trans_ID_3 = arg
            elif opt in ("-x", "--x-max"):
                self.m_xx = float(arg)
            elif opt in ("-y", "--y-max"):
                self.m_yy = float(arg)
            elif opt in ("-q", "--x-5p-max"):
                self.m_xx_5 = float(arg)
            elif opt in ("-r", "--y-5p-max"):
                self.m_yy_5 = float(arg)
            elif opt in ("-s", "--x-3p-max"):
                self.m_xx_3 = float(arg)
            elif opt in ("-t", "--y-3p-max"):
                self.m_yy_3 = float(arg)
            elif opt in ("-u", "--same-y"):
                self.is_same_y = True
            elif opt in ("-n", "--no-exon"):
                self.no_exon = True
            elif opt in ("-p", "--transcript_frames"):
                self.TR_frames = arg

        self.TR_frames = self.TR_frames.split(' ')

        transcript_lists = self.TR_frames[3].split(';')
        transcript_lists[:] = [TR_list.split('|') for TR_list in transcript_lists]

        self.in_frame_TRs_5 = [TR[:15] for TR in transcript_lists[0]]
        self.in_frame_TRs_3 = transcript_lists[1]
        self.out_frame_TRs_3 = transcript_lists[2]


    def get_gene_name(self):
        self.gene_name_5 = self.ID_to_name(self.five_prime_ID)
        self.gene_name_3 = self.ID_to_name(self.three_prime_ID)

    def get_transcript(self):
        #use gtf file offline
        if (self.five_prime_ID[:4] == 'ENSG'):
            self.TR_5 = self.find_best_transcript(self.five_prime_ID, self.strand_5, '5')
        else:
            self.TR_5 = self.five_prime_ID
        if (self.three_prime_ID[:4] == 'ENSG'):
            self.TR_3 = self.find_best_transcript(self.three_prime_ID, self.strand_3, '3')
        else:
            self.TR_3 = self.three_prime_ID
        #print("transcript ID: ", self.TR_5)
        #print("transcript ID: ", self.TR_3)

    def find_best_transcript(self, gene_ID, strand, side):
        #description = subprocess.check_output("grep %s %s | awk '$3 == \"transcript\"' | cut -f4,5,9" % (gene_ID, self.gene_model_directory), shell = True).rstrip('\n').split('\n')
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
            #exon_junctions = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f4,5" % (ID, self.gene_model_directory), shell = True).rstrip('\n').split('\n')
            exon_junctions = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f4,5" % (ID, self.gene_model_dir_2), shell = True).rstrip('\n').split('\n')
            exon_junctions[:] = [x.split('\t') for x in exon_junctions]
            TR_length = sum([(int(junction_pair[1]) - int(junction_pair[0])) for junction_pair in exon_junctions])
            length.append(TR_length)
            if ((side == '5') and (strand == '+')): # does fusion junction ever fall into the middle of an exon? consider exon junctions that occer before the fusion junction? use abs?
                distances = [abs(int(junction_pair[1]) - int(self.five_prime_fusion)) for junction_pair in exon_junctions]
            elif ((side == '3') and (strand == '-')):
                distances = [abs(int(junction_pair[1]) - int(self.three_prime_fusion)) for junction_pair in exon_junctions]
            elif ((side == '5') and (strand == '-')):
                distances = [abs(int(junction_pair[0]) - int(self.five_prime_fusion)) for junction_pair in exon_junctions]
            else:
                distances = [abs(int(junction_pair[0]) - int(self.three_prime_fusion)) for junction_pair in exon_junctions]
            min_dist_from_fusion.append(min(distances))

        in_frame = [int(ID in self.in_frame_TRs_5) for ID in TR_ID_list] if side == '5' else [int(ID in self.in_frame_TRs_3) for ID in TR_ID_list]

        TR_list = zip(TR_ID_list, biotype, in_frame, min_dist_from_fusion, tsl, length)
        TR_list.sort(key = lambda x: ((x[1] != 'protein_coding'), -x[2], x[3], x[4], -x[5]))
        #print("Transcript metrics: ", TR_list)
        return TR_list[0][0]

    def ID_to_name(self, gene_ID):
        #description = subprocess.check_output("grep %s %s |  awk '$3 == \"gene\"' | cut -f9" % (gene_ID, self.gene_model_directory), shell = True).rstrip('\n').split('; ')
        #print gene_ID,self.gene_model_dir_2
        description = subprocess.check_output("grep %s %s |  awk '$3 == \"gene\"' | cut -f9" % (gene_ID, self.gene_model_dir_2), shell = True).rstrip('\n').split('; ')
        description[:] = [x.split(' ') for x in description]
        return next(x for x in description if x[0] == 'gene_name')[1][1:-1]

    def get_exons(self):
        #self.exons_3 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f1,4,5" % (self.TR_3, self.gene_model_directory), shell = True).rstrip("\n").split("\n")
        #self.exons_5 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f1,4,5" % (self.TR_5, self.gene_model_directory), shell = True).rstrip("\n").split("\n")
        self.exons_3 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f1,4,5" % (self.TR_3, self.gene_model_dir_2), shell = True).rstrip("\n").split("\n")
        self.exons_5 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f1,4,5" % (self.TR_5, self.gene_model_dir_2), shell = True).rstrip("\n").split("\n")

        for i in range(0,len(self.exons_3)):
            self.exons_3[i] = self.exons_3[i].split("\t")

        for i in range(0,len(self.exons_5)):
            self.exons_5[i] = self.exons_5[i].split("\t")

        self.exons_3.sort(key = lambda x: x[1])
        self.exons_5.sort(key = lambda x: x[1])

        #print("\n" + "detected exons for " + self.gene_name_3 + ":" + '\n')
        #print(self.exons_3)
        #print("\n" + "detected exons for " + self.gene_name_5 + ":" + '\n')
        #print(self.exons_5)

        #not cast to int
        self.ch_3 = self.exons_3[0][0]
        self.ch_5 = self.exons_5[0][0]
        self.start_3 = int(self.exons_3[0][1])
        self.end_3 = int(self.exons_3[-1][2])
        self.start_5 = int(self.exons_5[0][1])
        self.end_5 = int(self.exons_5[-1][2])

    def get_directionality(self):
        #self.strand_3 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f7" % (self.three_prime_ID, self.gene_model_directory), shell = True).rstrip('\n')[0]
        #self.strand_5 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f7" % (self.five_prime_ID, self.gene_model_directory), shell = True).rstrip('\n')[0]
        #print(self.strand_3, self.strand_5)
        self.strand_3 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f7" % (self.three_prime_ID, self.gene_model_dir_2), shell = True).rstrip('\n')[0]
        self.strand_5 = subprocess.check_output("grep %s %s | awk '$3 == \"exon\"' | cut -f7" % (self.five_prime_ID, self.gene_model_dir_2), shell = True).rstrip('\n')[0]

    def get_reads(self):
        # if not os.isfile("%s.bam" % (gene_name_3)):
        subprocess.call("samtools view -b %s %s:%s-%s > %s.bam" % (self.bam_file_directory, self.ch_3, self.start_3, self.end_3, self.gene_name_3), shell = True)
        # if not os.isfile("%s.bam" % (gene_name_5)):
        subprocess.call("samtools view -b %s %s:%s-%s > %s.bam" % (self.bam_file_directory, self.ch_5, self.start_5, self.end_5, self.gene_name_5), shell = True)

        self.reads_3 = subprocess.check_output("bedtools genomecov -ibam %s.bam -bga" % (self.gene_name_3), shell = True)
        self.reads_5 = subprocess.check_output("bedtools genomecov -ibam %s.bam -bga" % (self.gene_name_5), shell = True)

        if self.reads_3: self.reads_3 = [read.split("\t") for read in self.reads_3.rstrip('\n').split('\n')]
        if self.reads_5: self.reads_5 = [read.split("\t") for read in self.reads_5.rstrip('\n').split('\n')]

        self.read_coords_3 = [0] * (self.end_3 - self.start_3 + 1)
        self.read_coords_5 = [0] * (self.end_5 - self.start_5 + 1)

        # read_coords contain intron
        if self.reads_3:
            for read in self.reads_3:
                count_start = max(self.start_3, int(read[1]))
                count_end = min(self.end_3, int(read[2]))
                count = int(read[3])
                for i in range(count_start - self.start_3, count_end - self.start_3 + 1):
                    self.read_coords_3[i] = count

        if self.reads_5:
            for read in self.reads_5:
                count_start = max(self.start_5, int(read[1]))
                count_end = min(self.end_5, int(read[2]))
                count = int(read[3])
                for i in range(count_start - self.start_5, count_end - self.start_5 + 1):
                    self.read_coords_5[i] = count

    def get_graph_coords(self):
        if self.coverage_mode == 'g': # make sure to highlight exons, optionally, shrink length of intron (not important for now)
            self.x_3 = range(0, len(self.read_coords_3))
            self.y_3 = [int(y) for y in self.read_coords_3]
            self.x_5 = range(0, len(self.read_coords_5))
            self.y_5 = [int(y) for y in self.read_coords_5]
        else:
            # graph_coords are exon only
            for exon in self.exons_3:
                exon_start = int(exon[1]) - int(self.start_3)
                exon_end = int(exon[2]) - int(self.start_3)
                self.graph_coords_3.extend(self.read_coords_3[exon_start:exon_end])
                # print("Exon start %s Exon end %s" % (exon_start,exon_end))
                # print(self.read_coords_3[exon_start:exon_end])
                self.graph_junctions_3.append(len(self.graph_coords_3))

            for exon in self.exons_5:
                exon_start = int(exon[1]) - int(self.start_5)
                exon_end = int(exon[2]) - int(self.start_5)
                self.graph_coords_5.extend(self.read_coords_5[exon_start:exon_end])
                self.graph_junctions_5.append(len(self.graph_coords_5))

            self.x_3 = range(0, len(self.graph_coords_3))
            self.y_3 = [int(y) for y in self.graph_coords_3]
            self.x_5 = range(0, len(self.graph_coords_5))
            self.y_5 = [int(y) for y in self.graph_coords_5]

            if (self.strand_3 == '-'):
                self.y_3.reverse()
                self.graph_junctions_3[:] = [(len(self.x_3) - junction) for junction in self.graph_junctions_3]
            if (self.strand_5 == '-'):
                self.y_5.reverse()
                self.graph_junctions_5[:] = [(len(self.x_5) - junction) for junction in self.graph_junctions_5]


    def get_fusion_coord(self):
        if self.coverage_mode == 'g':
            self.graph_fusion_3 = self.three_prime_fusion - self.start_3
            self.graph_fusion_5 = self.five_prime_fusion - self.start_5
        else:
            self.graph_fusion_3 = sum([(int(exon[2]) - int(exon[1])) for exon in self.exons_3 if int(exon[2]) <= self.three_prime_fusion])
            self.graph_fusion_5 = sum([(int(exon[2]) - int(exon[1])) for exon in self.exons_5 if int(exon[2]) <= self.five_prime_fusion])
        if (self.strand_3 == '-'):
            self.graph_fusion_3 = len(self.x_3) - self.graph_fusion_3
        if (self.strand_5 == '-'):
            self.graph_fusion_5 = len(self.x_5) - self.graph_fusion_5
        #print("fusion point on the graph: ", self.graph_fusion_3)
        #print("fusion point on the graph: ", self.graph_fusion_5)

    def draw(self):
        #pdb.set_trace()
        #tt=time.time()
        plt.rcParams["figure.figsize"] = self.fig_size
        fig = plt.figure()

        m_x_3=max(self.x_3)
        m_y_3=max(self.y_3)
        m_x_5=max(self.x_5)
        m_y_5=max(self.y_5)

        if self.m_xx>0:
            m_x_3=self.m_xx
            m_x_5=self.m_xx
        if self.m_yy>0:
            m_y_3=self.m_yy
            m_y_5=self.m_yy

        if self.m_xx_3>0:
            m_x_3=self.m_xx_3
        if self.m_yy_3>0:
            m_y_3=self.m_yy_3
        if self.m_xx_5>0:
            m_x_5=self.m_xx_5
        if self.m_xx_5>0:
            m_y_5=self.m_yy_5
        if self.is_same_y==True:
            m_y_3=max(m_y_3,m_y_5)
            m_y_5=m_y_3

        ax_3 = fig.add_subplot(122)
        #ax_3.scatter(self.x_3, self.y_3, s = 0.1, marker = ".")
        ax_3.plot(self.x_3, self.y_3, 'k', lw = 0.8)
        ax_3.set_xlim(-m_x_3*0.05, m_x_3*1.05)
        ax_3.set_ylim(-m_y_3*0.05, m_y_3*1.05)
        ax_3.set_title('%s' % (self.gene_name_3))
        #ax_3.set_xlabel('Nucleotide in %s' % (self.gene_name_3), fontsize = 14, color = 'black')
        ax_3.set_xlabel('Nucleotide', fontsize = 14, color = 'black')
        ax_3.set_ylabel('Read coverage', fontsize = 14, color = 'black')
        plt.axvline(x=self.graph_fusion_3, c = 'red', ls = '-', lw = 1, alpha=1)
        if self.no_exon==False:
            for xc in sorted(self.graph_junctions_3+[max(self.x_3)]):
                pp=len(self.x_3)*xc/max(self.x_3)-1
                if pp<0:
                    pp=0
                if pp>(len(self.x_3)-1):
                    pp=len(self.x_3)-1
                ymax=(self.y_3[pp]+m_y_3*0.05)/(m_y_3*1.1) if m_y_3 != 0 else 0
                #ymin=0.05/1.1
                ymin=0
                plt.axvline(x=xc, c = 'blue', ls = '-', lw = 0.4, alpha=1, ymin=ymin, ymax=ymax)
        #print "$$$$",time.time()-tt
        #tt=time.time()

        ax_5 = fig.add_subplot(121)
        #ax_5.scatter(self.x_5, self.y_5, s = 0.1, marker = ".")
        ax_5.plot(self.x_5, self.y_5, 'k', lw = 0.8)
        ax_5.set_xlim(-m_x_5*0.05, m_x_5*1.05)
        ax_5.set_ylim(-m_y_5*0.05, m_y_5*1.05)
        ax_5.set_title('%s' % (self.gene_name_5))
        #ax_5.set_xlabel('Position in %s' % (self.gene_name_5), fontsize = 14, color = 'black')
        ax_5.set_xlabel('Nucleotide', fontsize = 14, color = 'black')
        ax_5.set_ylabel('Read coverage', fontsize = 14, color = 'black')
        plt.axvline(x=self.graph_fusion_5, c = 'red', ls = '-', lw = 1, alpha=1)
        #print "$$$$",time.time()-tt
        #tt=time.time()
        #for xc in (self.graph_junctions_5+[max(self.x_5)]):
        #    plt.axvline(x=xc, c = 'blue', ls = '-', lw = 0.4, alpha=0.2)
        if self.no_exon==False:
            for xc in sorted(self.graph_junctions_5+[max(self.x_5)]):
                pp=len(self.x_5)*xc/max(self.x_5)-1
                if pp<0:
                    pp=0
                if pp>(len(self.x_5)-1):
                    pp=len(self.x_5)-1
                ymax=(self.y_5[pp]+m_y_5*0.05)/(m_y_5*1.1) if m_y_5 != 0 else 0
                #ymin=0.05/1.1
                ymin=0
                plt.axvline(x=xc, c = 'blue', ls = '-', lw = 0.4, alpha=1, ymin=ymin, ymax=ymax)
        #print "$$$$ 222",time.time()-tt
        #tt=time.time()

        fig.tight_layout()
        fig.savefig(self.output_dir,bbox_inches='tight')

        #print "$$$$",time.time()-tt


    def grep_gene(self, gene_ID_5, gene_ID_3, output_dir, suffix):
        if self.trans_ID_5!='None':
            subprocess.check_output("grep %s %s | awk '$3==\"gene\"{print}' > %s" % (gene_ID_5, self.gene_model_directory, output_dir+"."+suffix), shell = True)
            subprocess.check_output("grep %s %s >> %s" % (self.trans_ID_5, self.gene_model_directory, output_dir+"."+suffix), shell = True)
        else:
            subprocess.check_output("grep %s %s > %s" % (gene_ID_5, self.gene_model_directory, output_dir+"."+suffix), shell = True)
        if self.trans_ID_3!='None':
            subprocess.check_output("grep %s %s | awk '$3==\"gene\"{print}' >> %s" % (gene_ID_3, self.gene_model_directory, output_dir+"."+suffix), shell = True)
            subprocess.check_output("grep %s %s >> %s" % (self.trans_ID_3, self.gene_model_directory, output_dir+"."+suffix), shell = True)
        else:
            subprocess.check_output("grep %s %s >> %s" % (gene_ID_3, self.gene_model_directory, output_dir+"."+suffix), shell = True)
        self.gene_model_dir_2=output_dir+"."+suffix


    def rm_tmp(self, output_dir, suffix):
        os.remove(output_dir+"."+suffix)

def main(argv):
    t=time.time()
    panel = PANEL_C(argv[1:])
    #print "aaaa Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.grep_gene(panel.five_prime_ID, panel.three_prime_ID, panel.output_dir, "tmp.g.txt")
    #print "bbbb Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_gene_name()
    #print "cccc Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_directionality()
    #print "dddd Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_transcript()
    #print "eeee Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_exons()
    #print "eeee 111 Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_reads()
    #print "eeee 222 Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_graph_coords()
    #print "eeee 333Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.get_fusion_coord()
    #print "ffff Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.draw()
    #print "ffff 222 Time used: %0.2f" % (time.time()-t), "seconds."
    #t=time.time()
    panel.rm_tmp(panel.output_dir,"tmp.g.txt")
    panel.rm_tmp(panel.gene_name_3,"bam")
    panel.rm_tmp(panel.gene_name_5,"bam")
    print "Time used: %0.2f" % (time.time()-t), "seconds."

if __name__ == '__main__':
    sys.exit(main(sys.argv))

