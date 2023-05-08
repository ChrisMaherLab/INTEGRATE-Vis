import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import re
import sys
import os
import time
from numpy import sin, cos, pi, linspace
import readsNum_getter as rng
import math


class PANEL_E:
    fusion5gene = ''
    fusion3gene = ''
    fusion5strand = ''
    fusion3strand = ''
    fusion5chrom = ''
    fusion3chrom = ''
    fusion5pos = ''
    fusion3pos = ''
    fusion5annot = ''
    fusion3annot = ''
    back5pos = ''
    back3pos = ''
    back5annot = ''
    back3annot = ''
    ideogram = ''
    
    def __init__(self, argv):
        #This is not expected to be called manually - should only
        #be initialized after this script is called by vis_e_fcirc.py
        self.filename = argv[1]
        self.fusion5gene = argv[2]
        self.fusion3gene = argv[3]
        self.fusion5strand = argv[4]
        self.fusion3strand = argv[5]
        self.fusion5chrom = argv[6]
        self.fusion3chrom = argv[7]
        self.fusion5pos = argv[8]
        self.fusion3pos = argv[9]
        self.fusion5annot = argv[10]
        self.fusion3annot = argv[11]
        self.back5pos = argv[12]
        self.back3pos = argv[13]
        self.back5annot = argv[14]
        self.back3annot = argv[15]
        self.ideogram = argv[16]
        self.genePredFile = argv[17]
        self.bamfile = argv[18] # =="none" if no bam file was provided
        self.showlegend = argv[19] #==0 if no, 1 if yes
        
        self.initializePlot()
        self.setDefaults()
        self.parseAnnotations()
    
        #Calculate relative sizes within fcircRNA
        #self.calcProportions()
        
        self.writeGeneNames()
        self.writeTranscriptInfo()
        
        #Draw exons of circle
        #3' Fusion Gene (backsplice donor) 
        if(self.fusion3exon==self.back5exon): #Only want one exon
            self.drawExon(self.TOPBOUNDARY + self.NEAREND,
                          self.BOTTOMBOUNDARY - self.NEAREND,
                          self.fusion3exon, self.COLOR3)    
        elif(int(self.fusion3exon)==int(self.back5exon) - 1):
            self.drawExon(self.TOPBOUNDARY + self.NEAREND,
                          self.TOPBOUNDARY + (0.5*pi),
                          self.back5exon, self.COLOR3)
            self.drawExon(self.BOTTOMBOUNDARY - (0.5*pi),
                          self.BOTTOMBOUNDARY - self.NEAREND,
                          self.fusion3exon, self.COLOR3)
        else:
            self.drawExon(self.TOPBOUNDARY + self.NEAREND,
                          self.TOPBOUNDARY + self.FAREND,
                          self.back5exon, self.COLOR3)
            self.drawExon(self.BOTTOMBOUNDARY - self.FAREND,
                          self.BOTTOMBOUNDARY - self.NEAREND,
                          self.fusion3exon, self.COLOR3)
            self.drawExon(self.TOPBOUNDARY + self.FAREND,
                          self.BOTTOMBOUNDARY - self.FAREND,
                          "...", self.COLOR3)
            #self.drawConnection(self.TOPBOUNDARY + self.FAREND,
            #                    self.BOTTOMBOUNDARY - self.FAREND, self.COLOR3)
            #self.drawCircHash(False,self.COLOR3)
            #self.drawCircArrow(self.TOPBOUNDARY + self.FAREND,
            #                   self.BOTTOMBOUNDARY - self.FAREND,self.COLOR5)            

        #5' Fusion Gene (backsplice acceptor)
        if(self.fusion5exon==self.back3exon):
            self.drawExon(self.BOTTOMBOUNDARY + self.NEAREND,
                          self.TOPBOUNDARY - self.NEAREND + (2*pi),
                          self.fusion5exon, self.COLOR5)
        elif(int(self.fusion5exon) == (int(self.back3exon) + 1)):
            self.drawExon(self.BOTTOMBOUNDARY + self.NEAREND,
                          self.BOTTOMBOUNDARY + (0.5 * pi),
                          self.fusion5exon, self.COLOR5)
            self.drawExon(self.TOPBOUNDARY - self.NEAREND,
                          self.TOPBOUNDARY - (0.5 * pi),
                          self.back3exon, self.COLOR5)            
        else:
            self.drawExon(self.BOTTOMBOUNDARY + self.NEAREND,
                          self.BOTTOMBOUNDARY + self.FAREND,
                          self.fusion5exon, self.COLOR5)
            self.drawExon(self.TOPBOUNDARY - self.NEAREND,
                          self.TOPBOUNDARY - self.FAREND,
                          self.back3exon, self.COLOR5)
            self.drawExon(self.BOTTOMBOUNDARY + self.FAREND,
                          self.BOTTOMBOUNDARY + self.FAREND + (0.2 * pi),
                          "...", self.COLOR5)
            #self.drawConnection(self.BOTTOMBOUNDARY + self.FAREND,
            #                    self.TOPBOUNDARY - self.FAREND + (2*pi), self.COLOR5)
            #self.drawCircHash(True,self.COLOR5)
            #self.drawCircArrow(self.BOTTOMBOUNDARY + self.FAREND,
            #                   self.TOPBOUNDARY - self.FAREND + (2*pi), self.COLOR3)

            
        #Fusion/backsplice junctions in circle
        self.drawConnection(self.TOPBOUNDARY - self.NEAREND, 
                            self.TOPBOUNDARY + self.NEAREND, "black")
        self.drawConnection(self.BOTTOMBOUNDARY - self.NEAREND,
                            self.BOTTOMBOUNDARY + self.NEAREND, "black")
        ##Counterclockwise arrows
        #self.drawCircArrow(self.TOPBOUNDARY + self.NEAREND,
        #                   self.TOPBOUNDARY - self.NEAREND, "black")
        #self.drawCircArrow(self.BOTTOMBOUNDARY + self.NEAREND,
        #                   self.BOTTOMBOUNDARY - self.NEAREND, "black")
        #Clockwise arrows
        self.drawCircArrow(self.TOPBOUNDARY - self.NEAREND,
                           self.TOPBOUNDARY + self.NEAREND, "black")
        self.drawCircArrow(self.BOTTOMBOUNDARY - self.NEAREND,
                           self.BOTTOMBOUNDARY + self.NEAREND, "black")
        

        #Draw genes, pre-fusion
        # Add logic to deal with "exon 0"
        # Modify vis_e_circ to pass the total number of exons in the transcript
        self.drawLinearGene(self.COLOR5, self.fusion5exon, self.fusion5strand, True)
        self.drawLinearGene(self.COLOR3, self.fusion3exon, self.fusion3strand, False)
        self.labelBreakpos("black", self.fusion5chrom, self.fusion5pos, self.fusion5strand, True)
        self.labelBreakpos("black", self.fusion3chrom, self.fusion3pos, self.fusion3strand, False)
        
        
        #Draw linear fusion
        self.drawLinearFusion(self.COLOR5, self.COLOR3, self.fusion5exon,
                              self.back3exon, self.fusion3exon, self.back5exon)
        self.labelBacksplicepos("black")
        
        #Connect linear to circular
        back5x = 0
        back3x = 0
        if self.fusion3exon == self.back5exon:
            back5x = (self.FUSIONXOFFSET*self.INTERGENIC_L) + self.EXON_L + 1
        else:
            back5x = (self.FUSIONXOFFSET*self.INTERGENIC_L) + (2*self.EXON_L) + (3*self.INTERGENIC_L) + 1
        if self.fusion5exon == self.back3exon:
            back3x = -1 * ((self.FUSIONXOFFSET*self.INTERGENIC_L) + self.EXON_L + 1)
        else:
            back3x = -1 * ((self.FUSIONXOFFSET*self.INTERGENIC_L) + (2*self.EXON_L) + (3*self.INTERGENIC_L) + 1)
        self.connectLinearToCirc(back5x,back3x,
                                 self.TOPBOUNDARY-self.NEAREND,
                                 self.TOPBOUNDARY+self.NEAREND)
                                 #self.BOTTOMBOUNDARY-self.NEAREND,
                                 #self.BOTTOMBOUNDARY+self.NEAREND)
        
        
        # Cytoband
        self.getCytobandInfo()
        self.drawChromosomes()
        self.highlightRegions()
        
        # [Optional] Read support
        if self.bamfile != "none":
            self.addReadSupport()
            if self.showlegend == "1":
                self.drawLegend()
        
        self.draw()
    
    
    ## NOT CALLED
    def calcProportions(self):
        #As written, this assumes that the length of the fcircRNA is based
        # on breakpoints/backsplices (does not account for exonic vs intronic lengths).
        # Attempt to draw the circle such that if x% of the distance is from gene A,
        # then gene A makes up x% of the circle.
        # Current issues: 1) If one gene takes up to much space, things become
        # illegible and may overlap. 2) Doesn't account for intronic distance
        # 3) circle size is roughly proportional, but exon/intron sizes are not
        geneALength = int(self.back3pos) - int(self.fusion5pos)
        if(geneALength < 0):
            geneALength = geneALength * -1
        geneBLength = int(self.back5pos) - int(self.fusion3pos)
        if(geneBLength < 0):
            geneBLength = geneBLength * -1
        totalLength = geneALength + geneBLength
        propFromA = (geneALength * 1.0) / totalLength
        if (propFromA < 0.33):
            propFromA = 0.33
        if (propFromA > 0.66):
            propFromA = 0.66
        fusionBoundary = self.BOTTOMBOUNDARY - (propFromA * (2 * pi))
        self.TOPBOUNDARY = fusionBoundary
       
    def initializePlot(self):
        self.fig= plt.figure()
        self.ax = self.fig.add_subplot(111, aspect='equal')
    
    def setDefaults(self):
        self.OUTER_R = 9.5 # 6 # 9
        self.INNER_R = 7 # 4
        self.INTERGENIC_L = 2
        self.EXON_H = 2
        self.EXON_L = 4
        self.FUSIONXOFFSET = 0.4
        self.LINEARXOFFSET = 1.1
        self.LINEARYOFFSET = 27 #New
        self.FUSIONYOFFSET = 15 #New
        self.CHROMYOFFSET = 45
        self.TOPBOUNDARY = 2.5 * pi
        self.BOTTOMBOUNDARY = 3.5 * pi
        self.NEAREND = 0.15 * pi #originally 0.1
        self.FAREND = 0.40 * pi #originally 0.45
        self.COLOR5 = "red"
        self.COLOR3 = "blue"
        self.EXONTEXTSIZE = 8
        self.centromereRadius = 0.5 * self.EXON_H
        self.OUTERTRANSCRIPTEDGE = 31
        self.INNERTRANSCRIPTEDGE = 2
        self.TRANSCRIPTTOP = 34
        self.TRANSCRIPTBOTTOM = 35.75
        self.HIGHLIGHTBOXWIDTH = 1
        
        
        # self.OUTER_R = 6
        # self.INNER_R = 4
        # self.INTERGENIC_L = 2
        # self.EXON_H = 2
        # self.EXON_L = 4
        # self.FUSIONXOFFSET = 0.4
        # self.LINEARXOFFSET = 1.1
        # self.LINEARYOFFSET = -27
        # self.FUSIONYOFFSET = -15
        # self.CHROMYOFFSET = -45
        # self.TOPBOUNDARY = 2.5 * pi
        # self.BOTTOMBOUNDARY = 3.5 * pi
        # self.NEAREND = 0.1 * pi
        # self.FAREND = 0.45 * pi
        # self.COLOR5 = "red"
        # self.COLOR3 = "blue"
        # self.EXONTEXTSIZE = 8
        # self.centromereRadius = 0.5 * self.EXON_H
        # self.OUTERTRANSCRIPTEDGE = 31
        # self.INNERTRANSCRIPTEDGE = 2
        # self.TRANSCRIPTTOP = -34
        # self.TRANSCRIPTBOTTOM = -35.25
        # self.HIGHLIGHTBOXWIDTH = 1
    
    def parseAnnotations(self):
        self.fusion5transcript = self.fusion5annot.split(":")[0]
        self.fusion5exon = self.fusion5annot.split(":")[1]
        self.fusion3transcript = self.fusion3annot.split(":")[0]
        self.fusion3exon = self.fusion3annot.split(":")[1]
        self.back5transcript = self.back5annot.split(":")[0]
        self.back5exon = self.back5annot.split(":")[1]
        self.back3transcript = self.back3annot.split(":")[0]
        self.back3exon = self.back3annot.split(":")[1]
        
    def writeGeneNames(self):
        y = self.LINEARYOFFSET + (3*self.EXON_H)
        #5' gene
        direction = -1
        xoffset=((8*self.LINEARXOFFSET)*direction) - (1.5*self.EXON_L) - self.INTERGENIC_L
        #structure command uses fontsize 9 for this
        self.ax.text(xoffset,y,self.fusion5gene,fontsize=9,
                     ha='center',va='bottom',color='black')
        #3' gene
        xoffset=(8*self.LINEARXOFFSET) + (1.5*self.EXON_L) + self.INTERGENIC_L
        self.ax.text(xoffset,y,self.fusion3gene,fontsize=9,
                     ha='center',va='bottom',color='black')
    
    def writeTranscriptInfo(self):
        self.transcript5 = subprocess.check_output("grep -w %s %s" % (self.fusion5annot.split(":")[0], self.genePredFile), shell = True).rstrip().split('\t')
        self.transcript3 = subprocess.check_output("grep -w %s %s" % (self.fusion3annot.split(":")[0], self.genePredFile), shell = True).rstrip().split('\t')
        numExons5 = self.transcript5[7]
        numExons3 = self.transcript3[7]
        y = self.LINEARYOFFSET + (2.1*self.EXON_H)
        #structure command uses fontsize 9
        #5' gene
        direction = -1
        xoffset=((8*self.LINEARXOFFSET)*direction) - (1.5*self.EXON_L) - self.INTERGENIC_L
        self.ax.text(xoffset,y,"(" + self.fusion5annot.split(":")[0] + ", " + numExons5 + " exons)",
                     ha='center',va='bottom',color='black',fontsize=8)
        #Highlight the text
        xhighlight = [self.OUTERTRANSCRIPTEDGE,self.OUTERTRANSCRIPTEDGE,
             self.INNERTRANSCRIPTEDGE,self.INNERTRANSCRIPTEDGE] 
        yhighlight = [self.TRANSCRIPTTOP,self.TRANSCRIPTBOTTOM,
             self.TRANSCRIPTBOTTOM,self.TRANSCRIPTTOP]
        self.ax.plot([i * -1 for i in xhighlight],yhighlight,color=self.COLOR5,lw=1)
        #3' gene
        self.ax.text(-1 * xoffset,y,"(" + self.fusion3annot.split(":")[0] + ", " + numExons3 + " exons)",
                     ha='center',va='bottom',color='black',fontsize=8)
        #Highlight the text
        #Figure out how to remove these magic numbers
        self.ax.plot(xhighlight,yhighlight,color=self.COLOR3,lw=1)
        
    def drawExon(self,start, end, label, color):
        r1 = self.OUTER_R
        arc_angles = linspace(start, end, 20)
        arc_xs = r1 * cos(arc_angles)
        arc_ys = r1 * sin(arc_angles)
        outerx = arc_xs.tolist()
        outery = arc_ys.tolist()
        outerx.reverse()
        outery.reverse()
        #draw inner arc
        r2 = self.INNER_R
        arc_angles = linspace(start, end, 20)
        arc_xs = r2 * cos(arc_angles)
        arc_ys = r2 * sin(arc_angles)
        innerx = arc_xs.tolist()
        innery = arc_ys.tolist()
        # Draw polygon and fill it
        xs = innerx + outerx + [innerx[0]]
        ys = innery + outery + [innery[0]]
        self.ax.fill(xs,ys, color = color)
        self.ax.plot(xs,ys,color='black',lw=1)
        r3 = r1 -((r1-r2)/2)
        x = ((r3 * cos(arc_angles[9])) + (r3*cos(arc_angles[10]))) / 2.0
        y = ((r3 * sin(arc_angles[9])) + (r3*sin(arc_angles[10]))) / 2.0
        self.ax.text(x,y,label,fontsize=self.EXONTEXTSIZE, 
                     horizontalalignment = "center", 
                     verticalalignment = "center", 
                     color="white")
        
    def drawConnection(self, start, end, color):
        r = self.OUTER_R - ((self.OUTER_R - self.INNER_R)/2)
        arc_angles = linspace(start, end, 30)
        arc_xso = r * cos(arc_angles)
        arc_yso = r * sin(arc_angles)
        self.ax.plot(arc_xso[2:28], arc_yso[2:28], color=color, lw=2)
        
    def drawCircArrow(self, end, start, color):
        arrowsize = 0.6
        arc_angles = linspace(start, end, 20)
        r = self.OUTER_R - ((self.OUTER_R-self.INNER_R)/2) + arrowsize
        topstartx = r * cos(arc_angles[7])
        topstarty = r * sin(arc_angles[7])
        r = self.OUTER_R - ((self.OUTER_R-self.INNER_R)/2)
        endx = r * cos(arc_angles[13])
        endy = r * sin(arc_angles[13])
        r = r - arrowsize
        bottomstartx = r * cos(arc_angles[7])
        bottomstarty = r * sin(arc_angles[7])
        xs = [topstartx,endx,bottomstartx,topstartx]
        ys = [topstarty,endy,bottomstarty,topstarty]
        self.ax.fill(xs,ys,color=color)
        
    def drawCircHash(self,isSpliceAcceptor,color):
        x0 = self.INNER_R + 0.1 #4.1
        x1 = self.OUTER_R - 0.1 #5.9
        if not isSpliceAcceptor:
            x0 = x0 * -1
            x1 = x1 * -1
        y1 = 0.5
        y2 = -0.5
        self.ax.plot([x0,x1],[y1,y1],color=color,lw=1)
        self.ax.plot([x0,x1],[y2,y2],color=color,lw=1)

    def drawLinearGene(self, color, exon, strand, is5p):
        #Add logic for when <3 exons are needed
        exon = int(exon)
        direction = 1
        if(is5p):
            direction = -1
        xoffset=(8*self.LINEARXOFFSET)*direction
        exonLength=self.EXON_L * direction 
        intergenicLength=self.INTERGENIC_L*direction 
        exonHeight=self.EXON_H + 0.1
        yoffset=self.LINEARYOFFSET 
        #Draw closest exon to midline if an exon exists there
        drawMidExon = True
        if(is5p and strand=="-" and exon==1):
           drawMidExon = False
        if(is5p is False and strand=="+" and exon==1):
            drawMidExon = False
        xs=[xoffset,xoffset+exonLength,xoffset+exonLength,xoffset]
        ys=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight]
        xsl=[xoffset,xoffset+exonLength,xoffset+exonLength,xoffset,xoffset]
        ysl=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight,yoffset]
        if(drawMidExon):
            self.ax.fill(xs,ys,color=color)
        #Next distal exon
        xs2 = [x+intergenicLength+exonLength for x in xs]
        self.ax.fill(xs2,ys,color=color)
        drawDistalExon = True
        if(is5p and strand=="+" and exon==1):
            drawDistalExon = False
        if(is5p is False and strand=="-" and exon==1):
            drawDistalExon = False
        #Most distal exon
        xs3 = [x+intergenicLength+exonLength for x in xs2]
        
        if(drawDistalExon):
            self.ax.fill(xs3,ys,color=color)
        #connect exons
        self.ax.plot([xs[1],xs2[0]], [ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                     color=color, lw=2)
        if(drawMidExon):
            self.drawLinearArrow(xs[1],xs2[0],yoffset+(0.5*exonHeight),strand,color)
        if(drawDistalExon):
            self.drawLinearArrow(xs2[1],xs3[0],yoffset+(0.5*exonHeight),strand,color)
            self.ax.plot([xs2[1],xs3[0]],[yoffset+(0.5*exonHeight),yoffset+(0.5*exonHeight)],
                                          color=color,lw=2)

        
        #Extension
        if(drawDistalExon):
            self.ax.plot([xs3[1],xs3[1]+intergenicLength], [ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                         color=color, lw=2)
        else:
            self.ax.plot([xs2[1],xs2[1]+intergenicLength],[ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                         color=color, lw=2)
        if(drawMidExon):
            self.ax.plot([xs[0],xs[0]-intergenicLength], [ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                         color=color, lw=2)
        else:
            self.ax.plot([xs2[0],xs2[0]-intergenicLength],[ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                         color=color, lw=2)
        #Draw outlines
        if(True): #Not sure this feature is being kept
            xsl=[xoffset,xoffset+exonLength,xoffset+exonLength,xoffset,xoffset]
            ysl=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight,yoffset]
            xs2l = [x+intergenicLength+exonLength for x in xsl]
            xs3l = [x+intergenicLength+exonLength for x in xs2l]
            if(drawMidExon):
                self.ax.plot(xsl,ysl,color='black',lw=1)
            self.ax.plot(xs2l,ysl,color='black',lw=1)
            if(drawDistalExon):
                self.ax.plot(xs3l,ysl,color='black',lw=1)
        #Label exons
        proximalExon=0
        distalExon=0
        if(is5p):
            if(strand=="+"):
                proximalExon=exon+1
                distalExon=exon-1
            else:
                proximalExon=exon-1
                distalExon=exon+1
        else:
            if(strand=="+"):
                proximalExon=exon-1
                distalExon=exon+1
            else:
                proximalExon=exon+1
                distalExon=exon-1
        #Proximal exon
        if(drawMidExon):
            self.ax.text(xoffset+(0.5*exonLength),yoffset+(0.5*exonHeight),
                         proximalExon,
                         fontsize=self.EXONTEXTSIZE,
                         horizontalalignment="center",
                         verticalalignment="center",
                         color="white")
        #Middle exon
        self.ax.text(xs2[0]+(0.5*exonLength),yoffset+(0.5*exonHeight),
                exon,
                fontsize=self.EXONTEXTSIZE,
                horizontalalignment="center",
                verticalalignment="center",
                color="white")
        #Distal exon
        if(drawDistalExon):
            self.ax.text(xs3[0]+(0.5*exonLength),yoffset+(0.5*exonHeight),
                         distalExon,
                         fontsize=self.EXONTEXTSIZE,
                         horizontalalignment="center",
                         verticalalignment="center",
                         color="white")
        #Draw junction markers
        xstart=0
        if(is5p):
            if(strand=="+"):
                xstart=xs2[0]
            else:
                xstart=xs2[1]
        else:
            if(strand=="+"):
                xstart=xs2[0]
            else:
                xstart=xs2[1]
        xend=0
        if(is5p):
            xend=-((self.FUSIONXOFFSET*self.INTERGENIC_L) + 1) 
        else:
            xend=(self.FUSIONXOFFSET*self.INTERGENIC_L) + 1
        ystart=self.LINEARYOFFSET #+ self.EXON_H
        yend=self.FUSIONYOFFSET + self.EXON_H #- 0.2
        self.ax.plot([xstart,xend],[ystart,yend],
                     color="black",lw=1)
        
    def labelBreakpos(self, color, chrm, pos, strand, is5p):
        y = self.LINEARYOFFSET - self.EXON_H #Just below exons
        direction = 1
        if(is5p):
            direction = -1
        xoffset = (8*self.LINEARXOFFSET)*direction
        exonLength = self.EXON_L*direction
        intergenicLength = self.INTERGENIC_L*direction
        x = xoffset + exonLength + intergenicLength
        if (strand=="-" and is5p) or (strand=="-" and not is5p):
            x += exonLength
        # structure command uses font size of 6 here
        if is5p:
            self.ax.text(x-0.5,y,"chr" + str(chrm) + ":" + str(pos) + " ", fontsize=6, 
                         ha='right', va='bottom', color=color)
            self.ax.plot([x,x],[self.LINEARYOFFSET-0.3,y],color="black",lw=2)
        else:
            self.ax.text(x+0.5,y," chr" + str(chrm) + ":" + str(pos),fontsize=6,
                         ha='left', va='bottom', color=color)
            self.ax.plot([x,x],[self.LINEARYOFFSET-0.3,y],color='black',lw=2)
        
        
    def drawLinearArrow(self,x1,x2,y,strand,color):
        diff=x1-x2
        if(diff<0):
            diff=-1*diff
        left=x1
        right=x2
        if(x2<x1):
            left=x2
            right=x1
        begin=0
        end=0
        if(strand=="+"):
            begin=left
            end=right
        else:
            begin=right
            end=left
        xend=end
        xbegin=begin
        xbegin=begin+(0.5*diff)-(0.2*diff)
        xend=xbegin+(0.4*diff)
        if(begin>end):
            xbegin=begin-(0.5*diff)+(0.2*diff)
            xend=xbegin-(0.4*diff)
        xs=[xbegin,xend,xbegin]
        ys=[y+(0.3*self.EXON_H),y,y-(0.3*self.EXON_H)]
        self.ax.fill(xs,ys,color=color)
    
    def drawLinearFusion(self, color5, color3, exon5a, exon5b, exon3a, exon3b):
        yoffset=self.FUSIONYOFFSET 
        intergenicLength=self.INTERGENIC_L 
        exonHeight=self.EXON_H + 0.1
        exonLength=self.EXON_L 
        xoffset=(self.FUSIONXOFFSET*intergenicLength) + 1
        xs=[-xoffset,-xoffset-exonLength,-xoffset-exonLength,-xoffset]
        xsl=[-xoffset,-xoffset-exonLength,-xoffset-exonLength,-xoffset,-xoffset] #For outline
        ys=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight]
        ysl=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight,yoffset] #For outline
        ##Junction
        #Draw junction before exons, to ensure black doesn't overlap the exon
        self.ax.plot([-xoffset+0.28,xoffset-0.28],[ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
           color="black",lw=2)
        self.drawLinearArrow(-xoffset,xoffset,ys[0]+(0.5*exonHeight),"+","black")
        ##5' gene
        #Exons
        self.ax.fill(xs,ys,color=color5)
        
        mod = 3
        if(exon5a=="1"):
            mod = 1
        xs2 = [x-(mod*intergenicLength)-exonLength for x in xs]
        xs2l = [x-(mod*intergenicLength)-exonLength for x in xsl]
        self.ax.plot([xs[1],xs2[0]],[ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                         color=color5, lw=2)
        if(exon5a!="1"):
            self.ax.fill(xs2,ys,color=color5)
            
            self.ax.plot([xs2[1],xs2[1]-intergenicLength], [ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                         color=color5, lw=2)
            #hashs
            self.ax.plot([xs2[0]+(1.25*intergenicLength),xs2[0]+(1.40*intergenicLength)],
                         [ys[3]-0.1,ys[0]+0.1],
                         color=color5, lw=2)
            self.ax.plot([xs2[0]+(1.60*intergenicLength),xs2[0]+(1.75*intergenicLength)],
                         [ys[3]-0.1,ys[0]+0.1],
                         color=color5, lw=2)
            self.ax.plot(xs2l,ysl,color='black',lw=1) #outline distal exon
        self.ax.plot(xsl,ysl,color='black',lw=1) #outline proximal exon
        #Label 5' end
        #self.ax.text(xs2[1]-(1.5*intergenicLength),
        #             yoffset+(0.5*exonHeight),
        #             "5'",
        #             fontsize=self.EXONTEXTSIZE,
        #             horizontalalignment="center",
        #             verticalalignment="center",
        #             color="black")
        
        ##3' gene
        #Exons
        xs=[xoffset,xoffset+exonLength,xoffset+exonLength,xoffset]
        ys=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight]
        xsl=[xoffset,xoffset+exonLength,xoffset+exonLength,xoffset,xoffset] #outline
        ysl=[yoffset,yoffset,yoffset+exonHeight,yoffset+exonHeight,yoffset] #outline
        self.ax.fill(xs,ys,color=color3)
        xs2 = [x+(3*intergenicLength)+exonLength for x in xs]
        xs2l = [x+(3*intergenicLength)+exonLength for x in xsl]
        self.ax.fill(xs2,ys,color=color3)
        self.ax.plot([xs[1],xs2[0]], [ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                     color=color3, lw=2)
        self.ax.plot([xs2[1],xs2[1]+intergenicLength],[ys[0]+(0.5*exonHeight),ys[0]+(0.5*exonHeight)],
                     color=color3, lw=2)
        #hashes
        self.ax.plot([xs2[0]-(1.40*intergenicLength),xs2[0]-(1.25*intergenicLength)],
                     [ys[3]-0.1,ys[0]+0.1],
                     color=color3, lw=2)
        self.ax.plot([xs2[0]-(1.75*intergenicLength),xs2[0]-(1.60*intergenicLength)],
                     [ys[3]-0.1,ys[0]+0.1],
                     color=color3, lw=2)  
        self.ax.plot(xsl,ysl,color='black',lw=1) #outline
        self.ax.plot(xs2l,ysl,color='black',lw=1) #outline
    
        #Labels
        if exon5b == exon5a: #Only 1 exon in backsplice
            exon5b = 1
        self.ax.text(-xoffset-(0.5*exonLength),yoffset+(0.5*exonHeight),
                     exon5a, fontsize=self.EXONTEXTSIZE,
                     horizontalalignment="center",
                     verticalalignment="center",
                     color="white")
        self.ax.text(-xoffset-(0.5*exonLength)-(3*intergenicLength)-exonLength,yoffset+(0.5*exonHeight),
                     exon5b, fontsize=self.EXONTEXTSIZE,
                     horizontalalignment="center",
                     verticalalignment="center",
                     color="white")
        if exon3a == exon3b:
            exon3b = int(exon3b) + 1
        self.ax.text(xoffset+(0.5*exonLength),yoffset+(0.5*exonHeight),
                     exon3a,fontsize=self.EXONTEXTSIZE,
                     horizontalalignment="center",
                     verticalalignment="center",
                     color="white")
        self.ax.text(xoffset+(0.5*exonLength)+(3*intergenicLength)+exonLength,yoffset+(0.5*exonHeight),
                     exon3b,fontsize=self.EXONTEXTSIZE,
                     horizontalalignment="center",
                     verticalalignment="center",
                     color="white")
  
    def labelBacksplicepos(self,color):
        y = self.FUSIONYOFFSET - self.EXON_H
        # For 5' fusion gene
        direction = -1
        xoffset = direction * ((self.FUSIONXOFFSET*self.INTERGENIC_L) + 1)
        if(self.fusion5exon == self.back3exon):# or int(self.fusion5exon) == 1:
            xoffset -= self.EXON_L
        else:
            xoffset = xoffset - (2*self.EXON_L) - (3 * self.INTERGENIC_L)
        self.ax.text(xoffset-0.5, y, str(self.fusion5chrom) + ":" + str(self.back3pos) + " ",
                     ha='right', va='bottom', color=color, fontsize=6)
        self.ax.plot([xoffset,xoffset],[self.FUSIONYOFFSET-0.3,y],color="black",lw=2)
        # For 3' fusion gene
        direction = 1
        xoffset = (self.FUSIONXOFFSET * self.INTERGENIC_L) + 1
        right = True
        if(self.fusion3exon == self.back5exon): # or int(self.fusion3exon) == 1:
            xoffset += self.EXON_L
            right = False
        else:
            xoffset = xoffset + (2*self.EXON_L) + (3*self.INTERGENIC_L)
        #if right:
        #    self.ax.text(xoffset, y, str(self.fusion3chrom) + ":" + str(self.back5pos) + " ",
        #                 ha='right', va='bottom', color=color, fontsize = 6)
        #else:
        self.ax.text(xoffset+0.5, y, " " + str(self.fusion3chrom) + ":" + str(self.back5pos),
                     ha='left', va='bottom', color=color, fontsize = 6)
        self.ax.plot([xoffset,xoffset],[self.FUSIONYOFFSET-0.3,y],color='black',lw=2)

    
    
    def connectLinearToCirc(self,back5x,back3x,circExon5Pos,circExon3Pos):
        r1 = self.OUTER_R
        # Backsplice from 3' gene in fusion
        xb = r1*cos(circExon3Pos)
        ya=self.FUSIONYOFFSET#+self.EXON_H+0.1
        yb=r1*sin(circExon3Pos)
        self.ax.plot([back5x,xb],
                     [ya,yb],
                     color="black",lw=1)
        # Backsplice into 5' gene from fusion
        xb=r1*cos(circExon5Pos)
        ya=self.FUSIONYOFFSET#+self.EXON_H
        yb=r1*sin(circExon5Pos) - 0.2
        self.ax.plot([back3x,xb],
                     [ya,yb],
                     color="black",lw=1)
        
        
    def getCytobandInfo(self):
        #Read ideogram file for 5' gene
        self.fusion5chromAllbands = subprocess.check_output("grep -w chr%s %s" % (self.fusion5chrom, self.ideogram),
                                                         shell=True).rstrip('\n').split('\n')
        
        # Select band of interest and identify chromosome info for 5' gene
        prevRegion = "A"
        prevEnd = 1
        for i in range(0,len(self.fusion5chromAllbands)):
            region = self.fusion5chromAllbands[i].split('\t')
            if int(self.fusion5pos) >= int(region[1]) and int(self.fusion5pos) <= int(region[2]):
                self.fusion5Band = [int(region[1]),int(region[2]),region[3]] # Contains start,end and name
            if "p" in prevRegion and "q" in region[3]:
                self.fusion5shortArm = int(prevEnd)
            prevRegion = region[3]
            prevEnd = region[2]
            if i is len(self.fusion5chromAllbands) - 1:
                self.fusion5longArm = int(region[2])
        
            
        #Repeat above, for 3' gene
        self.fusion3chromAllbands = subprocess.check_output("grep -w chr%s %s" % (self.fusion3chrom, self.ideogram),
                                                            shell=True).rstrip('\n').split('\n')  
        prevRegion = "A"
        prevEnd = 1
        
        for i in range(0,len(self.fusion3chromAllbands)):
            region = self.fusion3chromAllbands[i].split('\t')
            if int(self.fusion3pos) >= int(region[1]) and int(self.fusion3pos) <= int(region[2]):
                self.fusion3Band = [int(region[1]),int(region[2]),region[3]]
            if "p" in prevRegion and "q" in region[3]:
                self.fusion3shortArm = int(prevEnd)
            prevRegion = region[3]
            prevEnd = region[2]
            if i is len(self.fusion3chromAllbands) - 1:
                self.fusion3longArm = int(region[2])
                
        #No longer needed
        self.fusion5chromAllbands = None
        self.fusion3chromAllbands = None
        
    
    def drawSingleChromosome(self, leftmostX, rightmostX, y, centromereRadius,
                             shortArm, longArm, roundingsize, buffer):
        shortArmProportionOfChrom = (1.0 * shortArm) / (longArm)
        fullLength = rightmostX - leftmostX
        shortArmLength = (shortArmProportionOfChrom * fullLength) #- centromereRadius
        longArmLength = ((1-shortArmProportionOfChrom) * fullLength) #- centromereRadius
        #Draw short arm
        self.ax.add_patch(matplotlib.patches.FancyBboxPatch(
            (leftmostX-centromereRadius,y-centromereRadius),
            shortArmLength,
            2*centromereRadius,
            facecolor='gray',ec='black',
            boxstyle=matplotlib.patches.BoxStyle("Round", rounding_size = roundingsize)))
        #Draw long arm
        self.ax.add_patch(matplotlib.patches.FancyBboxPatch(
            (rightmostX-longArmLength+centromereRadius,y-centromereRadius),
            longArmLength,
            2*centromereRadius,
            facecolor='gray',ec='black',
            boxstyle=matplotlib.patches.BoxStyle("Round", rounding_size = roundingsize)))
        #Draw centromere
        self.ax.add_patch(matplotlib.patches.Circle(
            (leftmostX+shortArmLength,y),
            centromereRadius - buffer,
            ec='none',
            color='black'))
        
        
    def drawChromosomes(self):
        yoffset = self.CHROMYOFFSET
        
        
        if self.fusion5chrom == self.fusion3chrom:
            leftmostX = -18 # ~Center of 5' gene name
            rightmostX = 18
            roundingsize = 1
            buffer = 0.1
            self.drawSingleChromosome(leftmostX, rightmostX, yoffset,
                                 self.centromereRadius, self.fusion5shortArm, self.fusion5longArm, 
                                 roundingsize, buffer)
            
        else:  # Draw 2 chromosomes
            leftmostX = ((8*self.LINEARXOFFSET)* - 1) - (3*self.INTERGENIC_L) - (3*self.EXON_L)
            rightmostX = ((8*self.LINEARXOFFSET) - self.INTERGENIC_L) * -1
            roundingsize = 1
            buffer = 0.1
            self.drawSingleChromosome(leftmostX, rightmostX, yoffset,
                                 self.centromereRadius, self.fusion5shortArm, self.fusion5longArm,
                                 roundingsize, buffer)
            leftmostX3 = -1 * rightmostX
            rightmostX3 = -1 * leftmostX
            self.drawSingleChromosome(leftmostX3, rightmostX3, yoffset,
                                 self.centromereRadius, self.fusion3shortArm, self.fusion3longArm,
                                 roundingsize, buffer)
     
            
    def highlightOneRegion(self,leftmostX,rightmostX,ybottom,
                           band,leftband,rightband,
                           shortArm,longArm,color,source):
        #Calculate cytoband position
        totalLength = rightmostX - leftmostX
        leftPosProp = (1.0 * leftband) / longArm
        buffDist = 0
        if "q" in band:
            buffDist = 2 * self.centromereRadius
        else:
            buffDist = -1 * self.centromereRadius
        leftPosX = (leftPosProp * totalLength) + leftmostX + buffDist
        #rightPosProp = (1.0 * rightband) / longArm
        #width = (rightPosProp * totalLength) + leftmostX + buffDist - leftPosX
        width = self.HIGHLIGHTBOXWIDTH # Use consistent width, despite band size
        height = 3.5 * self.centromereRadius
        leftPosX = leftPosX - (0.5 * width)
        #Draw connection lines from cytoband to transcript name
        botx = leftPosX + (0.5*width)
        #topy = ybottom+height
        topx = (self.OUTERTRANSCRIPTEDGE + self.INNERTRANSCRIPTEDGE) / 2.0
        if source is "five":
            #See writeTranscriptInfo
            self.ax.plot([botx,-1 * topx],
                         [ybottom,self.TRANSCRIPTBOTTOM],color=self.COLOR5,lw=1)
        elif source is "three":
            self.ax.plot([botx,topx],[ybottom,self.TRANSCRIPTBOTTOM],color=self.COLOR3,lw=1)
        else:
            self.ax.plot([botx,-1 * topx],[ybottom,self.TRANSCRIPTBOTTOM],color=self.COLOR5,lw=1)
            self.ax.plot([botx,topx],[ybottom,self.TRANSCRIPTBOTTOM],color=self.COLOR3,lw=1)
        #Draw box around band
        #do this last so it overlaps the lines drawn above
        self.ax.add_patch(matplotlib.patches.Rectangle(
            (leftPosX,ybottom),
            width,
            height,
            ec=color, alpha=1,fill=False,lw=1))
        return botx
    
    
    def labelCytoband(self,x,y,position,chrm,band,addBuffer):
        ostring = chrm+band
        bufferString = " "
        if addBuffer:
            bufferString += "  "
        if position is "left":
            ostring = ostring + bufferString
        elif position is "right":
            ostring = bufferString + ostring
        self.ax.text(x,y,ostring,fontsize = 8,
                     ha=position,va='bottom',color='black')
            
    
    def highlightRegions(self):
        ybottom = self.CHROMYOFFSET - (1.5 * self.centromereRadius)
        ytext = self.CHROMYOFFSET + (3 * self.centromereRadius) 
        if self.fusion5chrom == self.fusion3chrom:
            leftmostX = -18 # same as in drawChromosomes()
            rightmostX = 18
            if self.fusion5Band[2] == self.fusion3Band[2]:
                regionX = self.highlightOneRegion(leftmostX,rightmostX,ybottom,
                                        self.fusion5Band[2],self.fusion5Band[0],
                                        self.fusion5Band[1],self.fusion5shortArm,
                                        self.fusion5longArm,'black',"both")
                self.labelCytoband(regionX,ytext,"center",self.fusion5chrom,self.fusion5Band[2],False)
            else:
                regionX5 = self.highlightOneRegion(leftmostX,rightmostX,ybottom,
                                        self.fusion5Band[2],self.fusion5Band[0],
                                        self.fusion5Band[1],self.fusion5shortArm,
                                        self.fusion5longArm,self.COLOR5,"five")
                regionX3 = self.highlightOneRegion(leftmostX,rightmostX,ybottom,
                                        self.fusion3Band[2],self.fusion3Band[0],
                                        self.fusion3Band[1],self.fusion3shortArm,
                                        self.fusion3longArm,self.COLOR3,"three")
                addBuffer = False
                if regionX5 < regionX3:
                    if (regionX3 - regionX5) < 3:
                        addBuffer = True
                    self.labelCytoband(regionX5,ytext,"right",self.fusion5chrom,self.fusion5Band[2],addBuffer)
                    self.labelCytoband(regionX3,ytext,"left",self.fusion3chrom,self.fusion3Band[2],addBuffer)
                else:
                    if (regionX5 - regionX3) < 3:
                        addBuffer = True
                    self.labelCytoband(regionX5,ytext,"left",self.fusion5chrom,self.fusion5Band[2],addBuffer)
                    self.labelCytoband(regionX3,ytext,"right",self.fusion3chrom,self.fusion3Band[2],addBuffer)
        else:
            leftmostX = ((8*self.LINEARXOFFSET)* - 1) - (3*self.INTERGENIC_L) - (3*self.EXON_L)
            rightmostX = ((8*self.LINEARXOFFSET) - self.INTERGENIC_L) * -1
            regionX5 = self.highlightOneRegion(leftmostX,rightmostX,ybottom,
                                    self.fusion5Band[2],self.fusion5Band[0],
                                    self.fusion5Band[1],self.fusion5shortArm,
                                    self.fusion5longArm,self.COLOR5,"five")
            leftmostX3 = -1 * rightmostX
            rightmostX3 = -1 * leftmostX
            regionX3 = self.highlightOneRegion(leftmostX3,rightmostX3,ybottom,
                                    self.fusion3Band[2],self.fusion3Band[0],
                                    self.fusion3Band[1],self.fusion3shortArm,
                                    self.fusion3longArm,self.COLOR3,"three")
            self.labelCytoband(regionX5,ytext,"right",self.fusion5chrom,self.fusion5Band[2],False)
            self.labelCytoband(regionX3,ytext,"left",self.fusion3chrom,self.fusion3Band[2],False)
     
    
    def drawJunctionSupport(self,propScore,isFusion,isLegend,isMax):

        lengthBuffer = 0.2
        radiusBuffer = 3
        xscaleShift = 20
        yscaleShift = -4
        #outerR = self.OUTER_R + radiusBuffer + propScore
        #innerR = self.OUTER_R + radiusBuffer
        radius = self.OUTER_R + radiusBuffer #+ propScore
        outerR = radius + (0.25 * propScore)
        innerR = radius - (0.25 * propScore)
        boundary = self.BOTTOMBOUNDARY
        if not isFusion:
            radiusBuffer = -1.5
            #outerR = self.INNER_R + radiusBuffer 
            #innerR = self.INNER_R + radiusBuffer - propScore
            radius = self.INNER_R + radiusBuffer #- propScore
            outerR = radius + (0.25 * propScore)
            innerR = radius - (0.25 * propScore)
            boundary = self.TOPBOUNDARY
        
        if isMax:
                yscaleShift= yscaleShift - 5
        #Counterclockwise exon
        arc_angles = linspace(boundary + (lengthBuffer * pi),boundary,20)
        arc_xs = outerR * cos(arc_angles)
        arc_ys = outerR * sin(arc_angles)
        #arc_xs = radius * cos(arc_angles)
        #arc_ys = radius * sin(arc_angles)
        outerx = arc_xs.tolist()
        outery = arc_ys.tolist()
        outerx.reverse()
        outery.reverse()
        xs = outerx
        ys = outery
        arc_xs = innerR * cos(arc_angles)
        arc_ys = innerR * sin(arc_angles)
        innerx = arc_xs.tolist()
        innery = arc_ys.tolist()
        xs = innerx + outerx + [innerx[0]]
        ys = innery + outery + [innery[0]]
        if isFusion and not isLegend:
            self.ax.fill(xs,ys,color=self.COLOR5)
            #self.ax.plot(xs,ys,color=self.COLOR5,linewidth=propScore)
        elif not isFusion and not isLegend:
            self.ax.fill(xs,ys,color=self.COLOR3)
            #self.ax.plot(xs,ys,color=self.COLOR3,linewidth=propScore)
        elif isFusion and isLegend:
            shiftedXs = [x + xscaleShift for x in xs]
            shiftedYs = [y + yscaleShift for y in ys]
            self.ax.fill(shiftedXs,shiftedYs,color=self.COLOR5)
            #self.ax.plot(shiftedXs,shiftedYs,color=self.COLOR5,linewidth=propScore)
        else:
            shiftedXs = [x + xscaleShift for x in xs]
            shiftedYs = [y + yscaleShift for y in ys]
            self.ax.fill(shiftedXs,shiftedYs,color=self.COLOR3)
            #self.ax.plot(shiftedXs,shiftedYs,color=self.COLOR3,linewidth=propScore)
        #Colockwise exon
        arc_angles = linspace(boundary - (lengthBuffer * pi),boundary,20)
        arc_xs = outerR * cos(arc_angles)
        arc_ys = outerR * sin(arc_angles)
        #arc_xs = radius * cos(arc_angles)
        #arc_ys = radius * sin(arc_angles)
        outerx = arc_xs.tolist()
        outery = arc_ys.tolist()
        outerx.reverse()
        outery.reverse()
        xs = outerx
        ys = outery
        arc_xs = innerR * cos(arc_angles)
        arc_ys = innerR * sin(arc_angles)
        innerx = arc_xs.tolist()
        innery = arc_ys.tolist()
        xs = innerx + outerx + [innerx[0]]
        ys = innery + outery + [innery[0]]
        if isFusion and not isLegend:
            self.ax.fill(xs,ys,color=self.COLOR3)
            #self.ax.plot(xs,ys,color=self.COLOR3,linewidth=propScore)
        elif not isFusion and not isLegend:
            self.ax.fill(xs,ys,color=self.COLOR5)
            #self.ax.plot(xs,ys,color=self.COLOR5,linewidth=propScore)
        elif isFusion and isLegend:
            shiftedXs = [x + xscaleShift for x in xs]
            shiftedYs = [y + yscaleShift for y in ys]
            self.ax.fill(shiftedXs,shiftedYs,color=self.COLOR3)
            #self.ax.plot(shiftedXs,shiftedYs,color=self.COLOR3,linewidth=propScore)
        else:
            shiftedXs = [x + xscaleShift for x in xs]
            shiftedYs = [y + yscaleShift for y in ys]
            self.ax.fill(shiftedXs,shiftedYs,color=self.COLOR5)
            #self.ax.plot(shiftedXs,shiftedYs,color=self.COLOR5,linewidth=propScore)
        
    
            
    def addReadSupport(self):
        # Base the visual purely on spanning reads, as it can be difficult
        # to differentiate between encompassing reads for the fusion
        # and encompassing reads for the backsplice in some instances
        #Linear fusion
        en,sp2,sp1 = rng.get_number(self.bamfile, 
                                    "chr" + str(self.fusion5chrom),self.fusion5strand,
                                    int(self.transcript5[3]),
                                    int(self.transcript5[4]),int(self.fusion5pos),
                                    "chr" + str(self.fusion3chrom),self.fusion3strand,
                                    int(self.transcript3[3]),
                                    int(self.transcript3[4]),int(self.fusion3pos))
        #Backsplice
        ben,bsp2,bsp1 = rng.get_number(self.bamfile,
                                       "chr" + str(self.fusion3chrom),self.fusion3strand,
                                       int(self.transcript3[3]),
                                       int(self.transcript3[4]),int(self.back5pos),
                                       "chr" + str(self.fusion5chrom),self.fusion5strand,
                                       int(self.transcript5[3]),
                                       int(self.transcript5[4]),int(self.back3pos))
        
        fusionReads = int(sp1) + int(sp2) + 0.0
        backReads = int(bsp1) + int(bsp2) + 0.0
        totalReads = fusionReads + backReads
        fusionProp = fusionReads / totalReads
        backProp = backReads / totalReads
        maxScore = 5.0
        backScore = backProp * maxScore
        fusionScore = fusionProp * maxScore
        bs = " Read"
        if backReads > 1:
            bs += "s"
        fs = " Read"
        if fusionReads > 1:
            fs += "s"
        self.drawJunctionSupport(fusionScore,True,False,False)
        self.drawJunctionSupport(backScore,False,False,False)
        self.ax.text(0,0.2,str(int(backReads)) + bs,
                     fontsize=self.EXONTEXTSIZE, ha='center',va='bottom',color='black')
        self.ax.text(0,-17,str(int(fusionReads)) + fs,
                     fontsize=self.EXONTEXTSIZE, ha='center',va='bottom',color='black')
        
        self.minscaleValue = 1.0
        if fusionReads > 10 or backReads > 10:
            self.minscaleValue = 5.0
        self.minScaleReadScore = (self.minscaleValue / totalReads) * maxScore
        self.maxscaleValue = 10
        if fusionReads > 10 or backReads > 10:
            if fusionReads >= backReads:
                self.maxscaleValue = int(math.ceil(fusionReads / 10.0)) * 10
            else:
                self.maxscaleValue = int(math.ceil(backReads / 10.0)) * 10
        if fusionReads > 100 or backReads > 100:
            if fusionReads >= backReads:
                self.maxscaleValue = int(math.ceil(fusionReads / 100.0)) * 100
            else:
                self.maxscaleValue = int(math.ceil(backReads / 100.0)) * 100       
        self.maxScaleReadScore = (self.maxscaleValue / totalReads) * maxScore
        
     
    def drawLegend(self):
        #self.drawJunctionSupport(self.singleReadScore,False,True)
        self.drawJunctionSupport(self.minScaleReadScore,False,True,False)
        self.drawJunctionSupport(self.maxScaleReadScore,False,True,True)
        self.ax.text(20,3,"Scale",ha='center',va='bottom',color='black',fontsize=self.EXONTEXTSIZE)
        s = str(int(self.minscaleValue)) + " Read"
        if self.minscaleValue > 1:
            s += "s"
        self.ax.text(33,0,s,ha='right',va='bottom',color='black',fontsize=6)
        s = str(int(self.maxscaleValue)) + " Reads"
        self.ax.text(33,-5,s,ha='right',va='bottom',color='black',fontsize=6)
        xs = [14,14,35,35,14]
        ys = [-8,6,6,-8,-8]
        self.ax.plot(xs,ys,color='black',lw=1)
        
 
            
        
    def draw(self):
        plt.gca().set_aspect('equal')
        plt.xlim(-38,38)
        plt.ylim(-25,55)
        #plt.rcParams['figure.figsize'] = [20, 10]
        plt.tight_layout()
        plt.axis('off')
        plt.savefig(self.filename, format="pdf")
        
      
        


def main(argv):
    t=time.time()
    PANEL_E(argv)
    print "Time used: %0.2f" % (time.time()-t), "seconds."


if __name__ == '__main__':
    sys.exit(main(sys.argv))