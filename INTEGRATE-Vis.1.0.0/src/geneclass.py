#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 12:49:04 2022

@author: jacewebster
"""



class gene:
    def __init__(self,gene,strand,geneid):
        self.gene = gene
        self.strand = strand
        self.geneid = geneid
        self.exons = []
        
    
    def addExon(self, start, end):
        self.exons.append([start,end])
        
    
    def addExons(self, starts, ends):
        start = starts.split(",")[0:-1]
        end = ends.split(",")[0:-1]
        if self.strand == "-":
            start.reverse()
            end.reverse()
        for i in range(0,len(start)):
            self.addExon(start[i],end[i])
        