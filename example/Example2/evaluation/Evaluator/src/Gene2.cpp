/* 
 *  Created on: Mar 4, 2013
 *      Author: Jin Zhang
 */

#include "Gene2.h"



Gene::Gene() {
        // TODO Auto-generated constructor stub
}

Gene::~Gene() {
        // TODO Auto-generated destructor stub
}

bool myGeneSortFunc(gene_t i, gene_t j)
{
        if(i.chr.compare(j.chr)<0)
        {
                return true;
        }
        else if(i.chr.compare(j.chr)==0)
        {
                if(i.leftLimit<j.leftLimit)
                        return true;
                else
                        return false;
        }
        else
                return false;
}

bool myTransSortFunc(transcript_t i, transcript_t j)
{
        if(i.name2.compare(j.name2)<0)
                return true;
        else if(i.name2.compare(j.name2)==0)
        {
                if(i.chr.compare(j.chr)<0)
                {
                        return true;
                }
                else if(i.chr.compare(j.chr)==0)
                {
                        if(i.txStart<j.txStart)
                                return true;
                        else
                                return false;
                }
                else
                        return false;

        }
        else
                return false;
}

int Gene::loadGenesFromFile(char* file, int isTrunctOK) {

	uint32_t length=getFilelength(file);
	
        FILE *infile;

	infile = fopen(file, "r");
	if (!infile) {
		cerr<<"Couldn't open file for reading: "<<file<<"."<<endl;
		exit(1);
	}

	char * fileContent;

	try
	{
		fileContent= new char [length];
	}
	catch(exception& e)
	{
		cerr << "Trying to allocate Memory to load Genes. exception caught: " << e.what() << endl;
		return 1;
	}


	size_t result = fread (fileContent,1,length,infile);
	if (result != length)
	{
		cerr << "Fail to read the genes file"<<endl;
		exit (1);
	}

	char * line_buffer=fileContent;
	char* nextN=NULL;

	char* p=NULL;
	char* NC=NULL;
	char intChar [1024];

 	//int bin;
 	char nameC [1024];
 	char chromC [1024];
 	char strandC [1024];
 	uint32_t txStart;
 	uint32_t txEnd;
 	uint32_t cdsStart;//Dec 7, 2015
 	uint32_t cdsEnd;
 	int exonCount;
 	char exonStartsC [1000000]; //should be ok
 	char exonEndsC [1000000];	  //should be ok
 	//int score;
 	char name2C[1024];
 	char * chr;
 	string chrStr;

	int strand;

	nextN=strchr(line_buffer,'\n');
	line_buffer=nextN+1;



	int num=0;
	while (1) {
		nextN=strchr(line_buffer,'\n');

		//sscanf(line_buffer,"%d %s %s %s %u %u %u %u %d %s %s %d %s", &bin, nameC, chromC, strandC, &txStart, &txEnd,
		//	    			&cdsStart, &cdsEnd, &exonCount, exonStartsC, exonEndsC, &score, name2C);

		nextN[0]='\0';
		int numnum=sscanf(line_buffer,"%s\t%s\t%s\t%u\t%u\t%u\t%u\t%d\t%s\t%s\t%s", nameC, chromC, strandC, &txStart, &txEnd, &cdsStart, &cdsEnd, &exonCount, exonStartsC, exonEndsC, name2C);
		if(numnum!=11)
		{
			cerr<<"error loading genes at: "<<line_buffer<<endl;
			cerr<<"From 0.3.0, INTEGRATE also use cdsStart and cdsEnd, check your annotation file should have 11 columns."<<endl;
                	exit(1);
		}


	   	uint32_t * ps=new uint32_t [exonCount];
	   	uint32_t * pe=new uint32_t [exonCount];

	   	p=exonStartsC;
	   	for(int i=0;i<exonCount;i++)
	    {
	    	NC=strchr(p, ',');
	    	strncpy(intChar, p, NC-p);
	    	intChar[NC-p]='\0';
	    	ps[i]=atol(intChar);
	    	p=NC+1;
	    }

	   	p=exonEndsC;
	   	for(int i=0;i<exonCount;i++)
	   	{
	   		NC=strchr(p, ',');

	   		strncpy(intChar, p, NC-p);
	   		intChar[NC-p]='\0';
	   		pe[i]=atol(intChar);
	   		p=NC+1;
	    }

	   	if(strcmp(strandC,"+")==0)
	   		strand=0;
	    else
	    	strand=1;

	   	transcript_t tt;
	   	//tt.bin=bin;
	    tt.name=string(nameC);
	   	if(strstr(chromC,"chr"))
	   	{
	   		chr=chromC+3;
	    }
	   	else
	   		chr=chromC;
	    chrStr=string(chr);
	    tt.chr=chrStr;
            tt.strand=strand;
	    tt.txStart=txStart;
	    tt.txEnd=txEnd;
	    tt.cdsStart=cdsStart;//Dec 7, 2015
	    tt.cdsEnd=cdsEnd;
	    tt.exonCount=exonCount;
	    tt.exonStarts=ps;
	    tt.exonEnds=pe;
	    //tt.score=score;
	    tt.name2=string(name2C);

            //Dec 2015, let us not use the transcripts truncated in coding region
            int isRealAdd=1;
            if(tt.cdsStart!=tt.cdsEnd && ((tt.txStart==tt.cdsStart)||(tt.txEnd==tt.cdsEnd)))
                  isRealAdd=0;
            if(isRealAdd==1 || isTrunctOK==1)
	    {
                  transcripts.push_back(tt);
	          num++;
            }

	    if(nextN-fileContent >= length-1)
	    {
	    	cerr<<num<<" complete transcripts loaded."<<endl;
	    	break;
	    }
	    else
	    {
	    	line_buffer=nextN+1;
	    }

	}
	sort(transcripts.begin(),transcripts.end(),myTransSortFunc);
	return 0;

}

int Gene::setGene() {


	int fid=0;
	if(transcripts.size()<1)
	{
		cerr<<"No gene annotation"<<endl;
		exit(1);
	}
	else
	{
		gene_t gt;
		gt.name2=transcripts[0].name2;
		gt.strand=transcripts[0].strand;
		gt.transIds.push_back(0);
		gt.leftLimit=transcripts[0].txStart;
		gt.rightLimit=transcripts[0].txEnd;
		gt.chr=transcripts[0].chr;
		gt.fakeId=-1;
		genes.push_back(gt);
	}

	for(int i=1;i<transcripts.size();i++)
	{
		if(transcripts[i].name2.compare(transcripts[i-1].name2)!=0)
		{
			gene_t gt;
			gt.name2=transcripts[i].name2;
			gt.strand=transcripts[i].strand;
			gt.transIds.push_back(i);
			gt.leftLimit=transcripts[i].txStart;
			gt.rightLimit=transcripts[i].txEnd;
			gt.chr=transcripts[i].chr;
			gt.fakeId=-1;
			genes.push_back(gt);
		}
		else
		{

			if(transcripts[i].chr!=transcripts[i-1].chr)
			{
				gene_t gt;
				gt.name2=transcripts[i].name2;
				gt.strand=transcripts[i].strand;
				gt.transIds.push_back(i);
				gt.leftLimit=transcripts[i].txStart;
				gt.rightLimit=transcripts[i].txEnd;
				gt.chr=transcripts[i].chr;
				if(genes[genes.size()-1].fakeId!=-1)
					gt.fakeId=fid++;
				else
				{
					genes[genes.size()-1].fakeId=fid++;
					gt.fakeId=fid++;
				}
				genes.push_back(gt);
				continue;
			}
			else
			{
				if(transcripts[i].txStart > genes[genes.size()-1].rightLimit)
				{
					gene_t gt;
					gt.name2=transcripts[i].name2;
					gt.strand=transcripts[i].strand;
					gt.transIds.push_back(i);
					gt.leftLimit=transcripts[i].txStart;
					gt.rightLimit=transcripts[i].txEnd;
					gt.chr=transcripts[i].chr;
					if(genes[genes.size()-1].fakeId!=-1)
						gt.fakeId=fid++;
					else
					{
						genes[genes.size()-1].fakeId=fid++;
						gt.fakeId=fid++;
					}
					genes.push_back(gt);
					continue;
				}
			}

			//if(transcripts[i].txStart<genes[genes.size()-1].leftLimit)//this is not possible now
			//{
			//	genes[genes.size()-1].leftLimit=transcripts[i].txStart;
			//}
			if(transcripts[i].txEnd>genes[genes.size()-1].rightLimit)
			{
				genes[genes.size()-1].rightLimit=transcripts[i].txEnd;
			}
			genes[genes.size()-1].transIds.push_back(i);
		}
	}


	sort(genes.begin(),genes.end(),myGeneSortFunc);
	return 0;
}

int Gene::isInGene(string chr, uint32_t pos, vector<int>& geneIds) {
	gene_t dumbGene;
	dumbGene.chr=chr;
	dumbGene.leftLimit=pos;

	vector<gene_t>::iterator up=upper_bound(genes.begin(),genes.end(),dumbGene,myGeneSortFunc);
	if (up-genes.begin()==0)
		return 0;
	up--;
	while(up-genes.begin()>=0 && (*up).chr.compare(chr)==0 && pos > (*up).leftLimit)
	{
		if((*up).leftLimit < pos && (*up).rightLimit > pos)
		{
			geneIds.push_back(up-genes.begin());
		}
		up--;
	}
	if(geneIds.size()>0)
	{
		return 1;
	}
	else
		return 0;

}


gene_t * Gene::getGene(int index) {
	return &(genes[index]);
}


int Gene::getStrand(int geneId) {
	return genes[geneId].strand;
}

string Gene::getName2(int geneId) {
        return genes[geneId].name2;
}


uint32_t Gene::getLimitLeft(int geneId) {
	return genes[geneId].leftLimit;
}

uint32_t Gene::getLimitRight(int geneId) {
	return genes[geneId].rightLimit;
}

int Gene::getIndex(string name, vector<int> & ids) {
	for(int i=0;i<genes.size();i++)
	{
		if(genes[i].name2.compare(name)==0)
		{
			ids.push_back(i);
		}
	}
	return 0;
}