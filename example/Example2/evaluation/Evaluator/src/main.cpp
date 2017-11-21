/*
 * main.cpp
 *
 *  Created on: Mar 4, 2016
 *      Author: Jin Zhang
 */

#include "Bedpe.h"
#include "ExecuteByRule.h"
#include "Printer.h"
#include <iostream>
#include <getopt.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

map<int,char> intChar;
map<char,char> charChar;
map<string,char> tableAmino;

string version("0.1.4");

int usage()
{
    cerr<<endl;
    cerr<<"fusionToolEvaluator version "+version<<endl;
    cerr<<endl;
    cerr<<"Usage:"<<endl;
    cerr<<"fusionToolEvaluator --truth-file truth-file --result-file result-file --gene-annotation-file annot-file ";
    cerr<<"--subsetting-rule-file rule-file --output-file output-file"<<endl;
    cerr<<endl;
    cerr<<"Required parameters:"<<endl;
    cerr<<"    -t/--truth-file"<<endl;
    cerr<<"    -r/--result-file"<<endl;
    cerr<<"    -g/--gene-annotation-file"<<endl;
    cerr<<"    -s/--subsetting-rule-file"<<endl;
    cerr<<"    -o/--output-file"<<endl;
    cerr<<endl;
    cerr<<"Optional parameters:"<<endl;
    cerr<<"    -b/--base-resolution"<<endl;
    cerr<<"    -m/--max-diff"<<endl;
    cerr<<"    -p/--pseudo-counts"<<endl;
    cerr<<endl;
    cerr<<"Options:"<<endl;
    cerr<<"    -t/--truth-file              <string>    [ bedpe                  refer to SMC-RNA          ]"<<endl;
    cerr<<"    -r/--result-file             <string>    [ bedpe                  refer to SMC-RNA          ]"<<endl;
    cerr<<"    -g/--gene-annotation-file    <string>    [ 11 columns             refer to INTEGRATE 0.3.0  ]"<<endl;
    cerr<<"    -s/--subsetting-rule-file    <string>    [  2 columns             refer to SMC-RNA          ]"<<endl;
    cerr<<"    -o/--output-file             <string>    "<<endl;
    cerr<<"    -a/--all-transcrpt           <No value>  [ turn on all transcripts                          ]"<<endl;
    cerr<<"    -b/--base-resolution         <int>       [ default: 1                                       ]"<<endl;
    cerr<<"    -m/--max-diff                <int>       [ default: 20, tolerance for local homology        ]"<<endl;
    cerr<<"    -u/--use-number              <No value>  [ turn on theoretical metric for 0 denominator     ]"<<endl; 
    cerr<<"    -p/--pseudo-counts           <string>    [ default: 0,0,0,0,0,0                             ]"<<endl;
    cerr<<"                                             [ in order of tp_t,fp_t,truth_t,tp_g,fp_g,truth_g  ]"<<endl;
    cerr<<endl;
    return 0;
}


int main(int argc, char * argv[])
{

    int c;
    int option_index = 0;

    string truth_file="";
    string result_file="";
    string gene_file="";
    string output_file="";
    string rule_file="";
    int base_resolution=1;
    int max_diff=20;
    string pseudo_counts="0,0,0,0,0,0";
    int print_num=0;
    int isTruncOK=0;    

    static struct option long_options[] = {
        {"truth-file",            required_argument, 0,  't' },
        {"result-file",           required_argument, 0,  'r' },
        {"gene-annotation-file",  required_argument, 0,  'g' },
        {"output-file",           required_argument, 0,  'o' },
        {"subsetting-rule-file",  required_argument, 0,  's' },
        {"base-resolution",       required_argument, 0,  'b' },
        {"max-diff",              required_argument, 0,  'm' },
        {"pseudo-counts",         required_argument, 0,  'p' },
        {"all-transcript",        no_argument,       0,  'a' },
        {"use-number",            no_argument,       0,  'u' },
        {"help",                  no_argument,       0,  'h' },
        {0, 0, 0, 0}
    };

    while(1)
    {
        c = getopt_long(argc, argv, "t:r:g:o:s:b:m:p:auh",
                 long_options, &option_index);
        if (c==-1)
        {
            break;
        }
        switch(c)
        {
            case 'h':
                usage();
                exit(0);
            case 't':
                truth_file=optarg;
                break;
            case 'r':
                result_file=optarg;
                break;
            case 'g':
                gene_file=optarg;
                break;
            case 'o':
                output_file=optarg;
                break;
            case 's':
                rule_file=optarg;
                break;
            case 'b':
                base_resolution=atoi(optarg);
                break;
            case 'm':
                max_diff=atoi(optarg);
                break;
            case 'p':
                pseudo_counts=optarg;
                break;
            case 'a':
                isTruncOK=1;
                break;
            case 'u':
                print_num=1;
                break;
            default:
                break;
        }

    }


    if(truth_file.compare("")==0 || result_file.compare("")==0 || gene_file.compare("")==0 || output_file.compare("")==0|| rule_file.compare("")==0)
    {
        usage();
        exit(0);
    }

    pseudo_counts_t pct;
    get_pseudo_counts(pseudo_counts, pct);

    cerr<<"loading gene annotation file..."<<endl;    

    Gene g;
    g.loadGenesFromFile((char*)gene_file.c_str(),isTruncOK);
    g.setGene();    

    cerr<<"loading result file..."<<endl;

    Bedpe res;
    res.loadFromFile((char *)result_file.c_str());

    cerr<<"removing duplicate results..."<<endl;
    res.uniq();
    
    cerr<<"loading truth file..."<<endl;

    Bedpe truth;
    truth.loadFromFile((char*)truth_file.c_str());

    vector<evaluate_t> evaluates;

    cerr<<"comparing by subsetting rules..."<<endl;

    ExecuteByRule exe;
    exe.execute(res, truth, (char*)output_file.c_str(), (char*)rule_file.c_str(), evaluates, g, base_resolution, max_diff, pct, print_num);

    cerr<<"printing evaluation results..."<<endl;
   
    Printer per;
    per.print(evaluates, (char*)output_file.c_str(), gene_file, base_resolution, max_diff, pseudo_counts, version);

    return 0;
}