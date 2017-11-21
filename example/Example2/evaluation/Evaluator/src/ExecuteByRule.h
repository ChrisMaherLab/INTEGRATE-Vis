/*
 * ExecuteByRule.h
 * Mar 6 2016
 * Jin Zhang
 */

#ifndef ExecuteByRule_H_
#define ExecuteByRule_H_

#include "MyTypes.h"
#include "CutterByRule.h"
#include "BedCompare.h"

#include <iostream>
#include <string>
#include <cstring>
#include <vector>


using namespace std;

class ExecuteByRule
{
    public:
        ExecuteByRule(){};
        int execute(Bedpe & res, Bedpe & truth, char * outputfile, char * rulefile, vector<evaluate_t> & evaluates, Gene & g, int resolution, int max_diff, pseudo_counts_t & pct, int print_num);
};


#endif
