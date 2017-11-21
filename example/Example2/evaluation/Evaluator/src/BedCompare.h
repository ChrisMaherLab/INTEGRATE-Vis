/*
 * BedpeCompare.h
 * Mar 4 2016
 * Jin Zhang
 */

#ifndef BEDPECOMPARE_H_
#define BEDPECOMPARE_H_

#include "MyTypes.h"
#include "Bedpe.h"
#include "Gene2.h"
#include "Bed2SetPair.h"
#include "SetPair.h"
#include "Util.h"

#include <iostream>
#include <vector>

using namespace std;

class BedpeCompare
{
    public:
        BedpeCompare(){};
        int compare(Bedpe & res, Bedpe & truth, Gene & g, evaluate_t & et, int resolution, int max_diff, pseudo_counts_t & pct, int print_num);

};


#endif
