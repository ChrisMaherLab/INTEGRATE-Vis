/*
 * Bed2SetPair.h
 * Mar 4 2016
 * Jin Zhang
 */

#ifndef BED2SETPAIR_H_
#define BED2SETPAIR_H_

#include "MyTypes.h"
#include "Gene2.h"
#include "Bedpe.h"

#include <iostream>
#include <algorithm>
#include <array>

using namespace std;

class Bed2SetPair
{
    public:
        Bed2SetPair(){};
        int getSetPair(Bedpe & bp, vector<set_pair_t> & sp, Gene & g); 
};


#endif
