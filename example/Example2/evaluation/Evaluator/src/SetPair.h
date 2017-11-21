/*
 * SetPair.h
 * Mar 4 2016
 * Jin Zhang
 */

#ifndef SETPAIR_H_
#define SETPAIR_H_

#include "MyTypes.h"

#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

class SetPair
{
    public:
       SetPair(){};
       int merge(vector<set_pair_t> & spv);
       int numIntersectSet(vector<set_pair_t> & spv1,vector<set_pair_t> & spv2); 
};


#endif
