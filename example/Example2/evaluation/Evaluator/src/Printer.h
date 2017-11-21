/*
 * Printer.h
 * Mar 6 2016
 * Jin Zhang
 */

#ifndef PRINTER_H_
#define PRINTER_H_

#include "MyTypes.h"
#include "Util.h"

#include <iostream>
#include <fstream>

using namespace std;

class Printer
{
    public:
        Printer(){};
        int print(vector<evaluate_t> & evals, char * outputfile,  string gene_file, int base_resolution, int max_diff, string pseudo_counts, string version);
};


#endif
