#include "Printer.h"

int Printer::print(vector<evaluate_t> & evals, char * outputfile,  string gene_file, int base_resolution, int max_diff, string pseudo_counts, string version)
{
    ofstream ofs;
    ofs.open (outputfile, std::ofstream::out);

    ofs<<"#fusionToolEvaluator version "+version<<endl;
    ofs<<"#--gene-annotation-file="+gene_file+"; --base-resolution="+my_db2string(base_resolution)+"; --max-diff="+my_db2string(max_diff)+"; --pseudo-counts="+pseudo_counts<<endl;
    ofs<<"#Evaluation_Name\tNum_Res_Trans\tNum_Truth_Trans\tSensitivity_Trans\tPrecision_Trans\tF1_Trans\tNum_Res_Gene\tNum_Truth_Gene\tSensitivity_Gene\tPrecision_Gene\tF1_Gene"<<endl;
    for(int i=0;i<evals.size();i++)
    {
        evaluate_t et=evals[i];
        ofs<<et.name<<"\t";
        ofs<<et.num_res_trans<<"\t";
        ofs<<et.num_truth_trans<<"\t";
        ofs<<et.sensitivity_t<<"\t";
        ofs<<et.precision_t<<"\t";
        ofs<<et.f_t<<"\t";
        ofs<<et.num_res_gene<<"\t";
        ofs<<et.num_truth_gene<<"\t";
        ofs<<et.sensitivity_g<<"\t";
        ofs<<et.precision_g<<"\t";
        ofs<<et.f_g<<"\n";
    }
    ofs.close();

    return 0;
}
