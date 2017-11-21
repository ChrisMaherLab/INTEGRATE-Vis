#include "ExecuteByRule.h"

typedef struct
{
    string name;
    string rule;
} rule_t;

int ExecuteByRule::execute(Bedpe & res, Bedpe & truth, char * outputfile, char * rulefile, vector<evaluate_t> & evaluates, Gene & g, int resolution, int max_diff, pseudo_counts_t & pct, int print_num)
{
  string line;
  vector<rule_t> rtv;
  ifstream myfile (rulefile);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
        std::vector<std::string> tmp = my_split(line, '\t');
        rule_t rt;
        rt.name=tmp[0];
        rt.rule=tmp[1];
        rtv.push_back(rt);
    }
    myfile.close();
  }

  CutterByRule cbr;
  BedpeCompare bcomp;

  for(int i=0;i<rtv.size();i++)
  {
      string tmpfile1=string(outputfile)+"."+std::to_string((long long int)i)+"."+"result.bedpe";
      string tmpfile2=string(outputfile)+"."+std::to_string((long long int)i)+"."+"truth.bedpe";
      string ruleString=rtv[i].rule;
      Bedpe  outRes;
      Bedpe  outTruth;
      cbr.cut(res, outRes, ruleString, tmpfile1);
      cbr.cut(truth, outTruth, ruleString, tmpfile2);
     
      evaluate_t et; 
      et.name=rtv[i].name;
      bcomp.compare(outRes, outTruth, g, et, resolution, max_diff, pct, print_num);
      evaluates.push_back(et);
  }  


}

