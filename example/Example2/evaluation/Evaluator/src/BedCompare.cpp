#include "BedCompare.h"

bool isStrandKnown(string strand)
{
    if(strand.compare("+")==0 || strand.compare("-")==0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


int getRange(uint32_t tt_pos5, uint32_t tt_pos3, int diff, uint32_t & tt_pos5_l, uint32_t & tt_pos5_r, uint32_t & tt_pos3_l, uint32_t & tt_pos3_r)
{
    if(tt_pos5>diff)
        tt_pos5_l=tt_pos5-diff;
    else
        tt_pos5=1;
    if(tt_pos3>diff)
        tt_pos3_l=tt_pos3-diff;
    else
        tt_pos3=1;
    tt_pos5_r=tt_pos5+diff;
    tt_pos3_r=tt_pos3+diff;
    return 0;
}

int getRange1(uint32_t tt_pos5, int max_diff, uint32_t & tt_pos5_l, uint32_t & tt_pos5_r)
{
    if(tt_pos5>max_diff)
        tt_pos5_l=tt_pos5-max_diff;
    else
        tt_pos5=1;
    tt_pos5_r=tt_pos5+max_diff;
    return 0;

}


int getNewMaxDiff(const int & max_diff, const uint32_t & pos5, uint32_t & pos3)
{
    if(pos5<1 || pos3<1)
        cerr<<"check coordinates in truth file for values < 1."<<endl;
    
    int new_max_diff=max_diff;
    if(pos5<new_max_diff+1)   
        new_max_diff=pos5-1;    
    if(pos3<new_max_diff+1)
        new_max_diff=pos3-1;
    if(new_max_diff<0)
    {
        cerr<<"new_max_diff < 0, check truth file"<<endl;
        exit(1);
    }
    return new_max_diff;
}


int  getNewTruePos(const bedpe_t & tt, const bedpe_t & rr, const uint32_t & tt_pos5, const uint32_t & tt_pos5_l, const uint32_t & tt_pos5_r, const uint32_t & rr_pos5, const uint32_t & tt_pos3, uint32_t & newPos5, uint32_t & newPos3)
{
    bool extend;
    if(tt_pos5_l<=rr_pos5 && rr_pos5 <=tt_pos5_r)
    {
        if(tt.strand1.compare("+")==0)
        {
            if(rr_pos5>=tt_pos5)
                extend=true;
            else
                extend=false;
        }
        else
        {
            if(rr_pos5<=tt_pos5)
                extend=true;
            else
                extend=false;
        }        
    }

    int move=0;
    if(tt_pos5>=rr_pos5)
        move=tt_pos5-rr_pos5;
    else
        move=rr_pos5-tt_pos5;

    newPos5=rr_pos5;

    if(tt.strand2.compare("+")==0)
    {
        if(extend==true)//at 5p
            newPos3=tt_pos3+move;
        else
            newPos3=tt_pos3-move;
    }    
    else
    {
        if(extend==true)
            newPos3=tt_pos3-move;
        else
            newPos3=tt_pos3+move;
    }
    return 0;
}


int transcript_compare(Bedpe & res, Bedpe & truth, int resolution, int max_diff, double & sensitivity, double & precision, pseudo_counts_t & pct)
{
    int diff=resolution-1;
    if(diff<0)
    {
         cerr<<"resolution must >=1 base"<<endl;
         exit(1);
    }

    if(max_diff<0)
    {
         cerr<<"max-diff must >=0 base"<<endl;
         exit(1);
    }


    vector<int> found(truth.size(),0);
    vector<int> correct(res.size(),0);

    for(int i=0;i<truth.size();i++)
    {
        for(int j=0;j<res.size();j++)
        {
             bedpe_t tt=truth.getBedpe(i);
             bedpe_t rr=res.getBedpe(j);
           
             if(tt.chr1==rr.chr1 && tt.chr2==rr.chr2 && tt.strand1==rr.strand1 && tt.strand2==rr.strand2)
             {
                 uint32_t tt_pos5, tt_pos3, rr_pos5, rr_pos3;
                 truth.getPos(tt, tt_pos5, tt_pos3);
                 truth.getPos(rr, rr_pos5, rr_pos3);             
   
                 uint32_t tt_pos5_l, tt_pos5_r, tt_pos3_l, tt_pos3_r;

                 //for homology
                 if(isStrandKnown(tt.strand1) && isStrandKnown(tt.strand2) && isStrandKnown(rr.strand1) && isStrandKnown(rr.strand2))
                 {
                     int new_max_diff=getNewMaxDiff(max_diff, tt_pos5, tt_pos3);
                     getRange1(tt_pos5,new_max_diff,tt_pos5_l,tt_pos5_r);
                     if(tt_pos5_l<=rr_pos5 && rr_pos5 <=tt_pos5_r)
                     {
                         uint32_t newPos5,newPos3;
                         getNewTruePos(tt, rr, tt_pos5, tt_pos5_l, tt_pos5_r, rr_pos5, tt_pos3, newPos5, newPos3);
                         tt_pos5=newPos5;
                         tt_pos3=newPos3;
                     }
                 }
                 //end 

                 getRange(tt_pos5, tt_pos3, diff, tt_pos5_l, tt_pos5_r, tt_pos3_l, tt_pos3_r);

                 if(tt_pos5_l<=rr_pos5 && rr_pos5 <=tt_pos5_r && tt_pos3_l<=rr_pos3 && rr_pos3 <=tt_pos3_r) 
                 {
                     found[i]=1;
                     correct[j]=1;
                 }
             }
        }
    }
    
    int fd=0;
    
    for(int i=0;i<found.size();i++)
    {
        if(found[i]==1) fd++;
    }

    int ct=0; 
    for(int j=0;j<correct.size();j++)
    {
        if(correct[j]==1)
            ct++;
    }
    
    //code checking 0 is not here
     
 
    sensitivity = (fd+pct.t_t+0.0)/(found.size()+pct.truth_t);
    
    precision = (ct+pct.t_t+0.0)/(correct.size()+pct.t_t+pct.f_t);

    return 0;
}



int BedpeCompare::compare(Bedpe & res, Bedpe & truth, Gene & g, evaluate_t & et, int resolution, int max_diff, pseudo_counts_t & pct, int print_num)
{
//trans 1. 0 0; 2. compare;
    int num_res_trans=res.size()+pct.t_t+pct.f_t;;
    int num_truth_trans=truth.size()+pct.truth_t;
    
    et.num_res_trans=num_res_trans;
    et.num_truth_trans=num_truth_trans;   
  
    if((num_res_trans==0 || num_truth_trans==0) && print_num==1)
    {
        if(num_res_trans==0 && num_truth_trans!=0)
        {
            et.sensitivity_t="0";
            et.precision_t="1";
            et.f_t="0";
        }
 
        if(num_res_trans!=0 && num_truth_trans==0)
        {
            et.sensitivity_t="1";
            et.precision_t="0";
            et.f_t="0";
        }      

        if(num_res_trans==0 && num_truth_trans==0)
        {
            et.sensitivity_t="1";
            et.precision_t="1";
            et.f_t="1";
        }
    }
    else if((num_res_trans==0 || num_truth_trans==0) &&  print_num==0)
    {
        if(num_res_trans==0 && num_truth_trans!=0)
        {
            et.sensitivity_t="0";
            et.precision_t="NA";
            et.f_t="NA";
        }

        if(num_res_trans!=0 && num_truth_trans==0)
        {
            et.sensitivity_t="NA";
            et.precision_t="0";
            et.f_t="NA";
        }

        if(num_res_trans==0 && num_truth_trans==0)
        {
            et.sensitivity_t="NA";
            et.precision_t="NA";
            et.f_t="NA";
        }
    }
    else
    {
        double sens_t, prec_t;
        transcript_compare(res, truth, resolution, max_diff, sens_t, prec_t, pct);
        double f1_t=f_score(sens_t,prec_t);
        
        et.sensitivity_t=my_db2string(sens_t);
        et.precision_t=my_db2string(prec_t);
        if(sens_t==0 && prec_t==0)
        {
            if(print_num==1)
                et.f_t="0";
            else
                et.f_t="NA";
        } 
        else            
            et.f_t=my_db2string(f1_t);
        
    }

//gene 1. id set pair; 2. 0 0; 3. compare;
//1..
    Bed2SetPair bp;
    vector<set_pair_t> resSP;
    bp.getSetPair(res, resSP, g);
    vector<set_pair_t> truthSP;
    bp.getSetPair(truth, truthSP, g);
    SetPair sp;   
    sp.merge(resSP);
    sp.merge(truthSP);       
//2..
    int num_res_gene=resSP.size()+pct.t_g+pct.f_g;
    int num_truth_gene=truthSP.size()+pct.truth_g; 
    

    et.num_res_gene=num_res_gene;
    et.num_truth_gene=num_truth_gene;

    if((num_res_gene==0 || num_truth_gene==0) && print_num==1)
    {
        if(num_res_gene==0 && num_truth_gene!=0)
        {
            et.sensitivity_g="0";
            et.precision_g="1";
            et.f_g="0";
        }

        if(num_res_gene!=0 && num_truth_gene==0)
        {
            et.sensitivity_g="1";
            et.precision_g="0";
            et.f_g="0";
        }

        if(num_res_gene==0 && num_truth_gene==0)
        {
            et.sensitivity_g="1";
            et.precision_g="1";
            et.f_g="1";
        }
    }
    else if((num_res_gene==0 || num_truth_gene==0) && print_num==0)
    {
        if(num_res_gene==0 && num_truth_gene!=0)
        {
            et.sensitivity_g="0";
            et.precision_g="NA";
            et.f_g="NA";
        }

        if(num_res_gene!=0 && num_truth_gene==0)
        {
            et.sensitivity_g="NA";
            et.precision_g="0";
            et.f_g="NA";
        }

        if(num_res_gene==0 && num_truth_gene==0)
        {
            et.sensitivity_g="NA";
            et.precision_g="NA";
            et.f_g="NA";
        }
    }
    else //3
    {
        int num_inter_gene=sp.numIntersectSet(resSP, truthSP);
        num_inter_gene+=pct.t_g;
        //evaluate;
        double sens_g=(num_inter_gene+0.0)/num_truth_gene;
        double prec_g=(num_inter_gene+0.0)/num_res_gene;
        double f1_g=f_score(sens_g,prec_g);

        et.sensitivity_g=my_db2string(sens_g);
        et.precision_g=my_db2string(prec_g);
        et.f_g=my_db2string(f1_g);
        
        if(sens_g==0 && prec_g==0)
        {
            if(print_num==1)
                et.f_g="0";
            else
                et.f_g="NA";
        }
        else
            et.f_g=my_db2string(f1_g);

    }
    return 0;
}











