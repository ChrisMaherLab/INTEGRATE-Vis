#include "Bed2SetPair.h"


int Bed2SetPair::getSetPair(Bedpe & bp, vector<set_pair_t> & sp, Gene & g)
{
    for(int i=0;i<bp.size();i++)
    {
        bedpe_t bpt=bp.getBedpe(i);
        string chr1=bpt.chr1;
        string chr2=bpt.chr2;
        uint32_t pos5p,pos3p;

        bp.getPos(bpt,pos5p,pos3p);

        set_pair_t st;
 
        g.isInGene(chr1, pos5p, st.ids1);
        g.isInGene(chr2, pos3p, st.ids2);

        sp.push_back(st);     

    }
    return 0;    
}

