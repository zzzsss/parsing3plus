/*
 * Scores.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: zzs
 */

#include "Scores.h"

#include "Eisner.h"
#include "EisnerO2sib.h"
#include "EisnerO3g.h"


	const REAL_SCORES DEFAULT_SCORE = -10000000.0;
	//get the Scores
	Scores<REAL_SCORES>* get_the_scores(nn_input* input,REAL_SCORES* tab,int numc,int numl){
		WARNING_UNEQUAL(input->num_inst,numl,"Wrong numl.");	//numl*numc -- tab
		Scores<REAL_SCORES>* ret = new Scores<REAL_SCORES>(numc,numl,tab,DEFAULT_SCORE,input);
		int length = input->wordl->size();
		int numw = input->num_width;
		int the_index[4] = {0,0,0,0};	//h,m,s,g
		for(long i=0;i<numl*numw;i+=numw){
			for(int k=0;k<numw;k++)
				the_index[k] = input->inputs->at(i+k);
			long the_key = 0;
			int tmph = the_index[0];
			int tmpm = the_index[1];
			int tmps = the_index[2];
			int tmpg = the_index[3];
			if(tmps<0)
				tmps = tmph;
			if(tmpg<0)
				tmpg = 0;
			switch(numw){
			case 2:
				the_key = get_index2(length,tmph,tmpm);
				break;
			case 3:
				the_key = get_index2_o2sib(length,tmph,tmps,tmpm);
				break;
			case 4:
				the_key = get_index2_o3g(length,tmpg,tmph,tmps,tmpm);
				break;
			default:
				WARNING(numw,"Unknown order.");
				break;
			}
			ret->add_associate(the_key,i);
		}
		return ret;
	}

	//combine scores
	void combine_scores_o3g(Scores<REAL_SCORES>* so3,const Scores<REAL_SCORES>* so2,const Scores<REAL_SCORES>* so1)
	{
		if(!so1 && !so2)
			return;
		//combine some scores
		nn_input* the_inputs = static_cast<nn_input*>(so3->get_info());
		int length = the_inputs->wordl->size();
		int numw = the_inputs->num_width;
		for(int i=0;i<the_inputs->num_inst*numw;i+=numw){
			int tmph = the_inputs->inputs->at(i);
			int tmpm = the_inputs->inputs->at(i+1);
			int tmps = the_inputs->inputs->at(i+2);
			int tmpg = the_inputs->inputs->at(i+3);
			if(tmps<0)
				tmps = tmph;
			if(tmpg<0)
				tmpg = 0;
			long index_o3g = get_index2_o3g(length,tmpg,tmph,tmps,tmpm);
			long index_o2sib  = get_index2_o2sib(length,tmph,tmps,tmpm);
			long index_o1  = get_index2(length,tmph,tmpm);
			REAL_SCORES* tab_o3g = (*so3)[index_o3g];
			const REAL_SCORES* tab_o2sib = (so2) ? (*so2)[index_o2sib] : 0;
			const REAL_SCORES* tab_o1 = (so1) ? (*so1)[index_o1] : 0;
			for(int ll=0;ll<so3->get_numc();ll++){	//WARNING: tab_o2sib or tab_o1 might be zero, but maybe this could not happen
				tab_o3g[ll] += (so2) ? tab_o2sib[ll] : 0;
				tab_o3g[ll] += (so1) ? tab_o1[ll] : 0;
			}
		}
	}
	void combine_scores_o2sib(Scores<REAL_SCORES>* so2,const Scores<REAL_SCORES>* so1)
	{
		if(!so1)
			return;
		//combine some scores
		nn_input* the_inputs = static_cast<nn_input*>(so2->get_info());
		int length = the_inputs->wordl->size();
		int numw = the_inputs->num_width;
		for(int i=0;i<the_inputs->num_inst*numw;i+=numw){
			int tmph = the_inputs->inputs->at(i);
			int tmpm = the_inputs->inputs->at(i+1);
			int tmps = the_inputs->inputs->at(i+2);
			if(tmps<0)
				tmps = tmph;
			long index_o2sib  = get_index2_o2sib(length,tmph,tmps,tmpm);
			long index_o1  = get_index2(length,tmph,tmpm);
			REAL_SCORES* tab_o2sib = (*so2)[index_o2sib];
			const REAL_SCORES* tab_o1 = (*so1)[index_o1];
			for(int ll=0;ll<so2->get_numc();ll++){
				tab_o2sib[ll] += tab_o1[ll];
			}
		}
	}


