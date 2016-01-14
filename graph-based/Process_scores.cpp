/*
 * Process_scores.cpp
 *
 *  Created on: Oct 11, 2015
 *      Author: zzs
 */

#include "Process.h"
#include "../algorithms/Eisner.h"
#include "../algorithms/EisnerO2sib.h"
#include "../algorithms/EisnerO3g.h"
#include <exception>

//be careful about the magic numbers
static inline void TMP_push234(vector<int>* l,int h,int m,int c=IMPOSSIBLE_INDEX,int g=IMPOSSIBLE_INDEX)
{
	l->push_back(h);
	l->push_back(m);
	if(c > IMPOSSIBLE_INDEX){
		l->push_back(c);
		if(g > IMPOSSIBLE_INDEX)
			l->push_back(g);
	}
}
static inline void TMP_pop234(vector<int>* l,int h,int m,int c=IMPOSSIBLE_INDEX,int g=IMPOSSIBLE_INDEX)
{
	l->pop_back();
	l->pop_back();
	if(c > IMPOSSIBLE_INDEX){
		l->pop_back();
		if(g > IMPOSSIBLE_INDEX)
			l->pop_back();
	}
}
static inline bool TMP_check234(int h1,int h2,int c1=IMPOSSIBLE_INDEX,int c2=IMPOSSIBLE_INDEX,
		int g1=IMPOSSIBLE_INDEX,int g2=IMPOSSIBLE_INDEX)
{
	return (h1==h2)&&(c1==c2)&&(g1==g2);
}
static inline bool TMP_higho_sample(HypherParameters* HP,int h=IMPOSSIBLE_INDEX,int m=IMPOSSIBLE_INDEX,
		int s=IMPOSSIBLE_INDEX,int g=IMPOSSIBLE_INDEX)
{
	return drand48() < HP->CONF_higho_percent;
}

// the routine for getting scores
// -- m means the current machine(might not be the same as options)
// -- t for assigning inputs, need outside to delete it
// -- testing for the forward of nn
REAL* Process::forward_scores_o1(DependencyInstance* x,Csnn* mac,nn_input** t,nn_input_helper* helper,
		int testing,HypherParameters*hh)
{
	//default order1 parsing
	int odim = mac->get_odim();	//1 or 2 for no-labeled, otherwise...
	int length = x->length();
	//- for output goals
	bool is_labeled = (mac->get_classdim()>1);
	int nope_goal = mac->get_classdim();	//for no-rel
	//prepare scores
	int num_pair_togo = 0;
	int num_good=0,num_bad=0;
	vector<int>* the_inputs = new vector<int>();
	vector<int>* the_goals = new vector<int>();
	//loop --- (h,m)
	for(int m=1;m<length;m++){
		for(int h=0;h<length;h++){
			if(m != h){
				TMP_push234(the_inputs,h,m);
				if(!testing){	//when training, prepare the goals
					if(TMP_check234(x->heads->at(m),h)){
						the_goals->push_back(is_labeled?(x->index_deprels->at(m)):0);
						num_good++;
					}
					else{
						if(TMP_higho_sample(hh)){
							the_goals->push_back(nope_goal);
							num_bad++;
						}
						else{
							TMP_pop234(the_inputs,h,m);
							num_pair_togo -= 1;	//well, maybe the same
						}
					}
				}
				num_pair_togo ++;
			}
		}
	}
	(*t) = new nn_input(num_pair_togo,2,the_inputs,the_goals,x->index_forms,x->index_pos,helper,
			num_good,num_bad);
	//--fix the 0 bug
	Process::CHECK_EQUAL(num_pair_togo,(*t)->num_inst,"Forward bsize nope.");
	if(num_pair_togo==0)
		return 0;
	//--
	REAL* tmp_scores = mac->forward(*t,testing);
	return tmp_scores;
}
REAL* Process::forward_scores_o2sib(DependencyInstance* x,Csnn* mac,nn_input** t,nn_input_helper* h,
		int testing,bool* cut_o1,HypherParameters* hh)
{
	//o2sib
	int odim = mac->get_odim();	//1 or 2 for no-labeled, otherwise...
	int length = x->length();
	//- for output goals
	bool is_labeled = (mac->get_classdim()>1);
	int nope_goal = mac->get_classdim();	//for no-rel
	//prepare scores
	int num_togo = 0;
	int num_good=0,num_bad=0;
	vector<int>* the_inputs = new vector<int>();
	vector<int>* the_goals = new vector<int>();
	bool* score_o1 = cut_o1;
	//loop --- (h,m,s)
	for(int m=1;m<length;m++){
		//put the real one first if training
		// --- maybe doubled, but maybe does not matter (when testing: same value, when training: double positive)
		int real_center=IMPOSSIBLE_INDEX;
		int real_head=IMPOSSIBLE_INDEX;
		if(!testing){	//when training, prepare the goals
			real_center = x->siblings->at(m);
			real_head = x->heads->at(m);
			//must get the real one
			TMP_push234(the_inputs,real_head,m,real_center);
			the_goals->push_back(is_labeled?(x->index_deprels->at(m)):0);
			num_togo += 1;
			num_good++;
		}
		//others
		for(int h=0;h<length;h++){
			if(h==m)
				continue;
			bool norpob_hm = score_o1[get_index2(length,h,m)];
			//h,m,-1
			if(!norpob_hm){
				TMP_push234(the_inputs,h,m,-1);
				if(!testing){
					if(TMP_check234(real_head,h,real_center,-1))
					{
						TMP_pop234(the_inputs,h,m,-1);
						num_togo -= 1;
					}//don't add again here
					else{
						if(TMP_higho_sample(hh)){
							the_goals->push_back(nope_goal);
							num_bad++;
						}
						else{
							TMP_pop234(the_inputs,h,m,-1);
							num_togo -= 1;	//well, maybe the same
						}
					}
				}
				num_togo += 1;
			}
			//h,m,c
			int small = GET_MIN_ONE(m,h);
			int large = GET_MAX_ONE(m,h);
			if(!norpob_hm){
				for(int c=small+1;c<large;c++){
					if(!score_o1[get_index2(length,h,c)]){
						TMP_push234(the_inputs,h,m,c);
						if(!testing){
							if(TMP_check234(real_head,h,real_center,c))
							{
								TMP_pop234(the_inputs,h,m,c);
								num_togo -= 1;
							}//don't add again here
							else{
								if(TMP_higho_sample(hh)){
									the_goals->push_back(nope_goal);
									num_bad++;
								}
								else{
									TMP_pop234(the_inputs,h,m,c);
									num_togo -= 1;	//well, maybe the same
								}
							}
						}
						num_togo += 1;
					}
				}
			}
		}
	}
	(*t) = new nn_input(num_togo,3,the_inputs,the_goals,x->index_forms,x->index_pos,h,
			num_good,num_bad);	//!!DEBUG, goals and inputs reversed
	//--fix the 0 bug
	Process::CHECK_EQUAL(num_togo,(*t)->num_inst,"Forward bsize nope.");
	if(num_togo==0)
		return 0;
	//--
	REAL* tmp_scores = mac->forward(*t,testing);
	return tmp_scores;
}
REAL* Process::forward_scores_o3g(DependencyInstance* x,Csnn* mac,nn_input** t,nn_input_helper* h,
		int testing,bool* cut_o1,HypherParameters*hh)
{
	//o3g
	int odim = mac->get_odim();	//1 or 2 for no-labeled, otherwise...
	int length = x->length();
	//- for output goals
	bool is_labeled = (mac->get_classdim()>1);
	int nope_goal = mac->get_classdim();	//for no-rel
	//prepare scores
	int num_togo = 0;
	int num_good=0,num_bad=0;
	vector<int>* the_inputs = new vector<int>();
	vector<int>* the_goals = new vector<int>();
	bool* score_o1 = cut_o1;
	//loop --- (h,m,s,g)
	for(int m=1;m<length;m++){
		//put the real one first if training
		int real_center=IMPOSSIBLE_INDEX;
		int real_head=IMPOSSIBLE_INDEX;
		int real_grand=IMPOSSIBLE_INDEX;
		if(!testing){	//when training, prepare the goals
			real_center = x->siblings->at(m);
			real_head = x->heads->at(m);
			real_grand = x->heads->at(real_head);
			//must get the real one
			TMP_push234(the_inputs,real_head,m,real_center,real_grand);
			the_goals->push_back(is_labeled?(x->index_deprels->at(m)):0);
			num_togo += 1;
			num_good++;
		}
		//--others
		//1.when h==0
		{
			int h=0;
			bool norpob_hm = score_o1[get_index2(length,h,m)];
			if(!norpob_hm){
				//0,0,0,m as g,h,c,m
				TMP_push234(the_inputs,h,m,-1,-1);
				if(!testing){
					if(TMP_check234(real_head,h,real_center,-1,real_grand,-1))
					{
						TMP_pop234(the_inputs,h,m,-1,-1);
						num_togo -= 1;
					}//don't add again here
					else{
						if(TMP_higho_sample(hh)){
							the_goals->push_back(nope_goal);
							num_bad++;
						}
						else{
							TMP_pop234(the_inputs,h,m,-1,-1);
							num_togo -= 1;	//well, maybe the same
						}
					}
				}
				num_togo += 1;
				//0,0,c,m
				for(int c=1;c<m;c++){
					if(!score_o1[get_index2(length,h,c)]){
						TMP_push234(the_inputs,h,m,c,-1);
						if(!testing){
							if(TMP_check234(real_head,h,real_center,c,real_grand,-1))
							{
								TMP_pop234(the_inputs,h,m,c,-1);
								num_togo -= 1;
							}
							else{
								if(TMP_higho_sample(hh)){
									the_goals->push_back(nope_goal);
									num_bad++;
								}
								else{
									TMP_pop234(the_inputs,h,m,c,-1);
									num_togo -= 1;	//well, maybe the same
								}
							}
						}
						num_togo += 1;
					}
				}
			}
		}
		//2.others h>0
		for(int h=1;h<length;h++){
			if(h==m)
				continue;
			bool norpob_hm = score_o1[get_index2(length,h,m)];
			int small = GET_MIN_ONE(h,m);
			int large = GET_MAX_ONE(h,m);
			if(!norpob_hm){
				for(int g=0;g<length;g++){
					bool norpob_gh = score_o1[get_index2(length,g,h)];
					if((g>=small && g<=large) || norpob_gh)	//non-proj not allowed //!!DEBUG: yi-jing-wu-yu-le...
						continue;
					//for those c
					//g,h,-1,m
					TMP_push234(the_inputs,h,m,-1,g);
					if(!testing){
						if(TMP_check234(real_head,h,real_center,-1,real_grand,g))
						{
							TMP_pop234(the_inputs,h,m,-1,g);
							num_togo -= 1;
						}
						else{
							if(TMP_higho_sample(hh)){
								the_goals->push_back(nope_goal);
								num_bad++;
							}
							else{
								TMP_pop234(the_inputs,h,m,-1,g);
								num_togo -= 1;	//well, maybe the same
							}
						}
					}
					num_togo += 1;
					//g,h,c,m
					for(int c=small+1;c<large;c++){
						if(!score_o1[get_index2(length,h,c)]){
							TMP_push234(the_inputs,h,m,c,g);
							if(!testing){
								if(TMP_check234(real_head,h,real_center,c,real_grand,g))
								{
									TMP_pop234(the_inputs,h,m,c,g);
									num_togo -= 1;
								}
								else{
									if(TMP_higho_sample(hh)){
										the_goals->push_back(nope_goal);
										num_bad++;
									}
									else{
										TMP_pop234(the_inputs,h,m,c,g);
										num_togo -= 1;	//well, maybe the same
									}
								}
							}
							num_togo += 1;
						}
					}
				}
			}
		}
	}
	(*t) = new nn_input(num_togo,4,the_inputs,the_goals,x->index_forms,x->index_pos,h,num_good,num_bad);
	//--fix the 0 bug
	Process::CHECK_EQUAL(num_togo,(*t)->num_inst,"Forward bsize nope.");
	if(num_togo==0)
		return 0;
	//--
	REAL* tmp_scores = mac->forward(*t,testing);
	return tmp_scores;
}

