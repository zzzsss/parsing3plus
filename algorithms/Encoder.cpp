/*
 * Encoder.cpp
 *
 *  Created on:
 *      Author: zzs
 */

//inside-outside algorithm for order-1 model
// --- from MaxParser->DependendcyEncoder.cpp

#include "Eisner.h"
#include "Helper.h"

//not the same as get_index, but does not matter
static inline long getKey(int s, int t, int dir, int comp, long length){
	long key = s;
	key = key * length + t;
	key = key * 2 + dir;
	key = key * 2 + comp;
	return key;
}

//return z
static REAL_SCORES calc_inside(const int length, REAL_SCORES *beta,OneScores<REAL_SCORES>& probs)
{
	int key;
	for(int i = 0; i < length; i++){
		key = getKey(i, i, 0, 1, length);
		beta[key] = 0.0;
		key = getKey(i, i, 1, 1, length);
		beta[key] = 0.0;
	}
	for(int j = 1; j < length; j++){
		for(int s = 0; s + j < length; s++){
			int t = s + j;
			//double prodProb_st = probs[s][t][0];probs[get_index2(length,s,t)];
			//double prodProb_ts = probs[s][t][1];probs[get_index2(length,t,s)];
			REAL_SCORES prodProb_st = probs[get_index2(length,s,t)];
			REAL_SCORES prodProb_ts = probs[get_index2(length,t,s)];
			//init beta
			//incomplete spans
			int key_st_0 = getKey(s, t, 0, 0, length);
			beta[key_st_0] = 0.0;
			int key_ts_0 = getKey(s, t, 1, 0, length);
			beta[key_ts_0] = 0.0;
			//complete spans
			int key_st_1 = getKey(s, t, 0, 1, length);
			beta[key_st_1] = 0.0;
			int key_ts_1 = getKey(s, t, 1, 1, length);
			beta[key_ts_1] = 0.0;
			bool flg_st_0 = true, flg_ts_0 = true;
			bool flg_st_1 = true, flg_ts_1 = true;
			for(int r = s; r < t; r++){
				// first is direction, second is complete
				// _s means s is the parent
				int key1 = getKey(s, r, 0, 1, length);
				int key2 = getKey(r + 1, t, 1, 1, length);

				beta[key_st_0] = logsumexp(beta[key_st_0], beta[key1] + beta[key2] + prodProb_st, flg_st_0);
				flg_st_0 = false;

				beta[key_ts_0] = logsumexp(beta[key_ts_0], beta[key1] + beta[key2] + prodProb_ts, flg_ts_0);
				flg_ts_0 = false;
			}
			for(int r = s; r <= t; r++){
				if(r != s){
					int key1 = getKey(s, r, 0, 0, length);
					int key2 = getKey(r, t, 0, 1, length);
					beta[key_st_1] = logsumexp(beta[key_st_1], beta[key1] + beta[key2], flg_st_1);
					flg_st_1 = false;
				}
				if(r != t){
					int key1 = getKey(s, r, 1, 1, length);
					int key2 = getKey(r, t, 1, 0, length);
					beta[key_ts_1] = logsumexp(beta[key_ts_1], beta[key1] + beta[key2], flg_ts_1);
					flg_ts_1 = false;
				}
			}
		}
	}
	int key1 = getKey(0, length - 1, 0, 1, length);
	int key2 = getKey(0, length - 1, 1, 1, length);
	return logsumexp(beta[key1], beta[key2], false);
}

static void calc_outside(const int length,const REAL_SCORES *beta,OneScores<REAL_SCORES>& probs,REAL_SCORES *alpha)
{
	int key;
	int end = length - 1;
	for(int d = 0; d < 2; d++){
		for(int c = 0 ; c < 2; c++){
			key = getKey(0, end, d, c, length);
			alpha[key] = 0.0;
		}
	}
	for(int j = end; j >= 1; j--){
		for(int s = 0; s + j < length; s++){
			int t = s + j;
			//init alpha
			//incomplete spans
			int key_st_0 = getKey(s, t, 0, 0, length);
			alpha[key_st_0] = 0.0;
			int key_ts_0 = getKey(s, t, 1, 0, length);
			alpha[key_ts_0] = 0.0;
			//complete spans
			int key_st_1 = getKey(s, t, 0, 1, length);
			alpha[key_st_1] = 0.0;
			int key_ts_1 = getKey(s, t, 1, 1, length);
			alpha[key_ts_1] = 0.0;
			bool flg_st_0 = true, flg_ts_0 = true;
			bool flg_st_1 = true, flg_ts_1 = true;
			for(int r = 0; r < s; r++){
				//double prodProb_rt = probs[r][t][0];
				//double prodProb_tr = probs[r][t][1];
				REAL_SCORES prodProb_rt = probs[get_index2(length,r,t)];
				REAL_SCORES prodProb_tr = probs[get_index2(length,t,r)];
				//alpha[s][t][0][1]
				int key_b = getKey(r, s, 0, 0, length);
				int key_a = getKey(r, t, 0, 1, length);
				alpha[key_st_1] = logsumexp(alpha[key_st_1], beta[key_b] + alpha[key_a], flg_st_1);
				flg_st_1 = false;
				//alpha[s][t][1][1]
				key_b = getKey(r, s - 1, 0, 1, length);
				key_a = getKey(r, t, 0, 0, length);
				alpha[key_ts_1] = logsumexp(alpha[key_ts_1], beta[key_b] + alpha[key_a] + prodProb_rt, flg_ts_1);
				flg_ts_1 = false;
				key_a = getKey(r, t, 1, 0, length);
				alpha[key_ts_1] = logsumexp(alpha[key_ts_1], beta[key_b] + alpha[key_a] + prodProb_tr, flg_ts_1);
				flg_ts_1 = false;
			}
			for(int r = t + 1; r < length; r++){
				//double prodProb_sr = probs[s][r][0];
				//double prodProb_rs = probs[s][r][1];
				REAL_SCORES prodProb_sr = probs[get_index2(length,s,r)];
				REAL_SCORES prodProb_rs = probs[get_index2(length,r,s)];
				//alpha[s][t][0][1]
				int key_b = getKey(t + 1, r, 1, 1, length);
				int key_a = getKey(s, r, 0, 0, length);
				alpha[key_st_1] = logsumexp(alpha[key_st_1], beta[key_b] + alpha[key_a] + prodProb_sr, flg_st_1);
				flg_st_1 = false;
				key_a = getKey(s, r, 1, 0, length);
				alpha[key_st_1] = logsumexp(alpha[key_st_1], beta[key_b] + alpha[key_a] + prodProb_rs, flg_st_1);
				flg_st_1 = false;
				//alpha[s][t][1][1]
				key_b = getKey(t, r, 1, 0, length);
				key_a = getKey(s, r, 1, 1, length);
				alpha[key_ts_1] = logsumexp(alpha[key_ts_1], beta[key_b] + alpha[key_a], flg_ts_1);
				flg_ts_1 = false;
			}
			//alpha[s][t][0][0]
			for(int r = t; r < length; r++){
				int key_b = getKey(t, r, 0, 1, length);
				int key_a = getKey(s, r, 0, 1, length);
				alpha[key_st_0] = logsumexp(alpha[key_st_0], beta[key_b] + alpha[key_a], flg_st_0);
				flg_st_0 = false;
			}
			//alpha[s][t][1][0]
			for(int r = 0; r <= s; r++){
				int key_b = getKey(r, s, 1, 1, length);
				int key_a = getKey(r, t, 1, 1, length);
				alpha[key_ts_0] = logsumexp(alpha[key_ts_0], beta[key_b] + alpha[key_a], flg_ts_0);
				flg_ts_0 = false;
			}
		}
	}
}

REAL_SCORES* LencodeMarginals(const int length,Scores<REAL_SCORES>& scores)
{
	int ln = scores.get_numc();
	REAL_SCORES* marginals = new REAL_SCORES[scores.get_numl()*ln];	//the return table
	REAL_SCORES *beta = new REAL_SCORES[length * length * 2 * 2];
	REAL_SCORES *alpha = new REAL_SCORES[length * length * 2 * 2];
	//sumlabel score
	OneScores<REAL_SCORES> sum_scores = scores.get_scores_marginal();
	REAL_SCORES z = calc_inside(length, beta,sum_scores);
	calc_outside(length,beta,sum_scores,alpha);

	long all_num = 0;
	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			//i->j
			{
				long index_curr = get_index2(length,i,j);
				if(scores.has_value(index_curr)){
					all_num ++;
					int key_io = getKey(i, j, 0, 0, length);
					REAL_SCORES tmp_mscore = beta[key_io]+alpha[key_io]-z-sum_scores[get_index2(length,i,j)];

					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(tmp_mscore+from_assign[zl]);
				}
			}
			//j->i
			{
				long index_curr = get_index2(length,j,i);
				if(scores.has_value(index_curr)){
					all_num ++;
					int key_io = getKey(i, j, 1, 0, length);
					REAL_SCORES tmp_mscore = beta[key_io]+alpha[key_io]-z-sum_scores[get_index2(length,j,i)];

					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(tmp_mscore+from_assign[zl]);
				}
			}
		}
	}
	delete []beta;
	delete []alpha;

	WARNING_UNEQUAL(all_num,scores.get_numl(),"Unequal instances.");
	return marginals;
}

//-----------------------------
REAL_SCORES* encodeMarginals(const int length,REAL_SCORES* scores)
{
	OneScores<REAL_SCORES> onescores(length*length,0,scores,0);
	REAL_SCORES* marginals = new REAL_SCORES[length*length];	//use get_index2
	REAL_SCORES *beta = new REAL_SCORES[length * length * 2 * 2];
	REAL_SCORES *alpha = new REAL_SCORES[length * length * 2 * 2];
	REAL_SCORES z = calc_inside(length, beta,onescores);
	calc_outside(length,beta,onescores,alpha);

	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			//i->j
			int key_io = getKey(i, j, 0, 0, length);
			int key_assign = get_index2(length,i,j);
			marginals[key_assign] = exp(beta[key_io]+alpha[key_io]-z);
			//j->i
			key_io = getKey(i, j, 1, 0, length);
			key_assign = get_index2(length,j,i);
			marginals[key_assign] = exp(beta[key_io]+alpha[key_io]-z);
		}
	}
	delete []beta;
	delete []alpha;
	return marginals;
}
