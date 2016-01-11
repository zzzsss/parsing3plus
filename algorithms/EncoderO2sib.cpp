/*
 * EncoderO2sib.cpp
 *
 *  Created on: 2015.7.7
 *      Author: zzs
 */

// --- from MaxParser->DependendcyEncoder2OSibling.cpp

#include "EisnerO2sib.h"
#include "Helper.h"

inline static int getKey(int s, int t, int dir, int comp, int length){
	int key = s;
	key = key * length + t;
	key = key * 2 + dir;
	key = key * 3 + comp;
	return key;
}

//return z
static REAL_SCORES calc_inside(const int length, REAL_SCORES *beta,OneScores<REAL_SCORES>& probs)
{
	int key, key1, key2;

	for(int i = 0; i < length; i++){
		key = getKey(i, i, 0, 1, length);
		beta[key] = 0.0;
		key = getKey(i, i, 1, 1, length);
		beta[key] = 0.0;
	}

	for(int j = 1; j < length; j++){
		for(int s = 0; s + j < length; s++){
			int t = s + j;
			//double prodProb_st = probs[s][t][0];
			//double prodProb_ts = probs[s][t][1];

			//init beta
			//incomplete spans
			//r == s
			int key_st_0 = getKey(s, t, 0, 0, length);
			//double prodProb_sst = probs_trips[s][s][t] + probs_sibs[s][t][0] + prodProb_st;
			REAL_SCORES prodProb_sst = probs[get_index2_o2sib(length,s,s,t)];
			key1 = getKey(s, s, 0, 1, length);
			key2 = getKey(s + 1, t, 1, 1, length);
			beta[key_st_0] = logsumexp(beta[key_st_0], beta[key1] + beta[key2] + prodProb_sst, true);

			//r == t
			int key_ts_0 = getKey(s, t, 1, 0, length);
			//double prodProb_tts = probs_trips[t][t][s] + probs_sibs[t][s][0] + prodProb_ts;
			REAL_SCORES prodProb_tts = probs[get_index2_o2sib(length,t,t,s)];
			key1 = getKey(s, t - 1, 0, 1, length);
			key2 = getKey(t, t, 1, 1, length);
			beta[key_ts_0] = logsumexp(beta[key_ts_0], beta[key1] + beta[key2] + prodProb_tts, true);

			//sibling spans
			int key_st_2 = getKey(s, t, 0, 2, length);
			beta[key_st_2] = 0.0;
			int key_ts_2 = getKey(s, t, 1, 2, length);
			beta[key_ts_2] = 0.0;
			bool flg_st_2 = true, flg_ts_2 = true;

			//complete spans
			int key_st_1 = getKey(s, t, 0, 1, length);
			beta[key_st_1] = 0.0;
			int key_ts_1 = getKey(s, t, 1, 1, length);
			beta[key_ts_1] = 0.0;
			bool flg_st_1 = true, flg_ts_1 = true;

			//calc sibling spans
			for(int r = s; r < t; r++){
				key1 = getKey(s, r, 0 ,1, length);
				key2 = getKey(r + 1, t, 1, 1, length);

				beta[key_st_2] = logsumexp(beta[key_st_2], beta[key1] + beta[key2], flg_st_2);
				flg_st_2 = false;

				beta[key_ts_2] = logsumexp(beta[key_ts_2], beta[key1] + beta[key2], flg_ts_2);
				flg_ts_2 = false;
			}

			//calc incomplete spans
			for(int r = s + 1; r < t; r++){
				key1 = getKey(s, r, 0, 0, length);
				key2 = getKey(r, t, 0, 2, length);
				//double prodProb_srt = probs_trips[s][r][t] + probs_sibs[r][t][1] + prodProb_st;
				REAL_SCORES prodProb_srt = probs[get_index2_o2sib(length,s,r,t)];
				beta[key_st_0] = logsumexp(beta[key_st_0], beta[key1] + beta[key2] + prodProb_srt, false);

				key1 = getKey(s, r, 1, 2, length);
				key2 = getKey(r, t, 1, 0, length);
				//double prodProb_trs = probs_trips[t][r][s] + probs_sibs[r][s][1] + prodProb_ts;
				REAL_SCORES prodProb_trs = probs[get_index2_o2sib(length,t,r,s)];
				beta[key_ts_0] = logsumexp(beta[key_ts_0], beta[key1] + beta[key2] + prodProb_trs, false);
			}

			//calc complete spans
			for(int r = s; r <= t; r++){
				if(r != s){
					key1 = getKey(s, r, 0, 0, length);
					key2 = getKey(r, t, 0, 1, length);
					beta[key_st_1] = logsumexp(beta[key_st_1], beta[key1] + beta[key2], flg_st_1);
					flg_st_1 = false;
				}
				if(r != t){
					key1 = getKey(s, r, 1, 1, length);
					key2 = getKey(r, t, 1, 0, length);
					beta[key_ts_1] = logsumexp(beta[key_ts_1], beta[key1] + beta[key2], flg_ts_1);
					flg_ts_1 = false;
				}
			}
		}
	}

	key1 = getKey(0, length - 1, 0, 1, length);
	key2 = getKey(0, length - 1, 1, 1, length);
	return logsumexp(beta[key1], beta[key2], false);
}

static void calc_outside(const int length,const REAL_SCORES *beta,OneScores<REAL_SCORES>& probs,REAL_SCORES *alpha)
{
	int key;
	int end = length - 1;
	for(int d = 0; d < 2; d++){
		for(int c = 0 ; c < 3; c++){
			key = getKey(0, end, d, c, length);
			alpha[key] = 0.0;
		}
	}

	for(int j = end; j >= 1; j--){
		for(int s = 0; s + j < length; s++){
			int t = s + j;

			int key_a, key_b;

			//init alpha
			//sibling spans
			int key_st_2 = getKey(s, t, 0, 2, length);
			alpha[key_st_2] = 0.0;
			bool flg_st_2 = true;
			for(int r = 0; r < s; r++){
				//double prodProb_rst = probs_trips[r][s][t] + probs_sibs[s][t][1] + probs[r][t][0];
				REAL_SCORES prodProb_rst = probs[get_index2_o2sib(length,r,s,t)];
				key_b = getKey(r, s, 0, 0, length);
				key_a = getKey(r, t, 0, 0, length);
				alpha[key_st_2] = logsumexp(alpha[key_st_2], beta[key_b] + alpha[key_a] + prodProb_rst, flg_st_2);
				flg_st_2 = false;
			}
			for(int r = t + 1; r < length; r++){
				//double prodProb_rts = probs_trips[r][t][s] + probs_sibs[t][s][1] + probs[s][r][1];
				REAL_SCORES prodProb_rts = probs[get_index2_o2sib(length,r,t,s)];
				key_b = getKey(t, r, 1, 0, length);
				key_a = getKey(s, r, 1, 0, length);
				alpha[key_st_2] = logsumexp(alpha[key_st_2], beta[key_b] + alpha[key_a] + prodProb_rts, flg_st_2);
				flg_st_2 = false;
			}

			//complete spnas
			int key_st_1 = getKey(s, t, 0, 1, length);
			bool flg_st_1 = true;
			alpha[key_st_1] = 0.0;
			if(t + 1 < length){
				key_a = getKey(s, t + 1, 1, 0, length);
				//double prodProb = probs_trips[t + 1][t + 1][s] + probs_sibs[t + 1][s][0] + probs[s][t + 1][1];
				REAL_SCORES prodProb = probs[get_index2_o2sib(length,t+1,t+1,s)];
				alpha[key_st_1] = logsumexp(alpha[key_st_1], alpha[key_a] + prodProb, flg_st_1);
				flg_st_1 = false;
			}

			int key_ts_1 = getKey(s, t, 1, 1, length);
			bool flg_ts_1 = true;
			alpha[key_ts_1] = 0.0;
			if(s != 0){
				key_a = getKey(s - 1, t, 0, 0, length);
				//double prodProb = probs_trips[s - 1][s - 1][t] + probs_sibs[s - 1][t][0] + probs[s - 1][t][0];
				REAL_SCORES prodProb = probs[get_index2_o2sib(length,s-1,s-1,t)];
				alpha[key_ts_1] = logsumexp(alpha[key_ts_1], alpha[key_a] + prodProb, flg_ts_1);
				flg_ts_1 = false;
			}

			for(int r = 0; r < s; r++){
				key_b = getKey(r, s, 0, 0, length);
				key_a = getKey(r, t, 0, 1, length);
				alpha[key_st_1] = logsumexp(alpha[key_st_1], beta[key_b] + alpha[key_a], flg_st_1);
				flg_st_1 = false;

				if(!((r == 0) && (t == length -1))){
					key_b = getKey(r, s - 1, 0 ,1, length);
					key_a = getKey(r, t, 0, 2, length);
					alpha[key_ts_1] = logsumexp(alpha[key_ts_1], beta[key_b] + alpha[key_a], flg_ts_1);
					flg_ts_1 = false;
				}
			}
			for(int r = t + 1; r < length; r++){
				if(!((s == 0) && (r == length -1))){
					key_b = getKey(t + 1, r, 1, 1, length);
					key_a = getKey(s, r, 0 ,2, length);
					alpha[key_st_1] = logsumexp(alpha[key_st_1], beta[key_b] + alpha[key_a], flg_st_1);
					flg_st_1 = false;
				}

				key_b = getKey(t, r, 1, 0, length);
				key_a = getKey(s, r, 1, 1, length);
				alpha[key_ts_1] = logsumexp(alpha[key_ts_1], beta[key_b] + alpha[key_a], flg_ts_1);
				flg_ts_1 = false;
			}

			//incomplete spans
			int key_st_0 = getKey(s, t, 0, 0, length);
			alpha[key_st_0] = 0.0;
			bool flg_st_0 = true;

			int key_ts_0 = getKey(s, t, 1, 0, length);
			alpha[key_ts_0] = 0.0;
			bool flg_ts_0 = true;

			for(int r = t; r < length; r++){
				key_b = getKey(t, r, 0 ,1, length);
				key_a = getKey(s, r, 0 ,1, length);
				alpha[key_st_0] = logsumexp(alpha[key_st_0], beta[key_b] + alpha[key_a], flg_st_0);
				flg_st_0 = false;

				if(r != t){
					key_b = getKey(t, r, 0, 2, length);
					key_a = getKey(s, r, 0, 0, length);
					//double prodProb_str = probs_trips[s][t][r] + probs_sibs[t][r][1] + probs[s][r][0];
					REAL_SCORES prodProb_str = probs[get_index2_o2sib(length,s,t,r)];
					alpha[key_st_0] = logsumexp(alpha[key_st_0], beta[key_b] + alpha[key_a] + prodProb_str, flg_st_0);
					flg_st_0 = false;
				}
			}

			for(int r = 0; r <= s; r++){
				key_b = getKey(r, s, 1, 1, length);
				key_a = getKey(r, t, 1, 1, length);
				alpha[key_ts_0] = logsumexp(alpha[key_ts_0], beta[key_b] + alpha[key_a], flg_ts_0);
				flg_ts_0 = false;

				if(r != s){
					key_b = getKey(r, s, 0, 2, length);
					key_a = getKey(r, t, 1, 0, length);
					//double prodProb_tsr = probs_trips[t][s][r] + probs_sibs[s][r][1] + probs[r][t][1];
					REAL_SCORES prodProb_tsr = probs[get_index2_o2sib(length,t,s,r)];
					alpha[key_ts_0] = logsumexp(alpha[key_ts_0], beta[key_b] + alpha[key_a] + prodProb_tsr, flg_ts_0);
					flg_ts_0 = false;
				}
			}
		}
	}
}

//##MAGIC NUMBERS of MAX-Encoder##
//direction: 0,st(right);1,ts(left)
//spans: 0,incomplete;1,complete;2,sibling

REAL_SCORES* LencodeMarginals_o2sib(const long length,Scores<REAL_SCORES>& scores)
{
	int ln = scores.get_numc();
	REAL_SCORES* marginals = new REAL_SCORES[scores.get_numl()*ln];	//the return table
	REAL_SCORES *beta = new REAL_SCORES[length * length * 2 * 3];
	REAL_SCORES *alpha = new REAL_SCORES[length * length * 2 * 3];
	//sumlabel score
	OneScores<REAL_SCORES> sum_scores = scores.get_scores_marginal();
	REAL_SCORES z = calc_inside(length, beta,sum_scores);
	calc_outside(length,beta,sum_scores,alpha);

	//get them
	long all_num = 0;
	for(int s=0;s<length;s++){
		for(int t=s+1;t<length;t++){
			//sst
			{
				long index_curr = get_index2_o2sib(length,s,s,t);
				if(scores.has_value(index_curr)){
					all_num ++;
					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(beta[getKey(s+1,t,1,1,length)]+alpha[getKey(s,t,0,0,length)]+from_assign[zl]-z);
				}
			}
			for(int r=s+1;r<t;r++){
				//srt
				{
					long index_curr = get_index2_o2sib(length,s,r,t);
					if(scores.has_value(index_curr)){
						all_num ++;
						REAL_SCORES* from_assign = scores[index_curr];
						REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
						for(int zl=0;zl<ln;zl++)
							to_assign[zl] = exp(beta[getKey(s,r,0,0,length)]+beta[getKey(r,t,0,2,length)]
													+alpha[getKey(s,t,0,0,length)]+from_assign[zl]-z);
					}
				}
			}
			//tts
			{
				long index_curr = get_index2_o2sib(length,t,t,s);
				if(scores.has_value(index_curr)){
					all_num ++;
					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(beta[getKey(s,t-1,0,1,length)]+alpha[getKey(s,t,1,0,length)]+from_assign[zl]-z);
				}
			}
			for(int r=s+1;r<t;r++){
				//trs
				{
					long index_curr = get_index2_o2sib(length,t,r,s);
					if(scores.has_value(index_curr)){
						all_num ++;
						REAL_SCORES* from_assign = scores[index_curr];
						REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
						for(int zl=0;zl<ln;zl++)
							to_assign[zl] = exp(beta[getKey(r,t,1,0,length)]+beta[getKey(s,r,0,2,length)]
												+alpha[getKey(s,t,1,0,length)]+from_assign[zl]-z);
					}
				}
			}
		}
	}

	delete []beta;
	delete []alpha;
	WARNING_UNEQUAL(all_num,scores.get_numl(),"Unequal instances.");
	return marginals;
}


