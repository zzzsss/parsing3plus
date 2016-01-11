/*
 * EncoderO3g.cpp
 *
 *  Created on: 2015.7.15
 *      Author: zzs
 */

//from MaxParser.DependencyEncoder3OGSibling.cpp

#include "EisnerO3g.h"
#include "Helper.h"

inline static int getKey(int g, int s, int t, int dir, int comp, int length){
	int key = g;
	key = key * length + s;
	key = key * length + t;
	key = key * 2 + dir;
	key = key * 3 + comp;
	return key;
}

//return z
static REAL_SCORES calc_inside(const int length, REAL_SCORES *beta,OneScores<REAL_SCORES>& probs)
{
	int key, key1, key2;
	for(int s = 0; s < length; s++){
		for(int g = 0; g < length; g++){
			key = getKey(g, s, s, 0, 1, length);
			beta[key] = 0.0;
			key = getKey(g, s, s, 1, 1, length);
			beta[key] = 0.0;
		}
	}

	for(int j = 1; j < length; j++){
		for(int s = 0; s < length && s + j < length; s++){
			int t = s + j;

			//double prodProb_st = probs[s][t][0];
			//double prodProb_ts = probs[s][t][1];

			for(int g = 0; g < length; g++){
				if(g >= s && g <= t){
					continue;
				}

				//init beta
				//incomplete spans
				//r == s
				int key_gst_0 = getKey(g, s, t, 0, 0, length);
				REAL_SCORES prodProb_gsst = probs[get_index2_o3g(length,g,s,s,t)];
				key1 = getKey(g, s, s, 0, 1, length);
				key2 = getKey(s, s + 1, t, 1, 1, length);
				beta[key_gst_0] = logsumexp(beta[key_gst_0], beta[key1] + beta[key2] + prodProb_gsst, true);

				//r == t
				int key_gts_0 = getKey(g, s, t, 1, 0, length);
				REAL_SCORES prodProb_gtts = probs[get_index2_o3g(length,g,t,t,s)];
				key1 = getKey(t, s, t - 1, 0, 1, length);
				key2 = getKey(g, t, t, 1, 1, length);
				beta[key_gts_0] = logsumexp(beta[key_gts_0], beta[key1] + beta[key2] + prodProb_gtts, true);

				//sibling spans
				int key_gst_2 = getKey(g, s, t, 0, 2, length);
				beta[key_gst_2] = 0.0;
				int key_gts_2 = getKey(g, s, t, 1, 2, length);
				beta[key_gts_2] = 0.0;
				bool flg_gst_2 = true, flg_gts_2 = true;

				//complete spans
				int key_gst_1 = getKey(g, s, t, 0, 1, length);
				beta[key_gst_1] = 0.0;
				int key_gts_1 = getKey(g, s, t, 1, 1, length);
				beta[key_gts_1] = 0.0;
				bool flg_gst_1 = true, flg_gts_1 = true;

				//calc sibling spans
				for(int r = s; r < t; r++){
					key1 = getKey(g, s, r, 0, 1, length);
					key2 = getKey(g, r + 1, t, 1, 1, length);
					beta[key_gst_2] = logsumexp(beta[key_gst_2], beta[key1] + beta[key2], flg_gst_2);
					flg_gst_2 = false;

					beta[key_gts_2] = logsumexp(beta[key_gts_2], beta[key1] + beta[key2], flg_gts_2);
					flg_gts_2 = false;
				}

				//calc incomplete spans
				for(int r = s + 1; r < t; r++){
					// s -> (r,t)
					key1 = getKey(g, s, r, 0, 0, length);
					key2 = getKey(s, r, t, 0, 2, length);
					REAL_SCORES prodProb_gsrt = probs[get_index2_o3g(length,g,s,r,t)];
					beta[key_gst_0] = logsumexp(beta[key_gst_0], beta[key1] + beta[key2] + prodProb_gsrt, false);

					// t -> (r,s)
					key1 = getKey(t, s, r, 1, 2, length);
					key2 = getKey(g, r, t, 1, 0, length);
					REAL_SCORES prodProb_gtrs = probs[get_index2_o3g(length,g,t,r,s)];
					beta[key_gts_0] = logsumexp(beta[key_gts_0], beta[key1] + beta[key2] + prodProb_gtrs, false);
				}

				//calc complete spans
				for(int r = s; r <= t; r++){
					if(r != s){
						key1 = getKey(g, s, r, 0, 0, length);
						key2 = getKey(s, r, t, 0, 1, length);
						beta[key_gst_1] = logsumexp(beta[key_gst_1], beta[key1] + beta[key2], flg_gst_1);
						flg_gst_1 = false;
					}
					if(r != t){
						key1 = getKey(t, s, r, 1, 1, length);
						key2 = getKey(g, r, t, 1, 0, length);
						beta[key_gts_1] = logsumexp(beta[key_gts_1], beta[key1] + beta[key2], flg_gts_1);
						flg_gts_1 = false;
					}
				}
			}
		}
	}

	int end = length - 1;

	for(int j = 1; j < length; j++){
		int t = j, s = end - j;
		//init
		key = getKey(0, 0, t, 0, 0, length);
		REAL_SCORES prodProb_00t = probs[get_index2_o3g(length,0,0,0,t)];;
		key1 = getKey(0, 1, t, 1, 1, length);
		beta[key] = logsumexp(beta[key], beta[key1] + prodProb_00t, true);

		for(int r = 1; r < t; r++){
			REAL_SCORES prodProb_0rt = probs[get_index2_o3g(length,0,0,r,t)];;
			key1 = getKey(0, 0, r, 0, 0, length);
			key2 = getKey(0, r, t, 0, 2, length);
			beta[key] = logsumexp(beta[key], beta[key1] + beta[key2] + prodProb_0rt, false);
		}

		//init
		key = getKey(s, s, end, 1, 0, length);
		REAL_SCORES prodProb_nns = probs[get_index2_o3g(length,end,end,end,t)];;
		key1 = getKey(end, s, end - 1, 0, 1, length);
		beta[key] = logsumexp(beta[key], beta[key1] + prodProb_nns, true);

		for(int r = s + 1; r < end; r++){
			REAL_SCORES prodProb_nrs = probs[get_index2_o3g(length,end,end,r,t)];
			key1 = getKey(end, s, r, 0, 2, length);
			key2 = getKey(r, r, end, 1, 0, length);
			beta[key] = logsumexp(beta[key], beta[key1] + beta[key2] + prodProb_nrs, false);
		}
	}
	int key_0n_1 = getKey(0, 0, end, 0, 1, length);
	int key_n0_1 = getKey(0, 0, end, 1, 1, length);
	bool flg_0n_1 = true, flg_n0_1 = true;
	for(int r = 0; r < length; r++){
		if(r != 0){
			key1 = getKey(0, 0, r, 0, 0, length);
			key2 = getKey(0, r, end, 0, 1, length);
			beta[key_0n_1] = logsumexp(beta[key_0n_1], beta[key1] + beta[key2], flg_0n_1);
			flg_0n_1 = false;
		}
		if(r != end){
			key1 = getKey(end, 0, r, 1, 1, length);
			key2 = getKey(r, r, end, 1, 0, length);
			beta[key_n0_1] = logsumexp(beta[key_n0_1], beta[key1] + beta[key2], flg_n0_1);
			flg_n0_1 = false;
		}
	}
	return logsumexp(beta[key_0n_1], beta[key_n0_1], false);
}

static void calc_outside(const int length,const REAL_SCORES *beta,OneScores<REAL_SCORES>& probs,REAL_SCORES *alpha)
{
	int end = length - 1;
	int key;
	for(int c = 0; c < 3; c++){
		key = getKey(0, 0, end, 0, c, length);
		alpha[key] = 0.0;
		key = getKey(0, 0, end, 1, c, length);
		alpha[key] = 0.0;
	}

	for(int j = end; j >= 1; j--){
		int t = j, s = end - j;

		//init
		int key_t = getKey(0, 0, t, 0, 0, length);
		int key_b = getKey(0, t, end, 0, 1, length);
		int key_a = getKey(0, 0, end, 0, 1, length);
		alpha[key_t] = logsumexp(alpha[key_t], beta[key_b] + alpha[key_a], true);

		for(int r = t + 1; r < length; r++){
			REAL_SCORES prodProb_0tr = probs[get_index2_o3g(length,0,0,t,r)];
			key_b = getKey(0, t, r, 0, 2, length);
			key_a = getKey(0, 0, r, 0, 0, length);
			alpha[key_t] = logsumexp(alpha[key_t], beta[key_b] + alpha[key_a] + prodProb_0tr, false);
		}

		//init
		int key_s = getKey(s, s, end, 1, 0, length);
		key_b = getKey(end, 0, s, 1, 1, length);
		key_a = getKey(0, 0, end, 1, 1, length);
		alpha[key_s] = logsumexp(alpha[key_s], beta[key_b] + alpha[key_a], true);

		for(int r = 0; r < s; r++){
			REAL_SCORES prodProb_nsr = probs[get_index2_o3g(length,end,end,s,r)];;
			key_b = getKey(end, r, s, 0, 2, length);
			key_a = getKey(r, r, end, 1, 0, length);
			alpha[key_s] = logsumexp(alpha[key_s], beta[key_b] + alpha[key_a] + prodProb_nsr, false);
		}
	}

	for(int j = end; j >= 1; j--){
		for(int s = 0; s + j < length; s++){
			int t = s + j;
			int key_a, key_b;
			//bool coord;

			for(int g = 0; g < s; g++){
				//sibling spans
				int key_gst_2 = getKey(g, s, t, 0, 2, length);
				alpha[key_gst_2] = 0.0;
				bool flg_gst_2 = true;
				for(int r = 0; r < length; r++){
					if(r >= g && r <= t){
						continue;
					}
					REAL_SCORES prodProb_rgst = probs[get_index2_o3g(length,r,g,s,t)];
					key_b = getKey(r, g, s, 0, 0, length);
					key_a = getKey(r, g, t, 0, 0, length);
					alpha[key_gst_2] = logsumexp(alpha[key_gst_2], beta[key_b] + alpha[key_a] + prodProb_rgst, flg_gst_2);
					flg_gst_2 = false;
				}
				if(g == 0){
					REAL_SCORES prodProb_0st = probs[get_index2_o3g(length,0,0,s,t)];;
					key_b = getKey(0, 0, s, 0, 0, length);
					key_a = getKey(0, 0, t, 0, 0, length);
					alpha[key_gst_2] = logsumexp(alpha[key_gst_2], beta[key_b] + alpha[key_a] + prodProb_0st, flg_gst_2);
					flg_gst_2 = false;
				}

				//completed spans
				int key_gst_1 = getKey(g, s, t, 0, 1, length);
				alpha[key_gst_1] = 0.0;
				int key_gts_1 = getKey(g, s, t, 1, 1, length);
				alpha[key_gts_1] = 0.0;

				bool flg_gst_1 = true, flg_gts_1 = true;

				for(int r = 0; r < length; r++){
					if(r > t && r < length){
						key_b = getKey(g, t + 1, r, 1, 1, length);
						key_a = getKey(g, s, r, 0, 2, length);
						alpha[key_gst_1] = logsumexp(alpha[key_gst_1], beta[key_b] + alpha[key_a], flg_gst_1);
						flg_gst_1 = false;
					}
					if(r > g && r < s){
						key_b = getKey(g, r, s - 1, 0, 1, length);
						key_a = getKey(g, r, t, 0, 2, length);
						alpha[key_gts_1] = logsumexp(alpha[key_gts_1], beta[key_b] + alpha[key_a], flg_gts_1);
						flg_gts_1 = false;
					}
					if(r < g || r > t){
						key_b = getKey(r, g, s, 0, 0, length);
						key_a = getKey(r, g, t, 0, 1, length);
						alpha[key_gst_1] = logsumexp(alpha[key_gst_1], beta[key_b] + alpha[key_a], flg_gst_1);
						flg_gst_1 = false;

						if(g == s - 1){
							REAL_SCORES prodProb_rggt = probs[get_index2_o3g(length,r,g,g,t)];
							key_a = getKey(r, g, t, 0, 0, length);
							alpha[key_gts_1] = logsumexp(alpha[key_gts_1], alpha[key_a] + prodProb_rggt, flg_gts_1);
							flg_gts_1 = false;
						}
					}
				}
				if(g == 0 && t == end){
					key_b = getKey(0, 0, s, 0, 0, length);
					key_a = getKey(0, 0, end, 0, 1, length);
					alpha[key_gst_1] = logsumexp(alpha[key_gst_1], beta[key_b] + alpha[key_a], flg_gst_1);
					flg_gst_1 = false;
				}
				if(g == 0 && s == 1){
					REAL_SCORES prodProb_00t = probs[get_index2_o3g(length,0,0,0,t)];
					key_a = getKey(0, 0, t, 0, 0, length);
					alpha[key_gts_1] = logsumexp(alpha[key_gts_1], alpha[key_a] + prodProb_00t, flg_gts_1);
					flg_gts_1 = false;
				}
				//incompleted spans
				int key_gst_0 = getKey(g, s, t, 0, 0, length);
				alpha[key_gst_0] = 0.0;
				int key_gts_0 = getKey(g, s, t, 1, 0, length);
				alpha[key_gts_0] = 0.0;

				bool flg_gst_0 = true, flg_gts_0 = true;

				for(int r = t; r < length; r++){
					if(r != t){
						REAL_SCORES prodProb_gstr = probs[get_index2_o3g(length,g,s,t,r)];
						key_b = getKey(s, t, r, 0, 2, length);
						key_a = getKey(g, s, r, 0, 0, length);
						alpha[key_gst_0] = logsumexp(alpha[key_gst_0], beta[key_b] + alpha[key_a] + prodProb_gstr, flg_gst_0);
						flg_gst_0 = false;
					}
					key_b = getKey(s, t, r, 0, 1, length);
					key_a = getKey(g, s, r, 0, 1, length);
					alpha[key_gst_0] = logsumexp(alpha[key_gst_0], beta[key_b] + alpha[key_a], flg_gst_0);
					flg_gst_0 = false;
				}
				for(int r = g + 1; r <= s; r++){
					if(r != s){
						REAL_SCORES prodProb_gtsr = probs[get_index2_o3g(length,g,t,s,r)];
						key_b = getKey(t, r, s, 0, 2, length);
						key_a = getKey(g, r, t, 1, 0, length);
						alpha[key_gts_0] = logsumexp(alpha[key_gts_0], beta[key_b] + alpha[key_a] + prodProb_gtsr, flg_gts_0);
						flg_gts_0 = false;
					}
					key_b = getKey(t, r, s, 1, 1, length);
					key_a = getKey(g, r, t, 1, 1, length);
					alpha[key_gts_0] = logsumexp(alpha[key_gts_0], beta[key_b] + alpha[key_a], flg_gts_0);
					flg_gts_0 = false;
				}
			}
			for(int g = t + 1; g < length; g++){
				//sibling spans
				int key_gst_2 = getKey(g, s, t, 0, 2, length);
				alpha[key_gst_2] = 0.0;
				bool flg_gst_2 = true;
				for(int r = 0; r < length; r++){
					if(r >= s && r <= g){
						continue;
					}
					REAL_SCORES prodProb_rgts = probs[get_index2_o3g(length,r,g,t,s)];
					key_b = getKey(r, t, g, 1, 0, length);
					key_a = getKey(r, s, g, 1, 0, length);
					alpha[key_gst_2] = logsumexp(alpha[key_gst_2], beta[key_b] + alpha[key_a] + prodProb_rgts, flg_gst_2);
					flg_gst_2 = false;
				}
				if(g == end){
					REAL_SCORES prodProb_nts = probs[get_index2_o3g(length,end,end,t,s)];;
					key_b = getKey(t, t, end, 1, 0, length);
					key_a = getKey(s, s, end, 1, 0, length);
					alpha[key_gst_2] = logsumexp(alpha[key_gst_2], beta[key_b] + alpha[key_a] + prodProb_nts, flg_gst_2);
					flg_gst_2 = false;
				}
				//completed spans
				int key_gst_1 = getKey(g, s, t, 0, 1, length);
				alpha[key_gst_1] = 0.0;
				int key_gts_1 = getKey(g, s, t, 1, 1, length);
				alpha[key_gts_1] = 0.0;

				bool flg_gst_1 = true, flg_gts_1 = true;

				for(int r = 0; r < length; r++){
					if(r > t && r < g){
						key_b = getKey(g, t + 1, r, 1, 1, length);
						key_a = getKey(g, s, r, 0, 2, length);
						alpha[key_gst_1] = logsumexp(alpha[key_gst_1], beta[key_b] + alpha[key_a], flg_gst_1);
						flg_gst_1 = false;
					}
					if(r >= 0 && r < s){
						key_b = getKey(g, r, s - 1, 0, 1, length);
						key_a = getKey(g, r, t, 0, 2, length);
						alpha[key_gts_1] = logsumexp(alpha[key_gts_1], beta[key_b] + alpha[key_a], flg_gts_1);
						flg_gts_1 = false;
					}
					if(r < s || r > g){
						if(g == t + 1){
							REAL_SCORES prodProb_rggs = probs[get_index2_o3g(length,r,g,g,s)];
							key_a = getKey(r, s, g, 1, 0, length);
							alpha[key_gst_1] = logsumexp(alpha[key_gst_1], alpha[key_a] + prodProb_rggs, flg_gst_1);
							flg_gst_1 = false;
						}

						key_b = getKey(r, t, g, 1, 0, length);
						key_a = getKey(r, s, g, 1, 1, length);
						alpha[key_gts_1] = logsumexp(alpha[key_gts_1], beta[key_b] + alpha[key_a], flg_gts_1);
						flg_gts_1 = false;
					}
				}

				if(g == end && t == end - 1){
					REAL_SCORES prodProb_nns = probs[get_index2_o3g(length,end,end,end,s)];
					key_a = getKey(s, s, end, 1, 0, length);
					alpha[key_gst_1] = logsumexp(alpha[key_gst_1], alpha[key_a] + prodProb_nns, flg_gst_1);
					flg_gst_1 = false;
				}
				if(g == end && s == 0){
					key_b = getKey(t, t, end, 1, 0, length);
					key_a = getKey(0, 0, end, 1, 1, length);
					alpha[key_gts_1] = logsumexp(alpha[key_gts_1], beta[key_b] + alpha[key_a], flg_gts_1);
					flg_gts_1 = false;
				}

				//incompleted spans
				int key_gst_0 = getKey(g, s, t, 0, 0, length);
				alpha[key_gst_0] = 0.0;
				int key_gts_0 = getKey(g, s, t, 1, 0, length);
				alpha[key_gts_0] = 0.0;

				bool flg_gst_0 = true, flg_gts_0 = true;

				for(int r = t; r < g; r++){
					if(r != t){
						REAL_SCORES prodProb_gstr = probs[get_index2_o3g(length,g,s,t,r)];
						key_b = getKey(s, t, r, 0, 2, length);
						key_a = getKey(g, s, r, 0, 0, length);
						alpha[key_gst_0] = logsumexp(alpha[key_gst_0], beta[key_b] + alpha[key_a] + prodProb_gstr, flg_gst_0);
						flg_gst_0 = false;
					}
					key_b = getKey(s, t, r, 0, 1, length);
					key_a = getKey(g, s, r, 0, 1, length);
					alpha[key_gst_0] = logsumexp(alpha[key_gst_0], beta[key_b] + alpha[key_a], flg_gst_0);
					flg_gst_0 = false;
				}
				for(int r = 0; r <= s; r++){
					if(r != s){
						REAL_SCORES prodProb_gtsr = probs[get_index2_o3g(length,g,t,s,r)];
						key_b = getKey(t, r, s, 0, 2, length);
						key_a = getKey(g, r, t, 1, 0, length);
						alpha[key_gts_0] = logsumexp(alpha[key_gts_0], beta[key_b] + alpha[key_a] + prodProb_gtsr, flg_gts_0);
						flg_gts_0 = false;
					}
					key_b = getKey(t, r, s, 1, 1, length);
					key_a = getKey(g, r, t, 1, 1, length);
					alpha[key_gts_0] = logsumexp(alpha[key_gts_0], beta[key_b] + alpha[key_a], flg_gts_0);
					flg_gts_0 = false;
				}
			}
		}
	}
}

//##MAGIC NUMBERS of MAX-Encoder##
//direction: 0,st(right);1,ts(left)
//spans: 0,incomplete;1,complete;2,sibling


REAL_SCORES* LencodeMarginals_o3g(const int length,Scores<REAL_SCORES>& scores)
{
	int ln = scores.get_numc();
	REAL_SCORES* marginals = new REAL_SCORES[scores.get_numl()*ln];	//the return table
	REAL_SCORES *beta = new REAL_SCORES[length * length * length * 2 * 3];
	REAL_SCORES *alpha = new REAL_SCORES[length * length * length * 2 * 3];
	//sumlabel score
	OneScores<REAL_SCORES> sum_scores = scores.get_scores_marginal();
	REAL_SCORES z = calc_inside(length, beta,sum_scores);
	calc_outside(length,beta,sum_scores,alpha);

	//get them
	long all_num = 0;
	for(int s=0;s<length;s++){
		for(int t=s+1;t<length;t++){
			//loop on g
			for(int g=0;g<length;g++){
				if(g>=s && g<=t)
					continue;
			//gsst
			{
				long index_curr = get_index2_o3g(length,g,s,s,t);
				if(scores.has_value(index_curr)){
					all_num ++;
					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(beta[getKey(s,s+1,t,1,1,length)]+alpha[getKey(g,s,t,0,0,length)]+from_assign[zl]-z);
				}
			}
			for(int r=s+1;r<t;r++){
				//gsrt
				{
					long index_curr = get_index2_o3g(length,g,s,r,t);
					if(scores.has_value(index_curr)){
						all_num ++;
						REAL_SCORES* from_assign = scores[index_curr];
						REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
						for(int zl=0;zl<ln;zl++)
							to_assign[zl] = exp(beta[getKey(g,s,r,0,0,length)]+beta[getKey(s,r,t,0,2,length)]
															+alpha[getKey(g,s,t,0,0,length)]+from_assign[zl]-z);
					}
				}
			}
			//gtts
			{
				long index_curr = get_index2_o3g(length,g,t,t,s);
				if(scores.has_value(index_curr)){
					all_num ++;
					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(beta[getKey(t,s,t-1,0,1,length)]+alpha[getKey(g,s,t,1,0,length)]+from_assign[zl]-z);
				}
			}
			for(int r=s+1;r<t;r++){
				//gtrs
				{
					long index_curr = get_index2_o3g(length,g,t,r,s);
					if(scores.has_value(index_curr)){
						all_num ++;
						REAL_SCORES* from_assign = scores[index_curr];
						REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
						for(int zl=0;zl<ln;zl++)
							to_assign[zl] = exp(beta[getKey(g,r,t,1,0,length)]+beta[getKey(t,s,r,0,2,length)]
																+alpha[getKey(g,s,t,1,0,length)]+from_assign[zl]-z);
					}
				}
			}
			}
		}
	}
	for(int t=1;t<length;t++){
		//000t
		{
			long index_curr = get_index2_o3g(length,0,0,0,t);
			if(scores.has_value(index_curr)){
				all_num ++;
				REAL_SCORES* from_assign = scores[index_curr];
				REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
				for(int zl=0;zl<ln;zl++)
					to_assign[zl] = exp(beta[getKey(0,1,t,1,1,length)]+alpha[getKey(0,0,t,0,0,length)]+from_assign[zl]-z);
			}
		}
		//00rt
		for(int r=1;r<t;r++){
			{
				long index_curr = get_index2_o3g(length,0,0,r,t);
				if(scores.has_value(index_curr)){
					all_num ++;
					REAL_SCORES* from_assign = scores[index_curr];
					REAL_SCORES* to_assign = marginals + ln*scores.get_value(index_curr);
					for(int zl=0;zl<ln;zl++)
						to_assign[zl] = exp(beta[getKey(0,0,r,0,0,length)]+beta[getKey(0,r,t,0,2,length)]
														+alpha[getKey(0,0,t,0,0,length)]+from_assign[zl]-z);
				}
			}
		}
	}

	delete []beta;
	delete []alpha;
	WARNING_UNEQUAL(all_num,scores.get_numl(),"Unequal instances.");
	return marginals;
}



