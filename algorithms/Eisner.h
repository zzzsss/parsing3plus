/*
 * Eisner.h
 *
 *  Created on: Jul 10, 2014
 *      Author: zzs
 */

#ifndef EISNER_H_
#define EISNER_H_

#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;
#include "Scores.h"
using namespace the_scores;

#define Negative_Infinity -1e100

#define E_LEFT 0
#define E_RIGHT 1
#define E_INCOM 0
#define E_COM 1

//the index explanation --- C[len][len][2][2]
inline int get_index(int len,int s,int t,int lr,int c)
{
	int key = s;
	key = key * len + t;
	key = key * 2 + lr;
	key = key * 2 + c;
	return key;
}
//the index2 -- scores
// -- unified with S[h][m]
inline int get_index2(int len,int h,int m)
{
	int key = h;
	key = key * len + m;
	return key;
}

#define E_ERROR(str) \
	cout << str << endl;\
	exit(1)

extern vector<int>* decodeProjective(int length,OneScores<REAL_SCORES>& scores);	//only one class
extern REAL_SCORES* LencodeMarginals(const int length,Scores<REAL_SCORES>& scores);	//same as the ones in scores
extern REAL_SCORES* encodeMarginals(const int length,REAL_SCORES* scores);

#endif /* EISNER_H_ */
