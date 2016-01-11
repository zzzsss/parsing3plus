/*
 * Scores.h
 *
 *  Created on: Jan 11, 2016
 *      Author: zzs
 */

#ifndef ALGORITHMS_SCORES_H_
#define ALGORITHMS_SCORES_H_

#include<boost/unordered_map.hpp>
#include<boost/functional/hash.hpp>
#include "Helper.h"
using namespace std;
using namespace boost;

typedef float REAL_SCORES;
typedef boost::unordered_map<long, long> ScoreHashMap;	//index -> index in the table

namespace the_scores{

//warnings
template<class T>
static void WARNING_UNEQUAL(T a,T b,const char* x){
	if(a != b)
		cerr << "!! WARNING of " << x << ":" << a << " != " << b << endl;
}
template<class T>
static void WARNING(T a,const char* x){
	cerr << "!! WARNING directly of " << x << endl;
}

//the scores with num_classes == 1, used for algorithms
template<class T>
class OneScores{
private:
	//!!!! all from outside, so no need for copy,assign,destruction ...
	long the_num;	//number of keys
	ScoreHashMap* the_map;
	T*  the_table;
	T default_value;
public:
	OneScores(long n,ScoreHashMap* m,T* t,T dd):the_num(n),the_map(m),the_table(t),default_value(dd){}
	T operator[](long index){
		ScoreHashMap::iterator iter = the_map->find(index);
		if(iter == the_map->end())
			return default_value;
		else
			return the_table[iter->second];
	}
};

//The collection of the scores
template<class T>
class Scores{
private:
	static const long MAP_DEFAULT_LENGTH = 1000000;
	static const int MAP_ALPHA = 10;

	int  num_class;	//number of classes per link
	long num_links;	//number of keys
	ScoreHashMap* scores_map;
	T*  scores_table;	//the scores-table, origin one	!!OUTSIDE!!
	T default_value;

	T*  scores_max;	//special scores: max on
	int* index_max;	//special indexes: which one is the max
	T*  scores_mar;	//special scores: marginal one

public:
	int get_numc(){return num_class;}
	long get_numl()	{return num_links;}

	bool has_value(long ind){
		return scores_map->find(ind)!=scores_map->end();
	}
	long get_value(long ind){
		ScoreHashMap::iterator iter = scores_map->find(ind);
		if(iter == scores_map->end())
			return -1;
		else
			return iter->second;
	}

	Scores(int numc,long numl,T* tab,T dd):num_class(numc),num_links(numl),scores_table(tab),default_value(dd),
											scores_max(0),index_max(0),scores_mar(0){
		long msize = MAP_DEFAULT_LENGTH;
		if(msize < num_links*MAP_ALPHA)
			msize = num_links*MAP_ALPHA;
		scores_map = new ScoreHashMap(msize);	//haven't built associations yet
	}
	void add_associate(long key,long value){
		ScoreHashMap::iterator iter = scores_map->find(key);
		if(iter == scores_map->end())
			scores_map->insert(key,value);
		else
			WARNING(key,"Repeat key for Scores");
	}

	OneScores<T> get_scores_maxone(){	//here return object directly
		scores_max = new T[num_links];
		index_max = new int[num_links];
		T* to_assign = scores_table;
		for(int i=0;i<num_links;i++){
			T tmp_max = *to_assign;
			int tmp_max_index = 0;
			for(int ll=1;ll<num_class;ll++){
				if(to_assign[ll] > tmp_max){
					tmp_max = to_assign[ll];
					tmp_max_index = ll;
				}
			}
			scores_max[i] = tmp_max;
			index_max[i] = tmp_max_index;
			to_assign += num_class;
		}
		return OneScores<T>(num_links,scores_map,scores_max,default_value);
	}
	OneScores<T> get_scores_marginal(){	//here return object directly
		scores_mar = new T[num_links];
		T* to_assign = scores_table;
		for(int i=0;i<num_links;i++){
			T accu = *to_assign;
			for(int ll=1;ll<num_class;ll++)
				accu = logsumexp(accu,to_assign[ll],false);
			scores_mar[i] = accu;
			to_assign += num_class;
		}
		return OneScores<T>(num_links,scores_map,scores_mar,default_value);
	}

	T* operator[](long index){
		ScoreHashMap::iterator iter = scores_map->find(index);
		if(iter == scores_map->end())
			return 0;
		else
			return &scores_map[iter->second*num_class];
	}
};

}




#endif /* ALGORITHMS_SCORES_H_ */
