/*
 * Helper.h
 *
 *  Created on:
 *      Author: zzs
 */

#ifndef ALGORITHMS_HELPER_H_
#define ALGORITHMS_HELPER_H_

#include <cmath>
const double MINUS_LOG_EPSILON=50;

// log(exp(x) + exp(y));
template<class T>
inline T logsumexp(T x, T y, bool flg) {
	if (flg) return y;  // init mode
	const T vmin = (x>y)?y:x;
	const T vmax = (x<y)?y:x;
	if (vmax > vmin + MINUS_LOG_EPSILON) {
		return vmax;
	}else{
		return vmax + log(exp(vmin - vmax) + 1.0);
	}
}

#endif /* ALGORITHMS_HELPER_H_ */
