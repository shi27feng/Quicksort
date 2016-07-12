/* ************************************************************************* *\
               INTEL CORPORATION PROPRIETARY INFORMATION
     This software is supplied under the terms of a license agreement or 
     nondisclosure agreement with Intel Corporation and may not be copied 
     or disclosed except in accordance with the terms of that agreement. 
        Copyright (C) 2014 Intel Corporation. All Rights Reserved.
\* ************************************************************************* */

#ifndef QUICKSORT_H
#define QUICKSORT_H

#ifdef HOST
template <class T>
T median(T x1, T x2, T x3) {
	if (x1 < x2) {
		if (x2 < x3) {
			return x2;
		} else {
			if (x1 < x3) {
				return x3;
			} else {
				return x1;
			}
		}
	} else { // x1 >= x2
		if (x1 < x3) {
			return x1;
		} else { // x1 >= x3
			if (x2 < x3) {
				return x2;
			} else {
				return x3;
			}
		}
	}
}
#else // HOST
uint median(uint x1, uint x2, uint x3) {
	if (x1 < x2) {
		if (x2 < x3) {
			return x2;
		} else {
			if (x1 < x3) {
				return x3;
			} else {
				return x1;
			}
		}
	} else { // x1 >= x2
		if (x1 < x3) {
			return x1;
		} else { // x1 >= x3
			if (x2 < x3) {
				return x2;
			} else {
				return x3;
			}
		}
	}
}
#endif //HOST

//#define HASWELL 1
//#define GET_DETAILED_PERFORMANCE 1
#define TRUST_BUT_VERIFY 1
// Note that SORT_THRESHOLD should always be 2X LOCAL_THREADCOUNT due to the use of bitonic sort
// Always try LOCAL_THREADCOUNT to be 8X smaller than QUICKSORT_BLOCK_SIZE - then try everything else :)
#ifdef HASWELL
#define QUICKSORT_BLOCK_SIZE         1024 
#define GQSORT_LOCAL_WORKGROUP_SIZE   128 
#define LQSORT_LOCAL_WORKGROUP_SIZE   128 
#define SORT_THRESHOLD                256 
#define PRIVATE_SORT_THRESHOLD         12
#else
#define QUICKSORT_BLOCK_SIZE         1024 
#define GQSORT_LOCAL_WORKGROUP_SIZE   128 
#define LQSORT_LOCAL_WORKGROUP_SIZE   256 
#define SORT_THRESHOLD                512 
#define PRIVATE_SORT_THRESHOLD         10
#endif

#define EMPTY_RECORD             42

typedef struct work_record {
	uint start;
	uint end;
	uint pivot;
	uint direction;
#ifdef HOST
	work_record() : 
		start(0), end(0), pivot(0), direction(EMPTY_RECORD) {}
	work_record(uint s, uint e, uint p, uint d) : 
		start(s), end(e), pivot(p), direction(d) {}
#endif // HOST
} work_record;

typedef struct parent_record {
	uint sstart, send, oldstart, oldend, blockcount; 
#ifdef HOST
	parent_record(uint ss, uint se, uint os, uint oe, uint bc) : 
		sstart(ss), send(se), oldstart(os), oldend(oe), blockcount(bc) {}
#endif // HOST
} parent_record;

typedef struct block_record {
	uint start;
	uint end;
	uint pivot;
	uint direction;
	uint parent;
#ifdef HOST
	block_record(uint s, uint e, uint p, uint d, uint prnt) : 
		start(s), end(e), pivot(p), direction(d), parent(prnt) {}
#endif // HOST
} block_record;
#endif // QUICKSORT_H
