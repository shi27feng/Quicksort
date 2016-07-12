/* ************************************************************************* *\
               INTEL CORPORATION PROPRIETARY INFORMATION
     This software is supplied under the terms of a license agreement or 
     nondisclosure agreement with Intel Corporation and may not be copied 
     or disclosed except in accordance with the terms of that agreement. 
        Copyright (C) 2014 Intel Corporation. All Rights Reserved.
\* ************************************************************************* */

#include "Quicksort.h"

void plus_prescan(local uint *a, local uint *b) {
    uint av = *a;
	uint bv = *b;
    *a = bv;
    *b = bv + av;
}

/// bitonic_sort: sort 2*LOCAL_THREADCOUNT elements
void bitonic_sort(local uint* sh_data, const uint local_idx) {
	for (uint ulevel = 1; ulevel < LQSORT_LOCAL_WORKGROUP_SIZE; ulevel <<= 1) {
        for (uint j = ulevel; j > 0; j >>= 1) {
            uint pos = 2*local_idx - (local_idx & (j - 1));

			uint direction = local_idx & ulevel;
			uint av = sh_data[pos], bv = sh_data[pos + j];
			const bool sortThem = av > bv;
			const uint greater = select(av, bv, sortThem);
			const uint lesser  = select(bv, av, sortThem);

			sh_data[pos]     = select(greater, lesser, direction);
			sh_data[pos + j] = select(lesser, greater, direction);
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

	for (uint j = LQSORT_LOCAL_WORKGROUP_SIZE; j > 0; j >>= 1) {
        uint pos = 2*local_idx - (local_idx & (j - 1));

		uint av = sh_data[pos], bv = sh_data[pos + j];
		const bool sortThem = av > bv;
		sh_data[pos]      = select(av, bv, sortThem);
		sh_data[pos + j]  = select(bv, av, sortThem);

        barrier(CLK_LOCAL_MEM_FENCE);
    }
}

void sort_threshold(local uint* data_in, global uint* data_out,
					uint start, 
					uint end, local uint* temp, uint localid) {
	uint tsum = end - start;
	if (tsum == SORT_THRESHOLD) {
		bitonic_sort(data_in+start, localid);
		for (uint i = localid; i < SORT_THRESHOLD; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			data_out[start + i] = data_in[start + i];
		}
	} else if (tsum > 1) {
		for (uint i = localid; i < SORT_THRESHOLD; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			if (i < tsum) {
				temp[i] = data_in[start + i];
			} else {
				temp[i] = UINT_MAX;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		bitonic_sort(temp, localid);

		for (uint i = localid; i < tsum; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			data_out[start + i] = temp[i];
		}
	} else if (tsum == 1 && localid == 0) {
		data_out[start] = data_in[start];
	} 
}

//----------------------------------------------------------------------------
// Kernel implements gqsort_kernel
//----------------------------------------------------------------------------
kernel void gqsort_kernel(global uint* d, global uint*dn, global block_record* blocks, global parent_record* parents, global work_record* result) 
{
	const uint blockid     = get_group_id(0);
	const uint localid    = get_local_id(0);
	local uint lt[GQSORT_LOCAL_WORKGROUP_SIZE+1], gt[GQSORT_LOCAL_WORKGROUP_SIZE+1], ltsum, gtsum, lbeg, gbeg;
	uint i, lfrom, gfrom, lpivot, gpivot, tmp, ltp = 0, gtp = 0;

	// Get the sequence block assigned to this thread block
	block_record block = blocks[blockid];
	uint start = block.start, end = block.end, pivot = block.pivot, direction = block.direction;

	//printf("blockid %d, start %d, end %d, pivot %d, direction %d\n", blockid, start, end, pivot, direction);
	global parent_record* pparent = parents + block.parent; 
	global uint* psstart, *psend, *poldstart, *poldend, *pblockcount;
	global uint *s, *sn;

	if (direction == 1) {
		s = d;
		sn = dn;
	} else {
		s = dn;
		sn = d;
	}

	// Set thread local counters to zero
	lt[localid] = gt[localid] = 0;
	barrier(CLK_LOCAL_MEM_FENCE);

	// Align thread accesses for coalesced reads.
	i = start + localid;
	// Go through data...
	for(; i < end; i += GQSORT_LOCAL_WORKGROUP_SIZE) {
		tmp = s[i];
		// counting elements that are smaller ...
		if (tmp < pivot)
			ltp++;
		// or larger compared to the pivot.
		if (tmp > pivot) 
			gtp++;
	}
	lt[localid] = ltp;
	gt[localid] = gtp;
	barrier(CLK_LOCAL_MEM_FENCE);

	// calculate cumulative sums
	uint n;
	for(i = 1; i < GQSORT_LOCAL_WORKGROUP_SIZE; i <<= 1) {
		n = 2*i - 1;
		if ((localid & n) == n) {
			lt[localid] += lt[localid-i];
			gt[localid] += gt[localid-i];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if ((localid & n) == n) {
		lt[GQSORT_LOCAL_WORKGROUP_SIZE] = ltsum = lt[localid];
		gt[GQSORT_LOCAL_WORKGROUP_SIZE] = gtsum = gt[localid];
		lt[localid] = 0;
		gt[localid] = 0;
	}
		
	for(i = GQSORT_LOCAL_WORKGROUP_SIZE/2; i >= 1; i >>= 1) {
		n = 2*i - 1;
		if ((localid & n) == n) {
			plus_prescan(&lt[localid - i], &lt[localid]);
			plus_prescan(&gt[localid - i], &gt[localid]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// Allocate memory in the sequence this block is a part of
	if (localid == 0) {
		// get shared variables
		psstart = &pparent->sstart;
		psend = &pparent->send;
		poldstart = &pparent->oldstart;
		poldend = &pparent->oldend;
		pblockcount = &pparent->blockcount;
		// Atomic increment allocates memory to write to.
		lbeg = atomic_add(psstart, ltsum);
		// Atomic is necessary since multiple blocks access this
		gbeg = atomic_sub(psend, gtsum) - gtsum;
	}
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
	// Allocate locations for threads
	lfrom = lbeg + lt[localid];
	gfrom = gbeg + gt[localid];

	// go thru data again writing elements to their correct position
	i = start + localid;
	for(; i < end; i += GQSORT_LOCAL_WORKGROUP_SIZE) {
		tmp = s[i];
		// increment counts
		if (tmp < pivot) 
			sn[lfrom++] = tmp;

		if (tmp > pivot) 
			sn[gfrom++] = tmp;
	}
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	if (localid == 0) {
		if (atomic_dec(pblockcount) == 0) {
			uint sstart = *psstart;
			uint send = *psend;
			uint oldstart = *poldstart;
			uint oldend = *poldend;

			// Store the pivot value between the new sequences
			for(i = sstart; i < send; i ++) {
				d[i] = pivot;
			}
			lpivot = sn[oldstart];
			gpivot = sn[oldend-1];
			if (oldstart < sstart) {
				lpivot = median(lpivot,sn[(oldstart+sstart)/2], sn[sstart-1]);
			} 
			if (send < oldend) {
				gpivot = median(sn[send],sn[(oldend+send)/2], gpivot);
			}

			global work_record* result1 = result + 2*blockid;
			global work_record* result2 = result1 + 1;
			direction ^= 1;

			work_record r1 = {oldstart, sstart, lpivot, direction};
			*result1 = r1;

			work_record r2 = {send, oldend, gpivot, direction};
			*result2 = r2;
		}
	}
}

// record to push start of the sequence, end of the sequence and direction of sorting on internal stack
typedef struct workstack_record {
	uint start;
	uint end;
	uint direction;
} workstack_record;

#define PUSH(START, END) 			if (localid == 0) { \
										++workstack_pointer; \
                                        workstack_record wr = { (START), (END), direction ^ 1 }; \
										workstack[workstack_pointer] = wr; \
									} \
									barrier(CLK_LOCAL_MEM_FENCE);

//---------------------------------------------------------------------------------------
// Kernel implements the last stage of GPU-Quicksort, when all the subsequences are small
// enough to be processed in local memory. It uses similar algorithm to gqsort_kernel to 
// move items around the pivot and then switches to either bitonic sort for sequences in
// the range (PRIVATE_SORT_THRESHOLD, SORT_THRESHOLD] or shell sort for small sequences in
// the range [1, PRIVATE_SORT_THRESHOLD].
//
// d - input array
// dn - scratch array of the same size as the input array
// seqs - array of records to be sorted in a local memory, one sequence per work group.
//---------------------------------------------------------------------------------------
kernel void lqsort_kernel(global uint* d, global uint* dn, global work_record* seqs) 
{
	const uint blockid     = get_group_id(0);
	const uint localid     = get_local_id(0);

	// workstack: stores the start and end of the sequences, direction of sort
	// If the sequence is less that SORT_THRESHOLD, it gets sorted. 
	// It will only be pushed on the stack if it greater than the SORT_THRESHOLD. 
	// Note, that the sum of ltsum + gtsum is less than QUICKSORT_BLOCK_SIZE. 
	// The total sum of the length of records on the stack cannot exceed QUICKSORT_BLOCK_SIZE, 
	// but each individual record should be greater than SORT_THRESHOLD, so the maximum length 
	// of the stack is QUICKSORT_BLOCK_SIZE/SORT_THRESHOLD - in the case of BDW GT2 the length 
	// of the stack is 2 :)
	local workstack_record workstack[QUICKSORT_BLOCK_SIZE/SORT_THRESHOLD]; // 
	local int workstack_pointer;

	local uint mys[QUICKSORT_BLOCK_SIZE], mysn[QUICKSORT_BLOCK_SIZE];
	local uint *s, *sn;
	local uint lt[LQSORT_LOCAL_WORKGROUP_SIZE+1], gt[LQSORT_LOCAL_WORKGROUP_SIZE+1], ltsum, gtsum;
	local uint temp[SORT_THRESHOLD];
	uint i, tmp, ltp, gtp;
	
	const uint d_offset = seqs[blockid].start;
	uint start = 0; 
	uint end   = seqs[blockid].end - d_offset;

	uint direction = 1; // which direction to sort
	// initialize workstack and workstack_pointer: push the initial sequence on the stack
	if (localid == 0) {
		workstack_pointer = 0; // beginning of the stack
		workstack_record wr = { start, end, direction };
		workstack[0] = wr;
	}
	// copy block of data to be sorted by one workgroup into local memory
	// note that indeces of local data go from 0 to end-start-1
	if (seqs[blockid].direction == 1) {
		for (i = localid; i < end; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			mys[i] = d[i+d_offset];
		}
	} else {
		for (i = localid; i < end; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			mys[i] = dn[i+d_offset];
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	while (workstack_pointer >= 0) { 
		// pop up the stack
		workstack_record wr = workstack[workstack_pointer];
		start = wr.start;
		end = wr.end;
		direction = wr.direction;
		barrier(CLK_LOCAL_MEM_FENCE);
		if (localid == 0) {
			--workstack_pointer;

			ltsum = gtsum = 0;	
		}
		if (direction == 1) {
			s = mys;
			sn = mysn;
		} else {
			s = mysn;
			sn = mys;
		}
		// Set thread local counters to zero
		lt[localid] = gt[localid] = 0;
		ltp = gtp = 0;
		barrier(CLK_LOCAL_MEM_FENCE);

		// Pick a pivot
		uint pivot = s[start];
		if (start < end) {
			pivot = median(pivot, s[(start+end)/2], s[end-1]);
		}
		
		// Align thread accesses for coalesced reads.
		i = start + localid;
		// Go through data...
		for(; i < end; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			tmp = s[i];
			// counting elements that are smaller ...
			if (tmp < pivot)
				ltp++;
			// or larger compared to the pivot.
			if (tmp > pivot) 
				gtp++;
		}
		lt[localid] = ltp;
		gt[localid] = gtp;
		barrier(CLK_LOCAL_MEM_FENCE);
		
		// calculate cumulative sums
		uint n;
		for(i = 1; i < LQSORT_LOCAL_WORKGROUP_SIZE; i <<= 1) {
			n = 2*i - 1;
			if ((localid & n) == n) {
				lt[localid] += lt[localid-i];
				gt[localid] += gt[localid-i];
			}
			barrier(CLK_LOCAL_MEM_FENCE);
		}

		if ((localid & n) == n) {
			lt[LQSORT_LOCAL_WORKGROUP_SIZE] = ltsum = lt[localid];
			gt[LQSORT_LOCAL_WORKGROUP_SIZE] = gtsum = gt[localid];
			lt[localid] = 0;
			gt[localid] = 0;
		}
		
		for(i = LQSORT_LOCAL_WORKGROUP_SIZE/2; i >= 1; i >>= 1) {
			n = 2*i - 1;
			if ((localid & n) == n) {
				plus_prescan(&lt[localid - i], &lt[localid]);
				plus_prescan(&gt[localid - i], &gt[localid]);
			}
			barrier(CLK_LOCAL_MEM_FENCE);
		}

		// Allocate locations for threads
		uint lfrom = start + lt[localid];
		uint gfrom = end - gt[localid+1];

		// go thru data again writing elements to their correct position
		i = start + localid;
		for (; i < end; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			tmp = s[i];
			// increment counts
			if (tmp < pivot) 
				sn[lfrom++] = tmp;
			
			if (tmp > pivot) 
				sn[gfrom++] = tmp;
		}
		barrier(CLK_LOCAL_MEM_FENCE);

		// Store the pivot value between the new sequences
		i = start + ltsum + localid;
		for (;i < end - gtsum; i += LQSORT_LOCAL_WORKGROUP_SIZE) {
			d[i+d_offset] = pivot;
		}
		barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

		// fill in pivot, produce two new sequences
		if (ltsum <= SORT_THRESHOLD) {
			sort_threshold(sn, d+d_offset, start, start + ltsum, temp, localid);
		} else {
			PUSH(start, start + ltsum);
		}
		
		if (gtsum <= SORT_THRESHOLD) {
			sort_threshold(sn, d+d_offset, end - gtsum, end, temp, localid);
		} else {
			PUSH(end - gtsum, end);
		}
	}
}