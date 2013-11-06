#include "search.h"

bool BWExactFinalResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t start, pos;
	int16_t current_block;
	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		start = r_iterator->start;
		pos   = r_iterator->pos;

		k = r_iterator->k;
		l = r_iterator->l;

		current_block = pos / block_size;

		if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated
	
		} else {

			if (current_block > last_block) { //Not in last previsited block

				continue;

			} else { //I am in the last previsited block

				if ((pos + 1) % block_size) { //I am not in the first element of the block
				} else { //I am in the first element in the block (all the block must be processed)
					continue;
				}

			}

		}

		for(int16_t i=pos; i>=start; i--) {
			BWiteration(k, l, k, l, W[i], index);
			if (k > l) break;
		}

		if (k <= l) {
			change_result(r_iterator, k, l, start-1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWExactFinalResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t pos, end;
	int16_t current_block;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		pos   = r_iterator->pos;
		end   = r_iterator->end;

		k = r_iterator->k;
		l = r_iterator->l;

		current_block = pos / block_size;

		if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated

		} else {

			if (current_block < last_block) { //Not in last previsited block

				continue;

			} else { //I am in the last previsited block

				if (pos % block_size) { //I am not in the first element of the block
				} else { //I am in the first element in the block (all the block must be processed)
					continue;
				}

			}

		}

		for(int16_t i=pos; i<=end; i++) {
			BWiteration(k, l, k, l, W[i], index);
			if (k > l) break;
		}

		if (k <= l) {
			change_result(r_iterator, k, l, end+1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchFinalResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_final, int16_t block_size, int16_t last_block) {

	intmax_t k, l, k_aux, l_aux;
	int16_t start, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 0) continue;

		start = r_iterator->start;
		pos   = r_iterator->pos;

		if (pos < start) {
			if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos-1);
					if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					//Insertion
					if (b!=W[last_err_pos]) {
						change_result(r_iterator, k_aux, l_aux, pos);
						modify_last_mismatch2(r_iterator, INSERTION, b);
						if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

					//Mismatch
					if (b!=W[pos]) {
						change_result(r_iterator, k_aux, l_aux, pos-1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos-1);
				if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos-1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				//Mismatch
				if (b!=W[pos]) {

					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos-1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos-1);
			if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				//Mismatch
				if (b!=W[pos]) {
					r_iterator->pos = pos-1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchFinalResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_final, int16_t block_size, int16_t last_block) {

	intmax_t k, l, k_aux, l_aux;
	int16_t end, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 1) continue;

		end = r_iterator->end;
		pos = r_iterator->pos;

		if (pos > end) {
			add_result(r_iterator, rl_final);
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos+1);
					if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					//Insertion
					if (b!=W[last_err_pos]) {
						change_result(r_iterator, k_aux, l_aux, pos);
						modify_last_mismatch2(r_iterator, INSERTION, b);
						if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

					//Mismatch
					if (b!=W[pos]) {
						change_result(r_iterator, k_aux, l_aux, pos+1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos+1);
				if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos+1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				//Mismatch
				if (b!=W[pos]) {

					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos+1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos+1);
			if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				if (b!=W[pos]) { //Mismatch
					r_iterator->pos = pos+1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWExactPartialResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l, k_next, l_next;
	int16_t start, pos, current_block, last_block_pos;
	bool complete_search;

	result *r_iterator;
	intmax_t results, results_next;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		start  = r_iterator->start;
		pos    = r_iterator->pos;

		k_next = r_iterator->k;
		l_next = r_iterator->l;
		results_next = l_next - k_next;

		current_block = pos / block_size;

		if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated

			last_block_pos = start;
			complete_search = true;

		} else {

			if (current_block > last_block) { //Not in last previsited block

				if ((pos + 1) % block_size) { //Not in first element of the block
					last_block_pos = (current_block-1) * block_size;
				} else { //I am in the first element in the block (all the block must be processed)
					last_block_pos = current_block * block_size;
				}

				complete_search = false;

			} else { //I am in the last previsited block

				if ((pos + 1) % block_size) { //I am not in the first element of the block
					last_block_pos = start;
					complete_search = true;
				} else { //I am in the first element of the block (all the block must be processed)
					last_block_pos = current_block * block_size;
					complete_search = false;
				}

			}

		}

		for(int16_t i=pos; i>=last_block_pos; i--) {

			k = k_next;
			l = l_next;

			if (k > l) break;

			BWiteration(k, l, k_next, l_next, W[i], index);
			results      = results_next;
			results_next = l_next - k_next;
			if (results == results_next) continue;

			change_result(r_iterator, k, l, i);
			add_result(r_iterator, rl_next);

		}

		if (complete_search && k_next <= l_next) {
			change_result(r_iterator, k_next, l_next, start-1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWExactPartialResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l, k_next, l_next;
	int16_t pos, end, current_block, last_block_pos;
	bool complete_search;

	result *r_iterator;
	intmax_t results, results_next;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		pos   = r_iterator->pos;
		end   = r_iterator->end;

		k_next = r_iterator->k;
		l_next = r_iterator->l;
		results_next = l_next - k_next;

		current_block = pos / block_size;

		if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated
			last_block_pos = end;
			complete_search = true;

		} else {

			if (current_block < last_block) { //Not in last previsited block

				if (pos % block_size) { //Not in first element of the block
					last_block_pos = (current_block+2) * block_size - 1;
				} else { //I am in the first element in the block (all the block must be processed)
					last_block_pos = (current_block+1) * block_size - 1;
				}

				complete_search = false;

			} else { //I am in the last previsited block

				if (pos % block_size) { //I am not in the first element of the block
					last_block_pos = end;
					complete_search = true;
				} else { //I am in the first element in the block (all the block must be processed)
					last_block_pos = (current_block+1) * block_size - 1;
					complete_search = false;
				}

			}

		}

		for(int16_t i=pos; i<=last_block_pos; i++) {

			k = k_next;
			l = l_next;

			if (k > l) break;

			BWiteration(k, l, k_next, l_next, W[i], index);
			results      = results_next;
			results_next = l_next - k_next;

			if (results == results_next) continue;

			change_result(r_iterator, k, l, i);
			add_result(r_iterator, rl_next);

		}

		if (complete_search && k_next <= l_next) {
			change_result(r_iterator, k_next, l_next, end+1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchPartialResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next) {

	intmax_t k, l, k_aux, l_aux;
	int16_t start, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 0) continue;

		start = r_iterator->start;
		pos   = r_iterator->pos;

		if (pos < start) {
			add_result(r_iterator, rl_next);
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos-1);
					add_result(r_iterator, rl_next);
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					//Insertion
					if (b!=W[last_err_pos]) {
						change_result(r_iterator, k_aux, l_aux, pos);
						modify_last_mismatch2(r_iterator, INSERTION, b);
						add_result(r_iterator, rl_next);
					}

					//Mismatch
					if (b!=W[pos]) {
						change_result(r_iterator, k_aux, l_aux, pos-1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						add_result(r_iterator, rl_next);
					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos-1);
				add_result(r_iterator, rl_next);

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos-1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							add_result(r_iterator, rl_next);
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				add_result(r_iterator, rl_next);

				//Mismatch
				if (b!=W[pos]) {

					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos-1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						add_result(r_iterator, rl_next);
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos-1);
			add_result(r_iterator, rl_next);

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				add_result(r_iterator, rl_next);

				//Mismatch
				if (b!=W[pos]) {
					r_iterator->pos = pos-1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					add_result(r_iterator, rl_next);
				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchPartialResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next) {

	intmax_t k, l, k_aux, l_aux;
	int16_t end, pos=0, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 1) continue;

		end = r_iterator->end;
		pos = r_iterator->pos;

		if (pos > end) {
			add_result(r_iterator, rl_next);
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION


			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos+1);
					add_result(r_iterator, rl_next);
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					//Insertion
					if (b!=W[last_err_pos]) {
						change_result(r_iterator, k_aux, l_aux, pos);
						modify_last_mismatch2(r_iterator, INSERTION, b);
						add_result(r_iterator, rl_next);
					}

					//Mismatch
					if (b!=W[pos]) {
						change_result(r_iterator, k_aux, l_aux, pos+1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						add_result(r_iterator, rl_next);
					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos+1);
				add_result(r_iterator, rl_next);

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					//Mismatch
					if (b!=W[pos]) {

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos+1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							add_result(r_iterator, rl_next);
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				add_result(r_iterator, rl_next);

				//Mismatch
				if (b!=W[pos]) {

					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos+1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						add_result(r_iterator, rl_next);
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos+1);
			add_result(r_iterator, rl_next);

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				add_result(r_iterator, rl_next);

				if (b!=W[pos]) { //Mismatch
					r_iterator->pos = pos+1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					add_result(r_iterator, rl_next);
				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

void BWChangeDirectionBackward(bwt_index *backward, bwt_index *forward, results_list *r_list, int16_t end) {

	result *r_iterator;

	for (uintmax_t ii=r_list->num_results; ii > 0; ii--) {
		r_iterator = &r_list->list[ii-1];

		if (r_iterator->dir == 1) break;

		if (r_iterator->pos == (r_iterator->start-1)) {
			change_direction(backward, forward, r_iterator);
			r_iterator->pos = r_iterator->end+1;
			r_iterator->end = end;
			r_iterator->dir = 1;
		}

	}

}

void BWChangeDirectionForward(bwt_index *backward, bwt_index *forward, results_list *r_list, int16_t start) {

	result *r_iterator;

	for (uintmax_t ii=r_list->num_results; ii > 0; ii--) {
		r_iterator = &r_list->list[ii-1];

		if (r_iterator->dir == 0) break;

		if (r_iterator->pos == (r_iterator->end+1)) {

			change_direction(backward, forward, r_iterator);

			r_iterator->pos = r_iterator->start-1;
			r_iterator->start = start;
			r_iterator->dir = 0;

		}

	}

}

bool BWSearchCPU(uint8_t *W, uint64_t nW, bwt_index *backward, bwt_index *forward, results_list *rl_prev, results_list *rl_next, results_list *rl_prev_i, results_list *rl_next_i, results_list *rl_final, int16_t fragsize, bool type) {

	result r;

	int16_t fragments = nW / fragsize;
	int16_t half = fragments / 2;
	if (fragments % 2) half++;
	int err_count;

	bool flow;

	//printf("\n***** Size: %d, Fragments: %d Errors: %d\n", nW, fragments, fragments-1);

	if (fragments <= 1) {

		if (type) {
			init_result(&r, 0);
			change_result(&r, 0, size_SA(backward)-1, nW-1);
			bound_result(&r, 0, nW-1);
			BWExactSearchBackward(W, backward, &r);
			if (r.k<=r.l)
				add_result(&r, rl_final);
		} else {
			init_result(&r, 1);
			change_result(&r, 0, size_SA(forward)-1, 0);
			bound_result(&r, 0, nW-1);
			BWExactSearchForward(W, forward, &r);
			if (r.k<=r.l)
				add_result(&r, rl_final);
		}

		return false;

	}

	//////////////////////////////FORWARD///////////////////////////////////////////

	for (int16_t i = half-1; i > 0; i--) {

		flow = true;

		err_count = fragments-1;

		rl_prev->num_results = 0; rl_prev_i->num_results = 0;
		rl_next->num_results = 0; rl_next_i->num_results = 0;

		init_result(&r, 1);
		change_result(&r, 0, size_SA(forward)-1, fragsize*i);
		bound_result(&r, fragsize*i, fragsize*(i+1) - 1);
		BWExactSearchForward(W, forward, &r);
		r.end = nW-1;

		if (r.k <= r.l) {

			add_result(&r, rl_prev);

			while (err_count > 0) {

				if (BWExactPartialResultsForward(W, forward, rl_prev, rl_next, rl_prev_i, fragsize, half-1)) {flow = false; break;}
				BWChangeDirectionForward(forward, backward, rl_prev_i, 0);
				if (BWExactPartialResultsBackward(W, backward, rl_prev_i, rl_next_i, rl_final, fragsize, half-1)) {flow = false; break;}
				if (err_count==1) break;
				if (BWBranchPartialResultsForward(W, forward, rl_next, rl_prev)) {flow = false; break;}
				if (BWBranchPartialResultsBackward(W, backward, rl_next_i, rl_prev_i)) {flow = false; break;}
				err_count--;

			}

			if (flow) {
				BWBranchFinalResultsForward(W, forward, rl_next, rl_prev_i, fragsize, half-1);
				BWChangeDirectionForward(forward, backward, rl_prev_i, 0);
				BWExactFinalResultsBackward(W, backward, rl_prev_i, rl_final, fragsize, half-1);

				BWBranchFinalResultsBackward(W, backward, rl_next_i, rl_final, fragsize, half-1);
			}

		}

	}

	///////BLOCK 0/////////////////////////////////////
	//printf("\n****BLOCK %d****\n", 0);

	flow = true;

	err_count = fragments-1;

	rl_prev->num_results = 0; rl_prev_i->num_results = 0;
	rl_next->num_results = 0; rl_next_i->num_results = 0;

	init_result(&r, 1);
	change_result(&r, 0, size_SA(forward)-1, 0);
	bound_result(&r, 0, fragsize - 1);
	BWExactSearchForward(W, forward, &r);
	r.end = nW-1;

	if (r.k <= r.l) {

		add_result(&r, rl_prev);

		while (err_count > 0) {
			if (BWExactPartialResultsForward(W, forward, rl_prev, rl_next, rl_final, fragsize, half-1)) {flow = false; break;}
			if (err_count==1) break;
			if (BWBranchPartialResultsForward(W, forward, rl_next, rl_prev)) {flow = false; break;}
			err_count--;
		}

		if (flow) {
			BWBranchFinalResultsForward(W, forward, rl_next, rl_final, fragsize, half-1);
		}

	}

	//////////////////////////////BACKWARD///////////////////////////////////////////

	for (int16_t i = half; i<fragments-1; i++) {

		flow = true;

		/* printf("\n****BLOCK %d****\n", i); */
		err_count = fragments-1;

		rl_prev->num_results = 0; rl_prev_i->num_results = 0;
		rl_next->num_results = 0; rl_next_i->num_results = 0;

		init_result(&r, 0);
		change_result(&r, 0, size_SA(backward)-1, fragsize*(i+1) - 1);
		bound_result(&r, fragsize*i, fragsize*(i+1) - 1);
		BWExactSearchBackward(W, backward, &r);
		r.start = 0;

		if (r.k <= r.l) {

			add_result(&r, rl_prev);

			while (err_count > 0) {
				if (BWExactPartialResultsBackward(W, backward, rl_prev, rl_next, rl_prev_i, fragsize, 0)) {flow = false; break;}
				BWChangeDirectionBackward(backward, forward, rl_prev_i, nW-1);
				if (BWExactPartialResultsForward(W, forward, rl_prev_i, rl_next_i, rl_final, fragsize, 0)) {flow = false; break;}
				if (err_count==1) break;
				if (BWBranchPartialResultsBackward(W, backward, rl_next, rl_prev)) {flow = false; break;}
				if (BWBranchPartialResultsForward(W, forward, rl_next_i, rl_prev_i)) {flow = false; break;}
				err_count--;
			}

			if (flow) {
				BWBranchFinalResultsBackward(W, backward, rl_next, rl_prev_i, fragsize, 0);
				BWChangeDirectionBackward(backward, forward, rl_prev_i, nW-1);
				BWExactFinalResultsForward(W, forward, rl_prev_i, rl_final, fragsize, 0);

				BWBranchFinalResultsForward(W, forward, rl_next_i, rl_final, fragsize, 0);
			}

		}

	}

	///////BLOCK FRAGMENTS-1/////////////////////////////////////

	/* printf("\n****BLOCK %d****\n", fragments-1); */

	flow = true;

	err_count = fragments-1;

	rl_prev->num_results = 0; rl_prev_i->num_results = 0;
	rl_next->num_results = 0; rl_next_i->num_results = 0;

	init_result(&r, 0);
	change_result(&r, 0, size_SA(backward)-1, /*fragsize*fragments - 1 Last block is larger*/nW-1);
	bound_result(&r, fragsize*(fragments-1), /*fragsize*fragments - 1 Last block is larger*/nW-1);
	BWExactSearchBackward(W, backward, &r);
	r.start = 0;

	if (r.k <= r.l) {

		add_result(&r, rl_prev);

		while (err_count > 0) {
			if (BWExactPartialResultsBackward(W, backward, rl_prev, rl_next, rl_final, fragsize, 0)) {flow = false; break;}
			if(err_count==1) break;
			if (BWBranchPartialResultsBackward(W, backward, rl_next, rl_prev)) {flow = false; break;}
			err_count--;
		}

		if (flow) {
			BWBranchFinalResultsBackward(W, backward, rl_next, rl_final, fragsize, 0);
		}

	}

	return false;

}
