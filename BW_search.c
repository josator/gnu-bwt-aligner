#include "BW_search.h"

bool BWExactFinalResultsBackward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l;
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
			BWiteration(k, l, k, l, W[i], C, C1, O);
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

bool BWExactFinalResultsForward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l;
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
			BWiteration(k, l, k, l, W[i], C, C1, O);
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

bool BWBranchFinalResultsBackward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_final, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l, k_aux, l_aux;
	int16_t start, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	REF_TYPE last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 0) continue;

		start = r_iterator->start;
		pos   = r_iterator->pos;

		if (pos < start) {
			if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
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
			last_err_base = (REF_TYPE) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos-1);
					if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
				}

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

					if (k_aux > l_aux) continue;

					//Insertion
					if (b!=W[last_err_pos]) {
						change_result(r_iterator, k_aux, l_aux, pos);
						modify_last_mismatch2(r_iterator, INSERTION, b);
						if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
					}

					//Mismatch
					if (b!=W[pos]) {
						change_result(r_iterator, k_aux, l_aux, pos-1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos-1);
				if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos-1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

				//Mismatch
				if (b!=W[pos]) {

					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos-1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos-1);
			if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

				//Mismatch
				if (b!=W[pos]) {
					r_iterator->pos = pos-1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					if (BWExactFinalResultBackward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchFinalResultsForward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_final, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l, k_aux, l_aux;
	int16_t end, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	REF_TYPE last_err_base;
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
			last_err_base = (REF_TYPE) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos+1);
					if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
				}

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

					if (k_aux > l_aux) continue;

					//Insertion
					if (b!=W[last_err_pos]) {
						change_result(r_iterator, k_aux, l_aux, pos);
						modify_last_mismatch2(r_iterator, INSERTION, b);
						if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
					}

					//Mismatch
					if (b!=W[pos]) {
						change_result(r_iterator, k_aux, l_aux, pos+1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos+1);
				if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos+1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

				//Mismatch
				if (b!=W[pos]) {

					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos+1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos+1);
			if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

				if (k_aux > l_aux) continue;

				//Insertion
				change_result(r_iterator, k_aux, l_aux, pos);
				modify_last_mismatch2(r_iterator, INSERTION, b);
				if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;

				if (b!=W[pos]) { //Mismatch
					r_iterator->pos = pos+1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					if (BWExactFinalResultForward(W, C, C1, O, r_iterator, rl_final, block_size, last_block)) return true;
				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWExactPartialResultsBackward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l, k_next, l_next;
	int16_t start, pos, current_block, last_block_pos;
	bool complete_search;

	result *r_iterator;
	SA_TYPE results, results_next;

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

			BWiteration(k, l, k_next, l_next, W[i], C, C1, O);
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

bool BWExactPartialResultsForward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l, k_next, l_next;
	int16_t pos, end, current_block, last_block_pos;
	bool complete_search;

	result *r_iterator;
	SA_TYPE results, results_next;

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

			BWiteration(k, l, k_next, l_next, W[i], C, C1, O);
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

bool BWBranchPartialResultsBackward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

	SA_TYPE k, l, k_aux, l_aux;
	int16_t start, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	REF_TYPE last_err_base;
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
			last_err_base = (REF_TYPE) -1;
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

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

bool BWBranchPartialResultsForward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

	SA_TYPE k, l, k_aux, l_aux;
	int16_t end, pos=0, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	REF_TYPE last_err_base;
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
			last_err_base = (REF_TYPE) -1;
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

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

				for (REF_TYPE b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

			for (REF_TYPE b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

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

void BWChangeDirectionBackward(comp_vector *S, comp_vector *Ri, vector *C, comp_matrix *O, comp_matrix *Oi, results_list *r_list, int16_t end) {

	result *r_iterator;

	for (uintmax_t ii=r_list->num_results; ii > 0; ii--) {
		r_iterator = &r_list->list[ii-1];

		if (r_iterator->dir == 1) break;

		if (r_iterator->pos == (r_iterator->start-1)) {
			change_direction(S, Ri, C, O, Oi, r_iterator);
			r_iterator->pos = r_iterator->end+1;
			r_iterator->end = end;
			r_iterator->dir = 1;
		}

	}

}

void BWChangeDirectionForward(comp_vector *S, comp_vector *Ri, vector *C, comp_matrix *O, comp_matrix *Oi, results_list *r_list, int16_t start) {

	result *r_iterator;

	for (uintmax_t ii=r_list->num_results; ii > 0; ii--) {
		r_iterator = &r_list->list[ii-1];

		if (r_iterator->dir == 0) break;

		if (r_iterator->pos == (r_iterator->end+1)) {

			change_direction(S, Ri, C, O, Oi, r_iterator);

			r_iterator->pos = r_iterator->start-1;
			r_iterator->start = start;
			r_iterator->dir = 0;

		}

	}

}

bool BWSearchCPU(REF_TYPE *W, SA_TYPE nW, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, comp_vector *S, comp_vector *R, comp_vector *Si, comp_vector *Ri, results_list *rl_prev, results_list *rl_next, results_list *rl_prev_i, results_list *rl_next_i, results_list *rl_final, int16_t fragsize, bool type) {

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
			change_result(&r, 0, O->siz-2, nW-1);
			bound_result(&r, 0, nW-1);
			BWExactSearchBackward(W, C, C1, O, &r);
			if (r.k<=r.l)
				add_result(&r, rl_final);
		} else {
			init_result(&r, 1);
			change_result(&r, 0, O->siz-2, 0);
			bound_result(&r, 0, nW-1);
			BWExactSearchForward(W, C, C1, Oi, &r);
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
		change_result(&r, 0, O->siz-2, fragsize*i);
		bound_result(&r, fragsize*i, fragsize*(i+1) - 1);
		BWExactSearchForward(W, C, C1, Oi, &r);
		r.end = nW-1;

		if (r.k <= r.l) {

			add_result(&r, rl_prev);

			while (err_count > 0) {

				if (BWExactPartialResultsForward(W, C, C1, Oi, rl_prev, rl_next, rl_prev_i, fragsize, half-1)) {flow = false; break;}
				BWChangeDirectionForward(Si, R, C, Oi, O, rl_prev_i, 0);
				if (BWExactPartialResultsBackward(W, C, C1, O, rl_prev_i, rl_next_i, rl_final, fragsize, half-1)) {flow = false; break;}
				if (err_count==1) break;
				if (BWBranchPartialResultsForward(W, C, C1, Oi, rl_next, rl_prev)) {flow = false; break;}
				if (BWBranchPartialResultsBackward(W, C, C1, O, rl_next_i, rl_prev_i)) {flow = false; break;}
				err_count--;

			}

			if (flow) {
				BWBranchFinalResultsForward(W, C, C1, Oi, rl_next, rl_prev_i, fragsize, half-1);
				BWChangeDirectionForward(Si, R, C, Oi, O, rl_prev_i, 0);
				BWExactFinalResultsBackward(W, C, C1, O, rl_prev_i, rl_final, fragsize, half-1);

				BWBranchFinalResultsBackward(W, C, C1, O, rl_next_i, rl_final, fragsize, half-1);
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
	change_result(&r, 0, O->siz-2, 0);
	bound_result(&r, 0, fragsize - 1);
	BWExactSearchForward(W, C, C1, Oi, &r);
	r.end = nW-1;

	if (r.k <= r.l) {

		add_result(&r, rl_prev);

		while (err_count > 0) {
			if (BWExactPartialResultsForward(W, C, C1, Oi, rl_prev, rl_next, rl_final, fragsize, half-1)) {flow = false; break;}
			if (err_count==1) break;
			if (BWBranchPartialResultsForward(W, C, C1, Oi, rl_next, rl_prev)) {flow = false; break;}
			err_count--;
		}

		if (flow) {
			BWBranchFinalResultsForward(W, C, C1, Oi, rl_next, rl_final, fragsize, half-1);
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
		change_result(&r, 0, O->siz-2, fragsize*(i+1) - 1);
		bound_result(&r, fragsize*i, fragsize*(i+1) - 1);
		BWExactSearchBackward(W, C, C1, O, &r);
		r.start = 0;

		if (r.k <= r.l) {

			add_result(&r, rl_prev);

			while (err_count > 0) {
				if (BWExactPartialResultsBackward(W, C, C1, O, rl_prev, rl_next, rl_prev_i, fragsize, 0)) {flow = false; break;}
				BWChangeDirectionBackward(S, Ri, C, O, Oi, rl_prev_i, nW-1);
				if (BWExactPartialResultsForward(W, C, C1, Oi, rl_prev_i, rl_next_i, rl_final, fragsize, 0)) {flow = false; break;}
				if (err_count==1) break;
				if (BWBranchPartialResultsBackward(W, C, C1, O, rl_next, rl_prev)) {flow = false; break;}
				if (BWBranchPartialResultsForward(W, C, C1, Oi, rl_next_i, rl_prev_i)) {flow = false; break;}
				err_count--;
			}

			if (flow) {
				BWBranchFinalResultsBackward(W, C, C1, O, rl_next, rl_prev_i, fragsize, 0);
				BWChangeDirectionBackward(S, Ri, C, O, Oi, rl_prev_i, nW-1);
				BWExactFinalResultsForward(W, C, C1, Oi, rl_prev_i, rl_final, fragsize, 0);

				BWBranchFinalResultsForward(W, C, C1, Oi, rl_next_i, rl_final, fragsize, 0);
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
	change_result(&r, 0, O->siz-2, /*fragsize*fragments - 1 Last block is larger*/nW-1);
	bound_result(&r, fragsize*(fragments-1), /*fragsize*fragments - 1 Last block is larger*/nW-1);
	BWExactSearchBackward(W, C, C1, O, &r);
	r.start = 0;

	if (r.k <= r.l) {

		add_result(&r, rl_prev);

		while (err_count > 0) {
			if (BWExactPartialResultsBackward(W, C, C1, O, rl_prev, rl_next, rl_final, fragsize, 0)) {flow = false; break;}
			if(err_count==1) break;
			if (BWBranchPartialResultsBackward(W, C, C1, O, rl_next, rl_prev)) {flow = false; break;}
			err_count--;
		}

		if (flow) {
			BWBranchFinalResultsBackward(W, C, C1, O, rl_next, rl_final, fragsize, 0);
		}

	}

	return false;

}

bool BWSearch1GPUHelper(REF_TYPE *W, int16_t start, int16_t end, SA_TYPE *vec_k, SA_TYPE *vec_l, SA_TYPE *vec_ki, SA_TYPE *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list) {

	SA_TYPE _k, _l, _ki, _li, _k_aux, _l_aux, _ki_aux, _li_aux;
	SA_TYPE results, results_last;

	int16_t i, j, half, n;

	result r;

	n = end - start;
	half = n / 2;

	init_result(&r, 0);
	bound_result(&r, start, end);

	if (vec_k[0] <= vec_l[0]) {
		change_result(&r, vec_k[0], vec_l[0], -1);
		add_result(&r, r_list);
	}

	add_mismatch(&r, MATCH, -1, start);

	results = vec_l[0] - vec_k[0];

	results_last = results;
	_k  = vec_k[1];
	_l  = vec_l[1];
	results = _l  - _k;

	//printf("B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

	if (results != results_last) {

		//printf("*B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

		//printf("%d: %d -> %d\n", 0, results, results_last);

		//Deletion
		change_result(&r, _k, _l, -1);
		modify_last_mismatch3(&r, DELETION, -1, start);
		add_result(&r, r_list);

		for (REF_TYPE b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);
			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;
			//printf("*W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			REF_TYPE b_w = W[start];

			//Missmatch
			if (b!=b_w) {
				change_result(&r, _k_aux, _l_aux, -1);
				modify_last_mismatch2(&r, MISMATCH, b);
				add_result(&r, r_list);
			}

			//Insertion
			BWiteration(_k_aux, _l_aux, _k_aux, _l_aux, b_w, C, C1, O);

			if (_k_aux <= _l_aux) {
				change_result(&r, _k_aux, _l_aux, -1);
				modify_last_mismatch3(&r, 3, INSERTION, b);
				add_result(&r, r_list);
			}

		}

	}

	for (i=start+2, j=2; j<=half; i++, j++) {

		results_last = results;
		_k = vec_k[j];
		_l = vec_l[j];
		results = _l  - _k;

		//printf("B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

		if (results == results_last) continue;

		//printf("*B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

		//Deletion
		change_result(&r, _k, _l, i-2);
		modify_last_mismatch3(&r, DELETION, -1, i-1);
		BWExactSearchBackward(W, C, C1, O, &r);
		if (r.k<=r.l) add_result(&r, r_list);

		for (REF_TYPE b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i-1);
			modify_last_mismatch2(&r, INSERTION, b);
			BWExactSearchBackward(W, C, C1, O, &r);
			if (r.k<=r.l) add_result(&r, r_list);

			//Mismatch
			if (b!=W[i-1]) {
				change_result(&r, _k_aux, _l_aux, i-2);
				modify_last_mismatch1(&r, MISMATCH);
				BWExactSearchBackward(W, C, C1, O, &r);
				if (r.k<=r.l) add_result(&r, r_list);
			}

		}

	}

	//printf("\n");

	//TODO: Gestionar bien los errores de la busqueda al revés con Si y restas (precalcular Si con |X| - pos)
	half--;
	results = vec_li[n] - vec_ki[n];

	r.dir=1; //Change direction

	results_last = results;
	_ki  = vec_ki[n-1];
	_li  = vec_li[n-1];
	results = _li - _ki;

	//printf("F-> %d: %d -> %d, %u, %u\n", n, results, results_last, _ki, _li);

	if (results != results_last) {

		//printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

		//Deletion
		change_result(&r, _ki, _li, -1);
		modify_last_mismatch3(&r, DELETION, -1, end);
		add_result(&r, r_list);

		for (REF_TYPE b=0;b<nA;b++) {

			BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

			if (_ki_aux > _li_aux) continue;

			REF_TYPE b_w = W[end];

			//Mismatch
			if (b!=b_w) {
				change_result(&r, _ki_aux, _li_aux, -1);
				modify_last_mismatch2(&r, MISMATCH, b);
				add_result(&r, r_list);
			}

			//Insertion
			BWiteration(_ki_aux, _li_aux, _ki_aux, _li_aux, b_w, C, C1, Oi);

			//printf("\tI -> %d - %d\n", _ki_aux, _li_aux);

			if (_ki_aux <= _li_aux){
				change_result(&r, _ki_aux, _li_aux, -1);
				modify_last_mismatch2(&r, INSERTION, b);
				add_result(&r, r_list);
			}

		}

	}

	for(i=end-2,j=n-2; j>=half; i--, j--) {

		results_last = results;
		_ki  = vec_ki[j];
		_li  = vec_li[j];
		results = _li - _ki;

		//printf("F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

		if (results == results_last) continue;

		//printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

		//TODO: Anadir contador para podar cuando se que ya he encontrado las cadenas que permite la variabilidad en este simbolo.

		//Deletion
		change_result(&r, _ki, _li, i+2);
		modify_last_mismatch3(&r, DELETION, -1, i+1);
		BWExactSearchForward(W, C, C1, Oi, &r);
		if (r.k<=r.l) add_result(&r, r_list);

		for (REF_TYPE b=0;b<nA;b++) {

			BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

			//printf("W -> %d, %d - %d\n", b, _ki_aux, _li_aux);

			if (_ki_aux > _li_aux) continue;

			//Insertion
			change_result(&r, _ki_aux, _li_aux, i+1);
			modify_last_mismatch2(&r, INSERTION, b);
			BWExactSearchForward(W, C, C1, Oi, &r);
			if (r.k<=r.l) add_result(&r, r_list);

			//Mismatch
			if (b!= W[i+1]) {
				change_result(&r, _ki_aux, _li_aux, i+2);
				modify_last_mismatch1(&r, MISMATCH);
				BWExactSearchForward(W, C, C1, Oi, &r);
				if (r.k<=r.l) add_result(&r, r_list);
			}

		}

	}

	//printf("\n");

	return false;

}

///////////////////DEPRECATED FUNCTIONS/////////////////////

bool BWSimpleSearch1Backward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

	SA_TYPE _k,_l, _k_next, _l_next, _k_aux, _l_aux;
	SA_TYPE results, results_next;
	int16_t start, end, i;

	result r;

	start   = res->start;
	end     = res->pos;

	_k_next = res->k;
	_l_next = res->l;
	results_next = _l_next - _k_next;

	init_result(&r, 0);
	bound_result(&r, start, end);
	add_mismatch(&r, MATCH, -1, start);

	for(i=end; i>=start; i--) {

		_k = _k_next;
		_l = _l_next;

		//printf("%d:\n", i);

		if (_k > _l) {
			change_result(res, _k, _l, -1);
			return false;
		}

		BWiteration(_k, _l, _k_next, _l_next, W[i], C, C1, O);
		results      = results_next;
		results_next = _l_next - _k_next;
		//printf("(%lu, %lu, %lu)\t", results, _k, _l);

		if (results == results_next) continue;

		//Deletion
		change_result(&r, _k, _l, i-1);
		BWExactSearchBackward(W, C, C1, O, &r);
		if (r.k<=r.l) {
			modify_last_mismatch3(&r, DELETION, -1, i);
			add_result(&r, r_list);
		}

		for (REF_TYPE b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i);
			BWExactSearchBackward(W, C, C1, O, &r);
			if (r.k<=r.l) {
				modify_last_mismatch3(&r, INSERTION, b, i);
				add_result(&r, r_list);
			}

			//Mismatch
			if (b!=(int)W[i]) {
				change_result(&r, _k_aux, _l_aux, i-1);
				BWExactSearchBackward(W, C, C1, O, &r);
				if (r.k<=r.l) {
					modify_last_mismatch3(&r, MISMATCH, b, i);
					add_result(&r, r_list);
				}

			}

		}

	}

	//Match at exit in res
	change_result(res, _k_next, _l_next, -1);

	return false;

}

bool BWSimpleSearch1Forward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

	SA_TYPE _k, _l, _k_next, _l_next, _k_aux, _l_aux;
	SA_TYPE results, results_next;
	int16_t start, end, i;

	result r;

	start   = res->pos;
	end     = res->end;

	_k_next = res->k;
	_l_next = res->l;
	results_next = _l_next - _k_next;

	init_result(&r, 1);
	bound_result(&r, start, end);
	add_mismatch(&r, MATCH, -1, start);

	for(i=start; i<=end; i++) {

		_k = _k_next;
		_l = _l_next;

		//printf("%d:\n", i);

		if (_k > _l) {
			change_result(res, _k, _l, -1);
			return false;
		}

		BWiteration(_k, _l, _k_next, _l_next, W[i], C, C1, O);
		results      = results_next;
		results_next = _l_next - _k_next;
		if (results == results_next) continue;

		//Deletion
		change_result(&r, _k, _l, i+1);
		BWExactSearchForward(W, C, C1, O, &r);
		if (r.k<=r.l) {
			modify_last_mismatch3(&r, DELETION, -1, i);
			add_result(&r, r_list);
		}

		for (REF_TYPE b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i);
			BWExactSearchForward(W, C, C1, O, &r);
			if (r.k<=r.l) {
				modify_last_mismatch3(&r, INSERTION, b, i);
				add_result(&r, r_list);
			}

			//Mismatch
			if (b!=(int)W[i]) {
				change_result(&r, _k_aux, _l_aux, i+1);
				BWExactSearchForward(W, C, C1, O, &r);
				if (r.k<=r.l) {
					modify_last_mismatch3(&r, MISMATCH, b, i);
					add_result(&r, r_list);
				}
			}

		}

	}

	//Match at exit in res
	change_result(res, _k_next, _l_next, -1);

	return false;

}

bool BWSearch1CPU(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result *res, results_list *r_list) {

	int16_t start, end, half, n;;
	SA_TYPE _k, _l;

	_k = res->k;
	_l = res->l;

	start = res->start;
	end   = res->end;

	n = end - start + 1;
	half = n / 2;

	result r;

	init_result(&r, 0);
	bound_result(&r, half, end);
	change_result(&r, _k, _l, end);

	BWExactSearchBackward(W, C, C1, O, &r);

	if (r.k <= r.l) {
		r.start = start;
		r.pos = half-1;
		BWSimpleSearch1Backward(W, C, C1, O, &r, r_list);

		if (r.k <= r.l) add_result(&r, r_list); //Match
	}

	half--;

	init_result(&r, 1);
	bound_result(&r, start, half);
	change_result(&r, _k, _l, start);

	BWExactSearchForward(W, C, C1, Oi, &r);
	if (r.k <= r.l) {
		r.pos = half+1;
		r.end = end;
		BWSimpleSearch1Forward(W, C, C1, Oi, &r, r_list);
	}

	return false;

}

void BWExactSearchVectorBackward(REF_TYPE *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O) {

	if (k > l)       return;
	if (start > end) return;

	SA_TYPE k2, l2;
	int16_t last, i, j;

	last = end-start;

	k2 = k;
	l2 = l;

	for(i=end, j=last; i>=start; i--, j--) {

		BWiteration(k2, l2, k2, l2, W[i], C, C1, O);

		vec_k[j] = k2;
		vec_l[j] = l2;

		if (k2 > l2) {
			i--; j--;
			break;
		}

	}

	for(;i>=start; i--, j--) {
		vec_k[j] = k2;
		vec_l[j] = l2;
	}

}

void BWExactSearchVectorForward(REF_TYPE *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O) {

	if (k > l) return;
	if (start > end) return;

	SA_TYPE k2, l2;
	int16_t i, j;

	k2 = k;
	l2 = l;

	for(i=start, j=0; i<=end; i++, j++) {

		BWiteration(k2, l2, k2, l2, W[i], C, C1, O);

		vec_k[j] = k2;
		vec_l[j] = l2;

		if (k2 > l2) {
			i++; j++;
			break;
		}

	}

	for(; i<=end; i++, j++) {
		vec_k[j] = k2;
		vec_l[j] = l2;
	}

}

bool nextFASTAToken(FILE *queries_file, char *uncoded, REF_TYPE *coded, SA_TYPE *nquery, uint8_t *compressed, SA_TYPE* ncompress) {

	char line[MAXLINE];
	size_t length=0;

	*nquery=0;
	uncoded[0]='\0';

	while ( fgets(line, MAXLINE, queries_file) ) {

		if (line[0] == '>') {
			if (*nquery) break;
			else continue;
		}

		length=strlen(line);
		if (line[length-1]=='\n')
			length--;

		uncoded[*nquery] = '\0';
		strncpy(uncoded + *nquery, line , length);

		*nquery += length;

	}

	if (*nquery) {

		encode_bases(coded, uncoded, *nquery);

		if (compressed != NULL)
			*ncompress = comp4basesInByte(coded, *nquery, compressed);

		return true;

	} else {

		return false;

	}

}