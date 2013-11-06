#ifndef _SEARCH_SEARCH_
#define _SEARCH_SEARCH_

#include<inttypes.h>

#include "results.h"
#include "runtime.h"

inline void BWExactSearchBackward(uint8_t *W, bwt_index *index, result *r) {

	intmax_t k2, l2;
	int16_t i;

	k2 = r->k;
	l2 = r->l;

	printf("B1º -> %lu - %lu\n", k2, l2);

	for(i=r->pos; i>=r->start; i--) {

		BWiteration(k2, l2, k2, l2, W[i], index);
		printf("B -> %d -> %lu - %lu\n", i, k2, l2);
		if (k2 > l2) break;

	}

	r->k = k2;
	r->l = l2;
	r->pos = i;
}

inline void BWExactSearchForward(uint8_t *W, bwt_index *index, result *r) {

	intmax_t k2, l2;
	int16_t i;

	k2 = r->k;
	l2 = r->l;

	//printf("F1º -> %lu - %lu\n", k2, l2);

	for(i=r->pos; i<=r->end; i++) {

		BWiteration(k2, l2, k2, l2, W[i], index);
		//printf("F-> %d -> %lu - %lu\n", i, k2, l2);
		if (k2 > l2) break;

	}

	r->k = k2;
	r->l = l2;

	r->pos = i;

}

inline bool BWExactFinalResultBackward(uint8_t *W, bwt_index *index, result *r_iterator, results_list *rl_final, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t start, pos;
	int16_t current_block;

	start = r_iterator->start;
	pos   = r_iterator->pos;

	k = r_iterator->k;
	l = r_iterator->l;

	current_block = pos / block_size;

	if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated

	} else {

		if (current_block > last_block) { //Not in last previsited block

			return false;

		} else { //I am in the last previsited block

			if ((pos + 1) % block_size) { //I am not in the first element of the block
			} else { //I am in the first element in the block (all the block must be processed)
				return false;
			}

		}

	}

	for(int16_t i=pos; i>=start; i--) {
		BWiteration(k, l, k, l, W[i], index);
		if (k > l) break;
	}

	if (k <= l) {
		change_result(r_iterator, k, l, start-1);
		add_result(r_iterator, rl_final);
	}

	return false;

}

inline bool BWExactFinalResultForward(uint8_t *W, bwt_index *index, result *r_iterator, results_list *rl_final, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t pos, end;
	int16_t current_block;

	pos   = r_iterator->pos;
	end   = r_iterator->end;

	k = r_iterator->k;
	l = r_iterator->l;

	current_block = pos / block_size;

	if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated

	} else {

		if (current_block < last_block) { //Not in last previsited block

			return 0;

		} else { //I am in the last previsited block

			if (pos % block_size) { //I am not in the first element of the block
			} else { //I am in the first element in the block (all the block must be processed)
				return 0;
			}

		}

	}

	for(int16_t i=pos; i<=end; i++) {
		BWiteration(k, l, k, l, W[i], index);
		if (k > l) break;
	}

	if (k <= l) {
		change_result(r_iterator, k, l, end+1);
		add_result(r_iterator, rl_final);
	}

	return 0;

}

//inline void change_direction(comp_vector *S, comp_vector *Ri, vector *C, comp_matrix *O, comp_matrix *Oi, result *res) {
//
//	SA_TYPE k, l, ki, li, aux, aux2;
//	int16_t start, end, err_offset;
//
//	k  = res->k;
//	l  = res->l;
//	ki = O->siz-2;
//	li = 0;
//
//	start = res->start;
//	end   = res->end;
//
//	err_offset=0;
//
//	for (uint8_t rr=0; rr<res->num_mismatches; rr++) {
//
//		if      (res->err_kind[rr]==DELETION)
//			err_offset--;
//		else if (res->err_kind[rr]==INSERTION)
//			err_offset++;
//
//	}
//
//	if (S->ratio == 1) {
//
//		for (SA_TYPE i = k; i <= l; i++) {
//
//			aux  = S->siz - S->vector[i] - (end - start + 2) - err_offset;
//			aux2 = Ri->vector[aux];
//
//			if (aux2 < ki) ki = aux2;
//			if (aux2 > li) li = aux2;
//
//		}
//
//	} else {
//
//		for (SA_TYPE i = k; i <= l; i++) {
//
//			aux  = S->siz - getScompValue(i, S, C, O) - (end - start + 2) - err_offset;
//			aux2 = getRcompValue(aux, Ri, C, Oi);
//
//			if (aux2 < ki) ki = aux2;
//			if (aux2 > li) li = aux2;
//
//		}
//
//	}
//
//	res->k = ki;
//	res->l = li;
//
//}
//
//bool BWSearchCPU(uint8_t *W, SA_TYPE nW, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, comp_vector *S, comp_vector *R, comp_vector *Si, comp_vector *Ri, results_list *rl_prev, results_list *rl_next, results_list *rl_prev_i, results_list *rl_next_i, results_list *rl_final, int16_t fragsize, bool type);
//
//bool BWSearch1GPUHelper(uint8_t *W, int16_t start, int16_t end, SA_TYPE *vec_k, SA_TYPE *vec_l, SA_TYPE *vec_ki, SA_TYPE *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list);
//
///******DEPRECATED FUNCTIONS******/
//bool BWSearch1CPU(uint8_t *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result *res, results_list *r_list);
//bool BWSimpleSearch1Backward(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list);
//bool BWSimpleSearch1Forward(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list);
//
//void BWExactSearchVectorBackward(uint8_t *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O);
//void BWExactSearchVectorForward(uint8_t *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O);
///********************************/
//
//bool nextFASTAToken(FILE *queries_file, char *uncoded, uint8_t *coded, SA_TYPE *nquery); 

#endif
