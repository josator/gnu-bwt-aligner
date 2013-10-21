#ifndef _SEARCH_SEARCH_
#define _SEARCH_SEARCH_

#include "commons/commons.h"
#include "commons/string_utils.h"

#include "commons/BW_types.h"
#include "BW_results.h"

#if defined FM_COMP_32 || FM_COMP_64

#define BWiteration(k_in,l_in, k_out, l_out, b, C, C1, O)\
	do {\
		(k_out) = (C1)->vector[(b)] + getOcompValue((b), (k_in)  , (O));\
		(l_out) = (C)->vector[(b)]  + getOcompValue((b), (l_in)+1, (O));\
	} while (0);
//printf("k-> %lu, l-> %lu, O(k) -> %u, O(l) -> %u, C -> %u, C1 -> %u\n", (k_out), (l_out), getOcompValue((b), (k_in), (O)), getOcompValue((b), (l_in)+1, (O)), (C)->vector[(b)], (C1)->vector[(b)]);
#else

#define BWiteration(k_in, l_in, k_out, l_out, b, C, C1, O)\
	do {\
		(k_out) = (C1)->vector[(b)] + (O)->desp[(b)][(k_in)];\
		(l_out) = (C)->vector[(b)]  + (O)->desp[(b)][(l_in)+1];\
	} while(0);
//printf("k-> %lu, l-> %lu, O(k) -> %u, O(l) -> %u, C -> %u, C1 -> %u\n", k_out, l_out, (O)->desp[(b)][(k_in)], (O)->desp[(b)][(l_in)+1], (C)->vector[(b)], (C1)->vector[(b)]);
#endif

inline void BWExactSearchBackward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *r) {

	SA_TYPE k2, l2;
	int16_t i;

	k2 = r->k;
	l2 = r->l;

	//printf("B1º -> %lu - %lu\n", k2, l2);

	for(i=r->pos; i>=r->start; i--) {

		BWiteration(k2, l2, k2, l2, W[i], C, C1, O);
		//printf("B -> %d -> %lu - %lu\n", i, k2, l2);
		if (k2 > l2) break;

	}

	r->k = k2;
	r->l = l2;
	r->pos = i;
}

inline void BWExactSearchForward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *Oi, result *r) {

	SA_TYPE k2, l2;
	int16_t i;

	k2 = r->k;
	l2 = r->l;

	//printf("F1º -> %lu - %lu\n", k2, l2);

	for(i=r->pos; i<=r->end; i++) {

		BWiteration(k2, l2, k2, l2, W[i], C, C1, Oi);
		//printf("F-> %d -> %lu - %lu\n", i, k2, l2);
		if (k2 > l2) break;

	}

	r->k = k2;
	r->l = l2;

	r->pos = i;

}

inline bool BWExactFinalResultBackward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *r_iterator, results_list *rl_final, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l;
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
		BWiteration(k, l, k, l, W[i], C, C1, O);
		if (k > l) break;
	}

	if (k <= l) {
		change_result(r_iterator, k, l, start-1);
		add_result(r_iterator, rl_final);
	}

	return false;

}

inline bool BWExactFinalResultForward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *r_iterator, results_list *rl_final, int16_t block_size, int16_t last_block) {

	SA_TYPE k, l;
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
		BWiteration(k, l, k, l, W[i], C, C1, O);
		if (k > l) break;
	}

	if (k <= l) {
		change_result(r_iterator, k, l, end+1);
		add_result(r_iterator, rl_final);
	}

	return 0;

}

inline void change_direction(comp_vector *S, comp_vector *Ri, vector *C, comp_matrix *O, comp_matrix *Oi, result *res) {

	SA_TYPE k, l, ki, li, aux, aux2;
	int16_t start, end, err_offset;

	k  = res->k;
	l  = res->l;
	ki = O->siz-2;
	li = 0;

	start = res->start;
	end   = res->end;

	err_offset=0;

	for (uint8_t rr=0; rr<res->num_mismatches; rr++) {

		if      (res->err_kind[rr]==DELETION)
			err_offset--;
		else if (res->err_kind[rr]==INSERTION)
			err_offset++;

	}

	if (S->ratio == 1) {

		for (SA_TYPE i = k; i <= l; i++) {

			aux  = S->siz - S->vector[i] - (end - start + 2) - err_offset;
			aux2 = Ri->vector[aux];

			if (aux2 < ki) ki = aux2;
			if (aux2 > li) li = aux2;

		}

	} else {

		for (SA_TYPE i = k; i <= l; i++) {

			aux  = S->siz - getScompValue(i, S, C, O) - (end - start + 2) - err_offset;
			aux2 = getRcompValue(aux, Ri, C, Oi);

			if (aux2 < ki) ki = aux2;
			if (aux2 > li) li = aux2;

		}

	}

	res->k = ki;
	res->l = li;

}

bool BWSearchCPU(REF_TYPE *W, SA_TYPE nW, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, comp_vector *S, comp_vector *R, comp_vector *Si, comp_vector *Ri, results_list *rl_prev, results_list *rl_next, results_list *rl_prev_i, results_list *rl_next_i, results_list *rl_final, int16_t fragsize, bool type);

bool BWSearch1GPUHelper(REF_TYPE *W, int16_t start, int16_t end, SA_TYPE *vec_k, SA_TYPE *vec_l, SA_TYPE *vec_ki, SA_TYPE *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list);

/******DEPRECATED FUNCTIONS******/
bool BWSearch1CPU(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result *res, results_list *r_list);
bool BWSimpleSearch1Backward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list);
bool BWSimpleSearch1Forward(REF_TYPE *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list);

void BWExactSearchVectorBackward(REF_TYPE *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O);
void BWExactSearchVectorForward(REF_TYPE *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O);
/********************************/

bool nextFASTAToken(FILE *queries_file, char *uncoded, REF_TYPE *coded, SA_TYPE *nquery); 

#endif
