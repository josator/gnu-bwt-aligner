bool BWSearch1GPUHelper(uint8_t *W, int16_t start, int16_t end, SA_TYPE *vec_k, SA_TYPE *vec_l, SA_TYPE *vec_ki, SA_TYPE *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list) {

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

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);
			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;
			//printf("*W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			uint8_t b_w = W[start];

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

		for (uint8_t b=0;b<nA;b++) {

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

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

			if (_ki_aux > _li_aux) continue;

			uint8_t b_w = W[end];

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

		for (uint8_t b=0;b<nA;b++) {

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

bool BWSimpleSearch1Backward(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

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

		for (uint8_t b=0;b<nA;b++) {

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

bool BWSimpleSearch1Forward(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

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

		for (uint8_t b=0;b<nA;b++) {

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

bool BWSearch1CPU(uint8_t *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result *res, results_list *r_list) {

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

void BWExactSearchVectorBackward(uint8_t *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O) {

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

void BWExactSearchVectorForward(uint8_t *W, int16_t start, int16_t end, SA_TYPE k, SA_TYPE l, SA_TYPE *vec_k, SA_TYPE *vec_l, vector *C, vector *C1, comp_matrix *O) {

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
