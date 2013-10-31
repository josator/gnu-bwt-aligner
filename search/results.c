#include "results.h"

bool write_results(results_list *r_list, SA_TYPE *k, SA_TYPE *l, exome* ex, comp_vector *S, comp_vector *Si, vector *C, comp_matrix *O, comp_matrix *Oi, char *mapping, SA_TYPE nW, bool type, FILE *fp) {

  result *r;

  bool found = false;

  char search[MAXLINE+1];

  search[0] = '\0';
  strncat(search, mapping, nW);

	bool repeated;
  uintmax_t kl_count = 0;

  for (uintmax_t i=0;i<r_list->num_results; i++) {

		r = &r_list->list[i];

		repeated = false;
		uintmax_t kl;
		for (kl=0; kl < kl_count; kl++) {

			if(r->k == k[kl] && r->l == l[kl]) {
				repeated = true;
				break;
			}

		}

		if (!repeated) {
			k[kl] = r->k;
			l[kl] = r->l;
			kl_count++;
			manage_single_result(r, ex, S, Si, C, O, Oi, search, type, fp, r_list->read_index, &found);
		}

	}

	return found;

}
