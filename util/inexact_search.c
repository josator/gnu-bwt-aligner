#include <pthread.h>

#include "BW_search.h"

#define MAX_READ_THREAD 1024

char *h_Worig, *h_Worig2, *h_Worig3;
REF_TYPE *h_We, *h_We2, *h_We3;
SA_TYPE *h_nWe, *h_nWe2, *h_nWe3;

char *read_Worig, *gpu_Worig, *store_Worig, *swap_Worig;
REF_TYPE *read_We, *gpu_We, *store_We, *swap_We;
SA_TYPE *read_nWe, *gpu_nWe, *store_nWe, *swap_nWe;

comp_matrix h_O, h_rO, h_Oi, h_rOi;

vector h_C, h_rC, h_C1, h_rC1;
comp_vector S, Si, R, Ri;

results_list rl_prev, rl_next, rl_prev_i, rl_next_i;

results_list rl_final[MAX_READ_THREAD], rl_final2[MAX_READ_THREAD];
results_list rl_final_r[MAX_READ_THREAD], rl_final_r2[MAX_READ_THREAD];
results_list *gpu_rl_final,  *store_rl_final,  *swap_rl_final;
results_list *gpu_rl_final_r,  *store_rl_final_r,  *swap_rl_final_r;

SA_TYPE *k, *l;

timeval t1, t2, t1_write, t2_write;
float t_gpu=0, t_total=0, t_write=0;
uintmax_t descartadas=0;

pthread_mutex_t gpu_time_thread, print_results;
pthread_cond_t start_time, start_time_write, stop_time, stop_time_write;

int stop=0, start=0, end=0;
int start_write=0, stop_write=0, end_write=0;
uintmax_t tam_read_gpu=0, tam_read_gpu2=0;

FILE *output_file;

exome ex;

int num_errors, fragsize, RESULTS;

void *writeResults(void *threadid) {

	uintmax_t w=0;
	uintmax_t tam_write=0;

	unsigned long int contador=0;
	unsigned long int total=0;

	while(1) {

		pthread_mutex_lock(&print_results);

		//printf("W -> Para\n");

		while(start_write!=1) {
			pthread_cond_wait(&start_time_write, &print_results);
		}

		start_write=0;

		//printf("W -> Sigue\n");

		pthread_mutex_unlock(&print_results);

		if (end_write) {

			fflush(output_file);
			fclose(output_file);

			//printf("W -> Saliendo\n");
			//printf("%lu founds of %ju -> found %.2f, discarded %ju\n", contador, total, contador * 100.0 / total, descartadas * 100 / total);
			pthread_exit(NULL);

		}

		//printf("W -> Swap buffers\n");

		tam_write = tam_read_gpu2;

		swap_Worig = gpu_Worig;
		gpu_Worig = store_Worig;
		store_Worig = swap_Worig;

		swap_We = gpu_We;
		gpu_We = store_We;
		store_We = swap_We;

		swap_nWe = gpu_nWe;
		gpu_nWe = store_nWe;
		store_nWe = swap_nWe;

		swap_rl_final = gpu_rl_final;
		gpu_rl_final = store_rl_final;
		store_rl_final = swap_rl_final;

		swap_rl_final_r = gpu_rl_final_r;
		gpu_rl_final_r = store_rl_final_r;
		store_rl_final_r = swap_rl_final_r;

		pthread_mutex_lock(&print_results);

		//printf("W -> Ordena seguir a los hilos de búsqueda GPU\n");

		stop_write = 1;

		pthread_cond_broadcast(&stop_time_write);

		pthread_mutex_unlock(&print_results);

		//printf("W -> Empieza a escribir %d resultados\n", tam_write);

		gettimeofday(&t1_write, NULL);

		bool found, found2;
		for (uintmax_t i=0; i < tam_write; i++) {

			found = false; found2 = false;

			if (store_rl_final[i].num_results) {
				write_results(store_rl_final + i, k, l, &ex, &S, &Si, &h_C, &h_O, &h_Oi, store_Worig + i*MAXLINE, store_nWe[i], true, output_file);
				found = true;
			}

			if (store_rl_final_r[i].num_results) {
				write_results(store_rl_final_r + i, k, l, &ex, &S, &Si, &h_C, &h_O, &h_Oi, store_Worig + i*MAXLINE, store_nWe[i], false, output_file);
				found2 = true;
			}

			if ((!found) && (!found2)) {
			} else {
				contador++;
			}

			total++;

		} //read

		w++;

		gettimeofday(&t2_write, NULL);
		t_gpu = (t2_write.tv_sec-t1_write.tv_sec)*1e6+(t2_write.tv_usec-t1_write.tv_usec);
		t_write += t_gpu;

		//printf("W -> Termina de escribir hasta el resultado %d\n", w*MAX_READ_THREAD);

	}

}

void *cpuSearch(void *threadid) {

	uintmax_t num_read = 0;
	//int tid=0;
  //tid = (long)threadid;i

	REF_TYPE *g_We;
	SA_TYPE *g_nWe;

	g_We  = (REF_TYPE*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(REF_TYPE));
	check_malloc(g_We,  "cpuSearch");

	g_nWe  = (SA_TYPE*)malloc(MAX_READ_THREAD * sizeof(SA_TYPE));
	check_malloc(g_nWe, "cpuSearch");

	while(1) {

		pthread_mutex_lock(&gpu_time_thread);
		//printf("%d -> Espera lectura\n", tid);

		while(!start) {
			pthread_cond_wait(&start_time, &gpu_time_thread);
		}

		start=0;

		//printf("%d -> Sigue\n", tid);

		pthread_mutex_unlock(&gpu_time_thread);

		if (end) {

			pthread_mutex_lock(&print_results);

			//printf("%d -> Ordena terminar al hilo de escritura\n", tid);

			start_write++;
			end_write=1;

			pthread_cond_signal(&start_time_write);

			pthread_mutex_unlock(&print_results);

			//printf("%d -> Saliendo\n", tid);

			pthread_exit(NULL);

		}

		//printf("%d -> Copia vectores\n", tid);

		tam_read_gpu2 = tam_read_gpu;

		for (uintmax_t i=0; i < tam_read_gpu2 * MAXLINE; i++) {
			g_We[i] = read_We[i];
		}

		for (uintmax_t i=0; i < tam_read_gpu2; i++) {
			g_nWe[i] = read_nWe[i];
		}

		//printf("%d -> Ordena seguir al hilo de lectura\n", tid);

		pthread_mutex_lock(&gpu_time_thread);

		stop++;
		if (stop==1)
			pthread_cond_signal(&stop_time);

		pthread_mutex_unlock(&gpu_time_thread);

		for(uintmax_t i=0; i<tam_read_gpu2; i++) {

			gpu_rl_final[i].num_results = 0;
			gpu_rl_final_r[i].num_results = 0;
			gpu_rl_final[i].read_index = num_read;
			gpu_rl_final_r[i].read_index = num_read;

			int n_count = 0;
			for (SA_TYPE nn=0; nn<g_nWe[i]; nn++) {
				if (g_We[i * MAXLINE + nn] == 'N') n_count++;
				if (n_count>num_errors) {
					descartadas++;
					break;
				}
			}

			if (n_count<=num_errors) {

				int aux_fragsize = g_nWe[i] / (num_errors + 1);
				if (aux_fragsize < fragsize)
					aux_fragsize = fragsize;

				BWSearchCPU(
						g_We + i * MAXLINE,
						g_nWe[i],
						&h_C,
						&h_C1,
						&h_O,
						&h_Oi,
						&S,
						&R,
						&Si,
						&Ri,
						&rl_prev,
						&rl_next,
						&rl_prev_i,
						&rl_next_i,
						gpu_rl_final + i,
						aux_fragsize,
						1
						);

				BWSearchCPU(
						g_We + i * MAXLINE,
						g_nWe[i],
						&h_rC,
						&h_rC1,
						&h_rOi,
						&h_rO,
						&Si,
						&Ri,
						&S,
						&R,
						&rl_prev,
						&rl_next,
						&rl_prev_i,
						&rl_next_i,
						gpu_rl_final_r + i,
						aux_fragsize,
						0
						);

			}

			num_read++;

		} //read

		pthread_mutex_lock(&print_results);

		//printf("%d -> Ordena continuar al hilo de escritura\n", tid);
		start_write++;

		pthread_cond_signal(&start_time_write);

		//printf("%d -> Espera escritura\n", tid);

		while(!stop_write) {
			pthread_cond_wait(&stop_time_write, &print_results);
		}

		stop_write=0;

		//printf("%d -> Sigue\n", tid);

		pthread_mutex_unlock(&print_results);

	}

}

int main(int argc, char **argv) {

	pthread_t exec_thread;
	pthread_t write_thread;
	pthread_attr_t attr;
	int rc;

	gettimeofday(&t1, NULL);

	FILE *queries_file;

	if (argc!=8) {
		fprintf(stderr, "Syntax:\n\t%s f_mappings d_transform f_output search_tree_size num_errors min_fragment nucleotides\n", argv[0]);
		return 1;
	}

	RESULTS = atoi(argv[4]);
	num_errors = atoi(argv[5]);
	fragsize = atoi(argv[6]);
	init_replace_table(argv[7]);

	read_vector(&h_C,     argv[2], "C");
	read_vector(&h_C1,    argv[2], "C1");
	read_comp_matrix(&h_O,  argv[2], "O");
	read_comp_matrix(&h_Oi, argv[2], "Oi");

	read_comp_vector(&S,   argv[2], "Scomp");
	read_comp_vector(&Si,  argv[2], "Scompi");
	read_comp_vector(&R,   argv[2], "Rcomp");
	read_comp_vector(&Ri,  argv[2], "Rcompi");

	reverse_strand_C(&h_rC, &h_C, &h_rC1, &h_C1);
	reverse_strand_O(&h_rO, &h_O);
	reverse_strand_O(&h_rOi, &h_Oi);

	h_Worig  = (char*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(char));
	check_malloc(h_Worig,  "main");
	h_Worig2 = (char*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(char));
	check_malloc(h_Worig2, "main");
	h_Worig3 = (char*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(char));
	check_malloc(h_Worig3, "main");

	h_We  = (REF_TYPE*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(REF_TYPE));
	check_malloc(h_We,  "main");
	h_We2 = (REF_TYPE*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(REF_TYPE));
	check_malloc(h_We2, "main");
	h_We3 = (REF_TYPE*)malloc(MAX_READ_THREAD * MAXLINE * sizeof(REF_TYPE));
	check_malloc(h_We3, "main");

	h_nWe  = (SA_TYPE*)malloc(MAX_READ_THREAD * sizeof(SA_TYPE));
	check_malloc(h_nWe, "main");
	h_nWe2 = (SA_TYPE*)malloc(MAX_READ_THREAD * sizeof(SA_TYPE));
	check_malloc(h_nWe2, "main");
	h_nWe3 = (SA_TYPE*)malloc(MAX_READ_THREAD * sizeof(SA_TYPE));
	check_malloc(h_nWe3, "main");

	new_results_list(&rl_prev, RESULTS); new_results_list(&rl_prev_i, RESULTS);
	new_results_list(&rl_next, RESULTS); new_results_list(&rl_next_i, RESULTS);

	for(int i=0; i<MAX_READ_THREAD; i++) {
		new_results_list(rl_final  + i, RESULTS);
		new_results_list(rl_final2 + i, RESULTS);
		new_results_list(rl_final_r  + i, RESULTS);
		new_results_list(rl_final_r2 + i, RESULTS);
	}

	k = (SA_TYPE*)malloc(RESULTS * sizeof(SA_TYPE));
	l = (SA_TYPE*)malloc(RESULTS * sizeof(SA_TYPE));

	read_Worig = h_Worig; gpu_Worig = h_Worig2; store_Worig = h_Worig3;
	read_We = h_We; gpu_We = h_We2; store_We = h_We3;
	read_nWe = h_nWe; gpu_nWe = h_nWe2; store_nWe = h_nWe3;

	gpu_rl_final = rl_final; store_rl_final = rl_final2;
	gpu_rl_final_r = rl_final_r; store_rl_final_r = rl_final_r2;

	queries_file = fopen(argv[1], "r");
	check_file_open(queries_file, argv[1]);
	output_file = fopen(argv[3], "w");
	check_file_open(output_file, argv[3]);

	pthread_mutex_init(&gpu_time_thread,  NULL);
	pthread_mutex_init(&print_results,  NULL);

	pthread_cond_init(&start_time, NULL);
	pthread_cond_init(&stop_time,  NULL);

	pthread_cond_init(&start_time_write, NULL);
	pthread_cond_init(&stop_time_write,  NULL);

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	rc = pthread_create(&exec_thread, &attr, cpuSearch, (void *)0);
	if (rc){
		fprintf(stderr, "ERROR; return code from pthread_create() on cpuSearch thread %d is %d\n", 0, rc);
		return 1;
	}

	rc = pthread_create(&write_thread, &attr, writeResults, (void *)1);
	if (rc){                                                                  
		fprintf(stderr, "ERROR; return code from pthread_create() on writeResults thread %d is %d\n", 1, rc);
		return 1;
	}

	//printf("L -> Antes de cargar el fichero del exoma\n");
	load_exome_file(&ex, argv[2]);

	//printf("L -> Antes del bucle de lecturas\n");

	bool exit=true;

	while (exit) {

		//printf("L -> Empieza a leer\n");

		while (exit) {

			exit = nextFASTAToken(queries_file, read_Worig + tam_read_gpu * MAXLINE, read_We + tam_read_gpu * MAXLINE, read_nWe + tam_read_gpu);

			tam_read_gpu++;

			if (tam_read_gpu == MAX_READ_THREAD) break;

		}

		//printf("L -> Termina de leer tam_read_gpu = %d, status(%d)\n", tam_read_gpu, exit);

		if (!exit) {
			tam_read_gpu--;
			if (tam_read_gpu == 0) break;
		}

		pthread_mutex_lock(&gpu_time_thread);

		//printf("L -> Digo a los GPU que copien\n");

		start = 1;

		pthread_cond_broadcast(&start_time);

		//printf("L -> Antes de parar\n");

		while(stop!=1) {
			pthread_cond_wait(&stop_time, &gpu_time_thread);
		}

		stop = 0;

		pthread_mutex_unlock(&gpu_time_thread);

		//printf("L -> Swap buffers\n");

		swap_Worig = read_Worig;
		read_Worig = gpu_Worig;
		gpu_Worig = swap_Worig;

		swap_We = read_We;
		read_We = gpu_We;
		gpu_We = swap_We;

		swap_nWe = read_nWe;
		read_nWe = gpu_nWe;
		gpu_nWe = swap_nWe;

		tam_read_gpu = 0;

	}

	//Terminar los hilos
	pthread_mutex_lock(&gpu_time_thread);

	//printf("L -> Ordena salir a los hilos\n");
	start=1;
	end=1;

	pthread_cond_broadcast(&start_time);
	pthread_mutex_unlock(&gpu_time_thread);

	pthread_join(exec_thread, NULL);
	pthread_join(write_thread, NULL);

	pthread_mutex_destroy(&gpu_time_thread);
	pthread_mutex_destroy(&print_results);

	pthread_cond_destroy(&start_time);
	pthread_cond_destroy(&start_time_write);

	pthread_cond_destroy(&stop_time);
	pthread_cond_destroy(&stop_time_write);

	pthread_attr_destroy(&attr);

	free(h_C.vector);
	free(h_C1.vector);
	free(h_rC.vector);
	free(h_rC1.vector);

	free_comp_matrix(&h_rO, &h_O);
	free_comp_matrix(&h_rOi, &h_Oi);
	free(S.vector);
	free(Si.vector);
	free(R.vector);
	free(Ri.vector);

	free(h_We);

	free(h_nWe);
	free(h_nWe2);
	free(h_nWe3);

	free(k);
	free(l);

	free(rl_prev.list);
	free(rl_next.list);
	free(rl_prev_i.list);
	free(rl_next_i.list);

	for(int i=0; i<MAX_READ_THREAD; i++) {
		free(rl_final[i].list);
		free(rl_final2[i].list);
		free(rl_final_r[i].list);
		free(rl_final_r2[i].list);
	}

	fclose(queries_file);

	gettimeofday(&t2, NULL);
	t_total = (t2.tv_sec-t1.tv_sec)*1e6+(t2.tv_usec-t1.tv_usec);

	//printf("Write: %f\n", t_write / 1000000);
	//printf("Total: %f\n", t_total / 1000000);

	//printf("L -> Saliendo\n");

	return 0;

}
