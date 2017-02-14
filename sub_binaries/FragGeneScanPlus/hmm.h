#ifndef __HMM_H
#define __HMM_H

#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <semaphore.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "fasta.h" // STRINGLEN 

#define DONE_Q 0
#define EMPTY_Q 1
#define max_dbl 10000000000.0

#define LOG_53 -0.63487827243
#define LOG_16 -1.83258146375
#define LOG_30 -1.20397280433
#define LOG_25 -1.38629436112
#define LOG_83 -0.18632957819
#define LOG_10 -2.30258509299
#define LOG_07 -2.65926003693
#define LOG_95 -0.05129329438

#define A 0
#define C 1
#define G 2
#define T 3

#define NUM_STATE 29

#define NOSTATE -1
#define S_STATE 0
#define E_STATE 1
#define R_STATE 2
#define S_STATE_1 3
#define E_STATE_1 4
#define M1_STATE 5
#define M2_STATE 6
#define M3_STATE 7
#define M4_STATE 8
#define M5_STATE 9
#define M6_STATE 10
#define M1_STATE_1 11
#define M2_STATE_1 12
#define M3_STATE_1 13
#define M4_STATE_1 14
#define M5_STATE_1 15
#define M6_STATE_1 16
#define I1_STATE 17
#define I2_STATE 18
#define I3_STATE 19
#define I4_STATE 20
#define I5_STATE 21
#define I6_STATE 22
#define I1_STATE_1 23
#define I2_STATE_1 24
#define I3_STATE_1 25
#define I4_STATE_1 26
#define I5_STATE_1 27
#define I6_STATE_1 28

#define TR_MM 0
#define TR_MI 1
#define TR_MD 2
#define TR_II 3
#define TR_IM 4
#define TR_DD 5
#define TR_DM 6
#define TR_GE 7
#define TR_GG 8
#define TR_ER 9
#define TR_RS 10
#define TR_RR 11
#define TR_ES 12
#define TR_ES1 13

char hmm_file[STRINGLEN];
char aa_file[STRINGLEN];
char seq_file[STRINGLEN];
char out_file[STRINGLEN];
char dna_file[STRINGLEN];
char train_file[STRINGLEN];
char mstate_file[STRINGLEN];
char rstate_file[STRINGLEN];
char nstate_file[STRINGLEN];
char sstate_file[STRINGLEN];
char pstate_file[STRINGLEN];
char s1state_file[STRINGLEN];     /* stop codon of gene in - stand */
char p1state_file[STRINGLEN];
char dstate_file[STRINGLEN];
char train_dir[STRINGLEN];

// semaphores
#ifdef __APPLE__
typedef sem_t* SEM_T;
#elif __linux
typedef sem_t SEM_T;
#define sem_wait(x) sem_wait(&x)
#define sem_post(x) sem_post(&x)
#endif

SEM_T sema_Q;
SEM_T sema_R;
SEM_T sema_r;
SEM_T sema_w;


typedef struct {

    double  pi[29];    /* pi[1..N] pi[i] is the initial state distribution. */

    double tr[14];                 /* transition probability from a (delete/insert/match) state to a state */

    double e_M_1[6][16][4];      /* transition probability from a lowest-level state  to a  lowest-level state*/
    double e_M[6][16][4];

    double tr_R_R[4][4];
    double tr_I_I[4][4];
    double tr_M_I[4][4];

    double tr_S[61][64];
    double tr_E[61][64];
    double tr_S_1[61][64];
    double tr_E_1[61][64];

    double S_dist[6];  /*sigma, mu,alpha, sigma_r, mu_r, alpha_r */
    double E_dist[6];
    double S1_dist[6];
    double E1_dist[6];
} HMM;

typedef struct {

    double trans[44][6][16][4];
    double rtrans[44][6][16][4];
    double noncoding[44][4][4];
    double start[44][61][64];
    double stop[44][61][64];
    double start1[44][61][64];
    double stop1[44][61][64];

    double S_dist[44][6];
    double E_dist[44][6];
    double S1_dist[44][6];
    double E1_dist[44][6];

} TRAIN;

typedef struct thread_data {
    unsigned int wholegenome;
    unsigned int format;
    HMM *hmm;

    unsigned int *output_num_sequences;
    unsigned int *input_num_sequences;
    unsigned int  id;

    // all buffers allocated by master
    char*** input_buffer;
    char*** input_head_buffer;
    char*** output_buffer;
    char*** aa_buffer;
    char*** dna_buffer;

    // new stuff
    char* dna;
    char* dna1;
    char* dna_f;
    char* dna_f1;
    char* protein;
    int* insert;
    int* c_delete;
    char* temp_str;

    //unsigned int** acceptable_buffer;

    SEM_T sema_r;
    SEM_T sema_w;
} thread_data;

thread_data *thread_datas;

HMM hmm;
TRAIN train;

void init_thread_data(thread_data* td);
off_t read_file_into_buffer(FILE* fp, int fpos, thread_data* thread_data, unsigned int buf);
int read_seq_into_buffer(FASTAFILE* fp, thread_data* thread_data, unsigned int buf);

void get_prob_from_cg(HMM *hmm, TRAIN *train, char *O, int len_seq);
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename,
                         char *sfilename,char *pfilename,char *s1filename,char *p1filename, char *dfilename, TRAIN *train_ptr);
void viterbi(HMM *hmm_ptr, char *O, char* output_buffer, char* aa_buffer, char *dna_buffer,
             char *sequence_head, int whole_genome, int format, int len_seq,
             char* dna_ptr, char* dna1_ptr, char* dna_f_ptr, char* dna_f1_ptr, char* protein_ptr,
             int* insert_ptr, int* c_delete_ptr, char* temp_str_ptr);

void free_hmm(HMM *hmm);
void get_protein(char *dna, char *protein, int strand);
void get_rc_dna(char *dna, char *dna1);
void get_rc_dna_indel(char* dna_f, char* dna_f1);
void get_corrected_dna(char *dna, char *dna_f);
void print_usage();
void free_thread_data(thread_data* td);
void print_outputs(int codon_start, int start_t, int end_t, int frame, char* output_buffer, char* aa_buffer,
                   char* dna_buffer, char* sequence_head_short, char* dna, char* dna1, char* dna_f, char* dna_f1,
                   char* protein, int* insert, int* c_delete, int insert_id, int delete_id, int format, char* temp_str_ptr, unsigned int multiple);

// helper functions to cleanup the main function
void setTrainDirectory(char* train_path);
void conductWork();
void initializeThreads();
void destroySemaphores();
void initializeSemaphores();
void setupProgram(int argc, char** argv);
#endif
