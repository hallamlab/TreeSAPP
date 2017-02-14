#include "util_lib.h"

void viterbi(HMM *hmm_ptr, char *O, char* output_buffer, char* aa_buffer,
             char *dna_buffer, char *sequence_head, int whole_genome, int format,
             int len_seq, char* dna, char* dna1, char* dna_f, char* dna_f1,
             char* protein, int* insert, int* c_delete, char* temp_str_ptr) {

    int *vpath;                          // optimal path after backtracking
    int **path;                          // viterbi path array
    double **alpha;                      // viterbi prob array
    int i, j, t, kk;
    int from, from0, to;   /*from0: i-2 position, from: i-1 position */
    int from2;             /* from2: i-2, i-1 for condition in probability */
    int gene_len;
    int num_d;          		/* the number of delete */
    int freq_id;
    double h_kd, r_kd, p_kd;
    double temp_alpha, prob;
    double start_freq;
    double final_score;

    int codon_start = 0;
    int dna_id = 0;
    int dna_f_id = 0;
    int out_nt;
    int start_t = -1;
    int end_t;
    int prev_match;
    int start_orf;
    int frame;
    int insert_id, delete_id;
    int temp_i[6]   = {0,0,0,0,0,0};
    int temp_i_1[6] = {1,1,1,1,1,1};
    int num_N = 0;

    /***************************************************************/
    /* initialize                                                  */
    /***************************************************************/

    for(i =0; i < strlen(O); i++) {
        if( O[i]=='a' ) O[i]='A';
        if( O[i]=='t' ) O[i]='T';
        if( O[i]=='c' ) O[i]='C';
        if( O[i]=='g' ) O[i]='G';
    }

    if (whole_genome==1) {
        gene_len = 120;
    } else {
        gene_len = 60;
    }

    alpha = (double **)dmatrix(len_seq);
    path = (int **)imatrix(len_seq);
    vpath = (int *)ivector(len_seq);

    for (i = 0; i < NUM_STATE; i++) {
        alpha[i][0] = -1 * (hmm_ptr->pi[i]);
    }

    /* stop state */
    if ((O[0] == 'T') && (((O[1] == 'A') && (O[2] == 'A')) ||
                          ((O[1] == 'A') && (O[2] == 'G')) || ((O[1] == 'G') && (O[2] == 'A')))) {

        alpha[E_STATE][0] = max_dbl;
        alpha[E_STATE][1] = max_dbl;
        path[E_STATE][1] = E_STATE;
        path[E_STATE][2] = E_STATE;

        alpha[M6_STATE][2] = max_dbl;
        alpha[M5_STATE][1] = max_dbl;
        alpha[M4_STATE][0] = max_dbl;
        alpha[M3_STATE][2] = max_dbl;
        alpha[M2_STATE][1] = max_dbl;
        alpha[M1_STATE][0] = max_dbl;

        if ((O[1] == 'A') && (O[2] == 'A')) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - LOG_53;
        } else if ((O[1] == 'A') && (O[2] == 'G')) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - LOG_16;
        } else if ((O[1] == 'G') && (O[2] == 'A')) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - LOG_30;
        }
    }

    if ((O[2] == 'A') &&
            (((O[0] == 'T') && (O[1] == 'T')) ||
             ((O[0] == 'C') && (O[1] == 'T')) ||
             ((O[0] == 'T') && (O[1] == 'C')))) {
        alpha[S_STATE_1][0] = max_dbl;
        alpha[S_STATE_1][1] = max_dbl;
        alpha[S_STATE_1][2] = alpha[S_STATE][0];
        path[S_STATE_1][1] = S_STATE_1;
        path[S_STATE_1][2] = S_STATE_1;

        alpha[M3_STATE_1][2] = max_dbl;
        alpha[M6_STATE_1][2] = max_dbl;

        if ((O[0] == 'T') && (O[1] == 'T')) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - LOG_53;
        } else if ((O[0] == 'C') && (O[1] == 'T')) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - LOG_16;
        } else if ((O[0] == 'T') && (O[1] == 'C')) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - LOG_30;
        }
    }

    int multiple=0;
    /******************************************************************/
    /*  fill out the rest of the columns                              */
    /******************************************************************/
    for (t = 1; t < len_seq; t++) {
        from = nt2int(O[t-1]);
        if (t>1) {
            from0 = nt2int(O[t-2]);
        } else {
            from0 = 2;
        }
        to = nt2int(O[t]);

        /* if DNA is other than ACGT, do it later */
        if (from == 4) {
            from = 2;
        }
        if (from0 == 4) {
            from0 = 2;
        }
        if (to == 4) {
            to = 2;
            num_N += 1;
        } else {
            num_N = 0;
        }
        from2 = from0 * 4 + from;

        /******************/
        /* M state        */
        /******************/

        for (i=M1_STATE; i<=M6_STATE; i++)   {
            if (alpha[i][t]<max_dbl) {
                if (t==0) {
                } else {
                    if (i==M1_STATE) {
                        /* from M state */
                        j = M6_STATE;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) - (hmm_ptr->tr[TR_MM]) - (hmm_ptr->e_M[0][from2][to]);
                        path[i][t] = j;

                        /* from D state */
                        if (whole_genome==0) {
                            for (j=M5_STATE; j>=M1_STATE; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1<i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if(num_d>0) {
                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M[0][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }

                        /* from Start state */
                        temp_alpha = alpha[S_STATE][t-1] - (hmm_ptr->e_M[0][from2][to]);
                        if ( temp_alpha < alpha[i][t] ) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = S_STATE;
                        }

                    } else {  /*i ==M2-M6*/

                        /* from M state */
                        j = i - 1;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_MM]) -
                                      (hmm_ptr->e_M[i-M1_STATE][from2][to]);
                        path[i][t] = j;


                        /* from D state */
                        if (whole_genome == 0) {
                            for (j=M6_STATE; j>=M1_STATE; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1 < i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d > 0) {


                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M[i-M1_STATE][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }
                    }

                    /* from I state */
                    if (i == M1_STATE) {
                        j = I6_STATE;
                    } else {
                        j = I1_STATE + (i - M1_STATE -1);
                    }


                    /* to aviod stop codon */
                    if (t<2) {
                    } else if((i==M2_STATE || i==M5_STATE) && (O[temp_i[j-I1_STATE]] == 'T') &&
                              (((O[t] == 'A') && (O[t+1] =='A')) ||
                               ((O[t] == 'A') && (O[t+1] =='G')) ||
                               ((O[t] == 'G') && (O[t+1] =='A')))) {

                    } else if (((j-I1_STATE > 0) && (temp_i[j-I1_STATE]>0)) && ((i==M3_STATE || i==M6_STATE) && (O[temp_i[j-I1_STATE]-1] == 'T') &&
                               (((O[temp_i[j-I1_STATE]] == 'A') && (O[t] =='A')) ||
                                ((O[temp_i[j-I1_STATE]] == 'A') && (O[t] =='G')) ||
                                ((O[temp_i[j-I1_STATE]] == 'G') && (O[t] =='A'))))) {
                    } else {
                        temp_alpha = alpha[j][t-1]  - (hmm_ptr->tr[TR_IM]) - LOG_25;
                        if ( temp_alpha < alpha[i][t]) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;
                        }
                    }
                }
            }
        }

        /******************/
        /* I state        */
        /******************/
        for (i=I1_STATE; i<=I6_STATE; i++) {

            if (t==0) {
            } else {

                /* from I state */
                j = i;
                alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_II]) -
                              (hmm_ptr->tr_I_I[from][to]);
                path[i][t] = j;

                /* from M state */
                j = i - I1_STATE + M1_STATE ;
                if (i == I6_STATE) {
                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) -
                                 (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                } else {
                    temp_alpha = alpha[j][t-1]  -
                                 (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                }
                if (temp_alpha < alpha[i][t]) {
                    alpha[i][t] = temp_alpha;
                    path[i][t] = j;

                    temp_i[i-I1_STATE] = t-1;
                }
            }
        }

        /******************/
        /* M' state        */
        /******************/

        for (i = M1_STATE_1; i <= M6_STATE_1; i++)   {
            if  ((i==M1_STATE_1 || i==M4_STATE_1)&& t>=3 &&
                    (((O[t-3] == 'T') && (O[t-2] == 'T') && (O[t-1] == 'A')) ||
                     ((O[t-3] == 'C') && (O[t-2] == 'T') && (O[t-1] == 'A')) ||
                     ((O[t-3] == 'T') && (O[t-2] == 'C') && (O[t-1] == 'A')))) {

                /* from Start state  since this is actually stop codon in minus strand */
                alpha[i][t] = alpha[S_STATE_1][t-1] -
                              (hmm_ptr->e_M_1[i-M1_STATE_1][from2][to]);
                path[i][t] = S_STATE_1;

            } else {

                if (t==0) {
                } else {

                    if (i==M1_STATE_1 ) {

                        /* from M state */
                        j = M6_STATE_1;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) -
                                      (hmm_ptr->tr[TR_MM]) - (hmm_ptr->e_M_1[0][from2][to]);
                        path[i][t] = j;

                        /* from D state */
                        if (whole_genome==0) {
                            for (j=M5_STATE_1; j>=M1_STATE_1; j--) {
                                if (j >= i) {
                                    num_d = i-j+6;
                                } else if (j+1 <i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d > 0) {
                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M_1[0][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }

                    } else {

                        /* from M state */
                        j = i - 1;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_MM]) -
                                      (hmm_ptr->e_M_1[i-M1_STATE_1][from2][to]);
                        path[i][t] = j;

                        /* from D state */
                        if (whole_genome == 0) {
                            for (j=M6_STATE_1; j>=M1_STATE_1; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1 < i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d>0) {
                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M_1[i-M1_STATE_1][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }
                    }

                    /* from I state */
                    if (i == M1_STATE_1) {
                        j = I6_STATE_1;
                    } else {
                        j = I1_STATE_1 + (i - M1_STATE_1 -1);
                    }


                    //!! 11 --> 16
                    //!! What is the actual point of -1

                    /* to aviod stop codon */
                    if (t<2) {
                    } else  if((i==M2_STATE_1 ||
                                i==M5_STATE_1) &&
                               (O[t+1] == 'A' ) &&
                               (((O[temp_i_1[j-I1_STATE_1]] == 'T') && (O[t] =='T')) ||
                                ((O[temp_i_1[j-I1_STATE_1]] == 'C') && (O[t] =='T')) ||
                                ((O[temp_i_1[j-I1_STATE_1]] == 'T') && (O[t] =='C')))) {

                    } else if ((i==M3_STATE_1 ||
                                i==M6_STATE_1) &&
                               (O[t] == 'A' ) &&
                               (((O[temp_i_1[j-I1_STATE_1]-1] == 'T') && (O[temp_i_1[j-I1_STATE_1]] =='T')) ||
                                ((O[temp_i_1[j-I1_STATE_1]-1] == 'C') && (O[temp_i_1[j-I1_STATE_1]] =='T')) ||
                                ((O[temp_i_1[j-I1_STATE_1]-1] == 'T') && (O[temp_i_1[j-I1_STATE_1]] =='C')))) {
                    } else {

                        temp_alpha = alpha[j][t-1]  - (hmm_ptr->tr[TR_IM]) - LOG_25;
                        if ( temp_alpha < alpha[i][t]) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;
                        }
                    }
                }
            }
        }

        /******************/
        /* I' state        */
        /******************/
        for (i=I1_STATE_1; i<=I6_STATE_1; i++) {
            if (t == 0) {
            } else {
                /* from I state */
                j = i;
                alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_II]) -
                              (hmm_ptr->tr_I_I[from][to]);
                path[i][t] = j;

                if(t<5) continue;
                /* from M state */
                if (path[S_STATE_1][t-3]!= R_STATE && path[S_STATE_1][t-4] !=R_STATE &&
                        path[S_STATE_1][t-5] !=R_STATE) {
                    j = i - I1_STATE_1 + M1_STATE_1;
                    if (i==I6_STATE_1) {
                        temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) -
                                     (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                    } else {
                        temp_alpha = alpha[j][t-1]  -
                                     (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                    }
                    if (temp_alpha < alpha[i][t]) {
                        alpha[i][t] = temp_alpha;
                        path[i][t] = j;
                        //!! We are addressing the character array with this.
                        temp_i_1[i-I1_STATE_1] = t-1;
                    }
                }
            }
        }

        /***********************/
        /* Non_coding state    */
        /***********************/

        if (t==0) {
        } else {
            alpha[R_STATE][t] = alpha[R_STATE][t-1] - (hmm_ptr->tr_R_R[from][to]) -  (hmm_ptr->tr[TR_RR]);
            path[R_STATE][t] = R_STATE;

            temp_alpha = alpha[E_STATE][t-1]  - (hmm_ptr->tr[TR_ER])  ;
            if (temp_alpha < alpha[R_STATE][t] ) {
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE;
            }

            temp_alpha = alpha[E_STATE_1][t-1] - (hmm_ptr->tr[TR_ER]) ;
            if (temp_alpha < alpha[R_STATE][t] ) {
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE_1;
            }
            alpha[R_STATE][t] -= LOG_95;
        }

        /******************/
        /* END state      */
        /******************/
        if (alpha[E_STATE][t] == 0) {

            alpha[E_STATE][t] = max_dbl;
            path[E_STATE][t] = NOSTATE;

            if (t < len_seq -2 && (O[t] == 'T')  &&
                    (((O[t+1] == 'A') && (O[t+2] == 'A')) ||
                     ((O[t+1] == 'A') && (O[t+2] == 'G')) ||
                     ((O[t+1] == 'G') && (O[t+2] == 'A')))) {

                alpha[E_STATE][t+2] = max_dbl;
                /* transition from frame4,frame5,and frame6 */
                temp_alpha = alpha[M6_STATE][t-1] - (hmm_ptr->tr[TR_GE]);
                if (temp_alpha < alpha[E_STATE][t+2]) {
                    alpha[E_STATE][t+2] = temp_alpha;
                    path[E_STATE][t] = M6_STATE;
                }

                /* transition from frame1,frame2,and frame3 */
                temp_alpha  = alpha[M3_STATE][t-1] - (hmm_ptr->tr[TR_GE]);
                if (temp_alpha < alpha[E_STATE][t+2]) {
                    alpha[E_STATE][t+2] = temp_alpha;
                    path[E_STATE][t] = M3_STATE;
                }

                alpha[E_STATE][t] = max_dbl;
                alpha[E_STATE][t+1] = max_dbl;
                path[E_STATE][t+1] = E_STATE;
                path[E_STATE][t+2] = E_STATE;

                alpha[M6_STATE][t+2] = max_dbl;
                alpha[M5_STATE][t+1] = max_dbl;
                alpha[M4_STATE][t] = max_dbl;
                alpha[M3_STATE][t+2] = max_dbl;
                alpha[M2_STATE][t+1] = max_dbl;
                alpha[M1_STATE][t] = max_dbl;

                if ((O[t+1] == 'A') && (O[t+2] == 'A')) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - (0.54);
                } else if ((O[t+1] == 'A') && (O[t+2] == 'G')) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - LOG_16;
                } else if((O[t+1] == 'G') && (O[t+2] == 'A')) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - LOG_30;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;

                double sub_sum = 0;
                int sub_count = 0;

                if (t>=60) { /* bug reported by Yu-Wei */
                    for(i=-60; i<=-3; i++) {
                        if (t+i+2 < len_seq) {
                            start_freq -= (hmm_ptr->tr_E[i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                        }
                    }
                } else {
                    for(i=(-1*t); i<=-3; i++) {
                        if (t+i+2 < len_seq) {
                            sub_sum += (hmm_ptr->tr_E[i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                        }
                    }
                    sub_sum = sub_sum * 58 / (-3 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = hmm_ptr->E_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E_dist[1],2)/(2*pow(hmm_ptr->E_dist[0],2)));
                r_kd = hmm_ptr->E_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E_dist[4],2)/(2*pow(hmm_ptr->E_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log(p_kd);
            }
        }

        /*************************************************/
        /* START' state                                  */
        /* origianlly stop codon of genes in - strand    */
        /*************************************************/
        if (alpha[S_STATE_1][t] == 0) {

            alpha[S_STATE_1][t] = max_dbl;
            path[S_STATE_1][t] = NOSTATE;


            if (t<len_seq-2 && (O[t+2] == 'A') &&
                    (((O[t] == 'T') && (O[t+1] == 'T')) ||
                     ((O[t] == 'C') && (O[t+1] == 'T')) ||
                     ((O[t] == 'T') && (O[t+1] == 'C')))) {

                alpha[S_STATE_1][t] = max_dbl;
                path[S_STATE_1][t] = R_STATE;
                alpha[S_STATE_1][t+1] = max_dbl;
                alpha[S_STATE_1][t+2] = alpha[R_STATE][t-1] - (hmm_ptr->tr[TR_RS]);
                path[S_STATE_1][t+1] = S_STATE_1;
                path[S_STATE_1][t+2] = S_STATE_1;

                temp_alpha = alpha[E_STATE_1][t-1] - (hmm_ptr->tr[TR_ES]);
                if (temp_alpha < alpha[S_STATE_1][t+2]) {
                    alpha[S_STATE_1][t+2] = temp_alpha;
                    path[S_STATE_1][t] = E_STATE_1;
                }

                temp_alpha = alpha[E_STATE][t-1] - (hmm_ptr->tr[TR_ES1]);
                if (temp_alpha < alpha[S_STATE_1][t+2]) {
                    alpha[S_STATE_1][t+2] = temp_alpha;
                    path[S_STATE_1][t] = E_STATE;
                }

                alpha[M3_STATE_1][t+2] = max_dbl;
                alpha[M6_STATE_1][t+2] = max_dbl;

                if ((O[t] == 'T') && (O[t+1] == 'T')) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - (0.54);
                } else if ((O[t] == 'C') && (O[t+1] == 'T')) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - LOG_16;
                } else if((O[t] == 'T') && (O[t+1] == 'C')) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - LOG_30;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;
                for(i=3; i<=60; i++) {
                    if (t+i+2 < len_seq) {
                        start_freq -= (hmm_ptr->tr_S_1[i-3][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                    }
                }
                h_kd = hmm_ptr->S1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[1],2)/(2*pow(hmm_ptr->S1_dist[0],2)));
                r_kd = hmm_ptr->S1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[4],2)/(2*pow(hmm_ptr->S1_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log(p_kd);
            }
        }

        /************************/
        /* START state          */
        /************************/
        if (alpha[S_STATE][t] == 0) {

            alpha[S_STATE][t] = max_dbl;
            path[S_STATE][t] = NOSTATE;

            if (t<len_seq-2 &&  (O[t+1] == 'T') && (O[t+2] == 'G')&&
                    ((O[t] == 'A') || (O[t] == 'G') ||  (O[t] == 'T'))) {

                alpha[S_STATE][t] = max_dbl;
                alpha[S_STATE][t+1] = max_dbl;
                alpha[S_STATE][t+2] = alpha[R_STATE][t-1] - (hmm_ptr->tr[TR_RS]);
                path[S_STATE][t] = R_STATE;
                path[S_STATE][t+1] = S_STATE;
                path[S_STATE][t+2] = S_STATE;

                temp_alpha = alpha[E_STATE][t-1] - (hmm_ptr->tr[TR_ES]);
                if (temp_alpha < alpha[S_STATE][t+2]) {
                    alpha[S_STATE][t+2] = temp_alpha;
                    path[S_STATE][t] = E_STATE;
                }

                temp_alpha = alpha[E_STATE_1][t-1] - (hmm_ptr->tr[TR_ES1]);
                if (temp_alpha < alpha[S_STATE][t+2]) {
                    alpha[S_STATE][t+2] = temp_alpha;
                    path[S_STATE][t] = E_STATE_1;
                }


                if (O[t] == 'A') {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - LOG_83;
                } else if (O[t] == 'G') {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - LOG_10;
                } else if(O[t] == 'T') {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - LOG_07;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;

                double sub_sum = 0;
                int sub_count = 0;

                if (t>=30) {
                    for(i=-30; i<=30; i++) {
                        if (t+i+2 < len_seq) {
                            start_freq -= (hmm_ptr->tr_S[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                        }
                    }
                } else {
                    for(i=(-1*t); i<=30; i++) {
                        if (t+i+2 < len_seq) {
                            sub_sum += (hmm_ptr->tr_S[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                        }
                    }
                    sub_sum = sub_sum * 61 / (30 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = hmm_ptr->S_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S_dist[1],2)/(2*pow(hmm_ptr->S_dist[0],2)));
                r_kd = hmm_ptr->S_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S_dist[4],2)/(2*pow(hmm_ptr->S_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(p_kd);

            }
        }

        /**********************************************/
        /* END' state                                 */
        /* origianlly start codon of genes in - strand */
        /**********************************************/
        if (alpha[E_STATE_1][t] == 0) {

            alpha[E_STATE_1][t] = max_dbl;
            path[E_STATE_1][t] = NOSTATE;

            if (t < len_seq - 2 && (O[t] == 'C') && (O[t+1] == 'A') &&
                    ((O[t+2] == 'T') || (O[t+2] == 'C') || (O[t+2] == 'A'))) {

                /* transition from frame6 */
                alpha[E_STATE_1][t+2] = alpha[M6_STATE_1][t-1] - (hmm_ptr->tr[TR_GE]);
                path[E_STATE_1][t] = M6_STATE_1;
                alpha[E_STATE_1][t] = max_dbl;
                alpha[E_STATE_1][t+1] = max_dbl;
                path[E_STATE_1][t+1] = E_STATE_1;
                path[E_STATE_1][t+2] = E_STATE_1;

                if (O[t+2] == 'T') {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - LOG_83;
                } else if (O[t+2] == 'C' ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - LOG_10;
                } else if(O[t+2] == 'A' ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - LOG_07;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;

                double sub_sum = 0;
                int sub_count = 0;

                if (t>=30) {
                    for(i=-30; i<=30; i++) {
                        if (t+i+2 < len_seq) {
                            start_freq -= (hmm_ptr->tr_E_1[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                        }
                    }
                } else {
                    for(i=(-1*t); i<=30; i++) {
                        if (t+i+2 < len_seq) {
                            sub_sum += (hmm_ptr->tr_E_1[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])]);
                        }
                    }
                    sub_sum = sub_sum * 61 / (30 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = hmm_ptr->E1_dist[2] *
                       exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)))
                       ;
                r_kd = hmm_ptr->E1_dist[5] *
                       exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)))
                       ;
                p_kd = h_kd / (h_kd + r_kd);

                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(p_kd);
            }
        }
        if (num_N>9) {

            for (i=0; i<NUM_STATE; i++) {
                if (i!=R_STATE) {
                    alpha[i][t] = max_dbl;
                    path[i][t] = R_STATE;
                }
            }
        }
    }




    /***********************************************************/
    /* backtrack array to find the optimal path                */
    /***********************************************************/


    sprintf(output_buffer, "%s\n", sequence_head);

    /* find the state for O[N] with the highest probability */
    prob = max_dbl;
    for (i = 0; i < NUM_STATE; i++) {

        if (alpha[i][len_seq-1] < prob) {
            prob = alpha[i][len_seq-1];
            vpath[len_seq-1] = i;
        }
    }

    /* backtrack the opitmal path */
    for(t=len_seq-2; t>=0; t--) {
        if(t+1 < 0 || vpath[t+1] < 0) continue;
        vpath[t] = path[vpath[t+1]][t+1];
    }

    for (t=0; t<len_seq; t++) {

        if (codon_start == 0 && start_t < 0 &&
                ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
                 (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1) ||
                 vpath[t] == S_STATE || vpath[t] == S_STATE_1 )) {
            start_t=t+1;
        }

        if (codon_start == 0 &&
                (vpath[t]==M1_STATE || vpath[t]==M4_STATE ||
                 vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)) {

            stopMemset(dna, STRINGLEN);
            stopMemset(dna1, STRINGLEN);//
            stopMemset(dna_f, STRINGLEN);//
            stopMemset(dna_f1, STRINGLEN);//
            stopMemset(protein, STRINGLEN);
            stopMemset(insert, STRINGLEN);//
            stopMemset(c_delete, STRINGLEN);//

            insert_id = 0;
            delete_id = 0;
            dna_id = 0;
            dna_f_id = 0;
            dna[dna_id] = O[t];
            dna_f[dna_f_id] = O[t];
            start_orf = t+1;
            prev_match = vpath[t];

            if (vpath[t] < M6_STATE) {
                codon_start = 1;
            } else {
                codon_start = -1;
            }

        } else if (codon_start != 0 && (vpath[t] == E_STATE || vpath[t] == E_STATE_1 || t == len_seq-1)) {

            if (vpath[t] == E_STATE || vpath[t] == E_STATE_1) {
                end_t = t+3;
            } else {
                end_t = t+1;

                /* FGS1.12 start: remove incomplete codon */
                int temp_t = t;
                while(vpath[temp_t] != M1_STATE && vpath[temp_t] != M4_STATE  &&
                        vpath[temp_t] != M1_STATE_1  && vpath[temp_t] != M4_STATE_1) {
                    dna_f[dna_f_id] = '\0';
                    dna_f_id--;

                    dna[dna_id] = '\0';
                    dna_id--;

                    temp_t--;
                }
                /* FGS1.12 end: remove incomplete codon */
            }
            final_score = (alpha[vpath[end_t-4]][end_t-4] - alpha[vpath[start_t+2]][start_t+2] )/(end_t-start_t-5);
            frame = start_orf%3;
            if (frame==0) {
                frame=3;
            }

            //!! Transfer all of the output buffer writing code to another function. Modularize this.

            if (dna_id > gene_len) {
                print_outputs(codon_start, start_t, end_t, frame, output_buffer, aa_buffer, dna_buffer, sequence_head,
                              dna, dna1, dna_f, dna_f1, protein, insert, c_delete, insert_id, delete_id, format, temp_str_ptr,multiple);
                multiple++;
            }

            codon_start = 0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        } else if (codon_start != 0 &&
                   ((vpath[t] >= M1_STATE && vpath[t] <= M6_STATE) ||
                    (vpath[t] >= M1_STATE_1 && vpath[t] <= M6_STATE_1)) &&
                   vpath[t] - prev_match < 6) {

            if (vpath[t] < prev_match) {
                out_nt = vpath[t]+6-prev_match;
            } else {
                out_nt = vpath[t]-prev_match;
            }
            for (kk=0; kk<out_nt; kk++) {  /* for deleted nt in reads */
                dna_id ++;
                dna[dna_id] = 'N';
                dna_f_id++;
                dna_f[dna_f_id] = 'x';
                if (kk>0) {
                    c_delete[delete_id] = t+1;
                    delete_id++;
                }
            }
            dna[dna_id] = O[t];
            dna_f[dna_f_id] = O[t];
            prev_match = vpath[t];

        } else if (codon_start != 0 &&
                   ((vpath[t] >= I1_STATE && vpath[t] <= I6_STATE) ||
                    (vpath[t] >= I1_STATE_1 && vpath[t] <= I6_STATE_1))) {
            dna_f_id ++;
            dna_f[dna_f_id] = tolower(O[t]);
            insert[insert_id] = t+1;
            insert_id++;

        } else if (codon_start != 0 && vpath[t] == R_STATE) {
            /* for long NNNNNNNNN, pretend R state */
            codon_start = 0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        }
    }
    free_dmatrix(alpha);
    free_imatrix(path);
    free(vpath);

    vpath = 0;
    dna = 0;
    dna1 = 0;
    dna_f = 0;
    dna_f = 0;
    protein = 0;

}

void get_prob_from_cg(HMM *hmm_ptr, TRAIN *train_ptr, char *O, int len_seq) {

    int cg_id = -1;
    int cg_count=0;
    int i,j,k;

    for (i=0; i<len_seq; i++) {
        if ((O[i] == 'C'||O[i] =='c') || (O[i] == 'G'||O[i] == 'g') ) {
            cg_count++;
        }
    }

    cg_count = floor((cg_count*1.0/len_seq)*100)-26;

    if (cg_count < 0) {
        cg_count = 0;
    } else if (cg_count > 43) {
        cg_count = 43;
    }

    memcpy(hmm_ptr->e_M, train_ptr->trans[cg_count], sizeof(hmm_ptr->e_M));
    memcpy(hmm_ptr->e_M_1, train_ptr->rtrans[cg_count], sizeof(hmm_ptr->e_M_1));
    memcpy(hmm_ptr->tr_R_R, train_ptr->noncoding[cg_count], sizeof(hmm_ptr->tr_R_R));
    memcpy(hmm_ptr->tr_S, train_ptr->start[cg_count], sizeof(hmm_ptr->tr_S));
    memcpy(hmm_ptr->tr_E, train_ptr->stop[cg_count], sizeof(hmm_ptr->tr_E));
    memcpy(hmm_ptr->tr_S_1, train_ptr->start1[cg_count], sizeof(hmm_ptr->tr_S_1));
    memcpy(hmm_ptr->tr_E_1, train_ptr->stop1[cg_count], sizeof(hmm_ptr->tr_E_1));
    memcpy(hmm_ptr->S_dist, train_ptr->S_dist[cg_count], sizeof(hmm_ptr->S_dist));
    memcpy(hmm_ptr->E_dist, train_ptr->E_dist[cg_count], sizeof(hmm_ptr->E_dist));
    memcpy(hmm_ptr->S1_dist, train_ptr->S1_dist[cg_count], sizeof(hmm_ptr->S1_dist));
    memcpy(hmm_ptr->E1_dist, train_ptr->E1_dist[cg_count], sizeof(hmm_ptr->E1_dist));
}


//!! This function is for some strange reason definetly leaking memory. Need to find out why.
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename,
                         char *mfilename1, char *nfilename, char *sfilename, char *pfilename,
                         char *s1filename, char *p1filename, char *dfilename, TRAIN *train_ptr) {

    int i, j, k, p;
    double prob;
    FILE *fp, *fpm, *fpm1, *fpn, *fps, *fpp, *fps1, *fpp1, *fpd;

    char name[10];
    char head[20];
    char start[10];
    char end[10];

    /****************************************************/
    /* transition                                       */
    /****************************************************/
    fp = fopen (filename , "r");

    /* Transition */
    fscanf(fp, "%s", head);
    for (i=0; i<14; i++) {
        //!! This causes a memory leak for unknown reasons.
        fscanf(fp, "%s %lf", name, &prob);
        hmm_ptr->tr[tr2int(name)] = log(prob);
    }

    /* TransitionMI */
    fscanf(fp, "%s", head);
    for (i=0; i<16; i++) {
        fscanf(fp, "%s %s %lf\n", start, end, &prob);
        hmm_ptr->tr_M_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
    }

    /* TransitionII */
    fscanf(fp, "%s", head);
    for (i=0; i<16; i++) {
        fscanf(fp, "%s %s %lf", start, end, &prob);
        hmm_ptr->tr_I_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
    }

    /* PI */
    fscanf(fp, "%s", head);
    for (i=0; i<NUM_STATE; i++) {
        fscanf(fp, "%s %lf", name, &prob);
        hmm_ptr->pi[i] = log(prob);
    }
    fclose(fp);

    /****************************************************/
    /* M state transition                               */
    /****************************************************/
    fpm = fopen (mfilename , "r");
    for (p=0; p<44; p++) {                       /* cg */
        fscanf(fpm, "%s", head);
        for (i=0; i<6; i++) {                      /* period */
            for (j=0; j<16; j++) {                   /* condition */
                for (k=0; k<4; k++) {                  /* emission */
                    fscanf(fpm, "%lf", &prob);
                    train_ptr->trans[p][i][j][k] = log(prob);
                }
            }
        }
    }
    fclose(fpm);

    /****************************************************/
    /* M state_1 transition                             */
    /****************************************************/
    fpm1 = fopen (mfilename1 , "r");
    for (p=0; p<44; p++) {
        fscanf(fpm1, "%s", head);
        for (i=0; i<6; i++) {
            for (j=0; j<16; j++) {
                for (k=0; k<4; k++) {
                    fscanf(fpm1, "%lf", &prob);
                    train_ptr->rtrans[p][i][j][k] = log(prob);
                }
            }
        }
    }
    fclose(fpm1);

    /****************************************************/
    /* noncoding state  transition                      */
    /****************************************************/
    fpn = fopen (nfilename, "r");
    for (p=0; p<44; p++) {
        fscanf(fpn, "%s", head);
        for (j=0; j<4; j++) {
            for (k=0; k<4; k++) {
                fscanf(fpn, "%lf", &prob);
                train_ptr->noncoding[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpn);

    /****************************************************/
    /* start                                            */
    /****************************************************/
    fps = fopen (sfilename, "r");
    for (p=0; p<44; p++) {
        fscanf(fps, "%s", head);
        for (j=0; j<61; j++) {
            for (k=0; k<64; k++) {
                fscanf(fps, "%lf", &prob);
                train_ptr->start[p][j][k] = log(prob);
            }
        }
    }
    fclose(fps);

    /****************************************************/
    /* stop                                             */
    /****************************************************/
    fpp = fopen (sfilename, "r");
    for (p=0; p<44; p++) {
        fscanf(fpp, "%s", head);
        for (j=0; j<58; j++) {
            for (k=0; k<64; k++) {
                fscanf(fpp, "%lf", &prob);
                train_ptr->stop[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpp);

    /****************************************************/
    /* start1                                           */
    /****************************************************/
    fps1 = fopen (s1filename, "r");
    for (p=0; p<44; p++) {
        fscanf(fps1, "%s", head);
        for (j=0; j<58; j++) {
            for (k=0; k<64; k++) {
                fscanf(fps1, "%lf", &prob);
                train_ptr->start1[p][j][k] = log(prob);
            }
        }
    }
    fclose(fps1);

    /****************************************************/
    /* stop1                                            */
    /****************************************************/
    fpp1 = fopen (p1filename, "r");
    for (p=0; p<44; p++) {
        fscanf(fpp1, "%s", head);
        for (j=0; j<61; j++) {
            for (k=0; k<64; k++) {
                fscanf(fpp1, "%lf", &prob);
                train_ptr->stop1[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpp1);


    /****************************************************/
    /* pwm distribution                                 */
    /****************************************************/
    fpd = fopen (dfilename, "r");
    for (p=0; p<44; p++) {
        fscanf(fpd, "%s", head);
        for (k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->S_dist[p][k] = prob;
        }
        for (k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->E_dist[p][k] = prob;
        }
        for (k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->S1_dist[p][k] = prob;
        }
        for (k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->E1_dist[p][k] = prob;
        }
    }
    fclose(fpd);

}

void print_outputs(int codon_start, int start_t, int end_t, int frame, char* output_buffer, char* aa_buffer, char* dna_buffer,
                   char* sequence_head_short, char* dna, char* dna1, char* dna_f, char* dna_f1, char* protein,
                   int* insert, int* c_delete, int insert_id, int delete_id, int format, char* temp_str_ptr, unsigned int multiple) {

    int i;
    char tab[] = "\t";




    if (codon_start == 1) {

        //sprintf(temp_str_ptr, "%d\t%d\t+\t%d\t%lf\t", start_t, end_t, frame, final_score);
        sprintf(temp_str_ptr, "%d\t%d\t+\t%d\t", start_t, end_t, frame);
        strcat(output_buffer, temp_str_ptr);
        sprintf(temp_str_ptr, "I:");
        strcat(output_buffer, temp_str_ptr);


        for (i=0; i<insert_id; i++) {
            sprintf(temp_str_ptr, "%d,", insert[i]);
            strcat(output_buffer, temp_str_ptr);
        }

        sprintf(temp_str_ptr, "\tD:");
        strcat(output_buffer, temp_str_ptr);

        for (i=0; i<delete_id; i++) {
            sprintf(temp_str_ptr, "%d,", c_delete[i]);
            strcat(output_buffer, temp_str_ptr);
        }

        sprintf(temp_str_ptr, "\n");
        strcat(output_buffer, temp_str_ptr);

        sprintf(temp_str_ptr, "%s_%d_%d_+\n", sequence_head_short, start_t, end_t);

        if(multiple)  strcat(aa_buffer, tab);

        strcat(aa_buffer, temp_str_ptr);

        sprintf(temp_str_ptr, "%s_%d_%d_+\n", sequence_head_short, start_t, end_t);
        strcat(dna_buffer, temp_str_ptr);

        get_protein(dna,protein,1);
        sprintf(temp_str_ptr, "%s\n", protein);
        strcat(aa_buffer, temp_str_ptr);
        if (format == 0) {
            sprintf(temp_str_ptr, "%s\n", dna);
        } else if (format == 1) {
            sprintf(temp_str_ptr, "%s\n", dna_f);
        }
        strcat(dna_buffer, temp_str_ptr);

    } else if (codon_start == -1) {
        //sprintf(temp_str_ptr, "%d\t%d\t-\t%d\t%lf\t", start_t, end_t, frame, final_score);
        sprintf(temp_str_ptr, "%d\t%d\t-\t%d\t", start_t, end_t, frame);
        strcat(output_buffer, temp_str_ptr);
        sprintf(temp_str_ptr, "I:");
        strcat(output_buffer, temp_str_ptr);

        for (i = 0; i < insert_id; i++) {
            sprintf(temp_str_ptr, "%d,", insert[i]);
            strcat(output_buffer, temp_str_ptr);
        }

        sprintf(temp_str_ptr, "\tD:");
        strcat(output_buffer, temp_str_ptr);

        for (i=0; i<delete_id; i++) {
            sprintf(temp_str_ptr, "%d,", c_delete[i]);
            strcat(output_buffer, temp_str_ptr);
        }

        sprintf(temp_str_ptr, "\n");
        strcat(output_buffer, temp_str_ptr);

        sprintf(temp_str_ptr, "%s_%d_%d_-\n", sequence_head_short, start_t, end_t);

        if(multiple)  strcat(aa_buffer, tab);
        strcat(aa_buffer, temp_str_ptr);


        sprintf(temp_str_ptr, "%s_%d_%d_+\n", sequence_head_short, start_t, end_t);
        strcat(dna_buffer, temp_str_ptr);

        get_protein(dna,protein,1);
        sprintf(temp_str_ptr, "%s\n", protein);
        strcat(aa_buffer, temp_str_ptr);
        if (format == 0) {
            sprintf(temp_str_ptr, "%s\n", dna);
        } else if (format == 1) {
            sprintf(temp_str_ptr, "%s\n", dna_f);
        }
        strcat(dna_buffer, temp_str_ptr);

    } else if (codon_start == -1) {

        //sprintf(temp_str_ptr, "%d\t%d\t-\t%d\t%lf\t", start_t, end_t, frame, final_score);
        sprintf(temp_str_ptr, "%d\t%d\t-\t%d\t", start_t, end_t, frame);
        strcat(output_buffer, temp_str_ptr);
        sprintf(temp_str_ptr, "I:");
        strcat(output_buffer, temp_str_ptr);

        for (i = 0; i < insert_id; i++) {
            sprintf(temp_str_ptr, "%d,", insert[i]);
            strcat(output_buffer, temp_str_ptr);
        }

        sprintf(temp_str_ptr, "\tD:");
        strcat(aa_buffer, temp_str_ptr);
        sprintf(temp_str_ptr, "%s_%d_%d_-\n", sequence_head_short, start_t, end_t);
        strcat(dna_buffer, temp_str_ptr);

        get_protein(dna,protein,-1);
        get_rc_dna(dna, dna1);
        get_rc_dna_indel(dna_f, dna_f1);
        sprintf(temp_str_ptr, "%s\n", protein);
        if(multiple)  strcat(aa_buffer, tab);
        strcat(aa_buffer, temp_str_ptr);

        if (format == 0) {
            sprintf(temp_str_ptr, "%s\n", dna1);
        } else if (format == 1) {
            sprintf(temp_str_ptr, "%s\n", dna_f1);
        }
        strcat(dna_buffer, temp_str_ptr);
    }
}
