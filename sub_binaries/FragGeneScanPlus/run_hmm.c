#include "hmm.h"
#include "util_lib.h"
#include "fasta.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "run_hmm.h"

unsigned int ids; // thread id
unsigned int threadnum = 1;
unsigned int format;
unsigned int wholegenome;
unsigned int output_dna;
unsigned int output_meta;
unsigned int verbose;
unsigned int MAX_BYTES_PER_BUFFER;
unsigned int MAX_SEQS_PER_BUFFER;
unsigned int num_reads_flag = 0;

int max_mem = 0;

pthread_t writer_thread;
long long round_counter;
off_t stopped_at_fpos; // tracks how far we've read in the input

long writer_counter = 0;
long read_counter = 0;
long read_counter1 = 0;
long work_counter = 0;
long total_reads = -1;
long viterbi_counter = 0;
long num_writes = 0, num_reads=0;

SEM_T work_sema;
SEM_T stop_sema;
SEM_T counter_sema;

FILE* outfile_fp;
FILE* dna_outfile_fp;

char mystring[STRINGLEN];
char complete_sequence[STRINGLEN];

FASTAFILE* fp;

void parseArguments(int argc, char **argv) {
    /* read command line argument */
    //!! This argument reading should all be encapsulated in a single function, this will make reading the code much easier, right now we have to always move around it.
    if (argc <= 8) {
        fprintf(stderr, "ERROR: You missed some parameters for input\n");
        print_usage();
        exit(EXIT_FAILURE);
    }

    int c;

    while ((c=getopt(argc, argv, "fs:m:o:w:r:t:p:dev")) != -1) {
        switch (c) {
        case 's':
            strcpy(seq_file, optarg);

            if(strcmp(seq_file, "stdin") == 0) {
                break;
            } else if (access(seq_file, F_OK)==-1) {
                fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'm':
            max_mem = atoi(optarg);
            break;
        case 'w':
            wholegenome = atoi(optarg);
            if (wholegenome != 0 && wholegenome != 1) {
                fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'p':
            threadnum = atoi(optarg);
            if (threadnum < 1) {
                fprintf(stderr, "ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'o':
            strcpy(out_file, optarg);
            strcpy(aa_file, out_file);
            strcat(aa_file, ".faa");
            strcpy(dna_file, out_file);
            strcat(dna_file, ".ffn");
            break;
        case 'r':
            setTrainDirectory(optarg);
            break;
        case 't':
            strcpy(train_file, optarg);
            strcpy(hmm_file, train_dir);
            strcat(hmm_file, train_file);
            if (access(hmm_file, F_OK)==-1) {
                fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'f':
            format = 1;
            break;
        case 'd':
            output_dna = 1;
            break;
        case 'e':
            output_meta = 1;
            break;
        case 'v':
            verbose = 1;
            break;
        }
    }

    optind = 1;
    while ((c=getopt(argc, argv, "fs:m:o:w:r:t:p:dev")) != -1) {
        switch (c) {
        case 'r':
            setTrainDirectory(optarg);
            break;
        }
    }
}

void checkFiles() {
    /* check whether the specified files exist */
    if (access(mstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: Forward prob. file [%s] does not exist\n", mstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(rstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: Backward prob. file [%s] does not exist\n", rstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(nstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: noncoding prob. file [%s] does not exist\n", nstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(sstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: start prob. file [%s] does not exist\n", sstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(pstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: stop prob. file [%s] does not exist\n", pstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(s1state_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: start1 prob. file [%s] does not exist\n", s1state_file);
        exit(EXIT_FAILURE);
    }
    if (access(p1state_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: stop1 prob. file [%s] does not exist\n", p1state_file);
        exit(EXIT_FAILURE);
    }
    if (access(dstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: pwm dist. file [%s] does not exist\n", dstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(hmm_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: hmm file [%s] does not exist\n", hmm_file);
        exit(EXIT_FAILURE);
    }
}

void setMemoryLimits() {
    /* check for mem limit, allocate buffer */
    if (max_mem <= 0) {
        printf("Max memory limit specified invalid, defaulting to 1024MB\n");
        max_mem = 1024;
    }

    // 5 stands for the number of buffers we are currently using per thread
    MAX_BYTES_PER_BUFFER = max_mem*100000/(5*2*threadnum);
    MAX_SEQS_PER_BUFFER = MAX_BYTES_PER_BUFFER/STRINGLEN;
}

void checkOutputFiles() {

    // remove them, if they already exist
    remove(aa_file);
    if (output_meta) remove(out_file);
    if (output_dna) remove(dna_file);
}

void setTrainDirectory(char* train_path) {
    strcpy(train_dir, train_path);
    strcat(train_dir, "/");
    strcpy(mstate_file, train_dir);
    strcat(mstate_file, "gene");
    strcpy(rstate_file, train_dir);
    strcat(rstate_file, "rgene");
    strcpy(nstate_file, train_dir);
    strcat(nstate_file, "noncoding");
    strcpy(sstate_file, train_dir);
    strcat(sstate_file, "start");
    strcpy(pstate_file, train_dir);
    strcat(pstate_file, "stop");
    strcpy(s1state_file, train_dir);
    strcat(s1state_file, "stop1");
    strcpy(p1state_file, train_dir);
    strcat(p1state_file, "start1");
    strcpy(dstate_file, train_dir);
    strcat(dstate_file, "pwm");
}

void initializeSemaphores() {

#ifdef __APPLE__
    sem_unlink("/work_sema");
    if ( ( work_sema = sem_open("/work_sema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
        perror("ERROR: sem_open");
        exit(EXIT_FAILURE);
    }

    sem_unlink("/sema_Q");
    if (( sema_Q = sem_open("/sema_Q", O_CREAT, 0644, 1))== SEM_FAILED ) {
        perror("ERROR: sem_open");
        exit(EXIT_FAILURE);
    }

    sem_unlink("/sema_R");
    if (( sema_R = sem_open("/sema_R", O_CREAT, 0644, 1)) == SEM_FAILED ) {
        perror("ERROR: sem_open");
        exit(EXIT_FAILURE);
    }

    sem_unlink("/sema_r");
    if (( sema_r = sem_open("/sema_r", O_CREAT, 0644, 0)) == SEM_FAILED ) {
        perror("sem_open");
        exit(EXIT_FAILURE);
    }

    sem_unlink("/sema_w");
    if (( sema_w = sem_open("/sema_w", O_CREAT, 0644, 0)) == SEM_FAILED ) {
        perror("sem_open");
        exit(EXIT_FAILURE);
    }

    sem_unlink("/stop_sema");
    if (( stop_sema = sem_open("/stop_sema", O_CREAT, 0644, 0)) == SEM_FAILED ) {
        perror("ERROR: sem_open");
        exit(EXIT_FAILURE);
    }

    sem_unlink("/COUNTER_SEMA");
    if (( counter_sema = sem_open("/COUNTER_SEMA", O_CREAT, 0644, 1)) == SEM_FAILED ) {
        perror("ERROR: sem_open");
        exit(EXIT_FAILURE);
    }

#elif __linux
    sem_init(&work_sema, 0, 1);
    sem_init(&sema_Q, 0, 1);
    sem_init(&sema_R, 0, 1);
    sem_init(&sema_r, 0, 0);
    sem_init(&sema_w, 0, 0);
    sem_init(&stop_sema, 0, 0);
    sem_init(&counter_sema, 0, 1);
#endif
}

void destroySemaphores() {

#ifdef __APPLE__
    sem_unlink("/work_sema");
    sem_unlink("/sema_Q");
    sem_unlink("/sema_R");
    sem_unlink("/sema_r");
    sem_unlink("/sema_w");
    sem_unlink("/stop_sema");
    sem_unlink("/COUNTER_SEMA");

    char name[40];
    int j;
    for(j=0; j<threadnum; j++) {
        sprintf(name, "/sema_r%d", j);
        sem_unlink(name);

        sprintf(name, "/sema_w%d", j);
        sem_unlink(name);
    }

#elif __linux
    sem_destroy(&work_sema);
    sem_destroy(&sema_Q);
    sem_destroy(&sema_R);
    sem_destroy(&sema_r);
    sem_destroy(&sema_w);
    sem_destroy(&stop_sema);
    sem_destroy(&counter_sema);
#endif

}

void initializeThreads() {

    pthread_t *thread = calloc(threadnum, sizeof(pthread_t*));
    thread_datas = malloc(sizeof(thread_data) * threadnum);

    // allocate memory for each thread only once!
    if (verbose)
        printf("DEBUG: Allocating memory for all threads...\n");

    int k;
    for(k=0; k< threadnum; k++) {
        init_thread_data(thread_datas+k);
    }

    if (verbose) {
        printf("DEBUG: Allocated memory for all threads!\n");
        printf("DEBUG: Starting writer thread...\n");
    }

    pthread_create(&writer_thread, 0, writerThread, 0);

    fp = OpenFASTA(seq_file);

    if (!fp) {
        printf("ERROR! Could not open seq_file %s for reading...!\n", seq_file);
        exit(EXIT_FAILURE);
    }

    if (verbose)
        printf("INFO : Giving workers initial inputs...\n");

    int i,j;
    void *status;
    for(j=0; j<threadnum; j++)
        pthread_create(&thread[j], 0, workerThread, (void *)(thread_datas+j));

    for(j=0; j<threadnum; j++) {
        for(i=0; i<2; i++) {
            if( (stopped_at_fpos = read_seq_into_buffer(fp, thread_datas + j, i))!=0) {
                sem_post(thread_datas[j].sema_r);
            }
        }
    }

    if (verbose)
        printf("INFO : Initializing worker threads...\n");

}

void readerThread() {
    // master loop - while we haven't exhausted reading the file yet
    while (stopped_at_fpos!=0) {
        sem_wait(sema_r);

        sem_wait(sema_Q);
        QUEUE* temp;
        cutnpaste_q(&temp, EMPTY_Q);
        sem_post(sema_Q);

        while(temp) {
            sem_wait(sema_R);
            stopped_at_fpos = read_seq_into_buffer(fp,  temp->td, temp->buffer);

            if(stopped_at_fpos == 0) {
                num_reads_flag =1;
            }

            sem_post(sema_R);

            sem_post(temp->td->sema_r);
            temp = temp->next;
        }
    }
    CloseFASTA(fp);

    if (verbose)
        printf("INFO : Finished handing out all the work...\n");

    num_reads_flag =1;

    sem_wait(stop_sema);
}

int main (int argc, char **argv) {

    fp = 0;
    setTrainDirectory("train");

    parseArguments(argc, argv);

    checkFiles();

    setMemoryLimits();

    checkOutputFiles();

    initializeSemaphores();

    if (verbose)
        printf("DEBUG: Max number of sequences per thread : %d, max bytes per thread : %d\n", MAX_SEQS_PER_BUFFER, MAX_BYTES_PER_BUFFER*5);

    /* read all initial model */
    get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

    // prepare all of the worker threads as well as the writer thread
    initializeThreads();

    // master loop - while we haven't exhausted reading the file yet
    readerThread();

    // destroy the semaphores if we have a mac machine
    destroySemaphores();

    printf("Run finished with %d threads.\n", threadnum);
}

int read_seq_into_buffer(FASTAFILE* ffp, thread_data* thread_data, unsigned int buf) {

    char *seq;
    char *name;
    int   L;
    int status ;
    int count=0;

    while ( (count < MAX_SEQS_PER_BUFFER) && (status = ReadFASTA(ffp, &seq, &name, &L)) ==1 ) {
        strcpy(thread_data->input_head_buffer[buf][count], name);
        strcpy(thread_data->input_buffer[buf][count], seq);
        read_counter++;
        count++;
    }

    thread_data->input_num_sequences[buf] = count;
    read_counter1 += thread_data->input_num_sequences[buf];

    return count;

}

// read as much of the file as you can into the buffer, return the file position
// of the last carat character encountered before the memory limit was hit
off_t read_file_into_buffer(FILE* fp, int fpos, thread_data* thread_data, unsigned int buf) {

    int i = 0; // temp length of each sequence
    long count = 0; // index for the number of sequences
    off_t last_carat_position = 0;
    long long total_count = 0;
    int seq_length;

    memset(complete_sequence, 0, STRINGLEN);
    memset(mystring, 0, STRINGLEN);

    // loads the input
    while ( fgets(mystring, STRINGLEN, fp) && count < MAX_SEQS_PER_BUFFER-1 ) {
        if (mystring[strlen(mystring)-1] == '\n') mystring[strlen(mystring)-1] = 0;
        if (mystring[0] == '>') {
            if (i>0) {
                total_count += i;
                memset(thread_data->input_buffer[buf][count], 0, STRINGLEN);
                strcpy(thread_data->input_buffer[buf][count], complete_sequence);
                printf("%s\n",thread_data->input_buffer[buf][count]);
                memset(complete_sequence, 0, STRINGLEN);
            }
            count++;
            read_counter++;
            strcpy(thread_data->input_head_buffer[buf][count], mystring);
            last_carat_position = ftello(fp);
            i = 0;
        } else {
            strcat(complete_sequence, mystring);
            seq_length = strlen(mystring);
            while(mystring[seq_length-1] == 10 || mystring[seq_length-1]==13 || mystring[seq_length-1]==0) {
                seq_length --;
            }
            i += seq_length;
        }
    }

    strcpy(thread_data->input_buffer[buf][count], complete_sequence);
    thread_data->input_num_sequences[buf] = count;
    read_counter1 += thread_data->input_num_sequences[buf];

    if (ferror(fp)) {
        printf("Error with file encountered.\n");
    } else if(feof(fp)) {
        last_carat_position = -1; // terminating condition
    }
    num_reads++;

    return last_carat_position;
}


void init_thread_data(thread_data* td) {
    // Initialize thread data structure

    td->hmm = calloc(1, sizeof(HMM));
    memcpy(td->hmm, &hmm, sizeof(HMM));

    td->wholegenome = wholegenome;

    td->id = ids;
    ids++;

#ifdef __APPLE__
    char name[40];
    sprintf(name, "/sema_r%d", td->id);
    sem_unlink(name);

    if((td->sema_r = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) {
        perror("sem_open");
        exit(EXIT_FAILURE);
    }

    sprintf(name, "/sema_w%d", td->id);
    sem_unlink(name);
    if (( td->sema_w = sem_open(name, O_CREAT, 0644, 2)) == SEM_FAILED ) {
        perror("sem_open");
        exit(EXIT_FAILURE);
    }

#elif __linux
    sem_init(&td->sema_r, 0, 0);
    sem_init(&td->sema_w, 0, 2);
#endif

    // TODO : refactor to as many single large malloc calls as possible
    td->output_num_sequences = calloc(2, sizeof(int));
    td->input_num_sequences = calloc(2, sizeof(int));

    td->input_buffer = (char ***)malloc(sizeof(char**) * 2);
    td->input_head_buffer = (char ***)malloc(sizeof(char**) * 2);
    td->output_buffer = (char ***)malloc(sizeof(char**) * 2);
    td->aa_buffer = (char ***)malloc(sizeof(char**) * 2);
    td->dna_buffer = (char ***)malloc(sizeof(char**) * 2);

    td->dna	= calloc(STRINGLEN, sizeof(char));
    td->dna1 = calloc(STRINGLEN, sizeof(char));
    td->dna_f = calloc(STRINGLEN, sizeof(char));
    td->dna_f1 = calloc(STRINGLEN, sizeof(char));
    td->protein = calloc(STRINGLEN, sizeof(char));
    td->temp_str = calloc(STRINGLEN, sizeof(char));

    td->insert = malloc(sizeof(int) * STRINGLEN);
    td->c_delete = malloc(sizeof(int) * STRINGLEN);

    int i;
    for (i=0; i<2; i++) {
        td->input_buffer[i]	= (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
        td->input_head_buffer[i] = (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
        td->output_buffer[i]	= (char **) malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
        td->aa_buffer[i]	=   (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
        td->dna_buffer[i]	= (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);

        int j;
        for(j=0; j<MAX_SEQS_PER_BUFFER; j++) {
            td->input_buffer[i][j] = calloc(1, STRINGLEN);
            td->input_head_buffer[i][j] = calloc(1, STRINGLEN);
            td->aa_buffer[i][j] = calloc(1, STRINGLEN);
            td->dna_buffer[i][j] = calloc(1, STRINGLEN);
            td->output_buffer[i][j] = calloc(1, STRINGLEN);
        }
    }
}

void writeOutputFiles(FILE* aa_outfile_fp, thread_data* td, unsigned int buffer) {
    if (output_meta) {
        writeMeta();
    }
    if (output_dna) {
        writeDNA();
    }
    writeAminoAcids(aa_outfile_fp, td, buffer);
    num_writes++;
    if (verbose) printf("INFO: Wrote results for thread %d, buffer %d.\n", td->id, buffer );
}

void writeDNA() {
    dna_outfile_fp = fopen(dna_file, "a");
    if (!dna_outfile_fp) {
        printf("ERROR: Could not open dna output file %s for writing!\n", dna_file);
        exit(EXIT_FAILURE);
    }
}

void writeMeta() {
    outfile_fp = fopen(out_file, "a");
    if(!outfile_fp) {
        printf("ERROR: Could not open meta output file %s for writing!\n", out_file);
        exit(EXIT_FAILURE);
    }
}

void writeAminoAcids(FILE* aa_outfile_fp, thread_data* td, unsigned int buffer) {

    int j;
    for(j = 0; j < td->output_num_sequences[buffer]; j++) {
        writer_counter++;
        char *ptrc;
        if(td->aa_buffer[buffer][j][0]!=0) {
            ptrc=td->aa_buffer[buffer][j];
            while(*ptrc!='\0') {
                if(*ptrc=='\t') *ptrc='>';
                ptrc++;
            }
            fprintf(aa_outfile_fp, ">%s", td->aa_buffer[buffer][j]);
        }

        //!! Why are we clearing dna and output buff?
        stopMemset(td->output_buffer[buffer][j], STRINGLEN);
        stopMemset(td->aa_buffer[buffer][j], STRINGLEN);
        stopMemset(td->dna_buffer[buffer][j], STRINGLEN);
    }
}

FILE* openFilePointers() {

    FILE* aa_outfile_fp;
    if(strcmp(out_file, "stdout") == 0) {
        aa_outfile_fp = stdout;
    } else {
        aa_outfile_fp = fopen(aa_file, "a");
    }

    if(!aa_outfile_fp) {
        printf("ERROR: Could not open aa output file %s for writing!\n", aa_file);
        exit(EXIT_FAILURE);
    }

    return aa_outfile_fp;
}

void closeFilePointers( FILE** aa_outfile_fp, FILE** outfile_fp, FILE** dna_outfile_fp ) {

    fclose(*aa_outfile_fp);
    if (output_meta) fclose(*outfile_fp);
    if (output_dna) fclose(*dna_outfile_fp);

    *aa_outfile_fp = NULL;
    *outfile_fp = NULL;
    *dna_outfile_fp = NULL;
}

void* writerThread(void* args) {

    int j;
    FILE* aa_outfile_fp = openFilePointers();

    while(1) {

        sem_wait(sema_w);

        QUEUE* temp;

        sem_wait(sema_Q);
        cutnpaste_q(&temp, DONE_Q);
        sem_post(sema_Q);

        while(temp) {

            sem_wait(sema_R);

            thread_data* td = temp->td;
            unsigned int buffer = temp->buffer;

            writeOutputFiles(aa_outfile_fp, td, buffer);

            sem_post(sema_R);
            sem_post(td->sema_w);

            temp = temp->next;
        }

        if (num_reads_flag == 1 && writer_counter == read_counter) {
            sem_post(stop_sema);
            break;
        }

    }

    closeFilePointers(&aa_outfile_fp, &outfile_fp, &dna_outfile_fp );
}

void runViterbiOnBuffers(thread_data *td, unsigned int b) {

    unsigned int i;
    for (i=0; i < td->input_num_sequences[b]; i++) {
        unsigned int stringlength = strlen(td->input_buffer[b][i]);
        get_prob_from_cg(td->hmm, &train, td->input_buffer[b][i], stringlength);

        if (td->input_buffer[b][i] && td->input_head_buffer[b][i] ) {

            viterbi(td->hmm, td->input_buffer[b][i], td->output_buffer[b][i],
                    td->aa_buffer[b][i], td->dna_buffer[b][i],
                    td->input_head_buffer[b][i], td->wholegenome, td->format, stringlength,
                    td->dna, td->dna1, td->dna_f, td->dna_f1, td->protein,
                    td->insert, td->c_delete, td->temp_str);

            sem_wait(work_sema);
            work_counter++;
            sem_post(work_sema);
        }
    }

}

void* workerThread(void *_thread_datas) {

    thread_data *td = (thread_data*)_thread_datas;
    unsigned int b = 0;

    while(1) {
        sem_wait(td->sema_r);
        sem_wait(td->sema_w);

        sem_wait(counter_sema);
        viterbi_counter +=  td->input_num_sequences[b];
        sem_post(counter_sema);

        runViterbiOnBuffers(td, b);
        td->output_num_sequences[b] = td->input_num_sequences[b];

        if (verbose) {
            printf("INFO: Thread %d buffer %d done work on %d sequences!\n",
                   td->id, b, td->input_num_sequences[b]);
        }

        sem_wait(sema_Q);
        enqueue(td, b, EMPTY_Q);
        enqueue(td, b, DONE_Q);

        sem_post(sema_Q);
        sem_post(sema_r);
        sem_post(sema_w);

        b = (b + 1) % 2;
    }

    return (void*) 0;
}
