#ifndef __RUN_HMM_H
#define __RUN_HMM_H

void writeDNA();
void writeMeta();
void writeAminoAcids(FILE* aa_outfile_fp, thread_data* td, unsigned int buffer);

void parseArguments(int argc, char **argv);
void checkFiles();
void setMemoryLimits();
void checkOutputFiles();


void* writerThread(void* args);
void readerThread();
void* workerThread(void *_thread_datas);

#endif
