/*********

File        CHROMEISTER.c
Author      EPW <estebanpw@uma.es>
Description Computes hits and generates a dotplot

USAGE       Usage is described by calling ./CHROMEISTER --help



**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"

#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes
#define RANGE 2


uint64_t custom_kmer = 12; // Defined as external in structs.h

void init_args(int argc, char ** av, FILE ** database, FILE ** query);
long double get_random_num(const unsigned char * word, uint64_t ** w_array, uint64_t * w_sizes);
long double estimate_pi(uint64_t t_in_circle, uint64_t t_total);
int is_within_circle(long double x1, long double x2);

int main(int argc, char ** av){
    

    

    uint64_t i, j;

    //query to read kmers from, database to find seeds
    FILE * database = NULL, * query;
    long double MY_PI = (long double) 3.14159265359;
    
    init_args(argc, av, &database, &query);

    
    //Variables to account for positions
    //Print info
    fprintf(stdout, "[INFO] Generating K-mer tables.\n");
    //Variables to read kmers
    char c = 'N'; //Char to read character
    //Current length of array and variables for the buffer
    uint64_t idx = 0, r = 0;
    
    //Vector to read in batches
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }

    unsigned char curr_kmer[custom_kmer];
    curr_kmer[0] = '\0';
    uint64_t word_size = 0;

    //To hold all information related to database
    uint64_t current_len = 0;
    fseek(database, 0, SEEK_END);
    uint64_t aprox_len_query = ftell(database);
    rewind(database);
    uint64_t a_hundreth = (aprox_len_query/100);
    
    //To force reading from the buffer
    idx = READBUF + 1;

    //unsigned char aux_kmer[custom_kmer+1];
    
    //Vector to store query seq
    unsigned char * seq_vector_query = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_query == NULL) terror("Could not allocate memory for query vector");

    uint64_t w_sizes[custom_kmer];
    uint64_t ** w_array = (uint64_t **) calloc(custom_kmer, sizeof(uint64_t *));
    if(w_array == NULL) terror("Could not allocate wide array");

    for(i=0; i<custom_kmer; i++){
        w_sizes[i] = powl(4, i+1);
        w_array[i] = (uint64_t *) calloc(w_sizes[i], sizeof(uint64_t));
        if(w_array[i] == NULL){ fprintf(stdout, "Container %" PRIu64 ":\n", i); terror("Could not allocate container"); }
    }

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    while((!feof(database) || (feof(database) && idx < r))){
        if(c == '>'){
            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);  //Skip ID
            while(c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < custom_kmer) ++word_size;
                    ++current_len;
                    if(current_len % a_hundreth == 0){ 
                        fprintf(stdout, "\r%"PRIu64"%%...", 1+100*current_len/aprox_len_query); 
                        fflush(stdout);
                    }
                }else{ //It can be anything (including N, Y, X ...)
                    if(c != '\n' && c != '>'){
                        word_size = 0;
                        ++current_len;
                    } 
                }
                if(word_size == custom_kmer){
                    // Processing
                    for(j=0; j<custom_kmer; j++){
                        uint64_t hash = hashOfWord(curr_kmer, j+1, 0);
                        //fprintf(stdout, "at (%"PRIu64") generated hash %"PRIu64" for seq %s\n", j, hash, (char *)curr_kmer); getchar();
                        w_array[j][hash]++;
                    }
		            // Overlapping
                    memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                    --word_size;
                }
            }
            word_size = 0;   
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
    }

    fprintf(stdout, "[INFO] Rolling numbers.\n");

    //To force reading from the buffer
    idx = READBUF + 1;
    //rewind(database);

    uint64_t t_in_circle = 0, t_total = 0;
    long double x1, x2;
    unsigned char have_one_already = FALSE;
    current_len = 0;
    word_size = 0;

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    while((!feof(query) || (feof(query) && idx < r))){
        if(c == '>'){
            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);  //Skip ID
            while( c != '>' && (!feof(query) || (feof(query) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < custom_kmer) ++word_size;
                    ++current_len;
                    if(current_len % a_hundreth == 0 && t_total > 0){ 
                        //fprintf(stdout, "%"PRIu64"%%...\n", 1+100*current_len/aprox_len_query); 
                        fprintf(stdout, "PI estimation: %.17Le\n", estimate_pi(t_in_circle, t_total));
                        fprintf(stdout, "Real PI:       %.17Le\n", MY_PI);
                        fflush(stdout);
                    }
                }else{ //It can be anything (including N, Y, X ...)
                    if(c != '\n' && c != '>'){
                        word_size = 0;
                        ++current_len;
                    } 
                }
                if(word_size == custom_kmer){
                    // Processing
                    if(have_one_already == FALSE){
                        x1 = get_random_num(curr_kmer, w_array, w_sizes);
                        have_one_already = TRUE;
                    }else{
                        x2 = get_random_num(curr_kmer, w_array, w_sizes);
                        if(is_within_circle(x1, x2) == 1) ++t_in_circle;
                        ++t_total;
                        have_one_already = FALSE;
                    }
                    
                    //fprintf(stdout, "Word %s gives %e\n", (char *) curr_kmer, v_rand);
                    //getchar();
		            // Overlapping
                    memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                    --word_size;
                }
            }
            word_size = 0;   
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);    
        }
    }

    fprintf(stdout, "PI estimation: %.17Le\n", estimate_pi(t_in_circle, t_total));
    fprintf(stdout, "Real PI:       %.17Le\n", MY_PI);


    for(i=0; i<custom_kmer; i++){
        free(w_array[i]);
    }
    free(w_array);
    
    return 0;
}

long double get_random_num(const unsigned char * word, uint64_t ** w_array, uint64_t * w_sizes){

    uint64_t i, n = 0;
    uint64_t conv[256];
    conv[(unsigned char) 'A'] = 4771;
    conv[(unsigned char) 'C'] = 4505;
    conv[(unsigned char) 'G'] = 7352;
    conv[(unsigned char) 'T'] = 8032;

    for(i=0; i<custom_kmer; i++){
        // sum (Vi * Ti * 4^(k-i))
        //n = n + w_array[i][hashOfWord(word, i+1, 0)] * w_sizes[(custom_kmer-1)-i];
        n = n + ( ((uint64_t) conv[word[i]]) * (1+w_array[i][hashOfWord(word, i+1, 0)] * w_sizes[(custom_kmer-1)-i]));
    }
    //printf("before value n is %"PRIu64"\n", n);
    n = n % 22369620;
    //printf("after value n is %"PRIu64"\n", n);
    //getchar();
    return (long double) n / (long double) (22369620-1);
}

long double estimate_pi(uint64_t t_in_circle, uint64_t t_total){
    return ((long double) t_in_circle / (long double) t_total) / 0.25;
}

int is_within_circle(long double x1, long double x2){
    long double dsquared = ((long double)x1 - 0.5) * ((long double)x1 - 0.5) + ((long double)x2 - 0.5) * ((long double)x2 - 0.5);
    //printf(" dsquared %e\n", dsquared);
    if(dsquared <= 0.25) return 1;
    return 0;
}


void init_args(int argc, char ** av, FILE ** database, FILE ** query){


    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           RANDdna -in [database] -test [query]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        
        if(strcmp(av[pNum], "-in") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open input file");
        }
        if(strcmp(av[pNum], "-test") == 0){
            *query = fopen64(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open input file");
        }
        
        pNum++;
    }
    
    if(*database==NULL || *query==NULL) terror("An input file is required and a test is required");
}

