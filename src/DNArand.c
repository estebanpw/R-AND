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

void init_args(int argc, char ** av, FILE ** database);
long double get_random_num_from_bytes(const unsigned char * eight_bytes);
long double get_random_num_from_8_bytes(uint64_t eight_bytes);
long double get_random_num(const unsigned char * word, uint64_t ** w_array, uint64_t * w_sizes);
long double estimate_pi(uint64_t t_in_circle, uint64_t t_total);
int is_within_circle(long double x1, long double x2);

int main(int argc, char ** av){
    

    

    uint64_t i, j;

    //query to read kmers from, database to find seeds
    FILE * database = NULL;
    long double MY_PI = (long double) 3.14159265359;
    
    init_args(argc, av, &database);

    
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
    uint64_t a_thousandth = (aprox_len_query/1000);
    
    //To force reading from the buffer
    idx = READBUF + 1;

    //unsigned char aux_kmer[custom_kmer+1];
    
    //Vector to store query seq
    unsigned char * seq_vector_query = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_query == NULL) terror("Could not allocate memory for query vector");


    //To force reading from the buffer
    idx = READBUF + 1;
    //rewind(database);

    uint64_t t_in_circle = 0, t_total = 0;
    long double x1, x2;
    unsigned char have_one_already = FALSE, random_init_completed = FALSE;
    current_len = 0;
    word_size = 0;
    long double vAT = 0, vCG = 0;
    long double G = 0, C = 0, A = 0, T = 0;
    unsigned char curr_bytes[64];
    uint64_t bytes_taken = 0;
    uint64_t acum = 1;
    srand(time(NULL));
    uint64_t conv[256];
    conv[(unsigned char) 'A'] = 89465;
    conv[(unsigned char) 'C'] = 1265;
    conv[(unsigned char) 'G'] = 987987465;
    conv[(unsigned char) 'T'] = 2135364;
    long double prev_pi = 0;

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    while((!feof(database) || (feof(database) && idx < r))){
        if(c == '>'){
            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);  //Skip ID
            while( c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    //curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < custom_kmer) ++word_size;
                    ++current_len;

                    if(1 == 1){ 

                        //fprintf(stdout, "%d", (vAT >= 0));

                        //curr_bytes[bytes_taken++] = (unsigned char) (vAT >= 0);
                        
                        //printf("%d\n", (rand() % (2)));
                        //curr_bytes[bytes_taken++] = (rand() % (2));
                        acum = conv[(unsigned char)c] + (conv[(unsigned char)c] * acum);
                        //fprintf(stdout, "^AT = %d\n^CG = %d\n", (vAT >= 0), (vCG >= 0));
                        //getchar();

                        vAT = A - T;
                        vCG = C - G;

                        A = C = G = T = 0;
                    }

                    if(current_len % a_hundreth == 0 && current_len > 1){ 
                        //fprintf(stdout, "%"PRIu64"%%...\n", 1+100*current_len/aprox_len_query); 
                        
                        random_init_completed = TRUE;
                        //fprintf(stdout, "A %Le\nT %Le\nC %Le\nG %Le------\n", A, T, C, G);
                        
                        fprintf(stdout, "PI var: %.17Le\n", estimate_pi(t_in_circle, t_total)-prev_pi);
                        fprintf(stdout, "PI estimation: %.17Le\n", estimate_pi(t_in_circle, t_total));
                        fprintf(stdout, "Real PI:       %.17Le\n", MY_PI);
                        prev_pi = estimate_pi(t_in_circle, t_total);
                        //getchar();
                        fflush(stdout);
                    }
                }else{ //It can be anything (including N, Y, X ...)
                    if(c != '\n' && c != '>'){
                        word_size = 0;
                        ++current_len;
                    } 
                }
                // Treat DNA chars here
                switch(c){
                    case 'A': ++A; break;
                    case 'C': ++C; break;
                    case 'G': ++G; break;
                    case 'T': ++T; break;
                }

                


                if(current_len > 64){
                    // Processing
                    if(have_one_already == FALSE){
                        //x1 = get_random_num(curr_kmer, w_array, w_sizes);
                        x1 = get_random_num_from_8_bytes(acum);
                        //printf("%Le\n", x1);
                        have_one_already = TRUE;
                    }else{
                        //x2 = get_random_num(curr_kmer, w_array, w_sizes);
                        x2 = get_random_num_from_8_bytes(acum);
                        //printf("%Le\n", x2);
                        if(is_within_circle(x1, x2) == 1) ++t_in_circle;
                        ++t_total;
                        have_one_already = FALSE;
                    }

                    //printf(" t circle %"PRIu64" t total %"PRIu64"\n", t_in_circle, t_total);
                    
                    //fprintf(stdout, "Word %s gives %e\n", (char *) curr_kmer, v_rand);
                    //getchar();
		            // Overlapping

                    bytes_taken = 0;

                    memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                    --word_size;
                }
            }
            bytes_taken = 0;
            word_size = 0;   
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
    }

    fprintf(stdout, "PI estimation: %.17Le\n", estimate_pi(t_in_circle, t_total));
    fprintf(stdout, "Real PI:       %.17Le\n", MY_PI);

    fclose(database);
    
    return 0;
}

long double get_random_num_from_bytes(const unsigned char * eight_bytes){
    uint64_t i;
    uint64_t v = 0, c;
    for(i=0; i<64; i++){
        c = (uint64_t) eight_bytes[i];
        c = c << (63-i);
        v += c;
    }
    
    
    //printf("---------->%"PRIu64"\n", v);
    return (long double) ( (long double) v / (long double) ((uint64_t)0xFFFFFFFFFFFFFFFF));
}

long double get_random_num_from_8_bytes(uint64_t eight_bytes){
      
    
    //printf("---------->%"PRIu64"\n", v);
    return (long double) ( (long double) eight_bytes / (long double) ((uint64_t)0xFFFFFFFFFFFFFFFF));
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


void init_args(int argc, char ** av, FILE ** database){


    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           RANDdna -in [database]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        
        if(strcmp(av[pNum], "-in") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open input file");
        }
        
        pNum++;
    }
    
    if(*database==NULL) terror("An input file is required and a test is required");
}

