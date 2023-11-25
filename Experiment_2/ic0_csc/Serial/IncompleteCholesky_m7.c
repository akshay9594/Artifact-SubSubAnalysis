

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>


int* row_val; 
int* col_val; 
double* nnz_val;

int* row;
int* col;
double* val;
double* val_par;

int convertStrtoArr(char* str)
{
    // get length of string str
    int str_length = strlen(str);
  
    // create an array with size as string
    // length and initialize with 0
    int arr =0;
  
    int j = 0, i, sum = 0;
  
    // Traverse the string
    for (i = 0; str[i] != '\0'; i++) {
  
        // if str[i] is ', ' then split
        if (str[i] == ',')
            continue;
         if (str[i] == ' '){
            // Increment j to point to next
            // array location
            break;
        }
        else {
  
            // subtract str[i] by 48 to convert it to int
            // Generate number by multiplying 10 and adding
            // (int)(str[i])
            arr = arr * 10 + (str[i] - 48);
        }
    }
   
      return arr;
   
  
}

int main(int argc, char *argv[]) {


  if(argc < 2){
   printf("Input missing arguments, you need to specify input list file\n");
  }

   char* file_path = argv[1];
  //std::string inputMatrix;
  int i,j,k,m,l;

  struct timeval start,end;//,startTT, endTT;

    FILE * fp;
    char * line = NULL;
    int num_rows, num_cols, nonzeros, colCnt=0, nnzCnt=0,count=0;

    size_t len = 0;
    size_t read;
    char matrix_details[100];

    fp = fopen(file_path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1) {
        if(line[0] == '%'){
            continue;
        }
         else{
            break;
         }
    }

   char delim[] = " ";

	char *ptr = strtok(line, delim);
   char* rowstr = NULL;
   rowstr = (char *)malloc(sizeof(line)); 
   char* temp_str;

   while(ptr != NULL)
	{
      if(count ==0){  
         memcpy(rowstr, ptr, sizeof(ptr));
         num_rows = convertStrtoArr(ptr);
      }
      else if(count == 1){
        num_cols = convertStrtoArr(ptr);
      }
      else if(count == 2){
         nonzeros = convertStrtoArr(ptr);
      }
      count++;
		ptr = strtok(NULL, delim);
	}

   row_val = (int*)malloc(sizeof(int)*(nonzeros));
   col_val = (int*)malloc(sizeof(int)*nonzeros);
   nnz_val = (double*)malloc(sizeof(double)*nonzeros);


   i=0;
   while ((read = getline(&line, &len, fp)) != -1) {
       char *ptr = strtok(line, delim);
       count = 0;
       while(ptr != NULL){
         if(count ==0){  
            row_val[i]= convertStrtoArr(ptr);
         }
         else if(count == 1){
            col_val[i] = convertStrtoArr(ptr);
         }
         else if(count == 2){
            nnz_val[i] =strtod(ptr, &temp_str);
         }
          count++; 
          ptr = strtok(NULL, delim);
         
      }
       i++; 
   }

   row = malloc(sizeof(int)*nonzeros);
   col = malloc(sizeof(int)*nonzeros);
   val = malloc(sizeof(double)*nonzeros);

   col[0]=0;
   int x,y,col_num=0;

   for(j=0; j<nonzeros; j++){
      y = col_val[j];
      y=y-1;
      // x = row_val[j];
      // x=x-1;
      if(y==col_num){
         val[nnzCnt]=nnz_val[nnzCnt];
         //row[nnzCnt]=x;
         colCnt++; nnzCnt++;
      }
        else{//New col
         col[col_num+1]=col[col_num]+colCnt;
         col_num++;
         colCnt=1;
         val[nnzCnt]=nnz_val[nnzCnt];
         //row[nnzCnt]=x;
         nnzCnt++;
      }
      
   
   }
  
   for(i=0; i< nonzeros; i++){
      for(j=col[i]; j< col[i+1]; j++){
          x = row_val[j];
          x=x-1;
         row[j] = x;
         //printf("row[%d]=%d\n", j, row[j]);
      }
   }

   col[num_cols]= col[num_cols - 1] + colCnt;//last col

   int num_runs = 1;
   double seconds, total_time=0.0;

   for(k=0; k < num_runs; k++){

      gettimeofday(&start,NULL);

      //ic0_csr(num_cols,  val, col, row);

      for (i = 0; i < num_cols - 1; i++){
         val[col[i]] = sqrt(val[col[i]]);//S1

         for (m = col[i] + 1; m < col[i+1]; m++){
            val[m] = val[m] / val[col[i]];//S2
         }

         for (m = col[i] + 1; m < col[i+1]; m++) {
            for (k = col[row[m]] ; k < col[row[m]+1]; k++){
               for ( l = m; l < col[i+1] ; l++){
                  if (row[l] == row[k] && row[l+1] <= row[k]){
                     val[k] -= val[m]* val[l]; //S3
                  }
               }
            }
         }
      }

      gettimeofday(&end, NULL);

      double seconds = (end.tv_sec + (double)end.tv_usec/1000000) - (start.tv_sec + (double)start.tv_usec/1000000); 

      total_time+= seconds;

   }

   printf("Time taken by the Serial kernel= %f s\n",total_time/num_runs);


    fclose(fp);
    if (line)
      free(line);
    exit(EXIT_SUCCESS);

   return 0;
}

