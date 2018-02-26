#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>     /* assert */

int main ()
{
  FILE * pFile;
  long lSize;
  float * buffer;
  size_t result;

    pFile = fopen ( "HGG_TumorIC.bin" , "rb" );

  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);
    
    std::cout<<"lSize="<<lSize<<std::endl;

  // allocate memory to contain the whole file:
  buffer = (float*) malloc (sizeof(float)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

    
    
    
  /* the whole file is now loaded in the memory buffer. */
    for(int i=0; i< 3; ++i)
     std::cout<<"b["<<i<<"]="<<buffer[i]<<" (foat) buffer[i]="<< (float) buffer[i] <<std::endl;
    
    // terminate
    fclose (pFile);
    free (buffer);


  return 0;
}
