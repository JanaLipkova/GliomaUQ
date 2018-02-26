#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>     /* assert */

int main ()
{
    FILE * pFile;
    float buffer[3] = { 0.315 , 0.65 , 0.5 };
    pFile = fopen ("HGG_TumorIC.bin", "wb");
    fwrite (buffer , sizeof(float), sizeof(buffer), pFile);

    // terminate
    fclose (pFile);
    free (buffer);


  return 0;
}
