/*
 *  Copyright (c) 2020 Leonardo Araujo (leolca@gmail.com).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>
 * 
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include "FFTMP.h"

// compile: gcc FFTMP.c -std=c99 -lmpfr -lgmp -lmpc -o FFTMP
// test: echo -n "1 1 1 1 0 0 0 0" | ./FFTMP
// expected result: 
// Data: 1 1 1 1 0 0 0 0 
// FFT : 4 (1, -2.41421) 0 (1, -0.414214) 0 (1, 0.414214) 0 (1, 2.41421)

void PrintHelp()
{
    printf("Usage: FFTMP [OPTIONS]... [FILE]...\n");
    printf("Output the FFT from a given FILE to standard output.\n");
    printf("With no FILE, or when FILE is -, read standard input.\n");
    printf("  -i, --input          input filename (if not provided, read from stdin)\n");
    printf("  -h, --help           display this help and exit\n\n");
    printf("Example:\n");
    printf("echo -n '1 1 1 1 0 0 0 0' | ./FFTMP\n");
    printf("./FFTMP -i x_file\n");
    exit(1);
}

mpc_t *read_complex_data (FILE *stream, unsigned int *xLen)
{
  int mLen = DEFAULTLEN;
  mpc_t *x = malloc(mLen * sizeof(mpc_t));
  if (x == NULL) 
  { 
     printf("Memory not allocated.\n"); 
     exit(0); 
  } 

  int c, bLen = DEFAULTSTRLEN, bidx = 0;
  char *buffer = (char*) malloc(bLen * sizeof(char));
  while( (c = fgetc(stream)) != EOF ) 
  {
     if (bidx + 1 >= bLen )
     {
	bLen *= 2;
	buffer = realloc(buffer, bLen);
     }
     if ( !isspace(c) )
     {
	buffer[bidx++] = c;
     }
     else
     {
	if (bidx > 0)
	{
           mpc_init2 (x[*xLen], 256);
           mpc_set_str (x[(*xLen)++], buffer, 10, MPC_RNDDN);
           free(buffer);
           bLen = DEFAULTSTRLEN;
           buffer = (char*) malloc(bLen * sizeof(char));
           bidx = 0;
	}
     }
  }
  if (bidx > 0)
  {
    mpc_init2 (x[*xLen], 256);
    mpc_set_str (x[(*xLen)++], buffer, 10, MPC_RNDDN);
    free(buffer);
    bidx = 0;
  }

  return x; 
}

int main(int argc, char *argv[]) 
{
  int opt;
  FILE *stream = NULL;
  mpc_t *x, *X, *xr;
  unsigned int xLen = 0, XLen = 0;
  unsigned int mLen = DEFAULTLEN;

  const char* const short_opts = "n:wbi:h";
  static struct option long_options[] =
  {
       {"input", required_argument,  0, 'i'},
       {"help"  , no_argument,       0, 'h'},
       {0       , no_argument,       0,  0 }
  };

  while ((opt = getopt_long (argc, argv, short_opts, long_options, NULL)) != -1)
  {
    switch (opt) 
    {
       case 'i':
          if((stream = fopen(optarg,"r"))==NULL) { printf("Input file not found!\n"); exit(-1); }
          break;
       case 'h': // -h or --help
       case '?': // Unrecognized option
       default:
          PrintHelp();
          break;
    }
  }
  if (stream == NULL) stream = stdin;

  x = read_complex_data (stream, &xLen);

  printf("Data:\n");
  print_mpc_array (x, xLen, 2);

  printf("Data length: %d\n",xLen);
  if (!isPowerofTwo (xLen))
  {
    printf("zeropadding to next power of two\n");
    unsigned int newxLen = nextPowerofTwo(xLen);
    zeropad_mpc_array (x, &xLen, newxLen - xLen);
  }
  
  X = fft (x, xLen);
  printf("FFT:\n");
  print_mpc_array (X, xLen, 2);  
  xr = ifft (X, xLen);
  printf("iFFT:\n");
  print_mpc_array (xr, xLen, 2);

  clear_mpc_array (x, xLen);
  clear_mpc_array (X, xLen);
  clear_mpc_array (xr, xLen);
  fclose(stream);
  return 0;
}

