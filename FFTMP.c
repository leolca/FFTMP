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
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

#define DEFAULTLEN 256
#define DEFAULTSTRLEN 256

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

int isPowerofTwo(unsigned int x)
{
  if (x == 0) return 0;
  while (x > 1)
  {
     if (x % 2 != 0) return 0;
     x /= 2;
  }
  return 1;
}

unsigned int nextPowerofTwo(unsigned int x)
{
  // Bit Twiddling Hacks
  // By Sean Eron Anderson
  // https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
  x--;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x++;
  return x;
}

void clear_mpc_array(mpc_t *x, unsigned int xLen)
{
  for (int i=0; i<xLen; i++) mpc_clear(x[i]);
  free(x);
}

void _fft(mpc_t *x, mpc_t *X, unsigned int n, unsigned int step)
{
   if (step < n)
   {
     _fft(X, x, n, step*2);
     _fft(X + step, x + step, n, step*2);

     mpc_t r, I;
     mpfr_t PI;
     mpc_init2 (r, 256);
     mpc_init2 (I, 256);
     mpfr_inits2 (256, PI, (mpfr_ptr) 0);
     mpfr_const_pi (PI, MPFR_RNDN);
     mpc_set_ui_ui (I, 0, -1, MPC_RNDDN);
     for (int i = 0; i < n; i += 2 * step) 
     {
	// r = exp(-I * PI * i / n) * X[i + step];
	mpc_mul_fr (r, I, PI, MPC_RNDDN); // -I * pi
	mpc_mul_ui (r, r, i, MPC_RNDDN); // -I * pi * i
	mpc_div_ui (r, r, n, MPC_RNDDN); // -I * pi * i / n
	mpc_exp (r, r, MPC_RNDDN); // exp(-I * pi * i / n)
	mpc_mul (r, r, X[i + step], MPC_RNDDN);
	mpc_add (x[i/2], X[i], r, MPC_RNDDN); // x[i/2] = X[i] + r;
	mpc_sub (x[(i + n)/2], X[i], r, MPC_RNDDN); // x[(i + n)/2] = X[i] - r;
     }
     mpc_clear(r);
     mpc_clear(I);
     mpfr_clears (PI, (mpfr_ptr) 0);
   }
}

void fft(mpc_t *x, unsigned int n)
{
   mpc_t *X = malloc(n * sizeof(mpc_t));
   for (int i = 0; i < n; i++) 
   {
      mpc_init2 (X[i], 256);
      mpc_set(X[i], x[i], MPC_RNDDN);
   }
   _fft(x, X, n, 1);
   clear_mpc_array (X, n);
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

void zeropad_mpc_array (mpc_t * x, unsigned int *xLen, unsigned int n)
{
   for (int i = 0; i < n; i++)
   {
      mpc_init2 (x[(*xLen)+i], 256);
      mpc_set_ui (x[(*xLen)+i], 0, MPC_RNDDN);
   }
   *xLen += n;
}

void print_mpc_array (mpc_t *x, int xLen)
{
  char *buffer = (char*) malloc(DEFAULTSTRLEN * sizeof(char));
  for (int i=0; i<xLen; i++)
  {
     buffer = mpc_get_str (10, 0, x[i], MPC_RNDDN);
     printf("%s\n",buffer);
     mpc_free_str(buffer);
  }
}

int main(int argc, char *argv[]) 
{
  int opt;
  FILE *stream = NULL;
  mpc_t *x, *X;
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
  print_mpc_array (x, xLen);

  printf("Data length: %d\n",xLen);
  if (!isPowerofTwo (xLen))
  {
    printf("zeropadding to next power of two\n");
    unsigned int newxLen = nextPowerofTwo(xLen);
    zeropad_mpc_array (x, &xLen, newxLen - xLen);
  }
  
  fft (x, xLen);
  printf("FFT:\n");
  print_mpc_array (x, xLen);  

  clear_mpc_array (x, xLen);
  fclose(stream);
  return 0;
}

