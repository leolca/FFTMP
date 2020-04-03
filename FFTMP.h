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
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

#define DEFAULTLEN 256
#define DEFAULTSTRLEN 256
#define NUMBITS 1024

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

void zeropad_mpc_array (mpc_t * x, unsigned int *xLen, unsigned int n)
{
   for (int i = 0; i < n; i++)
   {
      mpc_init2 (x[(*xLen)+i], NUMBITS);
      mpc_set_ui (x[(*xLen)+i], 0, MPC_RNDDN);
   }
   *xLen += n;
}

void print_mpc (mpc_t x, size_t n)
{
  char *buffer = mpc_get_str (10, n, x, MPC_RNDDN);
  printf("%s\n",buffer);
  mpc_free_str(buffer); 
}

void print_mpc_array (mpc_t *x, int xLen, size_t n)
{
  for (int i=0; i<xLen; i++)
  {
     print_mpc (x[i], n);
  }
}

void print_mpz_array (mpz_t *x, unsigned int n)
{  
  for (int i = 0; i < n; i++)
    gmp_printf ("%Zd ", x[i]);
  printf("\n");
}

void clear_mpz_array (mpz_t *x, unsigned int n)
{
  for (int i = 0; i < n; i++)  mpz_clear (x[i]);
  free(x);
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
     mpc_init2 (r, NUMBITS);
     mpc_init2 (I, NUMBITS);
     mpfr_inits2 (NUMBITS, PI, (mpfr_ptr) 0);
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

mpc_t *fft(mpc_t *x, unsigned int n)
{
   mpc_t *X = malloc(n * sizeof(mpc_t));
   for (int i = 0; i < n; i++) 
   {
      mpc_init2 (X[i], NUMBITS);
      mpc_set(X[i], x[i], MPC_RNDDN);
   }
   _fft(X, x, n, 1);
   return X;
}

mpc_t * ifft(mpc_t *X, unsigned int n)
{
   mpc_t *x = malloc(n * sizeof(mpc_t));
   for (int i = 0; i < n; i++)
   {
      mpc_init2 (x[i], NUMBITS);
      mpc_conj(X[i], X[i], MPC_RNDDN);
      mpc_set(x[i], X[i], MPC_RNDDN);
   }
   _fft(x, X, n, 1);
   for (int i = 0; i < n; i++)
   {
      mpc_conj(x[i], x[i], MPC_RNDDN);
      mpc_div_ui(x[i], x[i], n, MPC_RNDDN);
   }
   return x;
}

