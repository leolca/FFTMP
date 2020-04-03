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
#include "FFTMP.h"

mpz_t *fftbincoeff (unsigned int n) 
{
  unsigned int n2 = nextPowerofTwo((n+1));
  mpz_t *c = malloc((n+1)*sizeof(mpz_t));
  mpc_t *x = malloc(n2*sizeof(mpc_t));
  mpfr_t cr; 
  mpfr_inits2 (NUMBITS, cr, (mpfr_ptr) 0);
  for (int i = 0; i < n2; i++)
  {
    mpc_init2 (x[i], NUMBITS);
    if (i < n+1) mpz_init (c[i]);
    if (i < 2) mpc_set_ui(x[i], 1, MPC_RNDDN);
    else mpc_set_ui(x[i], 0, MPC_RNDDN);
  }

  mpc_t *X = fft (x, n2);
  clear_mpc_array (x, n2);
  for (int i = 0; i < n2; i++)
  {
     mpc_pow_ui (X[i], X[i], n, MPC_RNDDN);
  }
  x = ifft (X, n2);
  clear_mpc_array (X, n2);
  for (int i = 0; i < n+1; i++)
  {
    mpc_real (cr , x[i], MPFR_RNDN);
    mpz_set_ui (c[i], mpfr_get_ui (cr, MPFR_RNDN));
  }
  clear_mpc_array (x, n2);
  mpfr_clear (cr);

  return c;
}

