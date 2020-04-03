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
#include "FFTMPbincoeff.h"

// compile: gcc FFTMPbincoeff.c -std=c99 -lmpfr -lgmp -lmpc -o FFTMPbincoeff

int main(int argc, char *argv[]) 
{
  for (unsigned int i = 0; i < 101; i++)
  {
    printf("n=%d : ",i);
    mpz_t *c = fftbincoeff (i);
    print_mpz_array (c, i+1);
    clear_mpz_array (c, i+1);
  }
  return 0;
}

