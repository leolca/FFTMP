# FFTMP
Multi-Precision FFT (implementation in C)

Implementation of FFT using Cooley–Tukey algorithm and GNU MPC (Multiple Precision Complex Library).

Compile:
```bash
$ gcc FFTMP.c -std=c99 -lmpfr -lgmp -lmpc -o FFTMP
```

Usage example:
```bash
$ echo -n "1 1 1 1 0 0 0 0" | ./FFTMP
$ ./FFTMP -i data_file
```
