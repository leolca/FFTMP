# FFTMP
Multi-Precision FFT (implementation in C)

Implementation of FFT using Cooley–Tukey algorithm and GNU MPC (Multiple Precision Complex Library).

Compile:
```bash
$ gcc FFTMP.c -std=c99 -lmpfr -lgmp -lmpc -o FFTMP
```

Help:
```bash
$ ./FFTMP -h
Usage: FFTMP [OPTIONS]... [FILE]...
Output the FFT from a given FILE to standard output.
With no FILE, or when FILE is -, read standard input.
  -i, --input          input filename (if not provided, read from stdin)
  -h, --help           display this help and exit

Example:
echo -n '1 1 1 1 0 0 0 0' | ./FFTMP
./FFTMP -i x_file
```

Usage example:
```bash
$ echo -n "1 1 1 1 0 0 0 0" | ./FFTMP
Data:
(1.0 +0)   
(1.0 +0)         
(1.0 +0)      
(1.0 +0)      
(+0 +0)     
(+0 +0)            
(+0 +0)     
(+0 +0)     
Data length: 8                            
FFT:                               
(4.0 +0)                                                                                        
(9.9e-1 -2.4)                                                                                   
(-0 +0)                                                                                         
(9.9e-1 -4.1e-1)                                                                                
(-0 +0)
(9.9e-1 4.1e-1)
(-0 +0)                            
(1.0 2.4)    
iFFT:
(9.9e-1 -0)
(9.9e-1 3.9e-309)
(1.0 -2.2e-78)
(1.0 -2.2e-78)
(2.1e-78 +0)
(2.1e-78 -3.9e-309)
(-0 2.2e-78)
(-0 2.2e-78)
```

## FFTMPbincoeff
FFTMPbincoeff uses FFTMP to compute binomial coefficients.


```bash
$ ./FFTMPbincoeff 11  
n=0 : 1                      
n=1 : 1 1                          
n=2 : 1 2 1                            
n=3 : 1 3 3 1                         
n=4 : 1 4 6 4 1                                                                                 
n=5 : 1 5 10 10 5 1                  
n=6 : 1 6 15 20 15 6 1                                                                          
n=7 : 1 7 21 35 35 21 7 1             
n=8 : 1 8 28 56 70 56 28 8 1         
n=9 : 1 9 36 84 126 126 84 36 9 1 
n=10 : 1 10 45 120 210 252 210 120 45 10 1 
```
