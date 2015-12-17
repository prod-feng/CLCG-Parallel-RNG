# CLCG-Parallel-RNG

This is a simple implemenration of multi-threads(OpenMP) parallel Combined Linear Congruential random number generator, MRG32K3A. It is OpenMP multiple threads safe.

>gcc -fopenmp pi_random_parallel.c -lgomp -o pi_random_parallel.exe

Running the program with different number of threads give the identical result:

>export OMP_NUM_THREADS=1
>
>time ./pi_clcg_random_parallel.exe 
>
> 10000000 trials, pi = 3.142059 
>
>real	0m11.107s
>
>user	0m11.103s
>
>sys	0m0.001s
>

>export OMP_NUM_THREADS=2
>
>time ./pi_clcg_random_parallel.exe 
>
> 10000000 trials, pi = 3.141130 
>
>real	0m6.812s
>
>user	0m13.601s
>
>sys	0m0.001s
>

>export OMP_NUM_THREADS=3
>
>time ./pi_clcg_random_parallel.exe 
>
> 10000000 trials, pi = 3.141223 
>
>real	0m4.713s
>
>user	0m14.106s
>
>sys	0m0.003s
>

>export OMP_NUM_THREADS=4
>
>time ./pi_clcg_random_parallel.exe 
>
> 10000000 trials, pi = 3.141506 
>
>real	0m3.512s
>
>user	0m14.008s
>
>sys	0m0.002s


Developed on CentOS 6.2, Kernel 2.6.32-220.2.1.el6.x86_64, gcc (GCC) 4.4.6 20110731
