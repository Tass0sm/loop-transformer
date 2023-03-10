/**
 * gemver.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 4000. */
#include "gemver.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE *alpha,
		 DATA_TYPE *beta,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(u1,N,n),
		 DATA_TYPE POLYBENCH_1D(v1,N,n),
		 DATA_TYPE POLYBENCH_1D(u2,N,n),
		 DATA_TYPE POLYBENCH_1D(v2,N,n),
		 DATA_TYPE POLYBENCH_1D(w,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

  *alpha = 43532;
  *beta = 12313;

  for (i = 0; i < n; i++)
    {
      u1[i] = i;
      u2[i] = (i+1)/n/2.0;
      v1[i] = (i+1)/n/4.0;
      v2[i] = (i+1)/n/6.0;
      y[i] = (i+1)/n/8.0;
      z[i] = (i+1)/n/9.0;
      x[i] = 0.0;
      w[i] = 0.0;
      for (j = 0; j < n; j++)
	A[i][j] = ((DATA_TYPE) i*j) / n;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(w,N,n))
{
  int i;

  for (i = 0; i < n; i++) {
    fprintf (stderr, DATA_PRINTF_MODIFIER, w[i]);
    if (i % 20 == 0) fprintf (stderr, "\n");
  }
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemver(int n,
		   DATA_TYPE alpha,
		   DATA_TYPE beta,
		   DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(u1,N,n),
		   DATA_TYPE POLYBENCH_1D(v1,N,n),
		   DATA_TYPE POLYBENCH_1D(u2,N,n),
		   DATA_TYPE POLYBENCH_1D(v2,N,n),
		   DATA_TYPE POLYBENCH_1D(w,N,n),
		   DATA_TYPE POLYBENCH_1D(x,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n),
		   DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

  /* ppcg generated CPU code */
  
  #define ppcg_min(x,y)    ({ __typeof__(x) _x = (x); __typeof__(y) _y = (y); _x < _y ? _x : _y; })
  {
    #pragma omp parallel for
    #pragma coalesce (c0,c1)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1)
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            A[c2][c3] = ((A[c2][c3] + (u1[c2] * v1[c3])) + (u2[c2] * v2[c3]));
    #pragma omp parallel for
    #pragma coalesce (c0)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 <= n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1) {
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            x[c2] = (x[c2] + ((beta * A[c3][c2]) * y[c3]));
          if (c1 + 31 >= n)
            x[c2] = (x[c2] + z[c2]);
        }
    #pragma omp parallel for
    #pragma coalesce (c0)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1)
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            w[c2] = (w[c2] + ((alpha * A[c2][c3]) * x[c3]));
  }
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemver_opt(int n,
                       DATA_TYPE alpha,
                       DATA_TYPE beta,
                       DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
                       DATA_TYPE POLYBENCH_1D(u1,N,n),
                       DATA_TYPE POLYBENCH_1D(v1,N,n),
                       DATA_TYPE POLYBENCH_1D(u2,N,n),
                       DATA_TYPE POLYBENCH_1D(v2,N,n),
                       DATA_TYPE POLYBENCH_1D(w,N,n),
                       DATA_TYPE POLYBENCH_1D(x,N,n),
                       DATA_TYPE POLYBENCH_1D(y,N,n),
                       DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

  /* ppcg generated CPU code */
{
  #pragma scop
  {
    #pragma omp parallel for
    #pragma coalesce (c0,c1)
    #pragma unroll (c2:4,c3:2), licm
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
      for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 4)
    {
      DATA_TYPE invariant_0 = u1[c2 + 0];
      DATA_TYPE invariant_1 = u2[c2 + 0];
      DATA_TYPE invariant_2 = u1[c2 + 1];
      DATA_TYPE invariant_3 = u2[c2 + 1];
      DATA_TYPE invariant_4 = u1[c2 + 2];
      DATA_TYPE invariant_5 = u2[c2 + 2];
      DATA_TYPE invariant_6 = u1[c2 + 3];
      DATA_TYPE invariant_7 = u2[c2 + 3];
      for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 2)
      {
        A[c2 + 0][c3 + 0] = (A[c2 + 0][c3 + 0] + (invariant_0 * v1[c3 + 0])) + (invariant_1 * v2[c3 + 0]);
        A[c2 + 0][c3 + 1] = (A[c2 + 0][c3 + 1] + (invariant_0 * v1[c3 + 1])) + (invariant_1 * v2[c3 + 1]);
        A[c2 + 1][c3 + 0] = (A[c2 + 1][c3 + 0] + (invariant_2 * v1[c3 + 0])) + (invariant_3 * v2[c3 + 0]);
        A[c2 + 1][c3 + 1] = (A[c2 + 1][c3 + 1] + (invariant_2 * v1[c3 + 1])) + (invariant_3 * v2[c3 + 1]);
        A[c2 + 2][c3 + 0] = (A[c2 + 2][c3 + 0] + (invariant_4 * v1[c3 + 0])) + (invariant_5 * v2[c3 + 0]);
        A[c2 + 2][c3 + 1] = (A[c2 + 2][c3 + 1] + (invariant_4 * v1[c3 + 1])) + (invariant_5 * v2[c3 + 1]);
        A[c2 + 3][c3 + 0] = (A[c2 + 3][c3 + 0] + (invariant_6 * v1[c3 + 0])) + (invariant_7 * v2[c3 + 0]);
        A[c2 + 3][c3 + 1] = (A[c2 + 3][c3 + 1] + (invariant_6 * v1[c3 + 1])) + (invariant_7 * v2[c3 + 1]);
      }

    }



    #pragma omp parallel for
    #pragma coalesce (c0)
    #pragma unroll (c2:4,c3:2), licm
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 <= n; c1 += 32)
      for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 4)
    {
      {
        DATA_TYPE invariant_0 = x[c2 + 0];
        for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 2)
        {
          x[c2 + 0] = invariant_0 + ((beta * A[c3 + 0][c2 + 0]) * y[c3 + 0]);
          x[c2 + 0] = invariant_0 + ((beta * A[c3 + 1][c2 + 0]) * y[c3 + 1]);
        }

      }
      if ((c1 + 31) >= n)
        x[c2 + 0] = x[c2 + 0] + z[c2 + 0];
      for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 2)
      {
        x[c2 + 1] = x[c2 + 1] + ((beta * A[c3 + 0][c2 + 1]) * y[c3 + 0]);
        x[c2 + 1] = x[c2 + 1] + ((beta * A[c3 + 1][c2 + 1]) * y[c3 + 1]);
      }

      if ((c1 + 31) >= n)
        x[c2 + 1] = x[c2 + 1] + z[c2 + 1];
      for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 2)
      {
        x[c2 + 2] = x[c2 + 2] + ((beta * A[c3 + 0][c2 + 2]) * y[c3 + 0]);
        x[c2 + 2] = x[c2 + 2] + ((beta * A[c3 + 1][c2 + 2]) * y[c3 + 1]);
      }

      if ((c1 + 31) >= n)
        x[c2 + 2] = x[c2 + 2] + z[c2 + 2];
      for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 2)
      {
        x[c2 + 3] = x[c2 + 3] + ((beta * A[c3 + 0][c2 + 3]) * y[c3 + 0]);
        x[c2 + 3] = x[c2 + 3] + ((beta * A[c3 + 1][c2 + 3]) * y[c3 + 1]);
      }

      if ((c1 + 31) >= n)
        x[c2 + 3] = x[c2 + 3] + z[c2 + 3];
    }



    #pragma omp parallel for
    #pragma coalesce (c0)
    #pragma unroll (c2:4,c3:2)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
      for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 4)
      for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 2)
    {
      w[c2 + 0] = w[c2 + 0] + ((alpha * A[c2 + 0][c3 + 0]) * x[c3 + 0]);
      w[c2 + 0] = w[c2 + 0] + ((alpha * A[c2 + 0][c3 + 1]) * x[c3 + 1]);
      w[c2 + 1] = w[c2 + 1] + ((alpha * A[c2 + 1][c3 + 0]) * x[c3 + 0]);
      w[c2 + 1] = w[c2 + 1] + ((alpha * A[c2 + 1][c3 + 1]) * x[c3 + 1]);
      w[c2 + 2] = w[c2 + 2] + ((alpha * A[c2 + 2][c3 + 0]) * x[c3 + 0]);
      w[c2 + 2] = w[c2 + 2] + ((alpha * A[c2 + 2][c3 + 1]) * x[c3 + 1]);
      w[c2 + 3] = w[c2 + 3] + ((alpha * A[c2 + 3][c3 + 0]) * x[c3 + 0]);
      w[c2 + 3] = w[c2 + 3] + ((alpha * A[c2 + 3][c3 + 1]) * x[c3 + 1]);
    }




  }
  #pragma endscop
}































}

static
void check (int n,
		 DATA_TYPE POLYBENCH_1D(q1,N,n),
		 DATA_TYPE POLYBENCH_1D(q2,N,n))
{
  int i, j;
  #define abs(x) ((x) >= 0 ? (x) : -(x))
  DATA_TYPE diff = 0.0;
  for (i = 0; i < _PB_N; i++)
    {
      diff += abs(q1[i] - q2[i]);
    }
  if (diff > 0.000001)
    printf ("CHECK FAIL\n");
  else
    printf ("CHECK PASS\n");
}

int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(u1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(v1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(u2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(v2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(w,     DATA_TYPE, N, n);  /* Output Array for Normal Kernel */
  POLYBENCH_1D_ARRAY_DECL(w_opt, DATA_TYPE, N, n);  /* Output Array for Opt Kernel */
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(z, DATA_TYPE, N, n);

#ifdef OPT
  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(u1),
	      POLYBENCH_ARRAY(v1),
	      POLYBENCH_ARRAY(u2),
	      POLYBENCH_ARRAY(v2),
	      POLYBENCH_ARRAY(w_opt),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y),
	      POLYBENCH_ARRAY(z));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemver_opt (n, alpha, beta,
                     POLYBENCH_ARRAY(A),
                     POLYBENCH_ARRAY(u1),
                     POLYBENCH_ARRAY(v1),
                     POLYBENCH_ARRAY(u2),
                     POLYBENCH_ARRAY(v2),
                     POLYBENCH_ARRAY(w_opt),
                     POLYBENCH_ARRAY(x),
                     POLYBENCH_ARRAY(y),
                     POLYBENCH_ARRAY(z));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;
#endif

  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(u1),
	      POLYBENCH_ARRAY(v1),
	      POLYBENCH_ARRAY(u2),
	      POLYBENCH_ARRAY(v2),
	      POLYBENCH_ARRAY(w),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y),
	      POLYBENCH_ARRAY(z));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemver (n, alpha, beta,
		 POLYBENCH_ARRAY(A),
		 POLYBENCH_ARRAY(u1),
		 POLYBENCH_ARRAY(v1),
		 POLYBENCH_ARRAY(u2),
		 POLYBENCH_ARRAY(v2),
		 POLYBENCH_ARRAY(w),
		 POLYBENCH_ARRAY(x),
		 POLYBENCH_ARRAY(y),
		 POLYBENCH_ARRAY(z));

#ifndef OPT
  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;
#endif

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w)));
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w_opt)));

#ifdef OPT
  check (n, POLYBENCH_ARRAY(w), POLYBENCH_ARRAY(w_opt));
#endif

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(u1);
  POLYBENCH_FREE_ARRAY(v1);
  POLYBENCH_FREE_ARRAY(u2);
  POLYBENCH_FREE_ARRAY(v2);
  POLYBENCH_FREE_ARRAY(w);
  POLYBENCH_FREE_ARRAY(w_opt);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(z);

  return 0;
}
