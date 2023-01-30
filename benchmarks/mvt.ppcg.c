/**
 * mvt.c: This file is part of the PolyBench/C 3.2 test suite.
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
#include "mvt.h"


/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE POLYBENCH_1D(x1,N,n),
		DATA_TYPE POLYBENCH_1D(x2,N,n),
		DATA_TYPE POLYBENCH_1D(y_1,N,n),
		DATA_TYPE POLYBENCH_1D(y_2,N,n),
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      x1[i] = ((DATA_TYPE) i) / n;
      x2[i] = ((DATA_TYPE) i + 1) / n;
      y_1[i] = ((DATA_TYPE) i + 3) / n;
      y_2[i] = ((DATA_TYPE) i + 4) / n;
      for (j = 0; j < n; j++)
	A[i][j] = ((DATA_TYPE) i*j) / N;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(x1,N,n),
		 DATA_TYPE POLYBENCH_1D(x2,N,n))

{
  int i;

  for (i = 0; i < n; i++) {
    fprintf (stderr, DATA_PRINTF_MODIFIER, x1[i]);
    fprintf (stderr, DATA_PRINTF_MODIFIER, x2[i]);
    if (i % 20 == 0) fprintf (stderr, "\n");
  }
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_mvt(int n,
		DATA_TYPE POLYBENCH_1D(x1,N,n),
		DATA_TYPE POLYBENCH_1D(x2,N,n),
		DATA_TYPE POLYBENCH_1D(y_1,N,n),
		DATA_TYPE POLYBENCH_1D(y_2,N,n),
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  /* ppcg generated CPU code */
  
  #define ppcg_min(x,y)    ({ __typeof__(x) _x = (x); __typeof__(y) _y = (y); _x < _y ? _x : _y; })
  {
    #pragma omp parallel for
    #pragma coalesce (c0)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1)
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            x1[c2] = (x1[c2] + (A[c2][c3] * y_1[c3]));
    #pragma omp parallel for
    #pragma coalesce (c0)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1)
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            x2[c2] = (x2[c2] + (A[c3][c2] * y_2[c3]));
  }

}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_mvt_opt(int n,
                    DATA_TYPE POLYBENCH_1D(x1,N,n),
                    DATA_TYPE POLYBENCH_1D(x2,N,n),
                    DATA_TYPE POLYBENCH_1D(y_1,N,n),
                    DATA_TYPE POLYBENCH_1D(y_2,N,n),
                    DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  /* ppcg generated CPU code */

  #define ppcg_min(x,y)    ({ __typeof__(x) _x = (x); __typeof__(y) _y = (y); _x < _y ? _x : _y; })
  {
    #pragma omp parallel for
    #pragma coalesce (c0)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1)
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            x1[c2] = (x1[c2] + (A[c2][c3] * y_1[c3]));
    #pragma omp parallel for
    #pragma coalesce (c0)
    for (int c0 = 0; c0 < n; c0 += 32)
      for (int c1 = 0; c1 < n; c1 += 32)
        for (int c2 = c0; c2 <= ppcg_min(n - 1, c0 + 31); c2 += 1)
          #pragma unroll (c2:4,c3:2), licm
          for (int c3 = c1; c3 <= ppcg_min(n - 1, c1 + 31); c3 += 1)
            x2[c2] = (x2[c2] + (A[c3][c2] * y_2[c3]));
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
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(x1,     DATA_TYPE, N, n); /* Output Arrays For Normal Kernel */
  POLYBENCH_1D_ARRAY_DECL(x1_opt, DATA_TYPE, N, n); /* Output Arrays For Opt Kernel */
  POLYBENCH_1D_ARRAY_DECL(x2,     DATA_TYPE, N, n); /* Output Arrays For Normal Kernel */
  POLYBENCH_1D_ARRAY_DECL(x2_opt, DATA_TYPE, N, n); /* Output Arrays For Opt Kernel */
  POLYBENCH_1D_ARRAY_DECL(y_1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y_2, DATA_TYPE, N, n);

#ifdef OPT
  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(x1_opt),
	      POLYBENCH_ARRAY(x2_opt),
	      POLYBENCH_ARRAY(y_1),
	      POLYBENCH_ARRAY(y_2),
	      POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_mvt_opt (n,
                  POLYBENCH_ARRAY(x1_opt),
                  POLYBENCH_ARRAY(x2_opt),
                  POLYBENCH_ARRAY(y_1),
                  POLYBENCH_ARRAY(y_2),
                  POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;
#endif

  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(x1),
	      POLYBENCH_ARRAY(x2),
	      POLYBENCH_ARRAY(y_1),
	      POLYBENCH_ARRAY(y_2),
	      POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_mvt_opt (n,
                  POLYBENCH_ARRAY(x1),
                  POLYBENCH_ARRAY(x2),
                  POLYBENCH_ARRAY(y_1),
                  POLYBENCH_ARRAY(y_2),
                  POLYBENCH_ARRAY(A));

#ifndef OPT
  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;
#endif

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x1), POLYBENCH_ARRAY(x2)));
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x1_opt), POLYBENCH_ARRAY(x2_opt)));

#ifdef OPT
  check (n, POLYBENCH_ARRAY(x1), POLYBENCH_ARRAY(x1_opt));
  check (n, POLYBENCH_ARRAY(x2), POLYBENCH_ARRAY(x2_opt));
#endif

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x1);
  POLYBENCH_FREE_ARRAY(x1_opt);
  POLYBENCH_FREE_ARRAY(x2);
  POLYBENCH_FREE_ARRAY(x2_opt);
  POLYBENCH_FREE_ARRAY(y_1);
  POLYBENCH_FREE_ARRAY(y_2);

  return 0;
}
