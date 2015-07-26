
 /**************************************************************/
 /*                                                            */
 /*  Monte Carlo Arithmetic Implementation for Linear Algebra  */
 /*                                                            */
 /*  Copyright (C) 2015 Valentin I. Labourdette :              */
 /*  valentin dot labourdette at gmail dot com                 */
 /*                                                            */
 /*  This file is part of the Monte Carlo Arithmetic           */
 /*  Implementation for Linear Algebra (MCAILA). MCAILA is     */
 /*  free software : you can redistribute it and/or modify it  */
 /*  under the terms of the GNU General Public License as      */
 /*  published by the Free Software Foundation, either         */
 /*  version 3 of the License, or (at your option) any later   */
 /*  version.                                                  */
 /*                                                            */
 /*  This program is distributed in the hope that it will be   */
 /*  useful, but WITHOUT ANY WARRANTY ; without even the       */
 /*  implied warranty of MERCHANTABILITY or FITNESS FOR A      */
 /*  PARTICULAR PURPOSE. See the GNU General Public License    */
 /*  for more details.                                         */
 /*                                                            */
 /*  You should have received a copy of the GNU General        */
 /*  Public License along with this program. If not, see       */
 /*  <http://www.gnu.org/licenses/>.                           */
 /*                                                            */
 /**************************************************************/




#ifndef __LAB_H__
#define __LAB_H__

#include <iostream>
#include <vector>
#include "mca.hpp"

template<typename T>
using matrix_t = std::vector<std::vector<T>>;

namespace mcaila{
  
  template<typename T> std::vector<std::vector<T>> make_matrix
  ( const size_t M, const size_t N)
  { return std::vector<std::vector<T>> (M, std::vector<T>(N)); }

  template<typename T/*, typename matrix_t*/>
  T norm (const matrix_t<T>& A)
  {
    T result = 0;
    T sum;
    for (size_t i = 0 ; i < A.size() ; i++)
      {
	sum = 0;
	for (size_t j = 0 ; j < A.size() ; j++)
	  {
	    sum += A[i][j];
	  }
	if (result < sum) result = sum;
      }
    return result;
  }
  
  template</*typename matrix_t, */typename T>
  void swap_rows (matrix_t<T>& A, size_t i, size_t j)
  {
    auto k = A[i];
    A[i] = A[j];
    A[j] = k;
  }

  template</*typename matrix_t, */typename T>
  matrix_t<T> matrix_product (matrix_t<T>& A, matrix_t<T>& B);
  
  template</*typename matrix_t, */typename T>
  std::vector<size_t> LU_factor (matrix_t<T>& A)
  {
    std::vector<size_t> pivot(A.size()) ;
    for (size_t i = 0 ; i < A.size()-1 ; i++)
      {
	//partial pivoting
	int j_max = i;
	for (size_t j = i+1 ; j < A.size() ; j++)
	  {
	    if (abs(A[j][i])>abs(A[j_max][i]))
	      j_max = j;
	  }
	if (j_max != i)
	  {
	    swap_rows<T>(A, i, j_max);
	  }
  
	pivot[i] = j_max;

	for (size_t j = i+1 ; j < A.size() ; j++)
	  {
	    A[j][i] /= A[i][i];
	  }

	
	for (size_t j = i+1 ; j < A.size() ; j++)
	  {
	    for (size_t k = i+1 ; k < A.size() ; k++)
	      {
		A[j][k] -= (A[j][i]*A[i][k]);
	      }
	  }
	
      }
 
    return pivot;
  }

  template</*typename matrix_t, */typename T>
  void LU_solve (matrix_t<T>& A, matrix_t<T>& b, std::vector<size_t> P)
  {
    for (size_t i = 0 ; i < A.size() ; i++)
      {
	swap_rows<T>(b, i, P[i]);
      }
    for (int i = A.size() - 1 ; i >= 0 ; i--)
      {
	for (int j = i+1 ; j < A.size() ; j++)
	  {
	    b[i][0] -= ( A[i][j] * b[j][0] );
	  }
	b[i][0] /= A[i][i];
      }
    
  }
  
  template</*typename matrix_t, */typename F, typename D>
  matrix_t<D> LU_mixte (matrix_t<D>& A, matrix_t<D>& b, D epsilon)
  {
    int iter = 0, ITER_MAX = 30;
    size_t n = A.size();
    matrix_t<F> lu, y, z;
    matrix_t<D> x, r;
    std::vector<int> pivot (n);
    for (size_t i = 0 ; i < n ; i++)
      {
	for (size_t j = 0 ; j < n ; j++)
	  lu[i][j] = static_cast<F>(A[i][j]);
	y[i][0] = static_cast<F>(b[i][0]);
      }
    pivot = LU_factor<F>(lu);

    LU_solve<F> (lu, y, pivot);
    
    for (size_t i = 0 ; i < n ; i++)
      { x[i][0] = static_cast<D>(y[i][0]); }
    r = matrix_product<D> (A, x);
    while ((norm(r) > epsilon) || (iter > ITER_MAX))
      {
	for (size_t i = 0 ; i < n ; i++)
	  { z[i][0] = static_cast<F>(b[i][0] - r[i][0]); }

	LU_solve<F>(lu, z, pivot);

	for (size_t i = 0 ; i < n ; i++)
	  { x[i][0] -= static_cast<D>(z[i][0]); }
	r = matrix_product<D>> (A, x);
      }
    
  }
  
  
  
}

#endif