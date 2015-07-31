
 /**************************************************************/
 /*                                                            */
 /*  Monte Carlo Arithmetic Implementation for Linear Algebra  */
 /*                                                            */
 /*  Copyright (C) 2015  
 /*  Christophe Denis :                                        */
 /*  christophe dot denis at cmla dot ens dash cachan dot fr   */
 /*  and Valentin Labourdette :                                */
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
	for (size_t j = 0 ; j < (A[0]).size() ; j++)
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
  matrix_t<T> matrix_product (matrix_t<T>& A, matrix_t<T>& B)
  { 
    size_t m = A.size();
    size_t p = (B[0]).size();
    matrix_t<T> AB = make_matrix<T> (m,p);
    size_t n = (A[0]).size(); 
    for (size_t i = 0 ; i < m ; i++)
      {
	for (size_t j = 0 ; j < p ; j++)
	  {
	    T temp = 0;
	    for (size_t k = 0 ; k < n ; k++)
	      {
		temp += A[i][k]*B[k][j];
	      }
	    AB[i][j] = temp;
	  }
      }
    return AB;
  }
  
  template</*typename matrix_t, */typename T>
  std::vector<size_t> LU_factor (matrix_t<T>& A)
  {
    std::vector<size_t> pivot (A.size());
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
	pivot[A.size() - 1] = A.size() - 1;
	
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

  template<typename T> matrix_t<T> LU (matrix_t<T>& A, matrix_t<T>& b)
  {
    size_t n = A.size();

    matrix_t<T> lu = make_matrix<T>(A.size(), A[0].size());
    matrix_t<T> x = make_matrix<T>(b.size(), b[0].size());	
    for (size_t i = 0 ; i < n ; i++)
      {
	for (size_t j = 0 ; j < n ; j++)
	  lu[i][j] = A[i][j];
	x[i][0] = b[i][0];
      }
    std::vector<size_t> pivot = LU_factor<T>(lu);
    LU_solve<T> (lu, x, pivot);
    return x;
  }
  
  template</*typename matrix_t, */typename F, typename D>
  matrix_t<D> LU_mixte (matrix_t<D>& A, matrix_t<D>& b, D epsilon)
  {
    int iter = 0, ITER_MAX = 30;
    size_t n = A.size();
    matrix_t<F> lu = make_matrix<F>(A.size(), A[0].size());
    matrix_t<F> y = make_matrix<F>(b.size(), b[0].size());
    matrix_t<F> z = make_matrix<F>(b.size(), b[0].size());
    matrix_t<D> x = make_matrix<D>(b.size(), b[0].size());
    matrix_t<D> r = make_matrix<D>(b.size(), b[0].size());

    for (size_t i = 0 ; i < n ; i++)
      {
	for (size_t j = 0 ; j < n ; j++)
	  lu[i][j] = static_cast<F>(A[i][j]);
	y[i][0] = static_cast<F>(b[i][0]);
      }
    std::vector<size_t> pivot = LU_factor<F>(lu);
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
	r = matrix_product<D> (A, x);
      }
    return x;
  }
  
  
  
}

#endif
