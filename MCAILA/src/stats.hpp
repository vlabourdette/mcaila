
 /**************************************************************/
 /*                                                            */
 /*  Monte Carlo Arithmetic Implementation for Linear Algebra  */
 /*                                                            */
 /*  Copyright (C) 2015                                        */
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



#ifndef __STATS_H__
#define __STATS_H__

#include <utility>
#include <functional>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "lab.hpp"

template<typename T>
using matrix_t = std::vector<std::vector<T>>;

namespace mcaila{

  template<typename T> std::pair<matrix_t<T>, matrix_t<T>> stat_analysis
  (std::function<matrix_t<T>(matrix_t<T>&, matrix_t<T>&)> method,
   matrix_t<T> A, matrix_t<T> b, int iter)
  {
    size_t n = b.size();
    matrix_t<T> average = make_matrix<T>(n, b[0].size());
    matrix_t<T> results = make_matrix<T>(n, iter);
    matrix_t<T> deviation = make_matrix<T>(n, b[0].size());
    matrix_t<T> temp = make_matrix<T>(n, b[0].size()); 

    
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (size_t i = 0 ; i < iter ; i++)
      {
	temp = method (A, b);
	for (size_t j = 0 ; j < n ; j++)
	  {
	    if (isnan(temp[j][0])) temp[j][0] = 0.;
	    results[j][i] = temp[j][0];
	  }
      }
    for (size_t i = 0 ; i < n ; i++)
      {
	    for (size_t j = 0 ; j < iter ; j++)
	      {
	    average[i][0] += results[i][j];
	  }
	    average[i][0] /= iter;
	    for (size_t j = 0 ; j < iter ; j++)
	      {
	    deviation[i][0] += (results[i][j] - average[i][0])
	      * (results[i][j] - average[i][0]);
	  }
	    deviation[i][0] /= (static_cast<T>(iter - 1));
	    deviation[i][0] = sqrt(deviation[i][0])/sqrt(iter);
      }
    return std::make_pair (average, deviation);
  }
  
    
}

#endif

