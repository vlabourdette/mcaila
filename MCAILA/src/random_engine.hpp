
 /**************************************************************/
 /*                                                            */
 /*  Monte Carlo Arithmetic Implementation for Linear Algebra  */
 /*                                                            */
 /*  Copyright (C) 2015                                        */
 /*  Marc Baboulin :                                           */
 /*  marc dot baboulin at gmail dot com                        */
 /*  Christophe Denis :                                        */
 /*  christophe dot denis at cmla dot ens dash cachan dot fr   */
 /*  Amal Khabou :                                             */
 /*  amal dot khabou at lri dot fr                             */
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





#ifndef __RANDOM_ENGINE_H__
#define __RANDOM_ENGINE_H__

#include <vector>
#include <chrono>
#include <random>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace mcaila{

  size_t thread_id ()
  {
#if defined(_OPENMP)
    return size_t(omp_get_thread_num());
#else
    return 0;
#endif
  }
  
  std::mt19937_64& random_engine(size_t thread_id = thread_id())
  {
    static std::vector<std::mt19937_64> r =
      []()->std::vector<std::mt19937_64>
      {
#if defined(_OPENMP)
	const size_t N = omp_get_max_threads();
#else
	const size_t N = 1;
#endif
	std::vector<std::mt19937_64> r0 (N);
	for (size_t i = 0 ; i < r0.size() ; ++i)
	  {
	    r0[i] = std::mt19937_64
	      (std::chrono::steady_clock::now().time_since_epoch().count() + i);
	  }
	return r0;
      }();
    return r[thread_id];
  }
  
}

#endif
