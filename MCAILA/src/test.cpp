
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




#include <functional>
#include "mca.hpp"
#include "lab.hpp"
#include "stats.hpp"
#include <iostream>

template<typename T>
using matrix_t = std::vector<std::vector<T>>;

int main ()
{
  const int PRECISION_F = 23;
  const int PRECISION_D = 52;
  
  /*
  auto A = mcaila::make_matrix<mcaila::wrapper<double, PRECISION_D>>(2,2);
  A[0][0] = mcaila::wrapper<double, PRECISION_D>
    (.2161, 4);
  A[0][1] = mcaila::wrapper<double, PRECISION_D>
    (.1441, 4);
  A[1][0] = mcaila::wrapper<double, PRECISION_D>
    (1.2969, 5);
  A[1][1] = mcaila::wrapper<double, PRECISION_D>
    (.8648, 4);

  auto b = mcaila::make_matrix<mcaila::wrapper<double, PRECISION_D>>(2,1);
  b[0][0] = mcaila::wrapper<double, PRECISION_D>
    (.1440, 4);
  b[1][0] = mcaila::wrapper<double, PRECISION_D>
    (.8642, 4);
  */
  
  
  auto A = mcaila::make_matrix<mcaila::wrapper<double, PRECISION_D>>(24,24);
  for (int i = 0 ; i < 24 ; i++)
    {
      for (int j = 0 ; j < 24 ; j++)
	{
	  if (i == j) A[i][j] = 1;
	  else if (i > j) A[i][j] = -1;
	  else A[i][j] = 0;
	}
    }
  auto b = mcaila::make_matrix<mcaila::wrapper<double, PRECISION_D>>(24,1);
  for (int i = 0 ; i < 24 ; i++)
    {
      b[i][0] = 1;
    }
  
  /*
  std::cout << "Etat initial de A : " << std::endl;
  std::cout << static_cast<double>(A[0][0]) << " "
	    << static_cast<double>(A[0][1]) << std::endl;
  std::cout << static_cast<double>(A[1][0]) << " "
	    << static_cast<double>(A[1][1]) << std::endl;
  std::cout << "\n";
  */

  /*
  std::vector<size_t> pivot =
    mcaila::LU_factor<mcaila::wrapper<double, PRECISION>>(A);
  
  mcaila::LU_solve<mcaila::wrapper<double, PRECISION>> (A, x, pivot);
  */
  /*  
  auto A = mcaila::make_matrix<double>(2, 2);
  A[0][0] = .2161;
  A[0][1] = .1441;
  A[1][0] = 1.2969;
  A[1][1] = .8648;

  auto b = mcaila::make_matrix<double>(2, 1);
  b[0][0] = .1440;
  b[1][0] = .8642;
  
  auto x = mcaila::LU_mixte<float, double>(A,b);
  

  
  mcaila::print<double>(x);
  */
  
  auto f = mcaila::LU_mixte<mcaila::wrapper<float, PRECISION_F>,
			    mcaila::wrapper<double, PRECISION_D>>;
  auto s = mcaila::LU<mcaila::wrapper<float, PRECISION_F>>;
  auto d = mcaila::LU<mcaila::wrapper<double, PRECISION_D>>;
  
  auto x = mcaila::stat_analysis<mcaila::wrapper<double, PRECISION_D>>
		(f, A, b, 100);

  std::cout << "Average : " << std::endl;
  mcaila::print<mcaila::wrapper<double, PRECISION_D>>(x.first);
  std::cout << "Standard error : " << std::endl;
  mcaila::print<mcaila::wrapper<double, PRECISION_D>>(x.second);
  
  /*
  std::cout << "Average : " << std::endl;
  for (int i = 0 ; i < 2 ; i++)
    {
      std::cout << static_cast<double>((x.first)[i][0]) << std::endl;
    }
  std::cout << "\n";
  std::cout << "Standard error : " << std::endl;
  for (int i = 0 ; i < 2 ; i++)
    {
      std::cout << static_cast<double>((x.second)[i][0]) << std::endl;
    }
  */
  return 0;
}
