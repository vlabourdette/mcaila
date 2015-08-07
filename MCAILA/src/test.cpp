
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




#include <functional>
#include "mca.hpp"
#include "lab.hpp"
#include "stats.hpp"
#include <iostream>

int main ()
{
 
  const int PRECISION_F = 23;
  const int PRECISION_D = 51;

  /* Première matrice de Kahan (dans le Parker) */
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
  
  /* Matrice de Wilkinson */
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

  /* Second membre égal à 1 */
  auto b = mcaila::make_matrix<mcaila::wrapper<double, PRECISION_D>>(24,1);
  for (int i = 0 ; i < 24 ; i++)
    {
      b[i][0] = 1;
    }
  


  /* LU mixte */
  auto f = mcaila::LU_mixte<mcaila::wrapper<float, PRECISION_F>,
			    mcaila::wrapper<double, PRECISION_D>>;
  /* LU simple */
  auto s = mcaila::LU<mcaila::wrapper<float, PRECISION_F>>;
  /* LU double */
  auto d = mcaila::LU<mcaila::wrapper<double, PRECISION_D>>;

  /* 
     Fonction d'analyse statistique avec MCA

     Prend en paramètre une méthode de résolution (ici LU mixte),
     un système linéaire et un nombre d'itérations MCA
     
     Pour lancer LU double, remplacer f par d
     Pour lancer LU simple, remplacer f par s, double par float et
     PRECISION_D par PRECISION_F

     Renvoie un std::pair qui contient les moyennes (x.first)
     et les déviations (x.second) pour chaque composante
  */
  auto x = mcaila::stat_analysis<mcaila::wrapper<double, PRECISION_D>>
		(f, A, b, 100);

  std::cout << "Average : " << std::endl;
  mcaila::print<mcaila::wrapper<double, PRECISION_D>>(x.first);
  std::cout << "Standard error : " << std::endl;
  mcaila::print<mcaila::wrapper<double, PRECISION_D>>(x.second);

  return 0;
}
