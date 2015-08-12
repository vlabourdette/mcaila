
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




#include <functional>
#include "mca.hpp"
#include "lab.hpp"
#include "stats.hpp"
#include "libmat.hpp"
#include <iostream>


int main ()
{
 
  const int PRECISION_F = 23;
  const int PRECISION_D = 51;

  /* Première matrice de Kahan (dans le Parker) */
  /*  
  auto A = mcaila::Parker1<mcaila::wrapper<double, PRECISION_D>>();

  auto b = mcaila::Parker1rhs<mcaila::wrapper<double, PRECISION_D>>();
  */
  
  
  /* Matrice de Wilkinson */
  /*
  auto A = mcaila::Wilkinson<mcaila::wrapper<double, PRECISION_D>>(24);
  */
  
  /* Second membre égal à 1 */
  /*
  auto b = mcaila::make_matrix<mcaila::wrapper<float, PRECISION_F>>(24,1);
  for (int i = 0 ; i < 24 ; i++)
    {
      b[i][0] = 1;
    }
  */


  /* Matrice aléatoire */
  auto A = mcaila::aleatoire<mcaila::wrapper<double, PRECISION_D>>
    (100, pow(10,10));

  /* Second membre égal à 1 */
  auto b = mcaila::make_matrix<mcaila::wrapper<double, PRECISION_D>>(100,1);
  for (int i = 0 ; i < 100 ; i++)
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
		(d, A, b, 100);

  std::cout << "Average : " << std::endl;
  mcaila::print<mcaila::wrapper<double, PRECISION_D>>(x.first);
  std::cout << "Standard error : " << std::endl;
  mcaila::print<mcaila::wrapper<double, PRECISION_D>>(x.second);

  return 0;
}
