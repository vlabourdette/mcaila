
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




#include "mca.hpp"
#include "lab.hpp"
#include <iostream>

#if defined(_OPENMP)
#include <omp.h>
#endif



int main ()
{
  const int PRECISION = 23;

  /*
    Génération de A : 
    On appelle explicitement le constructeur pour préciser le nombre de 
    chiffres significatifs.
    On pourrait aussi tout simplement écrire A[0][0] = .2161 puisque les
    opérateurs de conversion ont été surchargés.
   */
  auto A = mcaila::make_matrix<mcaila::wrapper<float, PRECISION>>(2,2);
  A[0][0] = mcaila::wrapper<float, PRECISION>
    (.2161, 4);
  A[0][1] = mcaila::wrapper<float, PRECISION>
    (.1441, 4);
  A[1][0] = mcaila::wrapper<float, PRECISION>
    (1.2969, 5);
  A[1][1] = mcaila::wrapper<float, PRECISION>
    (.8648, 4);

  /*
    Pour des raisons de généricité, même les vecteurs sont des matrices
   */
  auto x = mcaila::make_matrix<mcaila::wrapper<float, PRECISION>>(2,1);
  x[0][0] = mcaila::wrapper<float, PRECISION>
    (.1440, 4);
  x[1][0] = mcaila::wrapper<float, PRECISION>
    (.8642, 4);

  std::cout << "Etat initial de A : " << std::endl;
  std::cout << static_cast<float>(A[0][0]) << " "
	    << static_cast<float>(A[0][1]) << std::endl;
  std::cout << static_cast<float>(A[1][0]) << " "
	    << static_cast<float>(A[1][1]) << std::endl;
  std::cout << "\n";

  /*
    LU

    Grâce aux templates, les fonctions LU_factor et LU_solve marchent de la
    même façon pour les types flottants de base et pour leurs wrappers
   */
  std::vector<size_t> pivot =
    mcaila::LU_factor<mcaila::wrapper<float, PRECISION>>(A);
  
  mcaila::LU_solve<mcaila::wrapper<float, PRECISION>> (A, x, pivot);
  
  std::cout << static_cast<float>(x[0][0]) << " "
	    << static_cast<float>(x[1][0]) << std::endl;
  
  return 0;
}
