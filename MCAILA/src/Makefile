
 #**************************************************************/
 #*                                                            */
 #*  Monte Carlo Arithmetic Implementation for Linear Algebra  */
 #*                                                            */
 #*  Copyright (C) 2015                                        */
 #*  Marc Baboulin :                                           */
 #*  marc dot baboulin at gmail dot com                        */
 #*  Christophe Denis :                                        */
 #*  christophe dot denis at cmla dot ens dash cachan dot fr   */
 #*  Amal Khabou :                                             */
 #*  amal dot khabou at lri dot fr                             */
 #*  and Valentin Labourdette :                                */
 #*  valentin dot labourdette at gmail dot com                 */
 #*                                                            */
 #*  This file is part of the Monte Carlo Arithmetic           */
 #*  Implementation for Linear Algebra (MCAILA). MCAILA is     */
 #*  free software : you can redistribute it and/or modify it  */
 #*  under the terms of the GNU General Public License as      */
 #*  published by the Free Software Foundation, either         */
 #*  version 3 of the License, or (at your option) any later   */
 #*  version.                                                  */
 #*                                                            */
 #*  This program is distributed in the hope that it will be   */
 #*  useful, but WITHOUT ANY WARRANTY ; without even the       */
 #*  implied warranty of MERCHANTABILITY or FITNESS FOR A      */
 #*  PARTICULAR PURPOSE. See the GNU General Public License    */
 #*  for more details.                                         */
 #*                                                            */
 #*  You should have received a copy of the GNU General        */
 #*  Public License along with this program. If not, see       */
 #*  <http://www.gnu.org/licenses/>.                           */
 #*                                                            */
 #**************************************************************/




CC=g++ -std=c++11 -fopenmp

RELEASE=-Ofast

DEBUG=-Og -ggdb3

all: test

test: test.cpp
	$(CC) $(RELEASE) -o test test.cpp && mv test ../bin/


clear:
	rm *.o ../bin/test
