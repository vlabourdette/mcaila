
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





#ifndef __MCA_H__
#define __MCA_H__

#include <cmath>
#include <random>
#include <memory>
#include <functional>
#include "random_engine.hpp"
#include "iostream"

namespace mcaila{

  template<typename IEEE_type, int PRECISION>
  class wrapper
  { 
  public:
    
    wrapper ()
      : value(0),
	digits(std::numeric_limits<uint>::infinity())
    {}
    
    wrapper (const wrapper<IEEE_type, PRECISION>& x)
      : value(x.value),
	digits(x.digits)
    {}
    wrapper (IEEE_type x)
      : value(x),
	digits(0)
    {}
    wrapper (IEEE_type x, int s)
      : value(x),
	digits(s)
    {}
    template<typename castable> wrapper (castable x)
      : value(static_cast<IEEE_type>(x)),
	digits(0)
    {}
    template<typename castable> wrapper (castable x, int s)
      : value(static_cast<IEEE_type>(x)),
	digits(s)
    {}
    


    operator IEEE_type() const
    {
      return value;
    }

    template<typename castable> operator castable() const 
    {
      return static_cast<castable>(value);
    }

    template<typename castable>
    operator wrapper<castable, PRECISION>() const
    {
      wrapper<castable, PRECISION> y (static_cast<castable>(value),
					   digits);
      return y;
    }

    
    
    friend wrapper<IEEE_type, PRECISION> operator +
    (const wrapper<IEEE_type, PRECISION>& x,
     const wrapper<IEEE_type, PRECISION>& y)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      z += y;
      return z;
    }
    
    friend wrapper<IEEE_type, PRECISION> operator -
    (const wrapper<IEEE_type, PRECISION>& x,
     const wrapper<IEEE_type, PRECISION>& y)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      z -= y;
      return z;
    }
    friend wrapper<IEEE_type, PRECISION> operator +
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      wrapper<IEEE_type, PRECISION> un (1, z.digits+1);
      z *= un;
      return z;
    }
    friend wrapper<IEEE_type, PRECISION> operator -
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      wrapper<IEEE_type, PRECISION> moinsun (-1, z.digits+1);
      z *= moinsun;
      return z;
    }
    friend wrapper<IEEE_type, PRECISION> operator *
    (const wrapper<IEEE_type, PRECISION>& x,
     const wrapper<IEEE_type, PRECISION>& y)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      z *= y;
      return z;
    }
    friend wrapper<IEEE_type, PRECISION> operator /
    (const wrapper<IEEE_type, PRECISION>& x,
     const wrapper<IEEE_type, PRECISION>& y)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      z /= y;
      return z;
    }
    friend wrapper<IEEE_type, PRECISION> operator %
    (const wrapper<IEEE_type, PRECISION>& x,
     const wrapper<IEEE_type, PRECISION>& y)
    {
      wrapper<IEEE_type, PRECISION> z = x;
      z %= y;
      return z;
    }
    wrapper<IEEE_type, PRECISION>& operator ++ ()
    {
      wrapper<IEEE_type, PRECISION> z (1, z.digits+1);
      this->operator+=(z);
      return this;
    }
    friend wrapper<IEEE_type, PRECISION> operator ++
    (wrapper<IEEE_type, PRECISION>& x, int n)
    {
      return ++x;
    }
    wrapper<IEEE_type, PRECISION>& operator -- ()
    {
      wrapper<IEEE_type, PRECISION> z (1, z.digits+1);
      this->operator-=(z);
      return this;
    }
    friend wrapper<IEEE_type, PRECISION> operator --
    (wrapper<IEEE_type, PRECISION>& x, int n)
    {
      return --x;
    }
    
    friend bool operator == (const wrapper<IEEE_type, PRECISION>& x,
			     const wrapper<IEEE_type, PRECISION>& y)
    {
      if (x.digits != y.digits) return false;
      if (x.value != y.value) return false;
      return true;
    }
    friend bool operator != (const wrapper<IEEE_type, PRECISION>& x,
			     const wrapper<IEEE_type, PRECISION>& y)
    {
      return (!(x==y));
    }
    friend bool operator > (const wrapper<IEEE_type, PRECISION>& x,
			    const wrapper<IEEE_type, PRECISION>& y)
    {
      return (x.value > y.value);
    }
    friend bool operator < (const wrapper<IEEE_type, PRECISION>& x,
			    const wrapper<IEEE_type, PRECISION>& y)
    {
      return (x.value < y.value);
    }
    friend bool operator >= (const wrapper<IEEE_type, PRECISION>& x,
			     const wrapper<IEEE_type, PRECISION>& y)
    {
      return (x.value >= y.value);
    }
    friend bool operator <= (const wrapper<IEEE_type, PRECISION>& x,
			     const wrapper<IEEE_type, PRECISION>& y)
    {
      return (x.value <= y.value);
    }

    friend std::ostream& operator<< (std::ostream& os,
				     const wrapper<IEEE_type, PRECISION>& x)
    {
      os << x.value;
      return os;
    }
  
    
    wrapper<IEEE_type, PRECISION>& operator +=
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      this->t_digit_precision();
      wrapper<IEEE_type, PRECISION> x_copy = x;
      x_copy.t_digit_precision();
      value += x_copy.value;
      if (x.digits < digits) digits = x.digits;
      this->t_digit_precision();
      return *this;
    }
    wrapper<IEEE_type, PRECISION>& operator -=
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      this->t_digit_precision();
      wrapper<IEEE_type, PRECISION> x_copy = x;
      x_copy.t_digit_precision();
      value -= x_copy.value;
      if (x.digits < digits) digits = x.digits;
      this->t_digit_precision();
      return *this;
    }
    wrapper<IEEE_type, PRECISION>& operator *=
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      this->t_digit_precision();
      wrapper<IEEE_type, PRECISION> x_copy = x;
      x_copy.t_digit_precision();
      value *= x_copy.value;
      if (x.digits < digits) digits = x.digits;
      this->t_digit_precision();
      return *this;
    }  
    wrapper<IEEE_type, PRECISION>& operator /=
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      this->t_digit_precision();
      wrapper<IEEE_type, PRECISION> x_copy = x;
      x_copy.t_digit_precision();
      value /= x_copy.value;
      if (x.digits < digits) digits = x.digits;
      this->t_digit_precision();
      return *this;
    }
    wrapper<IEEE_type, PRECISION>& operator %=
    (const wrapper<IEEE_type, PRECISION>& x)
    {
      this->t_digit_precision();
      wrapper<IEEE_type, PRECISION> x_copy = x;
      x_copy.t_digit_precision();
      value %= x_copy.value;
      if (x.digits < digits) digits = x.digits;
      this->t_digit_precision();
      return *this;
    }
    
    static std::uniform_real_distribution<IEEE_type> DR;
    
  private:
    IEEE_type value;
    uint digits;
    
    void t_digit_precision ()
    {
      static int bias = 0;
      if (digits<PRECISION)
	{
	  IEEE_type xeta = DR(random_engine());
	  inexact(PRECISION, xeta);
	}
      bias++;
    }
    void inexact (int s, IEEE_type xeta)
    {
      value += pow(2,oom()-s) * xeta;
    }
    int oom ()
    {
      if (value != 0)
	{
	  return (int)(log(fabs(value))) + 1;
	}
      else return 0;
    }
    
  };

  template<typename IEEE_type, int PRECISION>
  std::uniform_real_distribution<IEEE_type> wrapper<IEEE_type, PRECISION>::DR
  = std::uniform_real_distribution<IEEE_type> (-.5, .5);
    
}

namespace std{

  template<typename IEEE_type, int PRECISION>
  mcaila::wrapper<IEEE_type, PRECISION> sqrt
  (mcaila::wrapper<IEEE_type, PRECISION> x)
  {
    mcaila::wrapper<IEEE_type, PRECISION> y
      (sqrt(static_cast<IEEE_type>(x)));
    return y;
  }

  template<typename IEEE_type, int PRECISION>
  mcaila::wrapper<IEEE_type, PRECISION> fabs
  (mcaila::wrapper<IEEE_type, PRECISION> x)
  {
    if (static_cast<IEEE_type>(x) < 0) return -x;
    else return x;
  }
  
}



#endif
