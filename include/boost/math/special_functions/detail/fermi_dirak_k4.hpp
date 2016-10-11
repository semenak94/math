//  Copyright (c) 2016 Semyon Kolganov
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FERMI_DIRAK_K4_HPP
#define BOOST_MATH_FERMI_DIRAK_K4_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/assert.hpp>

// New quick metod of calculation Fermi-Dirak function 
// with high precision

namespace boost { namespace math { namespace detail {

template <typename T>
T fermi_dirak_k4(const T x)
{
    const int k = 4;

    const T a[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0560148791230902149024568)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0351117957891800867706741)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0021834386943672331415760)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0002464861525522946634693)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000092228177886669241259))
	};

    const T b[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.0611726208769112866900252)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64,  0.0279968542816146833953639)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.0007512148294307540141223)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64,  0.0000860680747142919882956))
	};

    if (x < 0) {
        return get_negative_value(-x, k, a, b);
    }
    else {
        return pow(x, 5) / 5 + 8 * pow(x, 3) * I_1_0 + 8 * x * I_3_0 +
               get_negative_value(-x, k, a, b);
    }
}

}}} // namespaces

#endif // BOOST_MATH_FERMI_DIRAK_K4_HPP