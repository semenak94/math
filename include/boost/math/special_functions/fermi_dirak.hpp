
//  (C) Copyright Semyon Kolganov  2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FERMI_DIRAK_HPP
#define BOOST_MATH_SPECIAL_FERMI_DIRAK_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/policies/error_handling.hpp>

namespace boost {
namespace math {





template <typename T>
T fermi_dirak_k2(const T x)
{
    const int k = 2;
    const T a[] = { 
		0.2263816364340698560028783,
        0.0533684335574798857246766,
        0.0062904756340795211604491,
        0.0005023228274452983506998,
        0.0000189379675088061004880
	};

    const T b[] = { 
		0.0388816364340691133155655,
        0.0243043998742774445085992,
        0.0006290985326433190105734,
        0.0000657018161945458806177 
	};

    if (x < 0) {
        return get_negative_value(x, k, a, b);
    }
    else {
        return pow(x, 3) / 3 + 4 * x * I_1_0 + 
               get_negative_value(-x, k, a, b);
    }
}

template <typename T>
T fermi_dirak_k3(const T x)
{
    const int k = 3;
    const T a[] = { 
		0.1583482145380455955096383,
        0.0460645149909308107878344,
        0.0048861379108841469134267,
        0.0004336733305971515517559,
        0.0000173435613795895152436 
	};

    const T b[] = { 
		0.0125148812047107612191739,
		0.0266693407000929631393759,
		0.0003285431094547362504004,
		0.0000820910787890062715299 
	};

    if (x < 0) {
        return get_negative_value(x, k, a, b);
    }
    else {
        return pow(x, 4) / 4 + 6 * x * x * I_1_0 + 2 * I_3_0 - 
               get_negative_value(-x, k, a, b);
    }
    
}

template <typename T>
T fermi_dirak_k4(const T x)
{
    const int k = 4;

    const T a[] = { 
		0.0560148791230902149024568,
        0.0351117957891800867706741,
        0.0021834386943672331415760,
        0.0002464861525522946634693,
        0.0000092228177886669241259
	};

    const T b[] = { 
		-0.0611726208769112866900252,
         0.0279968542816146833953639,
        -0.0007512148294307540141223,
         0.0000860680747142919882956
	};

    if (x < 0) {
        return get_negative_value(-x, k, a, b);
    }
    else {
        return pow(x, 5) / 5 + 8 * pow(x, 3) * I_1_0 + 8 * x * I_3_0 +
               get_negative_value(-x, k, a, b);
    }
}


} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_FERMI_DIRAK_HPP