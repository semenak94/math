
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

// хочется сделать одну функцию, которая будет подтягивать коэффиуциенты в зависимости от индекса
template <typename T>
static T get_negative_value(const T x, const T k,
                            const T& a, const T& b)
{
    const size_t N_base = 4;
    T S1 = 0, S2 = 0;
    T y = log(1 + exp(x));

    for (size_t n = 0; n < N_base + 1; n++) {
        S1 = S1 + a[n] * pow(y, n + 1);
    }

    for (size_t m = 0; m < N_base; m++) {
        S2 = S2 + b[m] * pow(y, m + 1);
    }

    T I_negative = (1 + S1) / (1 + S2);

    // возвращаем значение аппроксимации согласно формуле (14) в статье
    return boost::math::tgamma(k + 1) * y * pow(I_negative, k);
}

template <typename T>
T fermi_dirak_k1(const T x)
{
    const int k = 1;
    static const T a[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.2715113138214362780964488)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0562661238060587633169245)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0067420740469345689743873)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0005169505155333205859985)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000194771836765773190602)) 
	};

    static const T b[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0215113138214352840651584)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0231105175729721417901084)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0003669081577365413477999)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000610424408732720110769)) 
	};

    if (x < 0) {
        return get_negative_value(x, k, a, b);
    }
    else {
        return x * x / 2 + 2 * I_1_0 - get_negative_value(-x, k, a, b);
    }
}

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