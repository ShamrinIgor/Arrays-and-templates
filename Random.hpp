#ifndef Random_hpp
#define Random_hpp

#include <iostream>
#include <random>
#include <cmath>
#include <complex>
namespace dotFunc
{
    template<typename C, typename R, typename URBG>//UniformRandomBitGenerator
    double RandomAtRadiusR1(C detachment, R radius, URBG& prng)
    {
        std::uniform_real_distribution<> randomX(-radius,radius);
        double NewX = randomX(prng);
        return detachment += NewX;
    }
    
    template <typename C, typename R, typename URBG>
    std::complex<C> RandomAtRadiusC(std::complex<C> z, R radius, URBG& prng)
    {
        using namespace std::complex_literals;
        std::uniform_real_distribution<> randomReal(-radius,radius);
        double NewReal = randomReal(prng);
        
        std::uniform_real_distribution<> randomImag(-radius,radius);
        double NewImag = randomImag(prng);
        while (sqrt(NewReal*NewReal + NewImag*NewImag) > radius)
        {
            NewImag = randomImag(prng);
        }
        std::complex<C> NewZ(NewReal, NewImag);
        z += NewZ;
        return z;
    }
}

#endif
