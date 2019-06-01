// VARIANT 8 v2.0
#include <iostream>
#include <random>
#include <complex>
#include "Random.hpp"
#include <array>
using namespace std::complex_literals;
const std::size_t QUANITY = 17; // Число отрядов. 17 по военной доктрине
template <typename C, std::size_t N>
C max_detachment_dist(const std::array<C, N>&);

template <typename C, std::size_t N>
C coord_lucky_detachment (const std::array<C, N>&, C );

template <typename C, std::size_t N>
std::complex<C> max_detachment_dist_Cmp (const std::array<std::complex<C>, N>&);

template <typename C, std::size_t N>
std::complex<C> coord_lucky_detachment_cmp (const std::array<std::complex<C>, N>&, C);

int main()
{
    int D; // Максимальное расстояние высадки отрядов от добрых прищельцев
    int d; // Радиус от точки высадки
    std::mt19937 prng {std::random_device{}()};
    std::cout << "Input D: ";
    if (!(std::cin >> D))
    {
        std::cerr << "Input error!\n";
        return 1;
    }
    if (D <= 0)
    {
        std::cerr << "Input error!\n";
        return 1;
    }
    std::cout << "Input d: ";
    if (!(std::cin >> d))
    {
        std::cerr << "Input error!\n";
        return 1;
    }
    if (d <= 0)
    {
        std::cerr << "Input error!\n";
        return 1;
    }
    /////////////// R1 //////////////
    std::array<double, QUANITY> detachments;
    //double max_detachment_dist(std::array<double, QUANITY>&);
    //double coord_lucky_detachment (std::array<double, QUANITY>&, double );
    double center_of_det = dotFunc::RandomAtRadiusR1(0, D, prng); // Нарандомили центр высадки отрядов
    
    for (std::size_t i = 0; i < QUANITY; i++)
    {
        detachments[i] = dotFunc::RandomAtRadiusR1(center_of_det, d, prng); // Нарандомили 17 отрядов в радиуса d от центра высадки
        std::cout<< "Detachment nubmer "<< i << " x=" << detachments[i] << ' ' << '\n' << '\n' << '\n';
    }
    double optimazed_distance = max_detachment_dist (detachments);
    std::cout << "Optimized puls Distance: " << optimazed_distance << '\n';
    std::cout << "Lucky detachment: x=" << coord_lucky_detachment(detachments, optimazed_distance) << '\n';
    
    /////////////////////////////////
    
    
    /////////////// C ///////////////
    std::complex<double> z(0.,0.);
    std::complex<double> complex_center_of_det = dotFunc::RandomAtRadiusC(z, D, prng);
    std::array<std::complex<double>, QUANITY> complex_detachments;
    //std::complex<double> max_detachment_dist_Cmp (std::array<std::complex<double>, QUANITY>&);
    //std::complex<double> coord_lucky_detachment_cmp (std::array<std::complex<double>, QUANITY>&, double);
    for (std::size_t i = 0; i < QUANITY; i++)
    {
        complex_detachments[i] = dotFunc::RandomAtRadiusC(complex_center_of_det, d, prng); // Нарандомили 17 отрядов в радиуса d от центра высадки
        std::cout << "Detachment nubmer "<< i << " Real(Z)=" << complex_detachments[i].real() << ' ' << " Imag(Z)=" << complex_detachments[i].imag() << '\n' << '\n' << '\n';
    }
    double optimazed_distance_cmp = abs(max_detachment_dist_Cmp (complex_detachments));
    std::cout << "Optimized puls Distance: " << optimazed_distance_cmp << '\n';
    std::cout << "Lucky detachment: " << coord_lucky_detachment_cmp (complex_detachments, optimazed_distance_cmp) << '\n';
    /////////////////////////////////
    return 0;
}

template <typename C, std::size_t N>
C max_detachment_dist (const std::array<C,N>& x)
{
    C max = abs(x[0]);
    C arr = x[0];
    for(std::size_t i = 1; i < QUANITY; i++)
    {
        if(abs(x[i]) > max)
        {
            max = abs(x[i]);
            arr = x[i];
        }
    }
    return arr;
}

template <typename C, std::size_t N>
C coord_lucky_detachment (const std::array<C, N>& x, C OptimazedDistance)
{
    std::size_t count = 0;
    std::array<C, QUANITY> DistanceDif;
    for (std::size_t i = 0; i < QUANITY; i++)
    {
        DistanceDif[i] = std::abs(OptimazedDistance/2 - x[i]);
    }
    C min = DistanceDif[0];
    for (std::size_t i = 1; i < QUANITY; i++)
    {
        if(DistanceDif[i] <= min)
        {
            min = DistanceDif[i];
            count = i;
        }
    }
    return x[count];
}
template <typename C, std::size_t N>
std::complex<C> max_detachment_dist_Cmp (const std::array<std::complex<C>, N>& z)
{
    C D[QUANITY];
    for(std::size_t i = 0; i < QUANITY; i++)
    {
        D[i] = abs(z[i]);
    }
    C max = D[0];
    std::complex<C> max_detachment_coord = z[0];
    for(std::size_t j = 1; j < QUANITY; j++)
    {
        if (D[j] > max)
        {
            max = D[j];
            max_detachment_coord = z[j];
        }
    }
    return max_detachment_coord;
}
template <typename C, std::size_t N>
std::complex<C> coord_lucky_detachment_cmp (const std::array<std::complex<C>, N>& z, C OptimazedDistance)
{
    std::complex<C> result = z[0];
    C D = OptimazedDistance/2;
    C DistanceDif[QUANITY];
    for (std::size_t i = 0; i < QUANITY; i++)
    {
        DistanceDif[i] = abs( abs(z[i]) - D);
    }
    C min = DistanceDif[0];
    for (std::size_t i = 1; i < QUANITY; i++)
    {
        if(DistanceDif[i] < min)
        {
            min = DistanceDif[i];
            result = z[i];
        }
    }
    return result;
}



