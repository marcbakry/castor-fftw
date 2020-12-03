#include <iostream>

#include "castor/matrix.hpp"
#include "castor/castor_fftw.hpp"

using namespace castor;

int main()
{
    std::cout << "+-----------------------------+" << std::endl;
    std::cout << "| TESTING CASTOR-FFTW: DOUBLE |" << std::endl;
    std::cout << "+-----------------------------+" << std::endl;

    //-------------------------------//
    // FORWARD/BACKWARD 1d-TRANSFORM //
    //-------------------------------//
    std::size_t N = 100;
    double tmin    = 0.;
    double tmax    = 1.;
    double dt      = (tmax - tmin)/N;

    // test the vector version
    matrix<double> t = zeros<double>(N,1);
    for(auto i=0; i<N; ++i) t(i) = i*dt;    // some time range

    matrix<double> f = sin(2*M_PI*t) + 7.5f*sin(3*2*M_PI*t); // some simple signal
    auto fhat       = fftw::fft(f)/static_cast<double>(N);   // do not forget normalization !
    disp(abs(norm(fftw::ifft(fhat) - f)));

    // test de matrix version
    std::size_t M = 50;
    matrix<double> A = rand<double>(M,N);
    auto Ahat_1 = fftw::fft(A)/static_cast<double>(M); // Fourier transform column by column
    disp(abs(norm(fftw::ifft(Ahat_1) - A)));
    auto Ahat_2 = fftw::fft(A,2)/static_cast<double>(N); // Fourier transform line by line
    disp(abs(norm(fftw::ifft(Ahat_2,2) - A)));

    //-------------------------------//
    // FORWARD/BACKWARD 2d-TRANSFORM //
    //-------------------------------//
    auto Ahat = fftw::fft2(A)/static_cast<double>(M*N);
    disp(abs(norm(fftw::ifft2(M,N,Ahat) - A)));

    // THE END
    return EXIT_SUCCESS;
}