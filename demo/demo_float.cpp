#include <iostream>

#include "castor/matrix.hpp"
#include "castor/castor_fftw.hpp"

using namespace castor;

int main()
{
    std::cout << "+----------------------------+" << std::endl;
    std::cout << "| TESTING CASTOR-FFTW: FLOAT |" << std::endl;
    std::cout << "+----------------------------+" << std::endl;

    //-------------------------------//
    // FORWARD/BACKWARD 1d-TRANSFORM //
    //-------------------------------//
    float N    = 100;
    float tmin = 0.;
    float tmax = 1.;
    float dt   = (tmax - tmin)/N;

    // test the vector version
    matrix<float> t = zeros<float>(N,1);
    for(auto i=0; i<N; ++i) t(i) = i*dt;    // some time range

    matrix<float> f = sin(2*(float)M_PI*t) + 7.5f*sin(3*2*(float)M_PI*t); // some simple signal
    auto fhat       = fftw::fft(f)/N;   // do not forget normalization !
    disp(abs(norm(fftw::ifft(fhat) - f)));

    // THE END
    return EXIT_SUCCESS;
}