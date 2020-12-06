#include <iostream>

#include "castor/matrix.hpp"
#include "castor/castor_fftwf.hpp"
#include "castor/castor_fftwm.hpp"

using namespace castor;

int main()
{
    //
    std::size_t N = 5;
    matrix<float> freq = fftw::fftfreq(N,1./N);
    disp(freq);
    disp(fftw::fftshift(freq));

    N = 4;
    freq = fftw::fftfreq(N,1./N);
    disp(freq);
    disp(fftw::fftshift(freq));
    return EXIT_SUCCESS;
}