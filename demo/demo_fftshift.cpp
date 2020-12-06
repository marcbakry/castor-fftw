#include <iostream>

#include "castor/matrix.hpp"
#include "castor/castor_fftwf.hpp"
#include "castor/castor_fftwm.hpp"

using namespace castor;

int main()
{
    std::cout << "+------------------------------------------+" << std::endl;
    std::cout << "| TESTING CASTOR-FFTW: FFTFREQ/(I)FFTSHIFT |" << std::endl;
    std::cout << "+------------------------------------------+" << std::endl;
    //
    std::size_t N = 7;
    matrix<float> freq = fftw::fftfreq(N,1./N);
    disp(freq,1);
    disp(fftw::fftshift(freq),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq)) - freq,"inf"),1);

    std::cout << std::endl;
    N = 6;
    freq = fftw::fftfreq(N,1./N);
    disp(freq,1);
    disp(fftw::fftshift(freq),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq)) - freq,"inf"),1);
    return EXIT_SUCCESS;
}