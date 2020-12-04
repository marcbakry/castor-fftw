
.. _label-howto:

How-to
======

The **Castor-FFTW** wrapper has been divided into two headers. The first one, ``castor_fftwf.hpp`` is for the ``float`` precision while ``castor_fftw.hpp`` is for the ``double`` precision. Both should be included if you plan to mix ``float`` and ``double``. However, the main interface is made such that the function call remains the same whatever the precision. We describe below a small example extracted from the ``demo/`` folder.

We recall that all functions of the **castor** project are defined within the namespace ``castor::``. However, the **castor** project already provides a FFT interface through the KISSFFT library. In order to avoid interferences, the functions of the **Castor-FFTW** wrapper are defined within the sub-namespace ``castor::fftw::``.

Basics
======

A minimal example
+++++++++++++++++

The minimal ``demo_float.cpp`` file should look like this

.. code:: c++

    #include "castor/matrix.hpp"
    #include "castor_fftwf.hpp"

    using namespace castor;

    int main()
    {
        // your code here
    }


We will compute Fourier transforms of a ``float`` matrix, compute the inverse transforms and compare them to the original ``matrix``. Everything can be easily translated to ``double`` precision.

First, let us create some data.

.. code:: c++

    std::size_t   N   = 100;
    matrix<float> A1d = rand<float>(N,1); // could have been rand<float>(1,N)

Then, we compute the forward Discrete Fourier Transform and we normalize it.

.. code:: c++

    auto Ahat = fftw::fft(A)/static_cast<float>(N);

Many things happened internally in the line above. First, the matrix ``A`` was *copied* into a ``matrix<std::complex<float>>``. **It would not have been the case if** ``A`` **had already been a** ``matrix<std::complex<float>>``. It is related to the current internals of the FFTW3 library. Then, the FFT was computed and returned as a ``matrix<std::complex<float>>``. Since the result of ``castor::fftw::fft()`` is not normalized, we divide it by the number of elements. Now, we compute the backward transform and compare the result with the original data. 

.. code:: c++

    disp(abs(norm(fftw::ifft(Ahat) - A)));

The result should be something like this depending on your computer

.. code:: text

    5.5672e-07

A value of the order of ``1e-06`` is totally acceptable.

**Remark:** It is *not* possible to write

.. code:: c++

    Ap = fftw::ifft(fftw::fft(A)/static_cast<float>(N)); // not possible

in our context. This is due to the wrapper which requires the use of a pointer to interact with the ``C`` functions. The issue will not be discussed further.

Now, let us perform a 1D transform along the dimensions of a matrix. The behavior of ``fftw::fft`` is the same as the corresponding Matlab function.

.. code:: c++

    std::size_t M = 50;
    A    = rand<float>(M,N);
    // default is column by column
    Ahat = fftw::fft(A)/static_cast<float>(M);
    disp(abs(norm(fftw::ifft(Ahat) - A)));
    // transform line by line
    Ahat = fftw::fft(A,2)/static_cast<float>(N);
    disp(abs(norm(fftw::ifft(Ahat,2) - A)));

The result should look like

.. code:: text

    4.06223e-06
    4.24932e-06

Finaly, we compute a 2D transform on a matrix. There are two ways to perform such a transform. By calling ``fft2(A)``, ``A`` is assumed to be two-dimensional, even if one dimension is equal to 1. By calling ``fft2(M,N,A)``, ``A`` can have any size as long as ``M*N == size(A,1)*size(A,2)``.

.. code:: c++

    Ahat = fftw::fft2(A)/static_cast<float>(M*N);
    disp(abs(norm(fftw::ifft2(M,N,Ahat) - A)));

.. code:: text

    6.5948e-06


**Note:** In future developments, the support for the three-dimensional FFT will be added through the ``fftw::fft3(M,N,K,A)`` interface where ``(M,N,K)`` are the dimensions.


Compilation
+++++++++++

Assuming that the **castor** project and **Castor-FFTW** have been installed in a standard location (meaning that the headers can be found automatically by the compiler), assuming that the compiler is ``g++``, the program above can be compiled easily with the following command line (**Ubuntu** and **MacOS**)

.. code:: text

    g++ demo_float.cpp -o test_float -lfftw3f

For the ``double`` version,

.. code:: text

    g++ demo_double.cpp -o test_double -lfftw3

Obviously, if ``float`` and ``double`` are mixed together, one can combine both 

.. code:: text

    g++ demo_double_float.cpp -o test_double -lfftw3 -lfftw3f

If, for one reason or the other, some headers cannot be found, it is possible to indicate their path to the compiler like 

.. code:: text

    g++ -I/path/to/missing/headers main.cpp -o myExecutable -l...

**Warning:** All the header files are assumed to be within a ``castor/`` subfolder. Consequently, the command line should be 

.. code:: text

    g++ -I/path/to/castor/folder -I/path/to/other/missing/headers main.cpp -o ...


Practical examples
==================

We give now two, somehow, practical examples in order to demonstrate the functionalities of **Castor-FFTW**. In the first example, we simply compute and plot the amplitude spectrum of a sum of sine functions. In the second example, we compute the transform of a rectangular function which we smooth using a gaussian filter.

In both examples, the plots will be made using the graphical functionalities of the **castor** project. We refer to the `corresponding documentation <http://leprojetcastor.gitlab.labos.polytechnique.fr/castor/graphics.html>`_.


Fourier transform of sine functions
+++++++++++++++++++++++++++++++++++

plop



Regularization using a gaussian filter
++++++++++++++++++++++++++++++++++++++

plop