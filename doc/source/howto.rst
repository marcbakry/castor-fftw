
.. _label-howto:

How-to
======

The **Castor-FFTW** wrapper has been divided into two headers. The first one, ``castor_fftwf.hpp`` is for the ``float`` precision while ``castor_fftw.hpp`` is for the ``double`` precision. Both should be included if you plan to mix ``float`` and ``double``. However, the main interface is made such that the function call remains the same whatever the precision. We describe below a small example extracted from the ``demo/`` folder.

The minimal ``main.cpp`` file should look like this

.. code:: c++

    #include "castor/matrix.hpp"
    #include "castor_fftwf.hpp"

    int main()
    {
        // your code here
    }


We will compute the one-dimensional Fourier transform of a ``float`` matrix, compute the inverse transform and compare it to the original ``matrix``. Everything can be easily translated to ``double`` precision.