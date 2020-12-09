.. _label-install:

Installation
++++++++++++

Obviously, **Castor-FFTW** requires the `castor project <http://leprojetcastor.gitlab.labos.polytechnique.fr/castor/>`_ to be installed on your computer. Then, installing the FFTW3 library is required.

Install the FFTW3
-----------------

Installing the FFTW3 depends on your operating system. If you have a package manager, it is probably the best way to install it. If you are running on **Ubuntu**, you can simply use the command

.. code:: text

    sudo apt install libfftw3-dev

If you are using **MacOS** and Homebrew, it can be installed using the following command

.. code:: text

    brew install fftw

It should be also available if you are using the Intel MKL library. Otherwise, you can compile it from source (the FFTW3 is standalone), see `this page <http://www.fftw.org/download.html>`_.


Install the headers
-------------------

In order to install **Castor-FFTW**, it is sufficient to copy the three header-files in the ``castor-fftw/include/castor/`` folder to the folder where the header-files of the **castor project** are installed. 


That's all! See the :ref:`label-howto` section for the compile instructions.