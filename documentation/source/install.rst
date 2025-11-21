.. _chap:install:

Installation
====================================

This chapter describes the installation of GRAMPC-S for Windows and Linux. 
GRAMPC-S depends on the Eigen library in version 3.4.90 or higher (tested in version 3.4.90) :footcite:p:`eigenweb`. 
If Eigen is not already installed on your system, clone the repository, create the folder **build**, and install Eigen. 
For **Linux**, execute the following commands:

::

    git clone https://gitlab.com/libeigen/eigen.git
    mkdir eigen/build
    cd eigen/build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    sudo make install

For Windows with MinGW compiler, execute the following commands:

::

    git clone https://gitlab.com/libeigen/eigen.git
    mkdir eigen\ build
    cd eigen\ build
    cmake -DCMAKE\BUILD\TYPE=Release .. -G "MinGW Makefiles"
    cmake --install .

After Eigen has been installed, GRAMPC-S can be built.
The source code of GRAMPC-S is available at https://github.com/grampc/grampc-s.
Clone the code including all submodules:

::

    git clone --recurse-submodules https://github.com/grampc/grampc-s.git

After cloning the repository, a new directory with the following subfolders is created:

- **documentation**: The PDF-file of the GRAMPC-S documentation can be found here.
- **examples**: Contains MPC-examples and a template for implementing your own code.
- **include**: The header files of GRAMPC-S are located in this folder.
- **libs**: Includes external libraries that are linked as git submodules.
- **matlab**: This folder contains Matlab files.
- **src**: The source code of GRAMPC-S can be found here.

As soon as the code has been downloaded, it can be built.
To build GRAMPC-S using Linux, run the CMakeLists.txt and call make afterwards:

::

    cmake -DCMAKE\BUILD\TYPE=Release CMakeLists.txt
    make

To build GRAMPC-S using Windows, the editor Visual Studio Code with the C++ extension is recommended. Open the folder **grampc-s** in Visual Studio Code, click the button **Build**, and select a compiler from the list.
After building the code, the folders **bin** and **build** are created in addition to the folders above.
The folder **bin** contains the executable files for the examples.

.. footbibliography::
