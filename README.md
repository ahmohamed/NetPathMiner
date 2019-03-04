# NetPathMiner
[![Travis-CI Build Status](https://travis-ci.org/ahmohamed/NetPathMiner.svg?branch=master)](https://travis-ci.org/ahmohamed/NetPathMiner)

NetPathMiner implements a flexible module-based process flow for network path mining and visualization,
which can be fully integrated with user-customized functions.
It supports construction of various types of genome scale networks from three different pathway
file formats, enabling its utility to most common pathway databases.
In addition, NetPathMiner provides different visualization techniques to facilitate the analysis of even
thousands of output paths.

To report bugs and arising issues, please visit https://github.com/ahmohamed/NetPathMiner

## Installation Instructions

### System Prerequisites
NetPathMiner depends on libxml2 and libsbml to process pathway files. Installation or running
certain functions MAY fail if these prerequisite libraries are
not available. Please read through the following instructions.

#### Prerequisites for Unix users (Linux and Mac OS)
##### Installing libxml2
Make sure your system has library libxml2 installed. In Mac OSX 10.6 or later, libxml2
are built in. For Linux users also, this is almost always the case, however, developing headers
may be missing. To install libxml2 and the headers:

```Shell
    sudo apt-get install libxml2   
    sudo apt-get install libxml2-dev   
```

##### Installing libSBML
Installing libSBML for Unix users is optional. However, NetPathMiner will not be able to process SBML
files. If you will not use SBML functions, you can skip this part.

From the website of libSBML http://sbml.org/Software/libSBML, you can directly download the
binaries suitable for your system from `Download libSBML` link. You can follow the installation instructions
on the website.

#### Prerequisites for Windows users
If you are installing the package through Bioconductor, you don't have to install external libraries. However, currently the Bioconductor version for Windows doesn't support SBML processing. Alternatively, we have prepared all dependencies in a tar file, downloadable from https://github.com/ahmohamed/NPM_dependencies . Please download the file and place in in the home directory of R (type <code>R RHOME</code> in command prompt to locate it), before installation.

Unless you want to use customized libraries, you can skip the rest of this section.To use customized libraries, you have to compile them and provide them to R at the time of installation. This is not a trivial task, please be sure you really need these custom libraries.

##### Installing libxml2
NetPathMiner expects an environment variable `LIB_XML` or `LIB_XML2` pointing to directory where
libxml2 is installed. This directory should have both the compiled library (DLL file) and the header files.

You can download libxml2.dll from http://sourceforge.net/projects/gnuwin32/files/libxml/ among other sources.
Please, place it in a `bin` folder under the installation directory.

You will need also the header files, which can be obtained from NPM_dependecies.tar file. After extracing it, copy
the include directory to the installation directory.

Finally, set the `LIB_XML2` variable to point to the installation directory, which should now contain dll files under `bin`
and header files under `include`.

##### Installing libSBML
Since libSBML is a C++ libraries, it needs to be compiled using GCC compiler. Unforturantely, there is no binary
version for Windows comipled with GCC. To use libSBML, you need to build it from source.

First, dowload source package from http://sourceforge.net/projects/sbml/files/libsbml/ , extract it. You will
need also MinGW http://www.mingw.org/ or the 64-bit version http://mingw-w64.sourceforge.net/ depending on your system.
Add `mingw/bin` to your PATH, by editing environment variables.

Second, you need CMake http://www.cmake.org/ . You can follow the instructions at http://sbml.org/Software/libSBML/docs/java-api/libsbml-installation.html#windows-configuring , however, choose "MinGW Makefiles" instead of "Visual Studio 10".

After finishing the CMake step, use the MinGW's `make.exe` to compile libSBML. Copy the dependencies you used
during the compilation to the `bin` directory. Set the environment variable `LIB_SBML` to point the installation
directory, which should now contain dll files under `bin` and header files under `include`


### R Package dependencies
NetPathMiner depends on package igraph to represent network objects. Installing igraph is required for the package
to work. You will also need devtools package to install directly from github.
NetPathMiner suggests package rBiopaxParser to process BioPAX files and RCurl to download annotations from the web. NetPathMiner can still work without installing the suggested packages, but you will not be able to use the aforementioned functionalities.

##### igraph
Package igraph is available at CRAN. To install it call:
```r
    install.packages("igraph")
```

##### devtools
Package devtools is available at CRAN. For Windows, this seems to depend on
having Rtools for Windows installed. You can download and install this from:
http://cran.r-project.org/bin/windows/Rtools/

To install R package devtools call:
```r
    install.packages("devtools")
```

##### RCurl
For Unix users, make sure your Linux has library libcurl installed. Check out:

```Shell
    locate libcurl   
    locate curl-config
```

If these are not found (usually the developer version is missing), most Linux
users will be able to fix this by running:
```Shell
    sudo apt-get install libcurl4-openssl-dev
```

You will now be able to install R package RCurl. In R console:
```r
    install.packages("RCurl")
```

If you encounter other problems check out http://www.omegahat.org/RCurl/FAQ.html

##### rBiopaxParser
Package rBiopaxParser is available on Bioconductor. For installation instructions check
out http://www.bioconductor.org/packages/release/bioc/html/rBiopaxParser.html or
call:

```r
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("rBiopaxParser")  
```

to install it right away.

### NetPathMiner Installation
If everything went well you will be able to install the NetPathMiner package.

#### From Bioconductor:
In R console, type:

```r
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("NetPathMiner")  
```


#### From GitHub using devtools:
In R console, type:

```r
    library(devtools)   
    install_github(repo="NetPathMiner", username="ahmohamed")
```
