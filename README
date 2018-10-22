## Package 'kinference': installation instructions

Package 'kinference' is presented as source code, and includes c++ (via Rcpp) that will require compilation before during install. The package makes use of the 'boost' c++ libraries and several R packages from the mvb ecosystem. 'kinference' install has been tested under Windows (Windows 10, 64-bit) and Linux (Ubuntu 16.04, 64-bit). Installation has not been tested for Mac. Before installing 'kinference', make sure you have installed:

* In R: 
  * from CRAN:
    * Oarray
    * mixtools
    * xtable
    * tools
    * devtools
    * Rcpp
    * roxygen2
  * from the MVB ecosystem:
    * mvbutils [NOT the out-of-date version from CRAN]
    * atease
    * Oarray.support
    * vecless
    * gbasics
    * genocalldart
* In your OS:
  * [Boost c++ libraries](https://www.boost.org/users/download/)

**GNU/Linux users**: you can install Boost using `sudo apt-get install libboost-dev-all`. If you install via apt-get, you *probably* will not have to update any 'path' variables.

**Windows users**: note the path to the root of the Boost - it will be a folder with a name similar to 'boost_1_68_0'. You *will* need to update some 'path' variables.

Additionally, Windows users will need to download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Installation process: GNU/Linux

In R, run the following:

    setwd( [path/to/kinference/] ) # make sure you have write permission for 'kinference'.
    library(devtools)
    install(kinference)
    
package will compile and install
    
    library(kinference)

you quickly check installation by:

    ?crapometer  ## should bring up an HTML help menu.

## Installation process: Windows

In R, run the following:

    Sys.setenv("PKG_CXXFLAGS"="-I/*path/to/boost/root/*") # fill in your boost path
    setwd( [path/to/kinference] ) # make sure you have write permission for 'kinference'.
    library(devtools)
    install(kinference)

package will compile and install

    library(kinference)

you can check installation by:

    ?crapometer  ## should bring up an HTML help menu.


