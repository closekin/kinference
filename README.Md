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

If installing from source, Windows users will need to download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Installation from source: GNU/Linux

In R, run the following:

    library(devtools)
    setwd( "path/to/kinference" ) # e.g., "~/R/local/kinference/"
                                  # make sure you have write permission for 'kinference'.
    install()
    
package will compile and install
    
    library(kinference)

you can check installation by:

    ?crapometer  ## should bring up an HTML help menu.

## Installation from mark's_semi-secret_repo: GNU/Linux or Windows

In R, run:

    install.packages(repos='https://markbravington.github.io/mark's_semi-secret_repo')

You should see a pop-up window with a list of packages. Double-click 'kinference'.

    library(kinference)

you can check installation by:

    ?crapometer  ## should bring up an HTML help menu.

## Installation from source: Windows

In R, run the following:

    library(devtools)
    setwd( "path/to/kinference" ) # e.g., "C:/Users/myID/R/local/kinference/"
                                  # make sure you have write permission for 'kinference'.
    install()

package will compile and install

    library(kinference)

you can check installation by:

    ?crapometer  ## should bring up an HTML help menu.

## Notes, in case you have bugs.

Bug 1: Errors like: "files in does not match files out".

 - R does not have write-permission in the root of 'kinference'.

Bug 2: Errors like: "file type not recognized: some_cpp.o" or "... RccpExports.o".

 - This bug is only known to occur under Windows. There is some compiled code where there shouldn't be, and your compiler is not happy about it. If you encounter this bug, just delete the following two files from kinference/src/, and run install() again:

  -- RcppExports.o
  -- some_cpp.o

Bug 3: Anything implying that c++11 isn't available.

 - 'kinference' inclues compiled c++ code, and that code includes some features that are only available in c++11 and later specifications. Some compilers cannot handle c++11, and it is possible that your compiler is one of those. You can either make sure you are using a compiler capable of c++11, or (if you are a Windows user) install a pre-compiled binary from Mark's semi-secret repo.

