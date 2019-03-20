## Package 'kinference': installation instructions

Package 'kinference' is presented as source code, and includes c++ (via Rcpp) that will require compilation before during install. The package requires the c++11 compilation standard, and makes use of several R packages from the mvb ecosystem. 'kinference' install has been tested under Windows (Windows 10, 64-bit) and Linux (Ubuntu 16.04, 64-bit). Installation has not been tested for Mac. Before installing 'kinference', make sure you have installed:

* In R: 
  * from CRAN:
    * Rcpp
    * RcppArmadillo
    * abind
    * BH
    * devtools
    * roxygen2
    * tensorA
  * from the MVB ecosystem:
    * mvbutils [NOT the out-of-date version from CRAN]
    * atease
    * vecless
    * gbasics
    * genocalldart
    * nicetime

If installing from source, Windows users will need to download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Installation from source: GNU/Linux

Download or clone this repository, then run the following in R: 

    library(devtools)
    setwd( "path/to/kinference" ) # e.g., "~/R/local/kinference/"
                                  # make sure you have write permission for 'kinference'.
    install()
    
package will compile and install.

## Installation from source: Windows

Download or clone this repository, then run the following in R: 

    library(devtools)
    setwd( "path/to/kinference" ) # e.g., "C:/Users/myID/R/local/kinference/"
                                  # make sure you have write permission for 'kinference'.
    install()

package will compile and install.

## Installation from local repositories: GNU/Linux or Windows (binary options for Windows)
 
The source code in this github repository is the development version of 'kinference'. A stable release
exists, and is managed by MVB. Installation is simple, but you will have to ask MVB for access to his
semi-secret repo.

There is also a development release hosted at CSIRO for internal use and testing, maintained by SMB. The CSIRO copy includes precompiled binaries for Windows users running R ver 3.3 - 3.5, and should be kept relatively up-to-date with the source version presented here. We intend to merge these two repositories in the short term.

To install from MVB's stable release, run:

    install.packages(repos='https://ask/mark/about/his/repo/')

You should see a pop-up window with a list of packages. Double-click 'kinference'.

To install from the CSIRO copy development release, run:

    install.packages(pkgs = "kinference", repos = 'file:///ask/shane/about/this/repo/')

package will install, without requiring compilation for Windows users.

## Checking the installation

If you have correctly installed the development version of 'kinference', 

    ?crapometer ## should bring up a help menu.

It will not return a help menu for the stable release.

## Notes and known bugs

Bug 1: Errors like: "files in does not match files out".

 - R does not have write-permission in the root of 'kinference'.

Bug 2: Errors like: "file type not recognized: some_cpp.o" or "... RccpExports.o".

 - This bug is only known to occur under Windows. There is some compiled code where there shouldn't be, and your compiler is not happy about it. If you encounter this bug, just delete the following two files from kinference/src/, and run install() again:

  -- RcppExports.o
  -- some_cpp.o

Bug 3: Anything implying that c++11 isn't available.

 - 'kinference' inclues compiled c++ code, and that code includes some features that are only available in c++11 and later specifications. Some compilers cannot handle c++11, and it is possible that your compiler is one of those. You can either make sure you are using a compiler capable of c++11, or (if you are a Windows user) install a pre-compiled binary from Mark's semi-secret repo.

