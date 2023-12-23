## Package 'kinference': installation instructions

This file that you're reading pertains only to the Github version of kinference. Windows users "just wanting the damn package" could read the "Alternative" bit at the end. OTOH, 

Installing the Github version  _from source_ is something that GNU/Linux OSes can do out-of-the-box, but Windows users need to add Rtools and Mac users need xcode. Those should be pretty standard these days.

### Rtools installation on Windows:
See [https://cran.r-project.org/bin/windows/Rtools/]

### xcode installation on Mac:
See [https://clanfear.github.io/CSSS508/docs/compiling.html]

### The CRAN packages: Rcpp and remotes

You will need to have `Rcpp` and `remotes` installed:

`install.packages("Rcpp")`
`install.packages("remotes")`

### The MVB support packages
`install.packages(pkgs = c("atease", "mvbutils", "vecless"), repos = "https://markbravington.github.io/Rmvb-repo")`

### kinference and gbasics

To install _from source_ the most-up-to-date (but possibly with new bugs--- as opposed to old bugs) version, it's something like this:

`remotes::install_github("closekin/gbasics")`
`remotes::install_github("closekin/kinference")`

but that doesn't seem to work, because the "path" to the github source isn't right, somehow. It might be because the github repo is private ATM. 


## Alternative installation for gbasics and kinference

This is the way most people might do it, because it corresponds closely to installation of most packages from CRAN and is very easy (not requiring compilation, if you're on Windows). But it relies on MVB having kept his own CRAN-alike repo up-to-date. (Note that Github's meaning of "repo" is different from CRAN's meaning. Not our fault.)

`install.packages(pkgs = c("gbasics", "kinference"), repos = "https://markbravington.github.io/Rmvb-repo")`
