## Package 'kinference': installation instructions

GNU/Linux OSes can do this out-of-the-box, but Windows users need to add
Rtools and Mac users need xcode.

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

Download `gbasics` and `kinference` as .zip files, and install using:

`remotes::install_local("/filepath/to/kinference-master.zip", subdir = "kinference")`
`remotes::install_local("/filepath/to/gbasics-master.zip", subdir = "gbasics")`

Note that you'll have to provide a `filepath/to/the/.zip/file` that is
appropriate to your OS.

## Alternative installation for gbasics and kinference
`install.packages(pkgs = c("gbasics", "kinference"), repos = "https://markbravington.github.io/Rmvb-repo")`
