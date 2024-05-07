## Package 'kinference': installation instructions

This file that you're reading pertains only to the Github version of kinference. Windows users "just wanting the damn package" could read the "Alternative" bit at the end.

Installing the Github version  _from source_ is something that GNU/Linux OSes can do out-of-the-box, whereas Windows users need to add Rtools and Mac users need xcode. Those should be pretty standard these days.

### Rtools installation on Windows:
See [https://cran.r-project.org/bin/windows/Rtools/]

### xcode installation on Mac:
See [https://clanfear.github.io/CSSS508/docs/compiling.html]

Then you need to make sure you have some CRAN packages, and some MVB packages, installed. Like so:

### The CRAN packages: Rcpp and remotes

`install.packages("Rcpp")`
`install.packages("remotes")`

### The MVB support packages

`install.packages(pkgs = c("atease", "mvbutils", "vecless"), repos = "https://markbravington.github.io/Rmvb-repo")`

### kinference and gbasics

ATM the Github repos for these two packages are private (though that's temporary; see below), so either we've specially given you access or more likely we've sent you zip files (which we ourselves obtained from Github, via Code->Download_ZIP). 

Then, if we've given you access to Github itself, you can do:

`remotes::install_github("closekin/gbasics")`
`remotes::install_github("closekin/kinference")`

If you only have the ZIPs, you need 

`remotes::install_local( "<path-to-gbasics-zip", subdir='gbasics')`
`remotes::install_local( "<path-to-kinference-zip", subdir='kinference')`


## Alternative installation for gbasics and kinference

This is the way most people might do it, because it corresponds closely to installation of most packages from CRAN and is very easy (not requiring compilation, if you're on Windows). It has worked well in the past, but it relies on MVB having kept his own CRAN-alike repo up-to-date, which is not always the case. (Note that Github's meaning of "repo" is different from CRAN's meaning. Not our fault.)

`install.packages(pkgs = c("gbasics", "kinference"), repos = "https://markbravington.github.io/Rmvb-repo")`


## Plans

These packages, plus the MVB support packages, will probably move soon to R-universe, which will become _the_ place to get them from (except for anyone needing bleeding-edge versions and actively developing gbasics or kinference, which is just SMB and MVB for now). R-universe is like CRAN from the user's PoV, but much less annoying from the developer's PoV. 

