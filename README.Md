## Package 'kinference': installation instructions

This file that you're reading pertains mainly to the Github version of kinference. However, you may also have dot-tar-dot-gz "tarballs" of 'gbasics' and 'kinference', which are good too. Users "just wanting the damn package" could read the "Alternative" bit at the end.

Installing the Github version, or the dot-tar-dot-gz versions,  _from source_ is something that GNU/Linux OSes can do out-of-the-box, whereas Windows users need to add Rtools and Mac users need xcode. Those should be pretty standard these days.

### Rtools installation on Windows:
See [https://cran.r-project.org/bin/windows/Rtools/]

### xcode installation on Mac:
See [https://clanfear.github.io/CSSS508/docs/compiling.html]

Then you need to make sure you have some CRAN packages, and some MVB packages, installed. Like so:

### The CRAN packages: Rcpp and remotes

`install.packages("Rcpp")`
`install.packages("remotes")`

### The MVB support packages

`install.packages(c("atease", "mvbutils", "vecless"), repos = "https://markbravington.r-universe.dev/")`

### kinference and gbasics

If you have dot-tar-dot-gz files, then 

`install.packages( "<path-to-gbasics-targz", repos=NULL)`
`install.packages( "<path-to-kinference-targz", repos=NULL)`

and you should be good to go.

If you are working off Github, ATM the Github repos for these two packages are private (though that's temporary; see below), so either we've specially given you access or more likely we've sent you zip files (which we ourselves obtained from Github, via Code->Download_ZIP). 

Then, if we've given you access to Github itself, you can do:

`remotes::install_github("closekin/gbasics")`
`remotes::install_github("closekin/kinference", build = TRUE, build_vignettes = TRUE)`


If you only have the ZIPs, you need 

`remotes::install_local( "<path-to-gbasics-zip", subdir='gbasics')`
`remotes::install_local( "<path-to-kinference-zip", subdir='kinference', build = TRUE, build_vignettes = TRUE)`

If you don't need the vignette, feel free to skip the `build` and `build_vignettes` arguments. Otherwise, do check it out - the kinference vignette is:

`vignette( "kinference-vignette")`

## Alternative installation for gbasics and kinference

This is the way most people might do it, because it corresponds closely to installation of most packages from CRAN and is very easy. It should work on Linux and Macs, as well as on Windows where it does not even require you to have Rtools. It has worked well in the past, but it relies on MVB having kept his own CRAN-alike repo up-to-date. (Note that Github's meaning of "repo" is different from CRAN's meaning. Not our fault.)

`install.packages(pkgs = c("gbasics", "kinference"), repos = "https://markbravington.github.io/Rmvb-repo")`

...except that we'll have to get back to you on the location of that repos, because Rmvb-repo is no longer with us.


## Plans

These packages, plus the MVB support packages, will probably move soon to R-universe, which will become _the_ place to get them from (except for anyone needing bleeding-edge versions and actively developing gbasics or kinference, which is just SMB and MVB for now). R-universe is like CRAN from the user's PoV, but much less annoying from the developer's PoV. 

