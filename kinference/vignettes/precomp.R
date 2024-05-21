local({ owd <- setwd( "D:/r2.0/kinference/vignettes");
try({
knitr::knit( "kinference-vignette.Rmd.orig", "kinference-vignette.Rmd");
});
setwd( owd);
})
## SMB's notes to self:
## If you want to precompile, setwd() to /kinference/kinference/vignettes/,
## then knitr::knit( "kinference-vignette.Rmd.orig", "kinference-vignette.Rmd").
## If you're not in /vignettes/, then figures-kinference/ will be generated
## in your WD, and R CMD check won't build the vignette properly
