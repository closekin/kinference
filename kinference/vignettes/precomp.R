local({ owd <- setwd( "D:/r2.0/kinference/vignettes");
try({
knitr::knit( "kinference-vignette.Rmd.orig", "kinference-vignette.Rmd");
});
setwd( owd);
})
