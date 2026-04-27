local({ owd <- setwd( "D:/r2.0/kinference/vignettes");
try({
knitr::knit( "kinference.Rmd.orig", "kinference.Rmd");
});
setwd( owd);
})
