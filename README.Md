# kinference will ultimately move to a different R-universe (not closekin), once the licence is sorted out. For now, install like this:

install.packages( 'kinference', repos = unique( c( closekin= 'https://closekin.r-universe.dev', markbravington= 'https://markbravington.r-universe.dev', getOption( 'repos')[ 'CRAN'], 'https://cloud.r-project.org')))
