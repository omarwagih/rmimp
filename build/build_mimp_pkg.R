setwd('~/Development/mimp/')

VERSION = '1.0'

detachPackage <- function(pkg){
  pkg = sprintf("package:%s", pkg)
  res = tryCatch({
    detach(pkg, unload=TRUE, character.only=T, force=T)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(res)
}

build_package <- function(){
  require(devtools)
  targz = sprintf('MIMP_%s.tar.gz', VERSION)
  # Move up one directory
  newf = file.path('./build', targz)
  # Compile things
  document('./')
  file.remove(newf)
  system('R CMD BUILD ./')
  re = file.rename(targz, newf)
  # Install package
  system(sprintf('R CMD INSTALL %s', newf))
  # Reload in current environment
  detachPackage('MIMP')
  require(MIMP)
}


# Build package
build_package()