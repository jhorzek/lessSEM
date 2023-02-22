# installing all dependencies one by one while also checking if they could be loaded

pkgs <- c('lavaan')

for(pkg in pkgs)
{
  
  print(paste0("Installing package ", pkg))
  
  # check if the package has already been installed
  installed_packages <- installed.packages()
  if(pkg %in% installed_packages)
    next
  
  install.packages(pkg, 
                   dependencies = c("Depends", 
                                    "Imports", 
                                    "LinkingTo"),
                   repos = 'http://cran.rstudio.com/')
  
  # if the package could not be installed: exit
  if(!require(pkg, character.only	= TRUE))
  {
    stop(paste0("Could not install package ", pkg, "."))
    quit(status = 1, save = "no")
  }
  
}
