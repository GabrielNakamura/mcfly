#####constructing McFly package#######


###packages
library(devtools)
library(usethis)

#initializing package
create_package(path = "~/OneDrive/Packages/mcfly")

#####creating package documentation####
document()

####checking for package documentation#####
check()
usethis::use_package("mcfly")
tools::Rdindex(RdFiles = here::here("man")) #gnerate index page


#####picking a licence for package######
use_mit_license("Gabriel Nakamura")

#####include Furnaridae data#####
load("Furnaridae.RData")
Furnariidae<- comm.subset
usethis::use_data(Furnariidae, mcfly, overwrite = T)
phylo_Furnariidae<- sub.phy
usethis::use_data(phylo_Furnariidae, overwrite = T)
envir<- env.full
usethis::use_data(envir, overwrite = T)

####creating a vignette#####
usethis::use_vignette("McFly_vignette")
