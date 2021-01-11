

## install R packages used in the CUT&RUNTools2.0
if(!require(BiocManager)){
    install.packages("BiocManager")
}
if(packageVersion("BiocManager") < "1.30.10"){
    install.packages("BiocManager")
}

message(paste("#### Install R packages #### "))

# CRAN packages
cran.pkgs <- c("devtools", "reticulate", "leiden", "data.table", "Matirx", "irlba", "ggplot2", "igraph",
  "uwot", "CENTIPEDE")
for(i in cran.pkgs){
    if(!require(i, character.only = T)) {
        message(paste("Install ", i, "\n"))
        install.packages(i, dependencies = TRUE, repos="http://R-Forge.R-project.org")
    }
}

# individual packages
indi.pkgs <- c("Rtsne", "RANN")
message(paste("Install Rtsne \n"))
devtools::install_github("jkrijthe/Rtsne")
message(paste("Install RANN \n"))
devtools::install_github("jefferis/RANN")

# biocManage packages
bioc.pkgs = c("rGREAT")
for(i in bioc.pkgs){
    if(!require(i, character.only = T)) {  
        message(paste("Install", i, "..."))
        BiocManager::install(i)
    }
}


# check if all the required packages were installed
# 
pkg <- c(cran.pkgs, indi.pkgs, bioc.pkgs)
if(any(! pkg %in% rownames(installed.packages()))){
    message(paste("#### All the R packages were installed successfully ! #### \n"))
} else {
    message(paste("#### Some packages were not installed successfully ! #### \n"))
    message(paste("#### Please install them manually #### \n"))
    pkg[! pkg %in% rownames(installed.packages())]
}


