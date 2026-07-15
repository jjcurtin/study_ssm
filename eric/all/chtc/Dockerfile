# Starting point
FROM --platform=linux/amd64 rocker/tidyverse:4.3.1

# Install R packages
RUN R -e "install.packages(c('data.table','KFAS','comprehenr','invgamma','truncnorm','fitdistrplus','MASS','MARSS','mvtnorm','tictoc','jsonlite'), repos = 'http://cran.us.r-project.org')"