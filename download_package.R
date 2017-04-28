# Set working directory to where zip file is downloaded 
unzip("gmwm2-master.zip")
file.rename("gmwm2-master", "gmwm2")
shell("R CMD build gmwm2")
# Install local gmwm2 package
install.packages("gmwm2_1.0.tar.gz", repos = NULL)
# Call gmwm2 Package
library(gmwm2)