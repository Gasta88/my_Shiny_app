# install.packages("devtools")
# devtools::install_github("bartongroup/proteusLabelFree")
# devtools::install_github("bartongroup/proteusTMT")
# devtools::install_github("bartongroup/proteusSILAC")
# devtools::install_github("bartongroup/Proteus", username="MarekGierlinski", auth_token = "7f457d5e442ac05d675c8de77ac6c7bea696d32e", build_vignettes = TRUE)
library(proteus)
library(shiny)
library(proteusLabelFree)

data(proteusLabelFree)
prodat.med <- normalizeData(prodat)
res <- limmaDE(prodat.med)

saveRDS(prodat.med,file="data/prodat.med.rds")
saveRDS(res,file="data/res.rds")

