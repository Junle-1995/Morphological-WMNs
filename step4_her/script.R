rm(list=ls())
b<- ("Z:/Shared/zhenzhen/VBN/handle")
savepath     <- "Z:/Shared/zhenzhen/VBN/result/"


functionpath <-"Z:/Shared/zhenzhen/VBN"
functionname <-'gretna_heritability.R'

setwd(functionpath)
source(functionname)
library("R.matlab")

{ 
  Files<-list.files("Z:/Shared/zhenzhen/VBN/handle")
  ACEresults <- list();
  
  {
    list_data   <- readMat(file.path(b,paste(Files, sep = "",collapse ="")))
    Data_var    <- list_data[[1]];
    result_var  <- gretna_heritability(Data_var);
    matfilename <- paste(savepath,Files,sep="")
    writeMat( matfilename,h2=result_var[,1])
  }
}