gretna_heritability<-function(Data){
  # ==================================================================================
  # This function is used to calculate the narrow sense of heritability using
  # ACE model based on OpenMx package (https://github.com/OpenMx/OpenMx).
  #
  # Syntax:  gretna_heritability <-function(Data)
  # 
  # Inputs:
  #          Data:
  #                A 2D (N*3) or 3D (N*3*M) data array, with N denoting the
  #                number of twin pairs, and M denoting the number of variables.
  #                For both cases, the first and second columns present observed
  #                data (each row for each pair), and the third column encodes
  #                the twin type (1 for monozygotic, and 2 for dizygotic). 
  #                
  # Outputs:
  #          a^2:
  #                Heritability.
  #          c^2:
  #
  #          e^2:
  # 
  # ZhenZhen Luo, IBRR, Guangzhou, 28/12/2020, luozhenzhen0419@gmail.com
  # Zhen Li,      IBRR, Guangzhou, 28/12/2020, Zhen.li.dada@gmail.com
  # Jinhui Wang,  IBRR, Guangzhou, 28/12/2020, jinhui.wang.1982@gmail.com
  # ================================================================================
  
  require(OpenMx)
  
  result_var_all  =  matrix(0,1,3)
  
  if (is.na(dim(Data)[3]))
  { num_var <-1
  }else{
    num_var <- dim(Data)[3]}
  
  for(ivar in 1:num_var)
  { if(is.na(dim(Data)[3]))
  { Data_var<-Data
  }else{
    Data_var= Data[,,ivar]}
    
    colnames(Data_var) <- c("twin1","twin2","zyg")
    
    selVars <- c('twin1','twin2')
    aceVars <- c("A1","C1","E1","A2","C2","E2")
    
    # Select data for analysis
    mzData <- subset(Data_var, Data_var[,3]==1, selVars)
    dzData <- subset(Data_var, Data_var[,3]==2, selVars)
    
    # Generate descriptive statistics
    Mean_mz  <- colMeans(mzData,na.rm=TRUE)
    Mean_dz  <- colMeans(dzData,na.rm=TRUE)
    Mean_two <- c(Mean_mz,Mean_dz)
    Col_Mean <- mean(Mean_two)
    
    # Path objects for multiple groups
    manifestVars =  selVars
    latentVars   =  aceVars
    
    # variances of latent variables
    latVariances <- mxPath( from=aceVars, arrows=2,free=FALSE, values=1 )
    # means of latent variables
    latMeans     <- mxPath( from="one", to=aceVars, arrows=1,free=FALSE, values=0 )
    # means of observed variables
    obsMeans     <- mxPath( from="one", to=selVars, arrows=1,free=TRUE, values=Col_Mean, labels="mean" )
    
    # path coefficients for twin 1
    pathAceT1  <- mxPath( from=c("A1","C1","E1"), to="twin1", arrows=1,free=TRUE, values=.5, label=c("a","c","e") )
    
    # path coefficients for twin 2
    pathAceT2  <- mxPath( from=c("A2","C2","E2"), to="twin2", arrows=1,free=TRUE, values=.5, label=c("a","c","e") )
    
    # Covariance between C1 & C2
    covC1C2    <- mxPath( from="C1", to="C2", arrows=2,free=FALSE, values=1 )
    
    # Covariance between A1 & A2 in MZ twins
    covA1A2_MZ <- mxPath( from="A1", to="A2", arrows=2, free=FALSE, values=1 )
    
    # Covariance between A1 & A2 in DZ twins
    covA1A2_DZ <- mxPath( from="A1", to="A2", arrows=2, free=FALSE, values=.5 )
    
    dataMZ   <- mxData( observed=mzData, type="raw" )
    dataDZ   <- mxData( observed=dzData, type="raw" )
    
    # Combine groups
    paths    <- list( latVariances, latMeans, obsMeans,pathAceT1, pathAceT2, covC1C2 )
    modelMZ  <- mxModel(model="MZ", type="RAM", manifestVars=selVars,latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
    modelDZ  <- mxModel(model="DZ", type="RAM", manifestVars=selVars,latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
    obj      <- mxFitFunctionMultigroup(c("MZ", "DZ"))
    modelACE <- mxModel(model="ACE", modelMZ, modelDZ, obj )
    
    # Run Model
    fitACE   <- mxTryHard(modelACE)
    
    sumACE <- summary(fitACE)
    A <- mxEval(a*a, fitACE)
    
    # Shared environmental variance, c^2
    C <- mxEval(c*c, fitACE)
    
    # Unique environmental variance, e^2
    E <- mxEval(e*e, fitACE)
    
    # Total variance
    V <- (A+C+E)
    
    # Standardized A
    a2 <- A/V
    
    # Standardized C
    c2 <- C/V
    
    # Standardized E
    e2 <- E/V
    
    estACE <- cbind(a2,c2,e2)
    result_var_all <- rbind(result_var_all,estACE)}
  
  result_var_all        <- result_var_all[-1,]
  names(result_var_all) <- c("a^2","c^2","e^2")
  
  estACE = result_var_all
  
  return(result_var_all)}
