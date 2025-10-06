library(stringr)
library(foreach)
library(doParallel)
library(doSNOW)
require(omics)
library(iterators)
library(lumi) # BiocManager::install("lumi")

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


# Create data with Phenotypes - Technical Variables and Omics Data
getfullOmicsPhenotype <- function( omicFile, metadataFile, phenotype, phenotypeID, phenotypeVariables, omicsLayer, population   ) 
{
    
    # Load omic Data
    if( tolower(tools::file_ext(omicFile)) == 'rds') {
        omicData <- readRDS(omicFile)
    } else if( tolower(tools::file_ext(omicFile)) == "rdata") {
        omicData <- get(load(omicFile))
    } else{
        stop("File extension unknown, allowed formats: RDS and .RData.")
    }
    
    if( !inherits(omicData, "ExpressionSet") & !inherits(omicData, "GenomicRatioSet") ) {
        stop("The omics data is not of the expected type, a GenomicRatioSet object or an ExpressionSet is expected")
    } 
    
    # Phenotypes data in ExpressinSet or GenomicRatioSet (technical data + common phenotypes)
    omicData_technicVars <- colnames(pData(omicData))
    
    # Transform str phenotypes to factor phenotypes
    phenotype <- as.data.frame(unclass(phenotype), stringsAsFactors=TRUE)
    
    # Load metadata with technical variables
    omicMetadata <- get(load(metadataFile)[which(load(metadataFile) %like% 'omics_metadata')])
    omicMetadataTechVars <- get(load(metadataFile)[which(load(metadataFile) %like% 'cdb_Full')])
    technicVars <- colnames(omicMetadata)[ which(colnames(omicMetadata) %like% omicsLayer) ]
    
    # Removing subjects withous variables from expressionSet or SummarizedExperiment
    omicData <- omicData[ ,sampleNames(omicData) %in% omicMetadata[, phenotypeID]] 
    
    # Filter by population (if not 'ALL' populations are selected)
    if( !'ALL' %in% population) {
        # Filter samples by population
        Sample.filter <- omicMetadata[, phenotypeID][which(omicMetadata$FINAL_ancestry %in% population)]
        if( length(Sample.filter)>0 ){
            omicData <- omicData[ ,sampleNames(omicData) %in% Sample.filter ] 
            omicMetadata <- omicMetadata[which(omicMetadata[, phenotypeID] %in% Sample.filter), ]
        } else {
            stop("No samples with ",population," ancestry  ")
        }
    }
    
    # Remove metadata from samples not present in GenomicRatioSet or ExpressionSet
    omicMetadata.f <- merge( data.frame( SampleID = pData(omicData)[, phenotypeID]) , omicMetadata, by = phenotypeID)
    omicMetadata.f.ord <- omicMetadata.f[ order( match(omicMetadata.f[, "SampleID"], pData(omicData)[, phenotypeID])), ] # Sort metadata SampleID == methylome SampleID
    
    # keep technical variables from omics expressionset or genomicratioset
    tokeep <- intersect(omicData_technicVars, c(sample_collect_cdb_Full[sample_collect_cdb_Full$Type %in% c("ID", "General", "Child", "Blood_cells", "Blood_collection", "DNA_extraction"),]$var))
    technicVars <- c( tokeep, technicVars)
    
    # Remove variables not related with Current Omics
    toRemove <- sample_collect_cdb_Full[ which( ( is.na( sample_collect_cdb_Full$RelatedOmics) | tolower( sample_collect_cdb_Full$RelatedOmics) != tolower(omicsLayer)) & 
                                                 !sample_collect_cdb_Full$var %in% technicVars ) ,]$var
    
    omicMetadata.f.ord <- omicMetadata.f.ord[ , -which(colnames(omicMetadata.f.ord) %in% toRemove)]    
    
    # Get common Vars in omics and phenotypes for merging purposes
    omicMetadataVars <- colnames(omicMetadata.f.ord)
    phenotypeVariables <- unique(phenotypeVariables)
    
    if( any(phenotypeVariables == 'ALL') || 'ALL' %in% phenotypeVariables){
        omicMetadata.f.ord <- merge( omicMetadata.f.ord, phenotype, by=phenotypeID)
        final_phenoVars <- colnames(phenotype)[-which(colnames(phenotype) == phenotypeID)]
    } else {
        phenotype <- phenotype[, unique(c(phenotypeID, phenotypeVariables))]
        omicMetadata.f.ord <- merge( omicMetadata.f.ord, phenotype, by = intersect(omicMetadataVars, colnames(phenotype)) )
        final_phenoVars <- phenotypeVariables[-which(phenotypeVariables == phenotypeID)]
    }
    
    omicData <- omicData[ ,sampleNames(omicData) %in% omicMetadata.f.ord$SampleID]
    rownames(omicMetadata.f.ord) <- omicMetadata.f.ord$SampleID
    
    if( inherits(omicData, "ExpressionSet") ) {
        phenoData(omicData) <- new("AnnotatedDataFrame", omicMetadata.f.ord)        
    } else if (inherits(omicData, "GenomicRatioSet")) {
        omicData <- minfi:::.pDataAdd(omicData, omicMetadata.f.ord)
    } 
    
    sampleNames(omicData) <- omicMetadata.f.ord$SampleID

    return(list( data = omicData, 
                 sampes = ncol(omicData),
                 nfeatures = nrow(omicData),
                 nvariables = ncol(pData(omicData)),
                 phenoVars = final_phenoVars,
                 technicVars = technicVars))
}


# Create formula to get SVA for each method: 
#       Interest: With all variables - only protect the phenotype of interest
#       Full: Protect all covariables - protect the phenotype of interest and all the covariates - apply SVA only in technical variables
#       Partial: Protect some covariables - protect the phenotype of interest and all come covariates - apply SVA in technical variables and some covariates

get_SVAFormulas <- function(svaMode, interest_var, covariates, to_protect_covariates, oType) {
    
    if(oType == 'ExpressionSet') {
        leftForm <- "t(exprs(omicSVA$data))"
    } else {
        # leftForm <- "assay(omicSVA$data)"
        leftForm <- "t(getBeta(omicSVA$data))" # Mandatory o deixar-ho com un opcional??
    }
    
    if('Interest' %in% svaMode) {
        if( interest_var !='' && !is.na(interest_var) ) {
            rightForm_inter <- interest_var
        }
    }
    
    if('Full' %in% svaMode) {
        if(length(covariates)>0 && covariates[1]!='' && !is.na(covariates[1]) && interest_var !='' && !is.na(interest_var)) {
            rightForm_full <- paste0( c(interest_var, covariates), collapse = " + ")
        }
    }
    
    if('Partial' %in% svaMode) {
        if(length(to_protect_covariates)>0 && to_protect_covariates[1]!='' && !is.na(to_protect_covariates[1]) ) {
            rightForm_partial <- paste0( c(interest_var, to_protect_covariates), collapse = " + ")
        }
    }
    return( list( f_Interest = ifelse( exists("rightForm_inter"), paste0(leftForm, " ~ ", rightForm_inter), NA ) , 
                  f_Full = ifelse( exists("rightForm_full"), paste0(leftForm, " ~ ", rightForm_full), NA ) , 
                  f_Partial = ifelse( exists( "rightForm_partial" ), paste0(leftForm, " ~ ", rightForm_partial), NA )
    )
    )
}



## ////////////////////////////////////////////////////////////////////////////////////////////////////////////
## ///                              COMPUTE SURROGATE VARIABLES (SVA)                                       ///
## ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Compute SVA variables
get_SVA <- function(svaMode, formulas, omicSVA, omicsLayer, interest_var, plot_vars = NULL, prefix="", respath = NULL, maxAssocPlots = 10, PCs = FALSE) {
    # Prepare paths
    if(is.null(respath) || respath==""){
        respath <- "./"
    } else {
        if(str_sub(respath, -1) != '/')
            respath <- paste0(respath, "/")
    }
    plot_path <- file.path(getwd(), respath, interest_var, "plots" )
    data_path <- file.path(getwd(), respath, interest_var )
    dir.create( data_path, showWarnings = FALSE)
    dir.create( plot_path, showWarnings = FALSE)
    
    # browser()
    SurrogateVars <- sapply(svaMode, function(mode){
        
        print(paste0("Computing ", mode, " model"))
        
        if(prefix!="") {
            bfilename = paste0( "SVs_", prefix, "_", omicsLayer, "_", mode)    
        } else {
            bfilename = paste0( "SVs_", omicsLayer, "_", mode)
        }
        
        useFormula <- names(formulas)[which(str_replace(names(formulas), "f_", "") == mode)]
        if( !is.na(formulas[[useFormula]]) ) {
            
            if(class(omicSVA$data)=="GenomicRatioSet") {
                assay(omicSVA$data)[assay(omicSVA$data)==0] <- 0.0001 #esto es importantisimo para trabajar luego con Mvalues (sino, aparecen -Inf)
            }
          
            options(na.action = "na.omit")
            Y.r <- t( resid( lm( as.formula(formulas[[useFormula]]), data = pData(omicSVA$data) ) ) )
            
            n.sv <- EstDimRMT( Y.r, TRUE )$dim + 1 # Adding one extra dimension to compensate potential loss of 1 degree of freedom
            options(na.action = "na.pass")
            mod <- model.matrix( as.formula(gsub('^.*(~)', '\\1',formulas[[useFormula]] )), data = pData(omicSVA$data) ) # Creating model
            nas_remove <- row.names(mod)[apply(mod, 1, function(x) any(is.na(x)))]
            options(na.action = "na.omit")
            mod <- mod[!rownames(mod) %in% nas_remove,]
            data_object <- gsub( "\\)$" ,"" , gsub( "^t\\("  ,"" , gsub(".~.*","", formulas[[useFormula]])))
            sv.obj <- tryCatch( {
                    smartsva.cpp( eval(parse(text = data_object))[,!colnames(eval(parse(text = data_object))) %in% nas_remove], mod, mod0 = NULL, n.sv = n.sv )   # creating SVs
                }, error=function(cond) {
                    message("Error processing smartsva.cpp - using sva() ")
                    sv.sva <- tryCatch( {
                            sva(  eval(parse(text = data_object))[,!colnames(eval(parse(text = data_object))) %in% nas_remove], mod, mod0 = NULL, n.sv = n.sv )
                        }, error=function(cond) {
                            message("Error processing smartsva.cpp and SVA - NO SVA COMPUTED ")
                            return(NA)
                        } )
                    return(sv.sva)
                } )
            if( all(!is.na(sv.obj))) {
                
                # Get data to be tested
                test_sv <- sv.obj$sv
                test_sv <- as.data.frame(test_sv)
                rownames(test_sv) <- rownames(pData(omicSVA$data))[!rownames(pData(omicSVA$data)) %in% nas_remove]
                colnames(test_sv) <- paste0("SV",1:ncol(test_sv))
                
               
                # Add phenotype data ¿Add more? 
                test_sv <- merge(test_sv, data.frame(pData(omicSVA$data)[!rownames(pData(omicSVA$data)) %in% nas_remove , c(interest_var, plot_vars)]), by = 0) %>% 
                    tibble::column_to_rownames(var = "Row.names")
                
                if( length(plot_vars) > 0 && all(plot_vars!="") && !is.null(plot_vars) && all(!is.na(plot_vars)) ) {
                    sapply(plot_vars, function(var) {
                        ggplot(test_sv, aes(SV1, SV2) )+
                            geom_point(aes(colour = var))
                        ggsave(paste0(plot_path, "/", bfilename, "_SV1_2_",var,".pdf"), width = 4, height = 4)
                        
                        ggplot(test_sv, aes(SV2, SV3) )+
                            geom_point(aes(colour = var))
                        ggsave(paste0(plot_path, "/", bfilename, "_SV2_3_",var,".pdf"), width = 4, height = 4)
                    })
                }
                
                pve <- variance_explained_sv(eval(parse(text = data_object))[,!colnames(eval(parse(text = data_object))) %in% nas_remove], sv.obj$pprob.gam, sv.obj$pprob.b, n.sv, bfilename, interest_var)
                if( length(omicSVA$phenoVars) > 0 ){
                    associat <- sapply( omicSVA$phenoVars, association, model = mode, df = pData(omicSVA$data)[!rownames(pData(omicSVA$data)) %in% nas_remove,], svs=sv.obj$sv, curvar = interest_var, maxsvs = maxAssocPlots, PCs = PCs, simplify = F)
                    
                    if( length(associat) > 0 ) {
                        write.csv2( do.call(rbind,associat), paste0( data_path, "/Summary_association_", bfilename, ".csv"))
                    }
                }
            }
        }
        return( SVs = sv.obj)
        
    }, simplify = F)
    
    # Create RObject with all the surrogate variables obtained with all methods
    save(SurrogateVars, file = paste0( data_path, "/SurrogateVars.RData") )    
    return(SurrogateVars)
}



## ///////////////////////////////////////////////////////////////////////////////////////////
## ///                              RESIDUALIZE DATA                                       ///
## ///////////////////////////////////////////////////////////////////////////////////////////

# Returns residualized data 
get_SVsResiduals <- function( resSVA, n_SVs, dd, SVs_list = NULL, respath = NULL) {
    
    # Prepare paths
    if(is.null(respath) || respath==""){
        respath <- "./"
    } else {
        if(str_sub(respath, -1) != '/')
            respath <- paste0(respath, "/")
    }
    
    
    residualized <- sapply( resSVA, function(data, SVs_list, dd) {
        
        if( length(names(data)) > 0 ) {
            
            if(length(SVs_list)>=1 && !is.na(SVs_list[1])){
                SVars <- SVs_list
            } else {
                if(is.na(n_SVs) || data$n.sv < n_SVs) 
                    n_SVs <- data$n.sv
                SVars <- c(1:n_SVs)
            }
            
            if(class(dd) == 'ExpressionSet') {
                
                model <- model.matrix( ~ data$sv[,SVars] )
                res <- residuals( lmFit( exprs(dd), model ), exprs(dd) )
                
                exprs(dd) <- res
                return(dd)
                
            } else {
                
                model <- model.matrix( ~ data$sv[,SVars] )
                assay(dd)[assay(dd)==0] <- 0.0001 
                res <- residuals( lmFit( getM(dd), model ), getM(dd) )
                
                assay(dd) <- res
                return(dd)
                    
            }    
        }else {
            return(NA)
        }
        
        
        
    }, SVs_list = SVs_list, dd = dd)
    
    # Create RObject with all the surrogate variables obtained with all methods
    # save(residualized, file = paste0( "Residualized.RData") )  
    save(residualized, file = paste0( respath, interest_var, "/Residualized.RData") )    
    return(residualized)
    
}


## ////////////////////////////////////////////////////////////////////////////////////
## ///                              DENOISING                                       ///
## ////////////////////////////////////////////////////////////////////////////////////


denoising_mlmer <- function( omicsData, omicLayer, noise_vars, ncores = 2, blocksize = 100) {
    data_path <- file.path(getwd(), "db", "descriptives" )
    dir.create( data_path, showWarnings = FALSE)

    # IMPUTE MISSING DATA:
    sapply( noise_vars, function(var) {
        if(is.factor( pData(omicsData)[,var])) {
            pData(omicsData)[ which(is.na(pData(omicsData)[,var])) , var] <- sample(pData(omicsData)[,var], size = length(sum( is.na(pData(omicsData)[,var]))))
        }else if(is.numeric( pData(omicsData)[,var])) {
            pData(omicsData)[ which(is.na(pData(omicsData)[,var])) , var] <- mean(pData(omicsData)[,var],na.rm=T)
        }else {
            message("Missing data not imputed")
        }
    } )
    
    noisedata <- as.data.frame(pData(omicsData)[ , noise_vars])
    sink(paste0(data_path,"/", omicLayer,"_denoising_initial_descriptive.log"))
        print(apply(noisedata,2,Hmisc::describe))
    sink(file = NULL)
    
    if(class(omicsData) == "ExpressionSet") {
        X <- t(exprs(omicsData))
    } else {
        assay(omicsData)[assay(omicsData)==0] <- 0.0001
        X <- t(getM(omicsData))
    }
        
    cols <- seq_len(ncol(X))

    if( ceiling(ncol(X)/blocksize) <= 1) {
        index <- split(cols, 1)
    } else {
        index <- split(cols, cut(cols, breaks=ceiling(ncol(X)/blocksize)))    
    }
    
    Xi <- lapply(index, function(v,index) v[,index], v=X) 
    
    rm(X)
    
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    t0 <- Sys.time()
    
    # Export mlmer_local to all cores
    clusterExport(cl,"mlmer_local")
    
    # Progress bar definition
    registerDoParallel(cl)
    registerDoSNOW(cl)
    iterations <- length(Xi)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Use %dopar% to parallelize all the process
    valDenoised <-  foreach(X_ = Xi, .combine=rbind, .options.snow= opts) %dopar% {
        strformulai <- paste0("X_ ~ ", paste0("( 1 | ", noise_vars, ")", collapse = " + ") )
        modeli <-  mlmer_local(formula( strformulai ), data = noisedata, vars = X_ , save.residuals = TRUE, save.ranks = FALSE, ncores = detectCores()-2)
        return(modeli)
    }
    # Printing time difference
    timec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    horas <- as.integer(timec / 3600)
    minutos <- as.integer((timec %% 3600) / 60)
    segundos <- as.integer(timec %% 60)
    tiempo_formateado <- sprintf("%02d:%02d:%02d", horas, minutos, segundos)
    print(paste0("Execution time: ", tiempo_formateado))
    
    close(pb)
    parallel::stopCluster(cl)
    
    omicsDenoised <- omicsData
    
    if(class(omicsData) == "ExpressionSet") {
        exprs(omicsDenoised) <- valDenoised
    } else {
        assay(omicsDenoised, withDimnames=FALSE) <- lumi::m2beta(valDenoised)    
    }
    
    dir.create(file.path(getwd(), "db" ), showWarnings = FALSE, recursive = TRUE)
    saveRDS(omicsDenoised, paste0(file.path(getwd(), "db" ),"/", omicLayer,"_denoised.rds"))
    
    return(omicsDenoised)
        
}


denoising_lm <- function( omicsData, omicLayer, noise_vars, ncores = 2, blocksize = 100) {
  
  data_path <- file.path(getwd(), "db", "descriptives" )
  dir.create( data_path, showWarnings = FALSE)
  
  # IMPUTE MISSING DATA:
  sapply( noise_vars, function(var) {
    if(is.factor( pData(omicsData)[,var])) {
      pData(omicsData)[ which(is.na(pData(omicsData)[,var])) , var] <- sample(pData(omicsData)[,var], size = length(sum( is.na(pData(omicsData)[,var]))))
    }else if(is.numeric( pData(omicsData)[,var])) {
      pData(omicsData)[ which(is.na(pData(omicsData)[,var])) , var] <- mean(pData(omicsData)[,var],na.rm=T)
    }else {
      message("Missing data not imputed")
    }
  } )
  
  noisedata <- as.data.frame(pData(omicsData)[ , noise_vars])
  sink(paste0(data_path,"/", omicLayer,"_denoising_initial_descriptive.log"))
  print(apply(noisedata,2,Hmisc::describe))
  sink(file = NULL)
  
  if(class(omicsData) == "ExpressionSet") {
    X <- t(exprs(omicsData))
  } else {
    assay(omicsData)[assay(omicsData)==0] <- 0.0001
    X <- t(getM(omicsData))
  }
  
  cols <- seq_len(ncol(X))
  
  if( ceiling(ncol(X)/blocksize) <= 1) {
    index <- split(cols, 1)
  } else {
    index <- split(cols, cut(cols, breaks=ceiling(ncol(X)/blocksize)))    
  }
  
  Xi <- lapply(index, function(v,index) v[,index], v=X) 
  
  rm(X)
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  t0 <- Sys.time()
  
  # Export mlmer_local to all cores
  clusterExport(cl,"lm_local")
  
  # Progress bar definition
  registerDoParallel(cl)
  registerDoSNOW(cl)
  iterations <- length(Xi)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Use %dopar% to parallelize all the process
  valDenoised <-  foreach(X_ = Xi, .combine=rbind, .options.snow= opts) %dopar% {
    strformulai <- paste0("X_ ~ ", paste0(noise_vars, collapse = " + ") )
    modeli <-  lm_local(formula( strformulai ), data = noisedata, vars = X_ , save.residuals = TRUE, save.ranks = FALSE, ncores = 6L)
    return(modeli)
  }
  timec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  horas <- as.integer(timec / 3600)
  minutos <- as.integer((timec %% 3600) / 60)
  segundos <- as.integer(timec %% 60)
  tiempo_formateado <- sprintf("%02d:%02d:%02d", horas, minutos, segundos)
  print(paste0("Execution time: ", tiempo_formateado))
  close(pb)
  parallel::stopCluster(cl)
  
  omicsDenoised <- omicsData
  
  if(class(omicsData) == "ExpressionSet") {
    exprs(omicsDenoised) <- valDenoised
  } else {
    assay(omicsDenoised, withDimnames=FALSE) <- lumi::m2beta(valDenoised)    
  }
  
  dir.create(file.path(getwd(), "db" ), showWarnings = FALSE, recursive = TRUE)
  saveRDS(omicsDenoised, paste0(file.path(getwd(), "db" ),"/", omicLayer,"_denoised.rds"))
  
  return(omicsDenoised)
  
}


check_noised_denoised <- function( original, denoised, Vars, prefix, topPCA = NULL, output_path) {
  data_path <- file.path(output_path)
  sapply(c('original', 'denoised'), function(data) {
    if( inherits(eval(parse(text=data)), "ExpressionSet")) {
      X <-  t(exprs(eval(parse(text=data))))
    } else if (inherits(eval(parse(text=data)), "GenomicRatioSet") ) {
      X <-  t(getBeta(eval(parse(text=data))))
    } else {
      stop("The omics data is not of the expected type, a GenomicRatioSet object or an ExpressionSet is expected")
    }
    
    # original data
    NoValidData <- which( rowSums(is.na(X[,-1]))>0 | rowSums(is.infinite(X[,-1]))>0)
    
    # Remove NA to get PCA
    if(length(NoValidData)>0) {
      pcaX <-  prcomp(X[-NoValidData,])
    } else {
      pcaX <-  prcomp(X)}
    
    ev <-  with(pcaX, sdev^2/sum(sdev^2))
    
    pdf( paste0(data_path,"/", prefix,"_",data,"Data.pdf"))
    plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
         ylab = "Proportion of EV", ylim = c(0, 1), cex = 0.3)
    points(cumsum(ev), pch = 19, col = "tomato", cex = 0.3)
    legend("topleft", pch = 19, col = c("navy", "tomato"),
           legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)
    dev.off();
    sum(!cumsum(ev) > 0.9)
    sink(paste0(data_path,"/", prefix,"_",data,"Data_VE.log"))
    cat( paste0("\nEV: \n", paste0(ev, collapse = ", "), "\n"))
    cat( paste0("\nNumber of PCs explaining 90% of variance: ", sum(!cumsum(ev) > 0.9)))
    cat( paste0("\nSamples excluded from PCA (NA or Infinite values): ", length(NoValidData) ))
    cat( paste0("\n\t", paste0( names(NoValidData), collapse = "\n\t")))
    sink(file = NULL)
    
    if( is.null(topPCA) || is.na(topPCA) || topPCA> ncol(pcaX$x)){
      topPCA == ncol(pcaX$x)
    }
    
    keycovar <- pData(eval(parse(text=data)))[, c(Vars)[which(!is.na(c(Vars)))]]
    
    rownames(keycovar) <- rownames( pData(eval(parse(text=data))) )
    rownames(pcaX$x) <- rownames( pData(eval(parse(text=data))) )
    
    toplot <- merge(pcaX$x[ , 1:topPCA], keycovar, by = 0 ) %>% 
      tibble::column_to_rownames(var = "Row.names") 
    
    res <- hetcor(toplot, parallel = TRUE, ncores=detectCores(logical=TRUE)-2)$correlations
    
    pdf( paste0(data_path,"/", prefix,"_",data,"_Correlation_Plot.pdf"), width = 8.3)
    corrplot.mixed(res, upper = "ellipse", lower = "number",
                   tl.pos = "lt", tl.col = "black", tl.offset=1, tl.srt = 40, 
                   tl.cex = 0.5, number.cex = 0.5)
    dev.off()
    
  })
}

