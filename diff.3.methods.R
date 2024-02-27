# Load required libraries
library(nlme)                    
library(lme4)                    
library(lattice)
library("stats")
library(dplyr)
library ("plyr")
library(broom)


# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
abundance <- args[1]  #"ecs_3.tsv"
metadata <- args[2]    #meta.txt

# Data reading and preparation
data <- read.csv(abundance, sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
data1 <- as.data.frame(t(data))
#data arcsine square-root transformation
#numeric_columns <- sapply(data1, is.numeric)
#data1[numeric_columns] <- asin(sqrt(data1[numeric_columns]))
data1$ID <- rownames(data1)
data1$ID<-sub("_level4ec", "", data1$ID)
data1 <- data1[, c("ID", setdiff(names(data1), "ID"))] 

meta <- read.csv(metadata, sep = "\t", check.names = FALSE, header = TRUE)
x <- merge(meta, data1, by = "ID")

# function for linear regression model  age+sex+antibiotics+CD
linear_regression_group <- function(x) {
  gene_name2<-c()
  gene_name<-c()
  pvalue_CD<-c()
  OR_CD<-c()
  LOW_CD<-c()
  HIGH_CD<-c()
  pvalue_UC<-c()
  OR_UC<-c()
  LOW_UC<-c()
  HIGH_UC<-c()
  pvalue_Male<-c()
  OR_Male<-c()
  LOW_Male<-c()
  HIGH_Male<-c()
  pvalue_Senior<-c()
  OR_Senior<-c()
  LOW_Senior<-c()
  HIGH_Senior<-c()
  number_1<-c()
  number_2<-c()
  for (i in 16:ncol(x))
  {
    g1=colnames(x)[i]
    tryCatch({ 
      gro <- x[,5]
      species<-x[,g1]
      aa <-na.omit(data.frame(gro,species))
      model <- lm(x[,i] ~ phenotype + age_category+gender+  antibiotics, data = x) #
      tidy_model <- tidy(model, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)
      gene_name2=c(gene_name2,g1)
      number<-ddply(aa,.(gro),nrow)
      number_1=c(number_1,number[1,2])
      number_2=c(number_2,number[2,2])
      pvalue_CD=c(pvalue_CD,as.numeric (tidy_model[2,5]))
      OR_CD=c(OR_CD,as.numeric (tidy_model[2,2]))
      LOW_CD=c(LOW_CD,as.numeric (tidy_model[2,6]))
      HIGH_CD=c(HIGH_CD,as.numeric (tidy_model[2,7]))
    }
    ,
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(data.frame(gene_name2,number_1,number_2,OR_CD,LOW_CD,HIGH_CD,pvalue_CD))
}

xIBD<-x[x$cohort!="Stinki",]
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD$phenotype <- factor(xIBD$phenotype)
xIBD$antibiotics <- factor(xIBD$antibiotics)
xIBD$age_category <- factor(xIBD$age_category)
xIBD$gender <- factor(xIBD$gender)
xIBD$phenotype <- relevel(xIBD$phenotype, ref = "HC")

myrsa_cor<-ddply(xIBD, .(cohort),linear_regression_group)
filename<-paste0(abundance,".lm.csv")
write.table(myrsa_cor,file=filename,sep = ",",row.names=FALSE)


# function for Linear Mixed-Effects -3 fixed cofactors
linear_mixed_group <- function(x) {
  results <- list()
  
  # Loop through the columns of interest
  for (i in 16:ncol(x)) {
    tryCatch({
      col_name <- colnames(x)[i]
      x$y <- x[, i]
      fit.m3 <- lme(y ~ phenotype   + age_category + gender+ antibiotics, #   
                    method = "ML", 
                    data = x, 
                    random = list(subject_accession = ~ 1, study_accession_database = ~ 1), # Corrected random effects formula
                    na.action = na.exclude)
      coef.fit.m3 <- summary(fit.m3)$tTable
      
      # Store results for each column in a list
      results[[col_name]] <- cbind(
        tax_name = col_name,
        diagnosisCD_value = coef.fit.m3[2, 1],
        diagnosisCD_t = coef.fit.m3[2, 4],
        diagnosisCD_p = coef.fit.m3[2, 5],
        #diagnosisUC_value = coef.fit.m3[3, 1],
        #diagnosisUC_t = coef.fit.m3[3, 4],
        #diagnosisUC_p = coef.fit.m3[3, 5],
        groupHC_mean = mean(x[x$phenotype == "HC", i], na.rm = TRUE),
        groupCD_mean = mean(x[x$phenotype == "CD", i], na.rm = TRUE) #,
        #groupUC_mean = mean(x[x$phenotype == "UC", i], na.rm = TRUE)
      )
    }, error = function(e){
      # Optionally, print the error message
      message("Error in column: ", col_name, "; Error: ", e$message)
    })
  }
  return (do.call(rbind, results))
}


myrsa_cor<-ddply(xIBD, .(cohort),linear_mixed_group)
filename<-paste0(abundance,".mix-effects.csv")
write.table(myrsa_cor,file=filename,sep = ",",row.names=FALSE)



#K-W diff#
diff<-function(x){
  gene_name2<-c()
  gene_name<-c()
  pvalue<-c()
  qvalue<-c()
  mean_1<-c()
  median_1<-c()
  mean_2<-c()
  median_2<-c()
  number_1<-c()
  number_2<-c()
  
  
  for (i in 16:ncol(x))
  {
    g1=colnames(x)[i]
    gene_name=c(gene_name,g1)
    val<-as.numeric(x[,i])
    gro <- x[,5]
    aa <-na.omit(data.frame(gro,val))
    tryCatch({y=kruskal.test(val~gro, data=aa)
    pvalue=c(pvalue,y$p.value)
    gene_name2=c(gene_name2,g1)
    mean <- tapply(aa[,2],aa[,1],mean)
    median <- tapply(aa[,2],aa[,1],median)
    number<-ddply(aa,.(gro),nrow)
    number_1=c(number_1,number[1,2])
    number_2=c(number_2,number[2,2])
    a<-as.data.frame(mean)
    b<-as.data.frame(median)
    mean_1=c(mean_1,a[1,])
    mean_2=c(mean_2,a[2,])
    median_1=c(median_1,b[1,])
    median_2=c(median_2,b[2,])
    }, error=function(e){
      print(i)})
  }  
  qvalue=p.adjust(pvalue,'BH')
  return (data.frame(gene_name2,mean_1, mean_2,number_1,median_1,median_2,number_2, pvalue,qvalue))
}


myrsa_cor<-ddply(xIBD, .(cohort),diff)
filename<-paste0(abundance,".kw.csv")
write.table(myrsa_cor,file=filename,sep = ",",row.names=FALSE)

