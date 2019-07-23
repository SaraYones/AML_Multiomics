mutationStatePerSample=read.maf('/Users/sarayones/Desktop/Projects/Integrative_AML_Project/Integrative_AML_Project/TARGET/WXS/WXS/L3/mutation/BCM/VerifiedSomatic/TARGET_AML_WXS_somatic_filtered_verified.maf⁩', useAll = TRUE)

#summarize mutations
maf =read.maf('TARGET_AML_WXS_somatic_filtered_verified1.maf', useAll = T, verbose = T)


Diagnosis<-read.xlsx("TARGET_AML_SampleMatrix_Discovery_20180118.xlsx",header = FALSE, sheetName = "Diagnosis")
Relapse<-read.xlsx("TARGET_AML_SampleMatrix_Discovery_20180118.xlsx",header = FALSE, sheetName = "Relapse")
Diagnosis=as.character(Diagnosis[,1])
Relapse=as.character(Relapse[,1])
#Work with CNV data
CNV<-read.table("eventsByPatientPerGene.target_aml.txt",header = FALSE)
colnames<-CNV[,1]
rownames<-unique(CNV[,3])
colnames(CNV)<-c("gene_id","cnv_status","samples")

#https://stackoverflow.com/questions/5890584/how-to-reshape-data-from-long-to-wide-format
View(reshape(CNV, idvar = "gene_id", timevar = "samples", direction = "wide"))

df_wide <- reshape(CNV, idvar="gene_id", timevar="samples", v.names="cnv_status", direction="wide")

#remove cnv.status from colnames https://stackoverflow.com/questions/37800704/r-remove-parts-of-column-name-after-certain-characters
colnames(df_wide)<-sub("cnv_status.", "", colnames(df_wide))

df_wide<- data.frame(df_wide, row.names = 1)

rownames=rownames(df_wide)
#Change the matrix from factors to charachters
df_wide <- data.frame(lapply(df_wide, as.character), stringsAsFactors=FALSE)

df_wide[is.na(df_wide)] <- 0

rownames(df_wide)<-rownames

hello=prepareDT(df_wide,c(1,2,3,4,5))

hello=t(hello)


mdata <- melt(CNV, id=c("gene_id","samples"))
subjmeans <- cast(mdata, gene_id+samples~value)

prepareDT=function(DT,values)
{
  print("hello")
 # temp=NULL
 # temp=DT[,1:dim(DT)[2]-1]
  temp=as.matrix(DT)
  print("i am here")
  temp[temp=="cnLOH"]<-values[1]
  temp[temp=="DELHom"]<-values[2]
  temp[temp=="GainSomatic"]<-values[3]
  temp[temp=="HETSomatic"]<-values[4]
  temp[temp=="LOH"]<-values[5]
  print("I am here at last")
  temp=apply(temp,2,as.numeric)
  print("here?")
  temp=as.data.frame(temp)
  #decisiontemp=unlist(lapply(DT$decisionSLE,as.character))
  #temp=cbind.data.frame(temp,decisiontemp)
  #temp=as.data.frame(apply(temp,2,as.factor))
  #names(temp)[names(temp)=="decisiontemp"]="decisionSLE"
  
  #temp$decision=as.factor(temp$decision)
  #temp$decisionSLE=as.integer(temp$decisionSLE)
  rownames(temp)<-rownames(DT)
  print("hi")
  return (temp)
}

library(stringr)

AddMutationToBinaryMatrix=function()
{
  
  AllGenesMutation= getGeneSummary(maf)$Hugo_Symbol
  BinaryMatrix =matrix(0L, nrow = dim(hello)[1], ncol = length(AllGenesMutation))
  rownames(BinaryMatrix)=rownames(hello)
  rownames(BinaryMatrix)<-str_replace_all(rownames(BinaryMatrix), "\\.", "-")
  colnames(BinaryMatrix)<-AllGenesMutation
 # row.names(BinaryMatrix)=SamplesMutationbarcode
#  row.names(detailedBinaryMatrix)=SamplesMutationbarcode
  mutationState=getSampleSummary(maf)
  #For SurvivalClassifier Script
  
  #significant_genes=unique(append(significant_genes_cnv,significant_genes_mutation))
  
  for (i in 1: length(AllGenesMutation))
  {
      print(i)
      MutatedSamplesPerGeneResults=genesToBarcodes(maf,genes=AllGenesMutation[i])
      #transform the result of Mutated Samples Per genes into matrix so that I can access it easily
      MutatedSamplesPerGeneResults=as.matrix(MutatedSamplesPerGeneResults[1])
      MutatedSamplesPerGene=as.character(MutatedSamplesPerGeneResults[[1]]$Tumor_Sample_Barcode)
      #Apply an Or operation with the already existing value inside the matrix
      if(length(MutatedSamplesPerGene)!=0)
      {
        for(j in 1:length(MutatedSamplesPerGene))
        {
          print("Hello ")
          #Because we are interested of samples which have mutation information and Structural variants information
          if(MutatedSamplesPerGene[j] %in% rownames(BinaryMatrix) )
          {
            #print(significant_genes[i])
            #CNVvalue=as.integer(BinaryMatrix[MutatedSamplesPerGene[j],significant_genes[i]]) 
         
           print("HI THERE PRINT")
            print(MutatedSamplesPerGene[j])
            print(AllGenesMutation[i])
              BinaryMatrix[MutatedSamplesPerGene[j],AllGenesMutation[i]]=1
  
            
          }
        }
      }
      
    }
  return(newList)  
  }
 
expression=read.csv("Normalized-Annotated.csv", header = TRUE)
expression=expression[, !(colnames(expression) %in% c("ENSMBL_ID","X"))]
expression<- data.frame(expression, row.names = 1)


  