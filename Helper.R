### helper 





local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.rstudio.com/"
  r["BioCsoft"] <- "https://bioconductor.org/packages/release/bioc"
  r["BioCann"] <- "https://bioconductor.org/packages/release/data/annotation"
  r["BioCexp"] <- "https://bioconductor.org/packages/release/data/experiment"
  options(repos = r)
})



library(shiny)
runGitHub("KJPMolgenLab/miRNA_Analyzer", ref="main")

setwd("C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/schubert/miRNApp/miRNA_App/")

input=c()
input$Readsfile="C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/schubert/Analysis_20220808_Hanna/QiagenSummary_Hanna IgA+-.xlsx"
input$Samplefile = "C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/schubert/Analysis_20220808_Hanna/Metadata_Hanna.xlsx"

#input$Readsfile="C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/schubert/Sputum20220510/20202804 8Proben 107316.all_samples.summary (1).xlsx"
#input$Samplefile = "C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/schubert/Sputum20220510/Metadata_BO Sputum FFM.xlsx"

#input$Samplefile=input$Readsfile

input$submit = TRUE
input$sampleID =   "SampleID"
input$target = "Target"
input$treatment = "IgApos"
input$covariates = "Gender"
input$reads.cutoff = 10e5
input$p.cut = 0.05
input$cor.cut = 0.8
input$run_now = T
input$mincount = 1
input$minInNetwork = 15
input$organism="hsapiens"
input$minsamples = 0
input$sdcutoff=1
input$cook="no"
input$geneunivers="no"

metaRaw <-read.xlsx(input$Samplefile, sheet=1)
ReadsRaw <-read.xlsx(input$Readsfile, sheet="miRNA_piRNA")

myvalues=c()
myvalues$treatment=input$treatment

getlog2FC=function(datavector, groupfactor){
  sums = tapply(datavector, groupfactor, sum, na.rm=T)
  idxRef = names(sums) == input$Reference
  log2FC=log2(sums[!idxRef]/sums[idxRef])
  return(log2FC)
}


log2FCs = apply(Reads, 1, function(x){getlog2FC(unlist(x),groupfactor = meta$group)})

datavector = unlist(NReads["hsa-miR-147b",])
groupfactor= dds$group
getlog2FC(datavector=unlist(Reads[1,]), groupfactor=meta$group)
plot(log2FCs[rownames(restab)], restab$log2FoldChange)
boxplot(datavector~groupfactor)

