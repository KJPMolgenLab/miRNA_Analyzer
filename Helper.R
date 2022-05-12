### helper 







library(shiny)
runGitHub("KJPMolgenLab/miRNA_Analyzer", ref="main")


input=c()
input$Readsfile="C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/Bär/Data20220511/AC_Baer_Summary.xlsx"
input$Samplefile = "C:/Users/andreas_chiocchetti/OneDrive/Documents/Frankfurt Uni/Kooperationen/Bär/Data20220511/AC_Baer_Meta_all.xlsx"
input$submit = TRUE
input$sampleID =   "SampleID"
input$target = "Target"
input$treatment = "PIBO"
input$covariates = NULL
input$reads.cutoff = 10e5
input$p.cut = 0.05
input$cor.cut = 0.8
input$run_now = T
input$mincount = 0
input$minInNetwork = 15
input$organism="hsapiens"
input$minsamples = 0
input$sdcutoff=3
input$cook="no"

metaRaw <-read.xlsx(input$Samplefile, sheet=1)
ReadsRaw <-read.xlsx(input$Readsfile, sheet="miRNA_piRNA")

myvalues=c()
myvalues$treatment="Hypox"

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

