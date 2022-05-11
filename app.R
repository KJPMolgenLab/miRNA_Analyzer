#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.rstudio.com/"
    r["BioCsoft"] <- "https://bioconductor.org/packages/release/bioc"
    r["BioCann"] <- "https://bioconductor.org/packages/release/data/annotation"
    r["BioCexp"] <- "https://bioconductor.org/packages/release/data/experiment"
    options(repos = r)
})

#update.packages(ask = F)


packages=c(
    "shiny",
    "openxlsx",
    "data.table",
    "DT",
    "DESeq2",
    "tidyverse",
    "compareGroups",
    "ggplot2",
    "gplots",
    "plotly",
    "pheatmap",
    "viridis",
    "RColorBrewer",
    "grid",
    "visNetwork",
    "igraph",
    "gprofiler2",
    "miRNAtap",
    "miRNAtap.db",
    "foreach")

for(p in packages){
    if(!require(p, character.only = TRUE)){
        install.packages(p, dependencies=T) 
    }
}


library(shiny)
library(openxlsx)
library(data.table)
library(DT)
library(DESeq2)
library(tidyverse)
library(compareGroups)
library(ggplot2)
library(gplots)
library(plotly)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(grid)
library(visNetwork)
library(igraph)
library(gprofiler2)
library(miRNAtap)
library(miRNAtap.db)
library(foreach)


getGOresults = function(geneset, genereference, organism = "mmusculus",
                        sources=c("GO:BP", "GO:MF", "GO:CC", "KEGG", "TF",
                                  #"MIRNA",
                                  "CORUM", "HP", "HPA")){
    
    resgo = gost(geneset, organism =organism,
                 correction_method = "gSCS",
                 domain_scope = "custom",
                 sources = sources,
                 evcodes = F,
                 custom_bg = genereference,
                 numeric_ns = "ENTREZGENE_ACC")
    if(length(resgo) != 0){
        return(resgo)
    } else {
        print("no significant GO terms identified")
        return(NULL)
    }
}

organismlist=list(hsapiens="hsa", 
                  mmusculus="mmu")

AGCcol=RColorBrewer::brewer.pal(8,"Dark2")
if( file.exists("data/genedict.RData")){
    load("data/genedict.RData")
} else {
    genedict=list()
}

# go profiler function


empty_plot <- function(title = NULL){
    p <- plotly_empty(type = "scatter", mode = "markers") %>%
        config(
            displayModeBar = FALSE
        ) %>%
        layout(
            title = list(
                text = title, 
                y=0.5
            )
        )
    return(p)
} 


# Define UI for application that draws a histogram
ui <- function(){
    fluidPage(
        
        # Application title
        titlePanel("miRNA Analyzer V2.1"),
        
        # Sidebar with a slider input for number of bins 
        sidebarLayout(
            
            sidebarPanel(width = 2,
                fileInput(inputId = 'Readsfile',
                          label = 'Upload Qiagen File',
                          accept = c('.xlsx')
                ),
                
                fileInput(inputId = 'Samplefile',
                          label = 'Upload Sample File',
                          accept = c('.xlsx')
                ),
                
                selectInput("sampleID", "SampleId", 
                            choices = NULL),
                
                selectInput("target", "Target", 
                            choices = "Target", selected="Target"),
                
                selectInput("treatment", "Treatment", 
                            choices = NULL),
                
                selectInput("covariates", "Covariates", 
                            multiple = TRUE, 
                            choices = NULL),
                selectInput("organism", "Organism",
                            choices = c("mmusculus", "hsapiens"), 
                            selected="hsapiens"),
                
                actionButton("run_now", "Run analysis"),
                hr(),
                h3("advanced settings"),
                numericInput("mincount", "QC: min Reads miRNA in at least 1 sample", 
                             min = 0, value = 10, step=1),
                numericInput("minsamples", "QC: % Sample with miRNA detected", 
                             min = 0, max=1, value = 0.5, step=0.05),
                numericInput("sdcutoff", "QC: exclude miRNA with a SD below cutoff", 
                             min = 0, value = 3, step=0.5),
                
                selectInput("reads.cutoff", "QC: Exclude samples with total Reads below", 
                            choices=c(0,
                                      10000,
                                      50000,
                                      100000, 
                                      500000,
                                      1000000), 
                            selected = 100000), 
                
                numericInput("p.cut", "DEG: p-value threshold for significance", 
                             0,1, value = 0.05, step = 0.001),
                selectInput("cook", "DEG: only display trustworthy p-values (cooksDistance corrected)", 
                             choices=c("yes", "no"), selected="no"),
                numericInput("minInNetwork", "Network: min Number of Cluster size", 
                             min=10, value=50, step=1),
                numericInput("cor.cut", "Network: Minimum correlation value between to miRNAs", 0,1, value=0.8, step = 0.01), 
                
            ),
            
            # Show a plot of the generated distribution
            mainPanel(width = 10,
                tabsetPanel(
                    type = "tabs",
                    tabPanel("Summary",
                             htmlOutput("summary")),
                    tabPanel("Differential Gene Expression", 
                             fluidRow("", 
                                      column(6, 
                                             fluidRow("", plotlyOutput("volcano_plot")),
                                             fluidRow("", hr()),
                                             fluidRow("", plotOutput("boxplot"))),
                                      column(6, DT::dataTableOutput("table_display_fc")))),
                    tabPanel("Heatmaps", 
                             fluidRow("",
                                      column(6, plotOutput("heatmap_samples")),
                                      column(6, plotOutput("DEG_heatmaps_Samples")),
                             ),
                             fluidRow("",
                                      column(6, plotOutput("DEG_heatmaps_miRNAs")),
                                      column(6, plotOutput("DEG_Corr_hist")))),
                    tabPanel("Network Plot", visNetworkOutput("Network", height="800px")),
                    
                    tabPanel("Target Enrichment", 
                             fluidRow("",tags$p("Analysis done with g:GOSt. 
                                                In addition to Gene Ontology (GO:MF = Molecular Function,GO:BP = Biological Processes, GO:CC = Cellular Compartment), 
                                                pathways from KEGG Reactome, 
                                                regulatory motif matches from TRANSFAC (TF), 
                                                tissue specificity from Human Protein Atlas (HPA); 
                                                protein complexes from CORUM and human disease phenotypes (HP) from Human Phenotype Ontology. 
                                                Website: https://biit.cs.ut.ee/gprofiler/gost")
                             ),
                             fluidRow("",
                                      tabsetPanel(
                                          tabPanel("All miRNAs",
                                                   fluidRow("",
                                                            column(12, plotlyOutput("GO_plot_all"))),
                                                   fluidRow("",
                                                            column(12, DT::dataTableOutput("table_GO_all")))),
                                          tabPanel("upregulated miRNAs",
                                                   fluidRow("",
                                                            column(12, plotlyOutput("GO_plot_up"))),
                                                   fluidRow("",
                                                            column(12, DT::dataTableOutput("table_GO_up")))),
                                          tabPanel("downregulated miRNAs",
                                                   fluidRow("",
                                                            column(12, plotlyOutput("GO_plot_dwn"))),
                                                   fluidRow("",
                                                            column(12, DT::dataTableOutput("table_GO_dwn"))))), 
                             )),
                    tabPanel("Samples Descriptives", 
                             fluidRow("",
                                      column(12, DT::dataTableOutput("table_display_samples"))),
                             fluidRow("", column(12, htmlOutput("covariate")))),
                    tabPanel("Reads Table", 
                             DT::dataTableOutput("table_display_reads")),
                    tabPanel("QC: Plots",
                             plotOutput("dds_plot"))
                    
                    
                )
            )
        )
    )
}


# Define server logic required to draw a histogram
server <- function(input, output, session){
    
    myvalues <-reactiveValues(
        ReadsRaw=NULL, 
        metaRaw=NULL, 
        dds=NULL, 
        meta=NULL, 
        Reads=NULL, 
        NReads=NULL, 
        exclude=NULL,
        sigReaddiff=NULL, 
        covariateComp=NULL, 
        Remove=NULL, 
        restab_DEseq=NULL)
    
    getSamples <- reactive({
        inFile <- input$Samplefile
        req(inFile)
        metaRaw <-read.xlsx(inFile$datapath)
        
        inFile <- input$Readsfile
        req(inFile)
        ReadsRaw <-read.xlsx(inFile$datapath, sheet="miRNA_piRNA")
        
        ReadsRaw <- remove_rownames(ReadsRaw)
        ReadsRaw <- column_to_rownames(ReadsRaw, "miRNA")
        ReadsRaw <- ReadsRaw[,grepl("READs",colnames(ReadsRaw))]
        colnames(ReadsRaw)=gsub("-READs","", colnames(ReadsRaw))
        ReadsRaw <- ReadsRaw[complete.cases(ReadsRaw),]
        myvalues$ReadsRaw=ReadsRaw
        myvalues$metaRaw=metaRaw
    })
    
    update_cols <- reactive({
        metaRaw <- myvalues$metaRaw
        vars <- names(metaRaw)
        # Update select input immediately after clicking on the action button.
        updateSelectInput(session, 
                          inputId = "target",
                          choices = vars)
        updateSelectInput(session, 
                          inputId = "sampleID",
                          choices = vars)
        updateSelectInput(session, 
                          inputId = "covariates",
                          choices = vars)
        
    })
    observeEvent(getSamples(),{
        update_cols()
    })
    
    
    
    observe({
        metaRaw <- myvalues[["metaRaw"]]
        CaCo = unique(metaRaw[,input$target])
        updateSelectInput(session, 
                          inputId = "treatment",
                          choices = CaCo, 
                          selected= CaCo[1])
    })
    
    observe({
        
        myvalues$treatment = input$treatment
        cat("treatment: ", myvalues$treatment, "\n")
    })
    
    observeEvent(input$run_now, {
        withProgress(message = "analyzing...", {
            print("getsamples")
            getSamples()
            incProgress(1/13)
            print("getsamples done")
            print("matchsamples")
            req(input$treatment)
            metaRaw <- myvalues$metaRaw
            Reads <- myvalues$ReadsRaw
            metaRaw <- remove_rownames(metaRaw)
            metaRaw <- column_to_rownames(metaRaw, input$sampleID)
            req(myvalues$treatment)
            meta <- metaRaw[!is.na(metaRaw[,input$target]),]
            Samples <- intersect(colnames(Reads), rownames(meta))
            if(length(Samples)==0){
                showNotification("Sampels between inputs files do not match (meta file must contain the same names as the Qiagen input wihtout the -READS flag)")
            }
            meta <-  meta[Samples, ]
            Reads <- Reads[,Samples]
            tmpnames=gsub("[^a-zA-Z0-9:]", "_", colnames(Reads))
            colnames(Reads) = tmpnames
            rownames(meta) = tmpnames
            incProgress(1/13)
            
            print("check metafile")
            meta$group <- as.factor(meta[,input$target]) 
            meta$group <- meta$group %>% relevel(ref=myvalues$treatment)
            meta$totalReads <- colSums(Reads, na.rm=T)
            
            idxRm <- meta$totalReads < as.numeric(input$reads.cutoff)
            if(sum(idxRm)>0.5*nrow(meta)){
                showNotification("More than half of the samples did not meet QC criteria; 
                                 all samples were thus included, please adjust the QC parameter reads cutoff")
                idxRm <- meta$totalReads < 0
            }
            
            
            Remove <- colnames(Reads)[idxRm]
            Reads <- Reads[,!idxRm]
            meta <- meta[!idxRm,]
            
            incProgress(1/13)
            covars = input$covariates
            exclude = names(which(apply(meta %>% dplyr::select(all_of(covars)), 2, function(x){length(unique(x))})==1))
            covars=covars[! covars %in% exclude]
            if(length(input$covariates)>0){
                modelform=as.formula(paste0("~1+",paste0(covars, sep="", collapse="+"), "+group"))
            }else {
                modelform=as.formula("~group")
            }
            
            Reads=Reads[,rownames(meta)]
            
            print("QC Reads")
            dds <- DESeqDataSetFromMatrix(Reads, colData = meta, design = modelform)
            
            keep <- rowSums(counts(dds)) >= input$mincount
            keep2 <- rowSums(counts(dds)> input$mincount) > input$minsamples*nrow(meta)
            keep3 <- apply(counts(dds), 1, sd)>input$sdcutoff
            
            myvalues$dropmirs <- sum(!keep | !keep2 | !keep3)
            
            dds <- dds[keep & keep2 & keep3,]
            incProgress(1/13)
            
            #dds <- estimateSizeFactors(dds)
            #dds <- estimateDispersions(dds)
            print("differential gene expression")
            dds <- DESeq(dds, fitType = "local", test = "Wald")
            incProgress(1/13)
            print("Descriptive stats meta file")
            NReads <- counts(dds, normalized=T)
            cs <- colSums(log(NReads+1))
            options(contrasts = c("contr.sum","contr.poly"))
            mod <- lm(cs~ meta$group)
            AOV <- drop1(mod, .~., test="F")
            p.val <- AOV$`Pr(>F)`[2]
            sig <- c()
            if(p.val>0.05){
                sig="not-significant"
            } else {
                sig="significant"
            }
            
            res <- compareGroups::compareGroups(group ~ ., meta[,c("group",input$target, input$covariates)])
            restab <- createTable(res, hide.no = "no")
            incProgress(1/13)
            print("generate DEG results tables")
            restab_DEseq <- results(dds, cooksCutoff=input$cook=="yes")
            restab_DEseq <- restab_DEseq[order(restab_DEseq$pvalue),]
            
            
            myvalues$sigNum=sum(restab_DEseq$padj<=input$p.cut, na.rm=T)
            myvalues$sigNumup=sum(restab_DEseq$padj<=input$p.cut & restab_DEseq$log2FoldChange>0, na.rm=T)
            myvalues$sigNumdwn=sum(restab_DEseq$padj<=input$p.cut  & restab_DEseq$log2FoldChange<0, na.rm=T)
            
            myvalues$nomsigNum=sum(restab_DEseq$pvalue<=input$p.cut, na.rm=T)
            myvalues$nomsigNumup=sum(restab_DEseq$pvalue<=input$p.cut& restab_DEseq$log2FoldChange>0, na.rm=T)
            myvalues$nomsigNumdwn=sum(restab_DEseq$pvalue<=input$p.cut& restab_DEseq$log2FoldChange<0, na.rm=T)
            incProgress(1/13)
            print("prepare Network files")
            if(myvalues$sigNum>=input$minInNetwork){
                myvalues$GraphAnal="adjusted p-value"
                myvalues$targetmirs=rownames(restab_DEseq)[which(restab_DEseq$padj<=input$p.cut)]
            } else {
                myvalues$GraphAnal="nominal p-value"
                myvalues$targetmirs=rownames(restab_DEseq)[which(restab_DEseq$pvalue<=input$p.cut)]}
            
            
            ## graph networks 
            myvalues$miRCo=cor(t(NReads[myvalues$targetmirs,]), method = "s")
            
            CorMAT = myvalues$miRCo
            
            resSig = restab_DEseq[myvalues$targetmirs,]
            
            CorMAT[ (CorMAT) < input$cor.cut] <- 0
            diag(CorMAT) <- 0
            print("calculate networks")
            graph <- graph.adjacency(CorMAT, weighted=TRUE, mode="lower")
            
            
            CEB=cluster_fast_greedy(graph) 
            
            rbPal <- colorRampPalette(c("red","green"))
            
            FCsig=resSig$log2FoldChange
            FCsig[is.infinite(FCsig)]=NA
            
            if(length(unique(meta$group))==2){
                LFCbnd=max(abs((FCsig)), na.rm=T)
                LFCsym=c(-LFCbnd, LFCbnd, FCsig)
                cols=rbPal(11)[as.numeric(cut(LFCsym, breaks=11))]
                cols=cols[-(1:2)]
                V(graph)$color=cols
                labsub="green=up-regulated, red=down-regulated"} else {
                    V(graph)$color=AGCcol[membership(CEB)]
                    labsub="colors = membership"
                }
            
            V(graph)$size=-log10(resSig$pvalue)
            plab=formatC(resSig$padj,format = "e", digits = 1)
            plab[resSig$padj>0.1]=""
            mirname=V(graph)$name
            V(graph)$name= paste(plab, "\n",V(graph)$name) 
            V(graph)$size=V(graph)$size
            
            V(graph)$label.cex=1
            V(graph)$label.color=AGCcol[membership(CEB)]
            V(graph)$label.dist=2
            
            E(graph)$width=(E(graph)$weight^5)*10
            
            layout=layout.fruchterman.reingold(graph, niter=10000)
            
            myvalues$graph = graph
            incProgress(1/13)
            print("write analysis outputs")
            
            myvalues$Rversion=R.version$version.string
            myvalues$NtotReads=round(sum(meta$totalReads)/10e5,digits = 1)
            myvalues$AVNreads=formatC(sum(meta$totalReads)/ncol(Reads), big.mark=",", format="fg")
            myvalues$Nmirna=nrow(Reads)
            
            myvalues$Nmirna_gt_5=sum(apply(Reads, 1, function(x){prod(x>=5)}))
            myvalues$Nmirna_gt_5_50=sum(apply(Reads, 1, function(x){sum(x>=5)>ncol(Reads)*0.5}))
            myvalues$Nmirna_gt_5=sum(apply(Reads, 1, function(x){prod(x>=5)}))
            
            myvalues$Confoundertext=paste0(covars, collapse=", ")
            
            myvalues$dds <- dds
            myvalues$meta <- meta 
            myvalues$Reads <- Reads 
            myvalues$NReads <- NReads 
            myvalues$exclude <- exclude
            myvalues$sigReaddiff <- sig 
            myvalues$covariateComp <- restab 
            myvalues$Remove <- Remove
            myvalues$restab_DEseq <- restab_DEseq
            
            
            if(myvalues$GraphAnal =="adjusted p-value"){
                index=which(restab_DEseq$padj<input$p.cut)
            } else 
            {index=which(restab_DEseq$pvalue<input$p.cut)}
            
            myvalues$restab_DEseq_sig  <-  restab_DEseq[index,]
            incProgress(1/13)
            miRNAunivers=rownames(restab_DEseq)
            miRNAunivers = grep("miR", miRNAunivers, value = T)
            
            #if the mir is not in the dictionary we check the proteins 
            print("enrichment analysis preparation")
            recalc = miRNAunivers[! miRNAunivers %in% names(genedict)]
            
            if(length(recalc)>0){
                print("reannotation of miRNAS")
                genedict_new=mapply(function(x){
                    targets <- getPredictedTargets(mirna = x, 
                                                   min_src = 4,
                                                   species = organismlist[[input$organism]])
                    #targets = head(targets, n=10)
                    return(targets)}, recalc)
                genedict <- c(genedict, genedict_new)
                genedict <- genedict[unique(names(genedict))]
                save(file = "data/genedict.RData",list = c("genedict"))
            }
            incProgress(1/13)
            
            filtgenedict= genedict[which(names(genedict)%in% miRNAunivers)]
            geneunivers = foreach(i=filtgenedict, .combine=c) %do% (rownames(i)) %>% unique()
            
            myvalues$geneunivers <- geneunivers
            myvalues$filtgenedict <- filtgenedict
            
            mirsreg <- grep("miR", rownames(myvalues$restab_DEseq_sig), 
                            value = T)
            
            genes_mirreg <- foreach(i=filtgenedict[mirsreg], 
                                    .combine=c) %do% (rownames(i)) %>% unique()
            mirsup <- grep("miR", rownames(myvalues$restab_DEseq_sig[myvalues$restab_DEseq_sig$log2FoldChange>0,]), 
                           value = T)
            genes_mirup <- foreach(i=filtgenedict[mirsup], 
                                   .combine=c) %do% (rownames(i)) %>% unique()
            
            mirsdwn <- grep("miR", rownames(myvalues$restab_DEseq_sig[myvalues$restab_DEseq_sig$log2FoldChange<0,]), 
                            value = T)
            genes_mirdwn <- foreach(i=filtgenedict[mirsdwn], 
                                    .combine=c) %do% (rownames(i)) %>% unique()
            
            myvalues$genes_mirreg <- genes_mirreg
            myvalues$genes_mirup <- genes_mirup
            myvalues$genes_mirdwn <- genes_mirdwn
            
            print("GOCalls")
            incProgress(1/13)
            
            
            myvalues$ResGO_all = getGOresults(myvalues$genes_mirreg,
                                     myvalues$geneunivers,
                                     input$organism)
            myvalues$ResGO_up = getGOresults(myvalues$genes_mirup,
                                    myvalues$geneunivers,
                                    input$organism)
            myvalues$ResGO_dwn = getGOresults(myvalues$genes_mirdwn,
                                     myvalues$geneunivers,
                                     input$organism)
            incProgress(1/13)
            showNotification("analysis complete")
        })
    })
    
    
    observe({
        req(myvalues$meta, myvalues$Reads, myvalues$dds, myvalues$graph)
        print("table_display_samples")
        output$table_display_samples <- renderDataTable(withProgress(message = "generating ...", 
                                                                     {
                                                                         f <- data.frame(myvalues[["meta"]] %>% as.data.frame() %>% 
                                                                                             dplyr::select(all_of(c("group", "totalReads",
                                                                                                                    input$covariates)))) #subsetting takes place here
                                                                         datatable(f, extensions = "Buttons",
                                                                                   options = list(
                                                                                       pageLength = 15,
                                                                                       info = FALSE,
                                                                                       lengthMenu = list(c(15,50, 100, -1), 
                                                                                                         c("15","50", "100" ,"All")
                                                                                       ), dom = 'Blfrtip',
                                                                                       buttons = c('copy', 'csv', 'excel', 'pdf')
                                                                                   )
                                                                         )
                                                                     }))
        print("table_display_reads")
        output$table_display_reads <- renderDataTable(withProgress(message = "generating ...", 
                                                                   {datatable(myvalues$Reads, 
                                                                              extensions = "Buttons",
                                                                              options = list(
                                                                                  autoWidth = TRUE,
                                                                                  pageLength = -1,
                                                                                  info = FALSE,
                                                                                  autoWidth=T, 
                                                                                  lengthMenu = list(c(-1, 15,50, 100), 
                                                                                                    c("All", "15","50", "100" )
                                                                                  ),dom = 'Blfrtip',
                                                                                  buttons = c('copy', 'csv', 'excel', 'pdf')
                                                                              )
                                                                   )
                                                                   }))
        
        print("summary")
        output$summary <- renderUI(withProgress(message = "generating ...", {
            HTML(paste0(
                "<h3>Summary</h3>",
                "<p>Total of ", ncol(myvalues$Reads), " overlapping Sample ID included from Sample file and Qiagen sheet</p>",
                "<p>Samples removed due to quality problems: ",
                paste0(myvalues$Remove, sep="", collapse=", "), "</p>",
                "<p>miRNAs removed due to quality problems: N=",myvalues$dropmirs, "</p>",
                "<p>Covariates removed due to no variance: ", paste0(myvalues$exclude, sep=" ", collapse="," ),"</p>",
                "<p>Total Number of Reads is ", formatC(sum(myvalues$Reads), format = "g"), " averaging to ", 
                formatC(sum(myvalues$Reads)/ncol(myvalues$Reads), format = "g", digits = 2), " Reads per sample</p>",
                "<p>Readcount after normlization is ", myvalues$sigReaddiff, " between the case and control cohorts</p>",
                "<p>To identify potential technical outliers hierachical clustering analysis is performed based on the Euclidean 
                distance between the sample based on all miRNAs across all Samples.</p>", 
                "<p> A total of ",myvalues$sigNum," miRNAs passed significance treshold for padj &#8804; ",input$p.cut, ", of these ", myvalues$sigNumup, 
                " were upregulated and ", myvalues$sigNumdwn, " were downregulated in ", input$treatment, " Samples. 
                At nominal threshold of uncorrected p-value &#8804; ",input$p.cut, ", a total of ",myvalues$nomsigNum," miRNas were significant 
                with ",myvalues$nomsigNumup," miRNAs higher expressed in ", input$treatment, " and ", 
                myvalues$nomsigNumdwn, " lower expressed in ", input$treatment, " samples.
                The Network analyses and figures are based on miRNAs differentially regulated with a significance of ", 
                myvalues$GraphAnal," &#8804; ", input$p.cut, ". </p>
                
                <h4>Methods</h4>
                
                <p>Differential expression analysis was performed in ", myvalues$Rversion, " (https://cran.r-project.org/).<br>
                Coverage files were converted into raw counts matrices using the Quiagen pipeline. <br>
                After NGS analysis we obtained ~ ", myvalues$NtotReads, " Mio reads 
                with an average read counts per sample of ", myvalues$AVNreads, " reads per sample. <br>
                Overall, we detected ",myvalues$Nmirna," miRNAs of which ",myvalues$Nmirna_gt_5, " had 5 or more reads in each sample and ", 
                myvalues$Nmirna_gt_5_50, " miRNAs were detected with 5 or more reads in at least 50% of the samples. A-priori 
                filtering for sparse read counts (min ",input$mincounts, " in one Sample, and min 1 read in ",round(input$minsamples*100,0),"% of all samples) has been applied. We performed hierarchical cluster analysis of raw count data based on the Euclidean distance between samples and the 
                Ward algorithm implemented in R using the normalized Read counts (function counts (Reads, normalize=T)) miRNAs.<br>
                For comparisons between two groups we loaded the respective raw read counts into DESeq2 using the 'DESeqDataSetFromMatrix' function. 
                Differential expression was estimated using the function 'DESeq' with fit-typ = 'local' based on the developer's recommendation. <br>
                Statistical models were corrected for ", myvalues$Confoundertext, "<br>
                Fdr correction (Benjamini Hochberg) was applied for each miRNA passing DESeq-quality thresholds. 
                miRNAs were considered to be differentially expressed with padj &#8804; ",input$p.cut, " 
                and not-differentially expressed with padj > 0.1. <br>
                Network analysis was performed using the igraph package in R assuming an undirected network based on the 
                correlation between miRNAs selected based on group differecens with ", myvalues$GraphAnal ," &#8804; ",input$p.cut, " 
                an aribtrary correlation treshold of ", input$cor.cut, " to identify subnetworks. 
                Sub-clusters were identified based on the edge betweenness using the cluster_edge_betweenness() function in igraph. </p>
                <p>Functional enrichment analysis was performed by testing the union of genes targeted by the differentially regulated miRNAs 
                using the gprofiler2 package, correcting for the multiple testing and the diacyclic stucture of the database (gSCS option). <br> 
                As targets we only considered genes predicted to be tageted by the respective miRNA in at least three out of five databases  (DIANA, Targetscan, PicTar, Miranda, and miRDB).
                The top 10 genes, ranked by the geometric mean of the normalized database scores were included. Target prediction and ranking was perfomed using the miRNAtap package in R and only miRNAs were tested.  
                </p>")
            )
        }))
        print("covariate descriptives")
        output$covariate <- renderUI(withProgress(message = "generating ...", {
            
            HTML(export2md(myvalues$covariateComp, format = "html"))
        }))
        
        print("ddsplot")
        output$dds_plot <- renderPlot(height = 700,
                                      withProgress(message = "generating ...",{
                                          req(myvalues$dds)
                                          par(mfrow=c(2,2), mar=c(5,5,3,5))
                                          dds_data = myvalues[["dds"]]
                                          plotDispEsts(dds_data)
                                          plotMA(dds_data)
                                          
                                          RawplotDat=log(myvalues$Reads+1)
                                          NormplotDat=log(myvalues$NReads+1)
                                          boxplot(RawplotDat, col=AGCcol[1+as.numeric(myvalues$meta$group)], main="raw counts", 
                                                  ylab="log(counts+1)", las=3, cex.axis=.7)
                                          legend(xpd=2, x=ncol(RawplotDat)+1, y=max(RawplotDat, na.rm=T), legend=levels(myvalues$meta$group), col=AGCcol[2:3], 
                                                 pch=15, bty="n")
                                          
                                          boxplot(NormplotDat, col=AGCcol[1+as.numeric(myvalues$meta$group)], main="normalized counts", 
                                                  ylab="log(counts+1)", las=3, cex.axis=.7)
                                          legend(xpd=2, x=ncol(NormplotDat)+1, y=max(NormplotDat, na.rm=T), legend=levels(myvalues$meta$group), col=AGCcol[2:3], 
                                                 pch=15, bty="n")
                                      }))
        
        print("DEG heatmaps")
        output$DEG_heatmaps_miRNAs <- renderPlot(
            withProgress(message = "generating ...",{
                
                gradient_base <- colorRampPalette(c("dodgerblue", "black"))(100)
                RStmp=rowSums(myvalues$NReads[colnames(myvalues$miRCo),],na.rm=T)
                RStmp_cal=100 * (RStmp-min(RStmp, na.rm=T))/max(RStmp,RStmp, na.rm=T)
                cc=gradient_base[floor(RStmp_cal)+1]
                
                heatmap.2(myvalues$miRCo, 
                          main="correlation of significant miRNAs",
                          col = viridis(100), 
                          trace="n", 
                          density.info ="n", 
                          ColSideColors = cc,
                          mar=c(8,8))
                
                legend("topright",title = "rel expr.", 
                       legend = c(round(min(RStmp,RStmp, na.rm=T),1),
                                  round(max(RStmp,RStmp, na.rm=T),1)), 
                       pch=15, col = c("dodgerblue", "black"),
                       cex = 0.7, xpd=T, bty="n")
                
            })) 
        
        print("DEG heatmaps_Sample")
        output$DEG_heatmaps_Samples <- renderPlot(
            withProgress(message = "generating ...", 
                         {
                             heatmap.2(myvalues$NReads[myvalues$targetmirs,], 
                                       main="Sample corr of miRNAs",
                                       col = viridis(100), 
                                       trace="n", 
                                       scale = "row",
                                       density.info ="n",
                                       ColSideColors = AGCcol[1+as.numeric(myvalues$meta$group)], 
                                       mar=c(8,8))
                             
                             legend("topright", 
                                    legend = sort(unique(myvalues$meta$group)), 
                                    pch=15, col = AGCcol[1+c(1:length(unique(myvalues$meta$group)))],
                                    cex = 0.7, xpd=T, bty="n")
                         })
        )
        
        print("DEG_Corr_hist")
        output$DEG_Corr_hist <- renderPlot(
            withProgress(message = "generating ...", 
                         hist(myvalues$miRCo, col=AGCcol[4], border = "#444545", 
                              main="corr between significant miRNAs")))
        print("Networkplot")
        output$Network <- renderVisNetwork(
            withProgress(message = "generating ...", 
                         visIgraph(myvalues$graph)))
        
        print("heatmap_samples")
        output$heatmap_samples <- renderPlot(
            withProgress(message ="generating...",
                         {
                             metaplot<- colData(myvalues[["dds"]])
                             metaplot <- metaplot[,c("group",input$covariates)] %>% as.data.frame()
                             metaplot <- metaplot %>% mutate_if(is.logical, as.character)%>% as.data.frame()
                             metaplot <- metaplot %>% mutate_if(is.character, as.factor)%>% as.data.frame()
                             DistM<-as.matrix(dist(t(log(myvalues$NReads+1)), 
                                                   method = "euclidean"))
                             colors<-viridis(255)
                             rownames(metaplot) <- colnames(DistM)
                             
                             heatmap.2(as.matrix(DistM), col = viridis(100),mar=c(8,8), 
                                       main="dissimilarity all miRNAs", density.info = "none", trace = "none", 
                                       hclustfun =function(x){hclust(x,method="ward.D2")},
                                       ColSideColors=AGCcol[1+as.numeric(myvalues$meta$group)])
                             
                             legend("topright", 
                                    legend = sort(unique(myvalues$meta$group)), 
                                    pch=15, col = AGCcol[1+c(1:length(unique(myvalues$meta$group)))],
                                    cex = 0.7, xpd=T, bty="n")
                             
                             
                         }))
        
        print("table_display_fc")
        output$table_display_fc <- renderDataTable(
            withProgress(message = "generating ...", 
                         {   s=input$volcano_plot_click
                         print(s)
                         dataset = as.data.frame(myvalues$restab_DEseq)
                         #v1 <- 1:nrow(row.names(dataset)
                         #cols1 <- ifelse(v1 %in% s,'orange','')
                         
                         datatable(dataset, 
                                   extensions = "Buttons",
                                   options = list(
                                       autoWidth = TRUE,
                                       pageLength = 100,
                                       info = FALSE,
                                       autoWidth=T, 
                                       lengthMenu = list(c(-1, 15,50, 100), 
                                                         c("All", "15","50", "100" )),
                                       dom = 'Blfrtip',
                                       buttons = c('copy', 'csv', 'excel', 'pdf'))
                         )
                         
                         }), server = F)
        
        
        print("volcano_plot")
        output$volcano_plot <- renderPlotly(
            withProgress(message ="generating...",
                         { s=input$table_display_fc_rows_selected
                         ResultsfinCaCo <- data.frame(myvalues$restab_DEseq)
                         par(mar=c(5,5,5,8))
                         
                         hline <- function(y = 0, color = "black") {
                             list(
                                 type = "line",
                                 x0 = 0,
                                 x1 = 1,
                                 xref = "paper",
                                 y0 = y,
                                 y1 = y,
                                 line = list(color = color)
                             )
                         }
                         fig <- plot_ly(ResultsfinCaCo)
                         
                         fig <- fig %>%
                             add_trace(               
                                 x=~log2FoldChange, 
                                 y=~-log10(pvalue),
                                 colors = c("black", "red"),
                                 mode="markers",
                                 type="scatter",
                                 alpha = 0.7,
                                 text = ~paste(rownames(ResultsfinCaCo),"<br>",
                                               "log2FC: ", log2FoldChange,"<br>",
                                               "padj: ", formatC(padj, format="e")),
                                 hoverinfo = 'text', showlegend = T, 
                                 color=~as.factor(padj<input$p.cut)) %>% 
                             
                             layout(legend = list(title=list(text=paste("<b> padj <",input$p.cut,"</b>"))),
                                    shapes = list(hline(-log10(input$p.cut)))) %>% 
                             
                             add_text(showlegend = FALSE, x=~min(log2FoldChange), 
                                      y=-log10(input$p.cut), textposition = "top right",
                                      text = paste("nominal p", input$p.cut), 
                                      textfont  = list(color = '#000'))
                         
                         fig <- fig %>% layout(legend = list(x=0, y=1.1, orientation = 'h'))
                         
                         if(length(s)){
                             fig <- fig %>% add_trace(x=~log2FoldChange[s], 
                                                      y=~-log10(pvalue[s]),
                                                      mode="markers", 
                                                      type="scatter", 
                                                      marker = list(
                                                          color = alpha("white", 0),
                                                          size = 10,
                                                          line = list(
                                                              color = 'rgb(231, 99, 250)',
                                                              width = 2
                                                          )), showlegend=F)
                         }
                         fig
                         
                         }))
        
        output$boxplot = renderPlot(
            withProgress(message ="generating...",
                         {
                             s=input$table_display_fc_rows_selected
                             if(length(s)>0){
                                 restab=data.frame(myvalues$restab_DEseq)
                                 NReads=data.frame(myvalues$NReads)
                                 meta=myvalues$meta
                                 mirtargest=rownames(restab)[s]
                                 plotdata <- NReads[mirtargest, ] %>% t() %>% as.data.frame()
                                 plotdata$Target = c(meta[rownames(plotdata),input$target])
                                 plotdata <- plotdata %>% reshape2::melt(value.name = "Nreads",variable.name="miRNA",
                                                                         measure.vars=mirtargest)
                                 ggplot(plotdata, aes(x=Target, y=log2(Nreads+1), fill=Target)) +
                                     geom_boxplot() + facet_wrap(~miRNA, scales="free") + theme_bw()
                                 }else{
                                 empty_plot()
                             }
                             
                         }
        ))

        
        print("GO_plot_all")
        output$GO_plot_all <- renderPlotly(
            withProgress(message ="generating...",
                         {
                             if(length(myvalues$ResGO_all)>0){
                                 p=gostplot(myvalues$ResGO_all)
                             } else{
                                 p = empty_plot("no significant enrichment identified")
                             }
                             p
                         }
            ))
        
        output$table_GO_all <- renderDataTable(
            withProgress(message = "generating ...", 
                         {
                             datatable(as.data.frame(myvalues$ResGO_all$result), 
                                       extensions = "Buttons",
                                       options = list(
                                           autoWidth = TRUE,
                                           pageLength = 100,
                                           info = FALSE,
                                           autoWidth=T, 
                                           lengthMenu = list(c(-1, 15,50, 100), 
                                                             c("All", "15","50", "100" )),
                                           dom = 'Blfrtip',
                                           buttons = c('copy', 'csv', 'excel', 'pdf')
                                       )
                             )
                         }))
        print("GO_plot_up")
        output$GO_plot_up <-renderPlotly(
            withProgress(message ="generating...",
                         {
                             if(length(myvalues$ResGO_up)>0){
                                 p=gostplot(myvalues$ResGO_up)
                             } else{
                                 p = empty_plot("no significant enrichment identified")
                             }
                             p
                         }
            ))
        
        
        
        output$table_GO_up <- renderDataTable(
            withProgress(message = "generating ...", 
                         {
                             datatable(as.data.frame(myvalues$ResGO_up$result), 
                                       extensions = "Buttons",
                                       options = list(
                                           autoWidth = TRUE,
                                           pageLength = 100,
                                           info = FALSE,
                                           autoWidth=T, 
                                           lengthMenu = list(c(-1, 15,50, 100), 
                                                             c("All", "15","50", "100" )),
                                           dom = 'Blfrtip',
                                           buttons = c('copy', 'csv', 'excel', 'pdf')
                                       )
                             )
                         }))
        print("GO_plot_dwn")
        output$GO_plot_dwn <-renderPlotly(
            withProgress(message ="generating...",
                         {
                             if(length(myvalues$ResGO_dwn)>0){
                                 p=gostplot(myvalues$ResGO_dwn)
                             } else{
                                 p = empty_plot("no significant enrichment identified")
                             }
                             p
                         }
            ))
        
        output$table_GO_dwn <- renderDataTable(
            withProgress(message = "generating ...", 
                         {
                             datatable(as.data.frame(myvalues$ResGO_dwn$result), 
                                       extensions = "Buttons",
                                       options = list(
                                           autoWidth = TRUE,
                                           pageLength = 100,
                                           info = FALSE,
                                           autoWidth=T, 
                                           lengthMenu = list(c(-1, 15,50, 100), 
                                                             c("All", "15","50", "100" )),
                                           dom = 'Blfrtip',
                                           buttons = c('copy', 'csv', 'excel', 'pdf')
                                       )
                             )
                         }))
        print("done")
        
        
        
        
    })
}


# Run the application 
shinyApp(ui = ui, server = server)
