#written in R-3.3.0, USE R-3.3.0!!!!
# Install these one at a time, shiny takes like a year
# Had issues installing shiny 08/06/18, added the binary and it worked just fine
#install.packages("shiny", type = "binary")
if (!"shiny" %in% rownames(installed.packages())){install.packages('shiny',dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"rhandsontable" %in% rownames(installed.packages())){install.packages('rhandsontable',dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"stringr" %in% rownames(installed.packages())){install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"ShortRead" %in% rownames(installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ShortRead",suppressUpdates=TRUE)}
if (!"shinyFiles"%in%rownames(installed.packages())){install.packages("shinyFiles",dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"shinythemes"%in%rownames(installed.packages())){install.packages("shinythemes",dependencies=TRUE, repos='http://cran.rstudio.com/')}
if (!"stringdist"%in%rownames(installed.packages())){install.packages("stringdist",dependencies=TRUE, repos='http://cran.rstudio.com/')}

library(shinythemes)
library(shiny)
library(ShortRead)
library(stringr)
library(rhandsontable)
library(shinyFiles)
library(stringdist)
rm(list=ls())

ui=navbarPage(theme=shinytheme("sandstone"),
              title=h1("Barcoded Virus Analysis Tool",style="line-height: 0;font-size:100;font-weight:bold"),
              windowTitle="Barcoded Virus Analysis Tool",
              tabPanel(h1("1: Choose Files",style="font-size:small"),  #establish tab 1
                       mainPanel(  
                         shinyFilesButton('runfile', 'Run File Select', 'Please select Run file', FALSE),
                         tableOutput("runfilepath"),
                         shinyFilesButton('reffile', "Stock Barcode Sequences Select(OPTIONAL)", 'Please select Reference file(OPTIONAL)', FALSE),
                         tableOutput("reffilepath")
                       )
              ),
              
              tabPanel(h1("2: Specify 5' Indexing Tag and Input Values",style="font-size:small"), #establish tab 2
                       mainPanel( 
                         rHandsontableOutput("barcodelist") #sets up the hands on table for barcode inputs
                       )
              ),
              tabPanel(h1("3: Specify Reference, Read Orientation, and Output Preferences",style="font-size:small"),
                       sidebarLayout(
                         sidebarPanel( #add all of the extra input information
                           textInput("Adapter","Add adapter sequence(leave blank if none, use Ns if random)",value="NNNN"), #adapter 
                           radioButtons("fwdrev","Is Sequence Forward or Reverse Complemented?",choices=c("Forward","Reversed"),selected="Reversed"),
                           textInput("ref","Add Reference Sequence",value="ATGGAAGAAAGACCTCCAGAAAATGAAG"),
                           numericInput("maxmismatch","Max # of Mismatches to Reference",value=2,min=0),
                           radioButtons("direction","Extract Upstream or Downstream of Reference?",choices=c("Upstream","Downstream")),
                           numericInput("numbases","Number of Bases to Extract",value=34,min=1),
                           radioButtons("cutoff","Cutoff Based on Input Value?",choices=c("Yes","No"),selected="Yes"),
                           numericInput("numparse","Maximum Number of Sequences to Process at a Time (Higher is Faster but Takes More RAM)",value=5*10^5,min=10^4),
                           numericInput("mindist","Flag If Hamming Distance to Another Sequence Less Than or Equal To (0 for no analysis)",value=1,
                                        min=0),
                           numericInput("contam","Percentage Threshold to be Considered Contaminated (100 for no analysis)",value=0.2,min=0,max=100)
                         ),
                         mainPanel(
                           plotOutput("refplot") #put the plot of text as a picture on the side
                         )   
                       )),
              tabPanel(h1("4: Run Analysis + Results",style="font-size:small"),
                       sidebarLayout(
                         sidebarPanel(
                           tableOutput("summary"),
                           actionButton("button","Start Analysis",style='background-color:#f44336'), #start button
                           span(textOutput("runfilemessage"),style="color:red"),
                           span(textOutput("reffilemessage"),style="color:red"),
                           span(textOutput("BarcodeMessage"),style="color:red"),
                           span(textOutput("ready"),style="color:white")
                         ),
                         mainPanel(
                           conditionalPanel(condition = "output.ready == 'ready'",
                                            downloadButton("downloadbarcodedata", "Download Barcode Analysis Data"), #download 1
                                            downloadButton("sequencestats", "Download Sequence Statistics"), #download 2
                                            tableOutput("check3"))
                         )
                       ))
)

server=function(input,output){
  output$ready=renderText({""})
  volumes <- c(User=str_sub(path.expand("~"),1,tail(unlist(gregexpr("/",path.expand("~"))),n=1)),getVolumes()())
  shinyFileChoose(input,update=2000, 'runfile', roots=volumes,filetypes=c('fastq','gz'))
  output$runfilepath=renderTable({parseFilePaths(volumes,input$runfile)[,1]})
  rememberpath<-eventReactive(input$runfile,{
    a=as.character(parseFilePaths(volumes,input$runfile)$datapath);
    occ=unlist(gregexpr("/",a))
    if (length(occ)==0){}
    else {occ=unlist(gregexpr("/",a))
    occ=occ[length(occ)]
    path=str_sub(a,1,occ-1)
    return(path)}
    
  })
  volumes2=reactive({c(LastPath=rememberpath(),volumes)})
  observe({shinyFileChoose(input, update=2000, 'reffile', roots=volumes2(),filetypes=c('fasta','gz'))})
  output$reffilepath=renderTable({parseFilePaths(volumes2(),input$reffile)[,1:2]})
  primer<-paste("P5.",1:400,sep="") #writes primer names, doesn't really matter if its VPX, INT, PBS, etc. 
  barcodeseqs<-c("TAAGGCGA","CGTACTAG","AGGCAGAA","TCCTGAGC","GGACTCCT","TAGGCATG","CTCTCTAC","CAGAGAGG","GCTACGCT","CGAGGCTG",
                 "AAGAGGCA","GTAGAGGA","TAGATCGC","CTCTCTAT","TATCCTCT","AGAGTAGA","GTAAGGAG","ACTGCATA","AAGGAGTA","CTAAGCCT",
                 "TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA","AGGCGAAG","TAATCTTA","CAGGACGT","GTACTGAC","ATTACTCG",'TCCGGAGA',
                 "CGCTCATT","GAGATTCC","ATTCAGAA","GAATTCGT","CTGAAGCT","TAATGCGC","CGGCTATG","TCCGCGAA",'TCTCGCGC','AGCGATAG',
                 rep("",360)) #40 multiplexed barcodes
  Include=rep(FALSE,400) #row 3
  Input=rep(0,400) #row 4
  barcode.list<-data.frame(primer,barcodeseqs,Include,Input) #this holds the primer number, sequence, inputs, and whether to keep 
  barcode.list[,1]=as.character(barcode.list[,1])
  barcode.list[,2]=as.character(barcode.list[,2])
  colnames(barcode.list)=c("5' Index Name","Index Sequence (blank if none)","Include?","Sequencer Template Input(0 if Unknown or Unused)") #give column names
  
  DF=reactive({if (is.null(input$barcodelist)) { 
    DF2 = barcode.list #if the barcode table is unchanged, the table remains what it was
  } else {
    DF2 = hot_to_r(input$barcodelist) #if the table changes, the table becomes updated
  }
    DF2
  })
  
  output$barcodelist <- renderRHandsontable({ #renders the actual table
    rhandsontable(DF())
  })
  output$summary=renderTable({data.frame(matrix(c("# of 5' Indexing Tags","Sequence Orientation",
                                                  "Cutoff?","# of Mismatches",
                                                  sum(DF()[,3]),input$fwdrev,input$cutoff,input$maxmismatch),
                                                nrow=4,ncol=2,dimnames=list(NULL,c("Option","Selected"))))})
  runfilemessage=reactive({if (nrow(parseFilePaths(volumes,input$runfile))==0){"Please Select Run File"}})
  reffilemessage=reactive({if (nrow(parseFilePaths(volumes,input$reffile))==0){"Please Select Reference File"}})
  BarcodeMessage=reactive({if (sum(DF()[,3])==0){"Please Include Barcodes"}})
  output$runfilemessage=renderText({runfilemessage()})
  output$reffilemessage=renderText({reffilemessage()})
  output$BarcodeMessage=renderText({BarcodeMessage()})
  
  ###################ignore this, this is to make the diagram on panel 3, totally irrelevant for functioning#####
  output$refplot=renderPlot({ #IGNORE THIS, THIS IS THE DIAGRAM ON PANEL 3#
    
    if (input$direction == "Upstream"){
      plot.new();
      text(0.5,0.8,input$ref,font=2,col="red",cex=1.8*15/nchar(input$ref));
      arrows(0.2,0.75,0.01,0.75,col="blue",font=2,lwd=4)
      text(0.1,0.7,paste("-",input$numbases," bases",sep=""),font=2,cex=1.5)
      rect(0.05,0.1,0.15,0.2,col="NA",border="red")
      text(0.1,0.15,input$Adapter,col="red",cex=1.3*4/nchar(input$Adapter),font=2)
      rect(0.16,0.1,0.4,0.2,col="NA",border="blue")
      text(0.28,0.15,"Index Sequence",col="blue",cex=1.4,font=2)
      rect(0.41,0.1,0.9,0.2,col="NA",border="darkgreen")
      text(0.65,0.15,"Sequence",col="darkgreen",cex=1.6,font=2)
      if (input$fwdrev=="Reversed"){text(0.9,0.07,"5'",font=2,cex=1.4)
        text(0.41,0.07,"3'",font=2,cex=1.4)
        arrows(0.8,0.05,0.51,0.05,lwd=4)}
      if (input$fwdrev=="Forward"){
        text(0.9,0.07,"3'",font=2,cex=1.4)
        text(0.41,0.07,"5'",font=2,cex=1.4)
        arrows(0.51,0.05,0.8,0.05,lwd=4)}
    }
    if (input$direction == "Downstream"){
      plot.new();
      text(0.5,0.8,input$ref,font=2,col="red",cex=1.8*15/nchar(input$ref));
      arrows(0.8,0.75,0.99,0.75,cex=1.7,col="blue",font=2,lwd=4)
      text(0.9,0.7,paste("+",input$numbases," bases",sep=""),font=2,cex=1.5)
      rect(0.05,0.1,0.15,0.2,col="NA",border="red")
      text(0.1,0.15,input$Adapter,col="red",cex=1.3*4/nchar(input$Adapter),font=2)
      rect(0.16,0.1,0.4,0.2,col="NA",border="blue")
      text(0.28,0.15,"Index Sequence",col="blue",cex=1.4,font=2)
      rect(0.41,0.1,0.9,0.2,col="NA",border="darkgreen")
      text(0.65,0.15,"Sequence",col="darkgreen",cex=1.6,font=2)
      if (input$fwdrev=="Reversed"){text(0.9,0.07,"5'",font=2,cex=1.4)
        text(0.41,0.07,"3'",font=2,cex=1.4)
        arrows(0.8,0.05,0.51,0.05,lwd=4)}
      if (input$fwdrev=="Forward"){
        text(0.9,0.07,"3'",font=2,cex=1.4)
        text(0.41,0.07,"5'",font=2,cex=1.4)
        arrows(0.51,0.05,0.8,0.05,lwd=4)}
      
    }
  })
  ####################################
  
  
  observeEvent(input$button,{ #if the start button is pressed:
    withProgress(message="Reading Run File",{
      barcode.list.cleaned=DF() #input list of barcodes
      path.to.file=as.character(parseFilePaths(volumes,input$runfile)$datapath) #input file path
      adapter=as.character(input$Adapter) #input adapter sequence
      reffilepath=as.character(parseFilePaths(volumes2(),input$reffile)$datapath) #input reference file
      fwdrev=input$fwdrev #specify if reads are fwd or rev
      ref=input$ref #specify the reference sequence
      mismatches=as.numeric(input$maxmismatch) #specify maximum mismatches to reference
      direction=input$direction #specify if to cut upstream or downstream
      numbases=input$numbases #specify size of cut
      cutoff=input$cutoff #specify whether to use the cutoff
      parse=input$numparse
      if (is.na(as.numeric(parse))|(as.numeric(parse)<10^4)){parse=5*10^5}
      reffile=readFasta(reffilepath) #read the reference fasta file
      reffileseq=as.character(sread(reffile)) #take out the sequences from the fasta file
      reffilename=as.character(ShortRead::id(reffile)) #take out names of the reference file
      mindist=input$mindist
      contam=input$contam
      
      stream=FastqStreamer(path.to.file,parse)
      
      if (sum(barcode.list.cleaned[,3])>0){
        barcode.list.cleaned=barcode.list.cleaned[barcode.list.cleaned[,3],]}#only use the specified barcodes
      if (sum(barcode.list.cleaned[,3])==0){barcode.list.cleaned=barcode.list.cleaned[1,];barcode.list.cleaned[1,2]=""} 
      #if no barcodes were specified, then 
      barcode.list.cleaned[,2]=as.character(barcode.list.cleaned[,2]) #convert barcode sequences into a character vector
      barcode.list.cleaned[barcode.list.cleaned[,4]==0,4]=Inf #if inputs are not given, we assume infinite input
      barcodeunique.hold<-as.data.frame(matrix(NA,nrow=1*10^4,ncol=(nrow(barcode.list.cleaned)*7))) #create a big holder matrix
      barcodeunique.hold[2,]=rep(c("Sequencer Input","# Barcode Sequences","Number of Barcodes","Proportion Uniques","Cross Contamination Percent Threshold",NA,NA),ncol(barcodeunique.hold)/7) 
      #second row of the big holder matrix has miseq input, seq output, # barcodes, proportion uniques
      barcodeunique.hold[5,]=rep(c("Barcode","Sequence","Counts","Proportion","Hamming Dist to Another Sequence","Contaminated?",NA),ncol(barcodeunique.hold)/7)
      #5th row of matrix has barcode, sequence, counts, and proportions, every 5 columns
      
      sequencestats=as.data.frame(matrix(0,nrow=nrow(barcode.list.cleaned),ncol=6)) #another matrix holding sequence stats
      colnames(sequencestats)=c("Primer","# of 5' Indexing Sequences","Total Sequences Extracted per Barcode",
                                "# of Sequences Matching Known Barcode","# of Barcodes","Cutoff Based on Input?(Inf means none)")
      #it tells you stats about the run. 
      sequencestats[,1]=barcode.list.cleaned[,1] #first column of sequence stats are the primers in the run
      #incProgress(amount=0.1,message=paste("Analyzing Barcode",barcode.list.cleaned[1,1])) #we begin barcode analysis:
      incProgress(amount=0.1,message=paste("Analyzing"))
      nthreads <- .Call(ShortRead:::.set_omp_threads, 12)
      on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))
      sum=0
      while(length(fq<-yield(stream))){
        file.seq=as.character(sread(fq))
        if (sum(str_detect(adapter,LETTERS[-14]))==0){ #if the adapter is only Ns
          file.seq=str_sub(file.seq,nchar(adapter)+1,nchar(file.seq)) #take out the adapter
        }
        if (sum(str_detect(adapter,LETTERS[-14]))>0){ #if the adapter has other letters 
          check.adapter=str_sub(file.seq,1,nchar(adapter)) #remove the adapter from all sequences
          file.seq=file.seq[check.adapter==adapter] #keep adapters only if the adapter matches the input adapter
          file.seq=str_sub(file.seq,nchar(adapter)+1,nchar(file.seq)) #remove the adapter
        }
        #read fastq file
        gc() #reallocate memory (fastq files have file name, sequence, and quality. We're giving back memory for name and quality)
        
        for (i in 1:nrow(barcode.list.cleaned)){ #loop over each barcode
          file.seq.4bpTrim=str_sub(file.seq,0,nchar(barcode.list.cleaned[i,2])) #cut out the length of the barcode
          index=file.seq.4bpTrim==barcode.list.cleaned[i,2] #see if the barcode matches the first reference 
          if ((sum(index)<(0.00001*length(file.seq.4bpTrim)/nrow(barcode.list.cleaned)))&(barcode.list.cleaned[i,4]==Inf)){
            ;next} #if there are fewer than .01% barcodes, move to the next barcode
          if (sum(index)==0 & barcode.list.cleaned[i,4]>0){next}
          rm(file.seq.4bpTrim);gc() #memory allocation
          file.seq.barcodesplit=str_sub(file.seq[index],nchar(barcode.list.cleaned[i,2])+1,nchar(file.seq[index]))
          #cut the sequences that have the barcode from right up the barcode to the end
          if (fwdrev=="Reversed"){ #if the file is reverse complemented, then reverse complement
            #incProgress(amount=0,detail="Reverse Complementing") #progress meter
            file.seq.barcodesplit=as.character(reverseComplement(DNAStringSet(file.seq.barcodesplit)))}
          
          if (direction=="Upstream"){index.downstream.num=aregexec(ref,file.seq.barcodesplit,max.distance=mismatches) 
          #uses aregexec to match to reference
          index.downstream.num=unlist(index.downstream.num)
          #lists where each match's start was
          file.seq.barcodesplit.noNA=file.seq.barcodesplit[index.downstream.num>=numbases+1]
          #keep only sequences where the match is found more than 34 bases from the beginning
          index.downstream.num.noNA=index.downstream.num[index.downstream.num>=numbases+1]
          #keep only indices that are found more than 34 bases from the beginning
          # 8/27/18, changing cutoff from 2 to 6 - Abbey Thorpe
          if (length(index.downstream.num.noNA)<6){next}
          seqs=str_sub(file.seq.barcodesplit.noNA,index.downstream.num.noNA-numbases,index.downstream.num.noNA-1) 
          #extracts sequences from 34 bp to 1 bp before the downstream region
          seqs=as.data.frame(DNAStringSet(seqs))}
          #convert the sequences to a data frame
          
          if (direction=="Downstream"){
            ref=as.character(reverseComplement(DNAString(ref)))
            file.seq.barcodesplit=as.character(reverseComplement(DNAStringSet(file.seq.barcodesplit)))
            index.downstream.num=aregexec(ref,file.seq.barcodesplit,max.distance=mismatches) 
            #uses aregexec to match to reference
            index.downstream.num=unlist(index.downstream.num)
            #lists where each match's start was
            file.seq.barcodesplit.noNA=file.seq.barcodesplit[index.downstream.num>=numbases+1]
            #keep only sequences where the match is found more than 34 bases from the beginning
            index.downstream.num.noNA=index.downstream.num[index.downstream.num>=numbases+1]
            #keep only indices that are found more than 34 bases from the beginning
            if (length(index.downstream.num.noNA)<2){next}
            seqs=str_sub(file.seq.barcodesplit.noNA,index.downstream.num.noNA-numbases,index.downstream.num.noNA-1) 
            #extracts sequences from 34 bp to 1 bp before the downstream region
            seqs=as.data.frame(reverseComplement(DNAStringSet(seqs)))}
          #convert the sequences to a data frame
          
          
          
          #sequence stats will record how many sequences there are
          seqs.Table=table(seqs) #essentially, finds duplicates
          rm(seqs); gc() #memory allocation
          seqs.Count=data.frame(seqs.Table) #converts duplicate finds to data frame format.
          rm(seqs.Table) #remove a useless object
          if (ncol(seqs.Count)<2){next}
          seqs.Count[,1]=as.character(seqs.Count[,1]) #convert first column of sequences to character
          seqs.Count[,2]=as.numeric(as.character(seqs.Count[,2])) #convert second column of numbers to numbers
          alreadyleft=barcodeunique.hold[6:(nrow(barcodeunique.hold)),(7*i-5):(7*i-4)]
          alreadyleft=alreadyleft[!is.na(alreadyleft[,1]),]
          alreadyleft[,2]=as.numeric(alreadyleft[,2])
          newSeqs.Count=data.frame(tapply(c(seqs.Count[,2],alreadyleft[,2]),c(seqs.Count[,1],alreadyleft[,1]),FUN = sum))
          newSeqs.Count[,2]=row.names(newSeqs.Count)
          newSeqs.Count=newSeqs.Count[,c(2,1)]
          seqs.Count=newSeqs.Count
          rm(newSeqs.Count)
          gc()
          #if the user said to only show barcode, keep those that're found in the reference file
          if (nrow(seqs.Count)<2){next} #if no sequences have vpr, move on
          
          ordering=seqs.Count[,1]%in%reffileseq
          sequencestats[i,4]=sum(ordering*seqs.Count[,2])
          seqs.Count[,3]=seqs.Count[,2]/(ordering*sum(seqs.Count[,2]))
          seqs.Count[,4]=""
          seqs.Count[ordering,4]=reffilename[match(seqs.Count[ordering,1],reffileseq)]
          seqs.Count=seqs.Count[,c(4,1,2,3)]
          seqs.Count=seqs.Count[order(-ordering,-seqs.Count[,3]),]
          seqs.Count[((sum(ordering)+1):(nrow(seqs.Count))),4]=0
          
          seqs.Count[((sum(ordering)+1):(length(seqs.Count[,1]))),1]=paste("Unique",1:sum(!ordering))
          
          if (sum(ordering)>0){
            seqs.Count[1:sum(ordering),4]=seqs.Count[1:sum(ordering),3]/sum(seqs.Count[1:sum(ordering),3])}
          colnames(seqs.Count)=c("Barcode","Sequence","Counts","Proportion")
          print(i)
          sequencestats[i,2]=sum(sum(index),sequencestats[i,2])
          if (nrow(seqs.Count)<2){next}
          if (length(seqs.Count)<2){next}
          barcodeunique.hold[6:(nrow(seqs.Count)+5),(7*i-6):(7*i-3)]=seqs.Count
          barcodeunique.hold[1,(7*i-6):(7*i-1)]=rep(barcode.list.cleaned[i,1],6)
        }
        sum=sum+length(fq)
        incProgress(amount=0,message=paste(sum/10^6,"Million Sequences Analyzed"))
      }
      
      incProgress(amount=0.7,message="Performing Post Processing Clean-Up")
      
      
      if ((contam<100)&(nrow(barcode.list.cleaned)>1)){data=as.data.frame(matrix(NA,nrow=10^4,ncol=length(barcode.list.cleaned[,2])*2))
      colnames(data)=rep(barcode.list.cleaned[,1],each=2)
      }
      
      for (i in 1:nrow(barcode.list.cleaned)){
        cur=barcodeunique.hold[6:nrow(barcodeunique.hold),(7*i-6):(7*i-4)]
        if (is.na(cur[1,1])){next}
        cur=cur[!is.na(cur[,3]),]
        cur[,3]=as.numeric(cur[,3])
        cur=cur[!is.na(cur[,3]),]
        cur=cur[cur[,3]>1,]
        if (length(cur[,3])==0){next}
        if (cutoff=="Yes"){
          if (all(str_sub(cur[,1],1,6)=="Unique")){cur[,4]=cur[,3]/sum(cur[,3])}
          if (length(cur[!((str_sub(cur[,1],1,6)=="Unique")),3])==0){}
          if (length(cur[!((str_sub(cur[,1],1,6)=="Unique")),3])>0){
            cur=cur[cur[!((str_sub(cur[,1],1,6)=="Unique")),3]/sum(cur[!((str_sub(cur[,1],1,6)=="Unique")),3])>=1/barcode.list.cleaned[i,4],]
            cur[!((str_sub(cur[,1],1,6)=="Unique")),4]=cur[!((str_sub(cur[,1],1,6)=="Unique")),3]/sum(cur[!((str_sub(cur[,1],1,6)=="Unique")),3])
            #if you choose to cutoff, then only sequences whose proportion is greater than 1/input is kept
          }}
        if ((cutoff=="No")){
          if (all(str_sub(cur[,1],1,6)=="Unique")){cur[,4]=cur[,3]/sum(cur[,3])}
          if (length(cur[!((str_sub(cur[,1],1,6)=="Unique")),3])>0){
            cur[!((str_sub(cur[,1],1,6)=="Unique")),4]=cur[!((str_sub(cur[,1],1,6)=="Unique")),3]/sum(cur[!((str_sub(cur[,1],1,6)=="Unique")),3])
          }
        }
        compare=cur[str_sub(cur[,1],1,6)!="Unique",]
        if ((nrow(compare)>2)&(mindist>0)){
          distmat=sapply(compare[,2],function(x){a=stringdist(x,compare[,2],method="hamming");a[a==0]=700;a[a>min(a)]==700;
          return(paste(a[which(a<=mindist)[1]],which(a<=mindist)[1],sep="@#"))})
          comp=as.data.frame(str_split_fixed(distmat,"@#",n=2))
          distmat[suppressWarnings(as.numeric(as.character(comp[,2]))>(1:length(comp[,2])))]="NA@#NA"
          comp=as.data.frame(str_split_fixed(distmat,"@#",n=2))
          comp[,2]=suppressWarnings(as.numeric(as.character(comp[,2])))
          
          newstuff=compare[comp[!is.na(comp[,2]),2],1]
          comp[,2]=as.character(comp[,2])
          comp[!is.na(comp[,2]),2]=newstuff
          distmat=paste(comp[,1],comp[,2],sep="-")
        }
        sequencestats[i,5]=nrow(cur)-sum(str_sub(cur[,1],1,6)=="Unique")
        sequencestats[i,4]=sum((!((str_sub(cur[,1],1,6)=="Unique")))*cur[,3])
        sequencestats[i,3]=sum(cur[,3])
        sequencestats[i,6]=barcode.list.cleaned[i,4]
        barcodeunique.hold[3,(7*i-6):(7*i-3)]=c(barcode.list.cleaned[i,4],sequencestats[i,4],sequencestats[i,5],min(1,1-sequencestats[i,4]/sequencestats[i,3]))
        barcodeunique.hold[6:nrow(barcodeunique.hold),(7*i-6):(7*i-3)]=NA
        barcodeunique.hold[6:(nrow(cur)+5),(7*i-6):(7*i-3)]=cur
        if ((nrow(compare)>2)&(mindist>0)){
          barcodeunique.hold[6:(length(distmat)+5),7*i-2]=distmat
          comp[,1]=suppressWarnings(as.numeric(as.character(comp[,1])))
          comp[is.na(comp[,1]),1]=0
          barcodeunique.hold[6:(length(distmat)+5),(7*i-6):(7*i-2)]=barcodeunique.hold[6:(length(distmat)+5),(7*i-6):(7*i-2)][order(suppressWarnings(as.numeric(as.character(comp[,1]))),
                                                                                                                                    -cur[1:length(distmat),3]),
                                                                                                                              ]
          if ((contam<100)&(nrow(barcode.list.cleaned)>1)){data[1:sum(comp[,1]==0),(2*i-1):(2*i)]=cur[1:length(distmat),2:3][comp[,1]==0,];
          barcodeunique.hold[3,7*i-2]=contam}
        }
        if ((contam<100)&(nrow(barcode.list.cleaned)>1)){if (((nrow(compare)<=2)|(mindist==0))){data[1:nrow(cur),(2*i-1):(2*i)]=cur[,2:3];barcodeunique.hold[3,7*i-2]=contam}}
      }
      
      datakeep=data
      trackbetween=function(x){ #track lineage between 2(or more) timepoints/samples
        #x=subset of samples
        mat=as.data.frame(matrix(0,nrow=length(reffileseq),ncol=ncol(x)/2+2))
        mat[,1]=reffileseq
        mat[,2]=reffilename
        for (i in seq(1,ncol(x),2)){
          current=x[,(i):(i+1)]
          current=current[!is.na(current[,2]),]
          current[,2]=as.numeric(as.character(current[,2]))
          find=match(current[,1],mat[,1])
          current=current[!is.na(find),]
          find=find[!is.na(find)]
          current[,2]
          current[,3]=current[,2]/sum(current[,2])
          mat[find,(i+1)/2+2]=0
          mat[find,(i+1)/2+2]=current[,2]
        }
        mat=mat[rowSums(mat[,3:ncol(mat)])>0,]
        mat=mat[order(-rowSums(mat[,3:ncol(mat)])),]
        return(mat)
      }
      if ((contam<100)&(length(reffileseq)>0)&(nrow(barcode.list.cleaned)>1)){
        lineage=trackbetween(data)
        colnames(lineage)=c("Sequence","Name",barcode.list.cleaned[,1])
        keep=lineage
        keep2=lineage
        meandiff=rep(NA,nrow(keep))
        newhold=as.data.frame(matrix(NA,nrow=10^4,ncol=(ncol(lineage)-2)*5))
        for (i in 1:nrow(keep)){
          dats=keep[i,3:ncol(keep)]
          max=max(dats)
          dats[max/dats>(100/contam)]=NA
          keep[i,3:ncol(keep)]=dats
        }
        for (j in 3:ncol(keep)){
          now=keep[,c(1,2,j)]
          now2=keep2[,c(1,2,j)]
          now=now[now2[,3]>0,]
          now2=now2[now2[,3]>0,]
          now=now[order(-now2[,3],now2[,1]),]
          now=now[,c(2,1,3)]
          change=barcodeunique.hold[6:(nrow(now)+5),(7*(j-2)-6):(7*(j-2)-1)]
          change[is.na(now[,3]),6]="YES"
          change[!is.na(now[,3]),6]="NO"
          change=change[order(change[,6],-as.numeric(change[,3])),]
          barcodeunique.hold[6:(nrow(now)+5),(7*(j-2)-6):(7*(j-2)-1)]=change
          barcodeunique.hold[3,(7*(j-2)-5)]=as.character(sum(as.numeric(change[change[,6]=="NO",3])))
          barcodeunique.hold[3,(7*(j-2)-4)]=as.character(length(as.numeric(change[change[,6]=="NO",3])))
          
        }}
      
      barcodeunique.hold[is.na(barcodeunique.hold)]=""
      colnames(barcodeunique.hold)=barcodeunique.hold[1,]
      barcodeunique.hold=barcodeunique.hold[-1,]
      
      output$check3=renderTable({barcodeunique.hold[1:10,1:5]})
      if (cutoff=="Yes"){cutofforna="Cutoff"} 
      if (cutoff=="No"){cutofforna="NOCutoff"}
      output$downloadbarcodedata=downloadHandler(
        filename = function(){paste(as.character(parseFilePaths(volumes,input$runfile)$name),"_Barcode Analysis_",cutofforna, ".csv", sep="")},
        content = function(file) {
          write.csv(barcodeunique.hold, file,row.names=F)
        }
      )
      output$sequencestats=downloadHandler(
        filename = function(){paste(as.character(parseFilePaths(volumes,input$runfile)$name),"_Sequence Statistics", ".csv", sep="")},
        content = function(file) {
          write.csv(sequencestats, file,row.names=F,col.names=F)
        }
      )
      incProgress(amount=1,message="Files Created. Ready to Download")
      output$ready=renderText({"ready"})
    })})}

options(shiny.launch.browser=TRUE)
shinyApp(ui,server)




