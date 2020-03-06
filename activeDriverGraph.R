library(ggplot2)
library(plyr)
library(tidyverse)

# helper function to collect all elements within the windowsize
# example 
# WSgenes = getElements(200000, 51097600, "chr7", mm9Ensembl)
getElements <- function(window_size, start, chr, elements) {
  #collect elements for the specified chromosome
  elements <- elements[elements$chr == chr,]
  
  #collect all elements within 
  within <- elements[elements$start >= start,]
  within <- within[within$end <= start + window_size,]
  
  #collect elements that have a ending location greater than the start of the window
  #change the elements starting position to the window start
  endEdge <- elements[elements$end >= start + window_size,]
  endEdge <- endEdge[endEdge$start >= start,]
  endEdge <- endEdge[endEdge$start <= start + window_size,]
  if (length(endEdge$start) > 0) {
    endEdge$end = start + window_size
  }
  
  #collect elements that have a starting location less than the end of the window
  #change the elements ending positiong to the window end
  startEdge <- elements[elements$start <= start,]
  startEdge <- startEdge[startEdge$end <= start + window_size,]
  startEdge <- startEdge[startEdge$end >= start,]
  if (length(startEdge$start) > 0){
    startEdge$start = start
  }
  
  #collect large elements that span across the window size
  #start before AND end after 
  spanAcross <- elements[elements$start <= start,]
  spanAcross <- spanAcross[spanAcross$end >= start + window_size,]
  if (length(spanAcross$start) > 0){
    spanAcross$start = start
    spanAcross$end = start + window_size
  }
  
  WSgenes = rbind(within, endEdge, startEdge, spanAcross)
  
  return(WSgenes)
}

# Function that plots a visual graph of the mutations within a gene element in a specified window size
# window_size = 50 000 by default

#work in progress
# activeDriverGraph(100000, "uc008ivl.2", "Notch1", mouseBased, mm9)
# activeDriverGraph(100000, "uc009tsz.1", "Zfx", mouseBased, mm9)
# activeDriverGraph(100000, "uc007lkz.1", "Jup", mouseBased, mm9)
# activeDriverGraph(100000, "uc008yoe.1", "Pcgf3", mouseBased, mm9)
# activeDriverGraph(100000, "uc008oao.2", "Pard6b", mouseBased, mm9)
# activeDriverGraph(100000, "uc009bxr.1", "Cbx3", mouseBased, mm9)

# activeDriverGraph(100000, "uc008ibx.1", "Eif3a", mouseBased, mm9)
# activeDriverGraph(100000, "uc008qsg.2", "Csde1", mouseBased, mm9)
# activeDriverGraph(100000, "uc007lvt.1", "Nsf", mouseBased, mm9)

activeDriverGraph <- function(window_size, gene, geneName, mutations, elements){
  #get the gene of interest 
  goi <- mm9[mm9$id == gene,]
  goiStart <- goi$start
  goiEnd <- goi$end
  goiChr <- goi$chr
  
  wsStart <- goiStart - window_size
  wsEnd <- goiEnd + window_size
  ws <- wsEnd - wsStart
  
  WSgenes = getElements(ws, wsStart, goiChr, elements)
  
  WSmut <- mutations[mutations$chr == goiChr,]
  WSmut <- WSmut[WSmut$pos1 >= wsStart,]
  WSmut <- WSmut[WSmut$pos1 <= wsEnd,]
  
  # return(WSmut)
  
  # round to nearest 100 bp
  WSmut$pos1 <- round_any(WSmut$pos1, 100)
  
  # make dataframe with: 
  ## integer of number of mutations at 100bp location
  ## boolean of whether this location is within an element
  ## rows = the bps starting at location argument in increments of 100bps
  
  # number of rows 
  bp = ws / 100
  bpStart = ws
  
  # sum up the number of mutations at each 00 bp
  mutCount = c()
  bp100 = c()
  count = c()
  element = c()
  
  curr = round_any(wsStart, 100)
  for (i in 1:bp) {
    mutCount[i] <- sum(WSmut$pos1 == curr)
    bp100[i] <- curr 
    # check if this position contains any elements
    for (j in 1:length(WSgenes$chr)){
      count[j] <- between(curr, WSgenes$start[j], WSgenes$end[j])
      #if element is the gene under investigation
      if (WSgenes$id[j] == gene &&  between(curr, WSgenes$start[j], WSgenes$end[j])) {
        element[i] <- "goi"
        break
      } else if (any(count == TRUE)) {
        element[i] <- "gene"
      } else {
        element[i] <- "noncoding"
      }
    }
    curr = curr + 100
  }
  
  # make into dataframe 
  df <- data.frame(mutCount, element, bp100)
  
  bp100 <- bp100[seq(1,length(bp100),100)]
  yaxis <- seq(0, max(mutCount),100)
  
  #plot as scatterplot
  result <- ggplot(data=df, aes(x=bp100, y=mutCount, fill=element)) +
    geom_bar(stat = "identity") + labs(colour="Insertion within Element",
                                       x="Location in Chromosome (100bps)", 
                                       y="Number of insertions") + 
    ggtitle(paste0(geneName, ": ", "Window Size: ", window_size)) +  
    theme_light() +
    scale_x_continuous(breaks=bp100) +
    scale_y_continuous(breaks = yaxis) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) + 
    geom_point(aes(color=element)) +
    scale_colour_manual(label = c("Gene", "Gene of Interest", "Non-coding region"),
                        values = c("goi" = "red2", "gene" = "grey46", "noncoding" = "grey77")) +
    scale_fill_manual(label = c("Gene of Interest", "Gene", "Non-coding region"),
                        values = c("goi" = "red2", "gene" = "grey46", "noncoding" = "grey77")) +
    guides(fill = FALSE)
  result
  # return(df)
}