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

#example 
# mut 
# graphWS(window_size = 200000, "chr17", 40971633, mouseBased, mm9)
# Trps
# graphWS(window_size = 200000, "chr15", 50386304, mouseBased, mm9)
# Met
# graphWS(window_size = 200000, "chr6", 17313956, mouseBased, mm9)
# ccdc41
# graphWS(window_size = 200000, "chr10", 94051534, mouseBased, mm9)
# Bag3
# graphWS(window_size = 200000, "chr7", 135567096, mouseBased, mm9)
# Gnb1 
# graphWS(window_size = 200000, "chr4", 154833378, mouseBased, mm9)
# Baiap2
# graphWS(window_size = 200000, "chr11", 119704405, mouseBased, mm9)
# Tspan12
# graphWS(window_size = 200000, "chr6", 21621394, mouseBased, mm9)
# Tab2 
# graphWS(window_size = 200000, "chr10", 7525445, 7625445, mouseBased, mm9)


graphWS <- function(window_size = 50000, chr, location, investGeneStart, mutations, elements){
  WSgenes = getElements(window_size, location, chr, elements)
  # collect all mutations that are within the window
  WSmut <- mutations[mutations$chr == chr,]
  WSmut <- WSmut[WSmut$pos1 >= location,]
  WSmut <- WSmut[WSmut$pos1 <= location + window_size,]
  
  # round to nearest 100 bp
  WSmut$pos1 <- round_any(WSmut$pos1, 100)
  
  # make dataframe with: 
  ## integer of number of mutations at 100bp location
  ## boolean of whether this location is within an element
  ## rows = the bps starting at location argument in increments of 100bps
  
  # number of rows 
  bp = window_size / 100
  bpStart = window_size 
  
  # sum up the number of mutations at each 00 bp
  mutCount = c()
  bp100 = c()
  count = c()
  eleBool = c()
  investGene = c()
  
  curr = round_any(location, 100)
  for (i in 1:bp) {
    mutCount[i] <- sum(WSmut$pos1 == curr)
    bp100[i] <- curr 
    # check if this position contains any elements
    for (j in 1:length(WSgenes$chr)){
      count[j] <- between(curr, WSgenes$start[j], WSgenes$end[j])
      if (any(count == TRUE)) {
        eleBool[i] <- any(count == TRUE) 
        #if element is the gene under investigation
        if (WSgenes$start[j] == investGeneStart) {
          investGene[i] = TRUE
        }
        investGene[i] = FALSE
      } else {
        eleBool[i] <- any(count == TRUE) 
        investGene[i] = FALSE
      }
    }
    curr = curr + 100
  }
  
  # make into dataframe 
  df <- data.frame(mutCount, eleBool, bp100, investGene)
  
  bp100 <- bp100[seq(1,length(bp100),20)]
  
  #plot as scatterplot
  result <- ggplot(data=df, aes(x=bp100, y=mutCount, fill=eleBool)) +
    geom_bar(stat = "identity") + labs(colour="Insertion within Element",
                                       x="Location in Chromosome by 100 bps", 
                                       y="Number of insertions") + 
    ggtitle(paste0(chr, " ", "Window Size: ", window_size)) +  
    theme_light() +
    # theme(axis.line = element_line(colour = "black")) +
    scale_x_continuous(breaks=bp100)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(aes(color=investGene)) + 
    
  
  result
  # return(df)
}