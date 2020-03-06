library(ggplot2)
library(plyr)
library(tidyverse)
library(GenomicRanges)
library(trackViewer)
library(IRanges)

# helper function to collect all elements within the windowsize
# example 
# WSgenes = getElements(50000, 51097600, "chr7", mm9Ensembl)
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

#example 
# uc009bxz.1	- Hoxa2
# df <- lolligraph(window_size = 200000, "chr6", 52012509, mouseBased, mm9, 2)
# uc008kqf.1 - Olfr119
# df <- lolligraph(window_size = 200933, "chr2", 88495897, mouseBased, mm9, 2)
# p53 
# df <- lolligraph(window_size = 100000, "chr11", 69343860, mouseBased, mm9, 2)

lolligraph <- function(window_size = 50000, chr, location, mutations, elements, threshold){
  WSgenes = getElements(window_size, location, chr, elements)
  
  # collect all mutations that are within the window
  WSmut <- mutations[mutations$chr == chr,]
  WSmut <- WSmut[WSmut$pos1 >= location,]
  WSmut <- WSmut[WSmut$pos1 <= location + window_size,]
  
  # round to nearest 100 bp
  WSmut$pos1 <- round_any(WSmut$pos1, 100)
  
  # number of rows 
  bp = window_size / 100
  bpStart = window_size 
  
  # sum up the number of mutations at each 00 bp
  mutCount = c()
  bp100 = c()
  curr = round_any(location, 100)
  for (i in 1:bp) {
    mutCount[i] <- sum(WSmut$pos1 == curr)
    bp100[i] <- curr 
    # check if this position contains any elements
    curr = curr + 100
  }
  
  # make into dataframe 
  df <- data.frame(mutCount, bp100)
  df <- df[df$mutCount > threshold,]
  
  #make lolliplot 
  insertions <- GRanges(chr, IRanges(c(df$bp100), width = 1, names = df$bp100))
  insertions$score <- df$mutCount
  insertions.rot <- insertions
  insertions.rot$label.parameter.rot <- 60
    
  elements <- GRanges(chr, IRanges(WSgenes$start, width=WSgenes$end - WSgenes$start + 1, names=WSgenes$id))
  elements$fill <- sample(colours(), length(WSgenes$chr))
  # 
  # insertions$color <- 

  lolliplot(insertions.rot, elements, ylab="Number of Insertions")
  
  title <- paste0(chr, " ", "Window Size: ", window_size)
  grid.text(title, x=.5, y=.98, just="top", 
            gp=gpar(cex=1.5, fontface="bold"))
  
  return(df)
}