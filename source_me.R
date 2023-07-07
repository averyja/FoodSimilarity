## A source file for food_loc analyses that contains some variables and helper functions
# Attach Packages
library(plyr)
library(tidyverse)    # data manipulation and visualization
library(ggsignif)
library(igraph)
library(reshape2)
library(CovTools)
library(RColorBrewer)

# you may not need these
library(psych)
library(sys)
library(R.utils)
library(cluster)
library(factoextra)
library(FactoMineR)
library(NbClust)

setwd("~/Documents/GitHub/FoodSimilarity/")
task = "food_loc"; 

## ========================================================================================== ##
## Convenience Functions =====
## ========================================================================================== ##
# Annotate p-values with significance labels
sig_label = function(p){
  ifelse(p>0.05,"NS",ifelse(p<0.001,"***",ifelse(p<0.01,"**","*")))
}
## ========================================================================================== ##
## Function to change the range of an input vector or matrix to between 0..1
rescale <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
## ========================================================================================== ##
## Make a vector from the upper or lower triangle of a matrix
ut <- function(x){x[upper.tri(x)]}
## ========================================================================================== ##
lt <- function(x){x[lower.tri(x)]}
## ========================================================================================== ##
graph_plot = function(graph,food_color){
  vertex_colors = color_community[food_color]
  edge_colors = sapply(1:length(E(graph)), function(x) color_community[V(graph)$comms[tail_of(graph,x)]])
  if(length(V(graph)$label)==0){V(graph)$label=V(graph)$food}
  labels = V(graph)$label
  plot.igraph(graph,
              # margin=-.4,
              #+++++ Vertex
              vertex.color=vertex_colors,
              vertex.frame.color = "white",
              vertex.size=V(graph)$strength*2,
              #+++++ Vertex Label
              vertex.label=labels,
              vertex.label.color="Black",
              vertex.label.family="Arial",
              vertex.label.dist=0,
              vertex.label.font=3,
              vertex.label.cex=1, #.5
              vertex.label.degree=0,
              #+++++ Edge
              edge.color=edge_colors,
              edge.width=(E(graph)$weight),
              edge.curved=0.45,
              layout= layout_nicely)
}
## ========================================================================================== ##
# Aggregate data to get summary stats table
summary_stats = function(dframe,factor,index){
  dframe = dframe[,c(factor,index)]; names(dframe) = c(factor,"idx")
  stats_table = ddply(dframe, c(factor), summarise, N = length(idx[!is.na(idx)]), mean = mean(idx, na.rm = T),
                      sd = sd(idx, na.rm = T), se = sd/sqrt(N), max_height = mean + se, T = mean/se, 
                      p.val = round(pt(-abs(T),df=N-1), digits = 3), label = sig_label(p.val))
  return(stats_table)
}
## ========================================================================================== ##
## New bar plotting function, using settings that I make use of a lot.
## Ideally, this table will have been melted to long format first, so it has only one numeric value for plotting
bar_plotting <- function(dframe,factor,index){
  # First, format input table 
  value = index
  variable = ifelse(length(factor) > 1,factor[2],factor[1]) #the name of the variable in your df, whether roi,task, or other factor
  sprintf("Setting variable as %s, value as %s",variable,value) #print out result for error checking
  
  df = summary_stats(dframe,factor,value) #create the data table using function above
  
  # set some df variables used for plotting
  nval = ifelse(length(factor) > 1,length(levels(df[,2])),length(levels(df[,1]))) #if present, number of levels of 2nd factor
  xval = ifelse(length(factor) > 1,length(levels(df[,1])),1) #if present, number of levels of 2nd factor
  # some extra df variables used for plotting
  df[,value] = 0
  df$x = rep(c(1:nval),xval)
  
  # set these variables as parameters for ggplot
  start_height=max(df$max_height); bar_size = 0; text_size = 3; delta_y = start_height/10; ypos = start_height + 2*delta_y; 
  
  myplot <- ggplot(dframe, aes_string(variable, value, fill=variable)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2) +
    geom_hline(yintercept=0) +
    geom_text(y = ypos, aes(label = label), size = text_size, data = df, vjust = "inward") + 
    annotate("text", label = "", size = text_size, x = c(1:nval), y = c(rep((ypos+delta_y),nval))) +
    labs(fill = variable)
  
  if (length(factor) > 1){
    myplot = myplot + facet_wrap(as.formula(paste("~", factor[1]))) + theme(strip.text.x = element_text(size = 8, colour = "black")) 
  } 
  
  # set names and output df and plot
  output = setNames(list(df,myplot),c("df","myplot"))
  return(output)
}

## ========================================================================================== ##
## Multivariate Distance Metrics =====
## ========================================================================================== ##
## ========================================================================================== ##
## Function for calculating an RDM using regular mahalanobis distance
mahal.dist = function(your_data,your_labels=colnames(your_data)){
  X = array(0,c(length(your_labels),length(your_labels)),dimnames = list(your_labels,your_labels)) #make array for all stims
  your_data = your_data[,your_labels]
  oas_cov = (CovEst.2010OAS(t(your_data)))$S
  for (i in 1:(length(your_labels)-1)){
    x = your_data[,i]
    for (j in (i+1):length(your_labels)){
      y = your_data[,j]
      X[i,j] = mahalanobis(x, y, oas_cov)
    }
  }
  return(X+t(X))
}
## ========================================================================================== ##
## Correlation distance function, implemented in order to match format of other distance functions
corr.dist <- function(your_data,your_labels=colnames(your_data)){
  your_data = your_data[,your_labels]
  X = 1 - cor(your_data, method = "pearson")
  row.names(X) = your_labels; colnames(X) = your_labels
  return(X)
}
## ========================================================================================== ###
## simplified function for euclidean distance. Makes the function structure similar to cor or cos.dist 
## input and output
euc.dist = function(your_data,your_labels=colnames(your_data)){
  your_data = your_data[,your_labels]
  X = as.matrix(dist(t(your_data), diag = T, upper = T))
  row.names(X) = your_labels; colnames(X) = your_labels
  return(X)
}

## ========================================================================================== ##
## Useful variables =====
## ========================================================================================== ##
r_colors = c(brewer.pal(9,"Reds")[6],"gold",brewer.pal(9,"Blues")[5],brewer.pal(9,"Oranges")[5],brewer.pal(9,"Purples")[5],brewer.pal(9,"Greens")[6])
r_adj_colors = c(brewer.pal(9,"Reds")[8],"goldenrod4",brewer.pal(9,"Blues")[8],brewer.pal(9,"Oranges")[7],brewer.pal(9,"Purples")[8],brewer.pal(9,"Greens")[8])
## ========================================================================================== ##
food_frame = read.csv("FLS_DATA/FLS_foods.csv", row.names = 1)
food_vars = read.csv("FoodTriplet_DATA/food_pics_data.csv", row.names = 1)
pval_mat = read.csv("FoodTriplet_DATA/food_similarity_matrix.csv", row.names = 1)
foods = row.names(food_vars)
categories = c("HFLS","HFHS","LFLS","LFHS")
clusters = c("veggies","fats","fruits","starches","sweets") #assign names to food clusters 
## ========================================================================================== ##
## ========================================================================================== ##
limbic = c("Dorsal_ACC","Left_Amyg","Left_dMI","Left_vAI","Left_OFC","Left_VS","Right_MI","Right_OFC","Right_Amyg")
prefrontal = c("Left_IFG","Left_MFG","Pre_SMA","Right_IFG","Right_IPS","Right_MFG","Right_PCG")
roi_order = c(prefrontal,limbic)
visual = c("V1","Left_mOG","Right_mOG","Right_FG","Left_iTG","Right_PHG")

## ========================================================================================== ##
## Flattened upper triangle vectors used for RSA Analyses
PC1_vec = ut(as.matrix(dist(food_vars[,1])))
PC2_vec = ut(as.matrix(dist(food_vars[,2])))

cluster_vec = ut(read.csv("FLS_DATA/cluster_rdm.csv",row.names = 1))
category_vec = ut(read.csv("FLS_DATA/category_rdm.csv",row.names = 1))

subj_data = read.csv("FLS_DATA/RSA_subject_data_food_loc_zscores.csv",row.names = 1,check.names = F)
subjects = as.character(levels(subj_data$sublabel)); rois = as.character(levels(subj_data$roi))
## ========================================================================================== ##


