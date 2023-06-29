## Analyses and Figure prep for Avery et al., 2023 Food Similarity Paper
setwd("~/Documents/GitHub/FoodSimilarity/")
source("source_me.R")

task = "food_loc"; 
out_dir = "DATA";

## ========================================================================================== ##
## ========================================================================================== ##
## Network Graph of Food ROIs
## ========================================================================================== ##
## read in food mat data, if starting from scratch
food_roi_mat = read.csv(sprintf("%s/food_roi_mat.csv",out_dir), row.names = 1)
food_roi_corr_mat = read.csv(sprintf("%s/food_roi_corr_mat.csv",out_dir), row.names = 1)

## ========================================================================================== ##
## Test for outliers in ROI
corr_dist_mat = 1 - food_roi_corr_mat
euc_dist_mat = euc.dist(food_roi_mat)
## ========================================================================================== ##
## Calculate Mahalanobis distance for of each ROI from mean of distribution

oas_cov = (CovEst.2010OAS(t(food_roi_mat)))$S
roi_mahal_dist = sapply(1:length(rois), function(x) mahalanobis(food_roi_mat[,x],rowMeans(food_roi_mat[,-x]),oas_cov))
names(roi_mahal_dist) = rois
# V1 has greatest MD of all ROIs
# Dorsal_ACC  Left_Amyg   Left_dMI   Left_IFG   Left_MFG   Left_OFC   Left_vAI    Left_VS    Pre_SMA Right_Amyg 
# 19.40317   19.63901   19.69311   19.29921   19.75567   19.36727   19.85763   19.78370   19.56824   19.77713 
# Right_IFG  Right_IPS  Right_MFG   Right_MI  Right_OFC  Right_PCG         V1 
# 19.31467   19.51503   19.46706   19.22710   19.86877   19.89786   20.16307 
scale(roi_mahal_dist) # V1 = 2.09 SD greater than mean
## V1 is clear outlier
scale(rowMeans(corr.dist(food_roi_mat))) # V1 = 1.94 SD greater than mean
scale(rowMeans(euc.dist(food_roi_mat))) # V1 = 2.47 SD greater than mean
scale(rowMeans(mahal.dist(food_roi_mat))) # V1 = 2.09 SD greater than mean
## ========================================================================================== ##
## For following analyses, V1 ROI is removed
food_roi_corr_mat = cor(food_roi_mat[,1:16])

## ========================================================================================== ##
# Figure 2C: Make ROI heatmap
png(sprintf("%s/Food_ROI_heatmap.png",out_dir),width=6,height=6,units="in",res=300)
  heatmap(as.matrix(food_roi_corr_mat),labCol = NA,cexRow=1)
dev.off()
## ========================================================================================== ##
## Test similarity between network ROIs
df = food_roi_corr_mat; diag(df) = NA; df[lower.tri(df)] = NA
df = melt(as.matrix(df), na.rm = T, value.name = "weight" )

df$type = "between"
df[df$Var1 %in% hi_road&df$Var2 %in% hi_road,"type"] = "within"
df[df$Var1 %in% lo_road&df$Var2 %in% lo_road,"type"] = "within"
df$type = as.factor(df$type)

t.test(weight ~ type, data = df)
# data:  weight by type
# t = -10.207, df = 121.12, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1758694 -0.1187291
# sample estimates:
#   mean in group between  mean in group within 
# 0.3647266             0.5120259 
## ========================================================================================== ##
# Figure 2D: Make network graph
## ========================================================================================== ##
df = food_roi_corr_mat; diag(df) = 0; df[lower.tri(df)] = 0
df = melt(as.matrix(df), na.rm = T, value.name = "weight" )

net <- graph.data.frame(df, directed = FALSE)
net = delete_edges(net, which(E(net)$weight == 0))

V(net)$"label" = rownames(food_roi_corr_mat)
V(net)$"strength" = rowSums(food_roi_corr_mat, na.rm=TRUE)

comms <- cluster_optimal(net, weights = E(net)$weight) # find optimal network clusters
V(net)$comms = membership(comms); # write clusters back to graph

roi_color = V(net)$comms

vertex_colors = r_colors[roi_color]

find_edge_color = function(net,x){
  edge_color = "black"
  head_color = r_colors[V(net)$comms[head_of(net,x)]]; head_rgb = col2rgb(head_color)
  tail_color = r_colors[V(net)$comms[tail_of(net,x)]]; tail_rgb = col2rgb(tail_color)
  mean_color = rgb(t((head_rgb+tail_rgb)/2)/256)
  edge_color = ifelse(head_color == tail_color, tail_color, mean_color)
}

edge_colors = sapply(1:length(E(net)), function(x) find_edge_color(net,x))

labels = V(net)$label

png(sprintf("%s/Food_network_graph.png",out_dir),width=8,height=8,units="in",res=300)
plot.igraph(net,
            # margin=-.4,
            #+++++ Vertex
            vertex.color=vertex_colors,
            vertex.frame.color = "black",
            vertex.size=V(net)$strength*3,
            #+++++ Vertex Label
            vertex.label=labels,
            vertex.label.color="Black",
            vertex.label.family="Helvetica",
            vertex.label.dist=0,
            vertex.label.font=2,
            vertex.label.cex=1, #.5
            vertex.label.degree=0,
            #+++++ Edge
            edge.color=edge_colors,
            edge.width=(E(net)$weight),
            edge.curved=0.45,
            layout= layout_nicely)
dev.off()
## ========================================================================================== ##
## ========================================================================================== ##
## Do analyses of category and PC1 data
# Create a data table for all your RSA outputs
# Loop through subjects then ROIs
working_dir = "rsa_files_2022-10-30"
subj_data = list()

for (roi in rois){
  # make large list for all subj data
  roi_mat = sapply(subjects, function(x) ut(1-cor(read.csv(sprintf("%s/%s.%s.%s.food_mat.csv",working_dir,x,task,roi)))[,foods]))
  PC1_dat = apply(roi_mat,2, function(x) fisherz(cor(x,PC1_vec, method = "spearman", use = "complete")))
  PC2_dat = apply(roi_mat,2, function(x) fisherz(cor(x,PC2_vec, method = "spearman", use = "complete")))
  cluster_dat = apply(roi_mat,2, function(x) fisherz(cor(x,cluster_vec, method = "spearman", use = "complete")))
  category_dat = apply(roi_mat,2, function(x) fisherz(cor(x,category_vec, method = "spearman", use = "complete")))
  
  subj_data[[roi]] = cbind(subjects,roi,PC1_dat,PC2_dat,cluster_dat,category_dat)
}
subj_data = data.frame(do.call("rbind", subj_data)) 
names(subj_data) = c("subj","roi","PC1","PC2","cluster","category")

for (x in 3:ncol(subj_data)){
  subj_data[,x]=as.double(as.character(subj_data[,x]));
}
subj_data$diff = subj_data[,"cluster"] - subj_data[,"category"]

subj_data$roi = factor(subj_data$roi, levels = rois)

subj_data$network = NA
subj_data[subj_data$roi %in% hi_road,"network"] = "Prefrontal"
subj_data[subj_data$roi %in% lo_road,"network"] = "Limbic"

subj_data$network = factor(subj_data$network, levels = c("Prefrontal","Limbic"))

write.csv(subj_data, file = sprintf("%s/RSA_subject_data_%s_zscores.csv",out_dir,task), quote = TRUE, eol = "\n", na = "NA", row.names = TRUE)
## ========================================================================================== ##
subj_data = read.csv(sprintf("%s/RSA_subject_data_%s_zscores.csv",out_dir,task),row.names = 1,check.names = F)

grouping_data = melt(subj_data[,c("subj","roi","network","cluster","category")])
grouping_data$network = factor(grouping_data$network, levels = c("Prefrontal","Limbic"))

D = summary_stats(grouping_data,c('network','variable'),"value")

summary(aov(value ~ variable*network+Error((subj)) , data = grouping_data))

gtests = paired_ttests("variable",grouping_data,"value")

t.test(subj_data$cluster,subj_data$category,paired = T)
# data:  subj_data$cluster and subj_data$category
# t = 4.3728, df = 1219, p-value = 1.331e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.002887704 0.007587552
# sample estimates:
#   mean of the differences 
# 0.005237628
t.test(diff ~ network, data = subj_data)
# data:  diff by network
# t = 3.0931, df = 855.79, p-value = 0.002045
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.003114433 0.013930119
# sample estimates:
#   mean in group Prefrontal     mean in group Limbic 
# 0.011159458              0.002637181 

## ========================================================================================== ##
## Figure 3D: Compare cluster and category models in Both networks
## ========================================================================================== ##
plotting_data = bar_plotting(grouping_data,c('network','variable'),'value')
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

plot(myplot)
text_size = 8
bar_plot = myplot + 
  annotate("text", label = "", size = text_size, x = c(1:2), y = c(rep((0.03),2))) +
  geom_signif(comparisons=list(c("cluster","category")), test = "t.test",
              test.args=list(alternative = "greater", var.equal = TRUE, paired=TRUE),
              map_signif_level=TRUE,
              y_position = c(0.025), tip_length = 0, vjust=-0.3, textsize=3) +
  theme(axis.title.y = element_blank()) + 
  theme(strip.text = element_text(hjust = 0.5, color="black", size=24, face="bold"))+
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  labs(fill = "Model") + theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("%s/RSA_model_comparison.png",out_dir), width = 4, height = 3,dpi=300)
## ========================================================================================== ##
## Figure 3E. Compare diff b/t cluster and category b/t networks
## ========================================================================================== ##
subj_data$network = factor(subj_data$network, levels = c("Prefrontal","Limbic"))

plotting_data = bar_plotting(subj_data,"network","diff")
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

plot(myplot)
text_size = 8
bar_plot = myplot + 
  annotate("text", label = "", size = 3, x = c(1:2), y = c(rep((0.03),2))) +
  geom_signif(comparisons=list(c("Prefrontal","Limbic")), test = "t.test",
              test.args=list(alternative = "greater", var.equal = TRUE, paired=FALSE),
              map_signif_level=TRUE,
              y_position = c(0.02), tip_length = 0, vjust=-0.3, textsize=3) +
  theme(axis.title.y = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  scale_fill_manual(values=c(r_colors[2],r_colors[1],"green")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("%s/RSA_network_comparison_model_fits.png",out_dir), width = 4, height = 3,dpi=300)
## ========================================================================================== ##
## Figure 3F: Compare Both PC models across both networks
## ========================================================================================== ##
pc_data = melt(subj_data[,c("subj","roi","network","PC1","PC2")])
pc_data$network = factor(pc_data$network, levels = c("Prefrontal","Limbic"))

plotting_data = bar_plotting(pc_data,c('variable','network'),'value')
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

plot(myplot)

text_size = 8
bar_plot = myplot + 
  annotate("text", label = "", size = 3, x = c(1:2), y = c(rep((0.03),2))) +
  geom_signif(comparisons=list(c("Prefrontal","Limbic")), test = "t.test",
              test.args=list(alternative = "greater", var.equal = TRUE, paired=FALSE),
              map_signif_level=TRUE,
              y_position = c(0.025), tip_length = 0, vjust=-0.3, textsize=3) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(hjust = 0.5, color="black", size=24, face="bold"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  scale_fill_manual(values=c(r_colors[2],r_colors[1],"green")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("%s/RSA_network_comparison_2PCS.png",out_dir), width = 4, height = 3,dpi=300)

t.test(PC1 ~ network, data = subj_data)
# t = 3.1474, df = 813.4, p-value = 0.001707
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.004065009 0.017537684
# sample estimates:
#   mean in group Prefrontal     mean in group Limbic 
# 0.015335950              0.004534603 
## ========================================================================================== ##
## Figure 5: PSC Analyses

# read from csv files
ples_ratings = read.csv("DATA/PSC_ples_ratings_averages.csv", row.names = 1)
self_ratings = read.csv("DATA/PSC_selfctrl_ratings_averages.csv", row.names = 1)

# Now look at the correlations b/t ples and selfctrl ratings...
psc_corrs = sapply(levels(psc_data$subj), function(x) cor(ples_ratings[,x],self_ratings[,x],use = "complete"))
mean(abs(psc_corrs))
# [1] 0.9140807 # So the average correlation is pretty high

# correlation of mean ratings, by item
ples_means = rowMeans(ples_ratings,na.rm = T)
self_means = rowMeans(self_ratings,na.rm = T)
cor(ples_means,self_means)
# [1]  0.9793447 # And the correlation of the averages is also pretty high
## ========================================================================================== ##
## Now make a scatterplot with the mean ratings
setwd("~/Documents/OneDrive - National Institutes of Health/Studies/FoodLocSim/DATA")
library(ggimage)
df = data.frame(ples_means,self_means)
names(df) = c("ples","selfctrl")

df$food_file = sapply(as.character(row.names(df)), function(x) grep(sprintf("x_%s",x),list.files("psc_images"), value = T))
df$image = sprintf("psc_images/%s",df$food_file)

options(ggrepel.max.overlaps = Inf) #set max overlaps for ggrepel, for cluster plot

ratings_model = lm(selfctrl ~ ples, data = df)
rsqr = signif(summary(ratings_model)$adj.r.squared, 2)

ggplot(df, aes(x = ples, y = selfctrl)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 12, x = 3, y = 4.5) +
  geom_image(aes(image=image), size=.08) +
  geom_smooth(method="lm",formula=y~x, col = "blue") +
  theme(axis.text.x = element_text(size = 12,  face="bold")) +
  theme(axis.text.y = element_text(size = 12,  face="bold")) +
  xlab("Pleasantness Ratings - Average") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  ylab("Self-Control Ratings - Average") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold"))+
  theme(plot.title = element_text(size=22, face="bold", hjust = 0.5)) + #edit plot title size, font, and center
  labs(title = "Pleasantness vs. Self-Control Ratings")
ggsave(sprintf("FLS_analyses_2023-04-04/PSC_ples_vs_selfctrl_scatterplot_with_images.png",out_dir), width = 8, height = 6,dpi=300)


## ========================================================================================== ##
psc_data = read.csv("3dROIstats.ples.roi_data.04-27-2023.csv")
# psc_data = read.csv("~/Documents/GitHub/FoodSimilarity/DATA/3dROIstats.psc.roi_data.04-21-2023.csv")
names(psc_data) = c("subj","ROI","ples_UM","ples_AM")

psc_data = psc_data[!psc_data$ROI %in% visual,]
psc_data$ROI = factor(psc_data$ROI, levels = c(hi_road,lo_road))

psc_data$network = NA
psc_data[psc_data$ROI %in% hi_road,"network"] = "Prefrontal"
psc_data[psc_data$ROI %in% lo_road,"network"] = "Limbic"
psc_data$network = factor(psc_data$network, levels = c("Prefrontal","Limbic"))

# psc_um_data = melt(psc_data[,c(1,2,7,3,4)])
# psc_am_data = melt(psc_data[,c(1,2,7,5,6)])
psc_um_data = melt(psc_data[,c(1,2,5,3)])
psc_am_data = melt(psc_data[,c(1,2,5,4)])

psc_am_data$network = factor(psc_am_data$network, levels = c("Prefrontal","Limbic"))

## =================
## Plot Unmodulated data

plotting_data = bar_plotting(psc_um_data,c('network'),'value')
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

## =================
## =================
## Plot Amplitude Modulated data
plotting_data = bar_plotting(psc_am_data,c('network'),'value')
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

plot(myplot)
text_size = 8
bar_plot = myplot + 
  annotate("text", label = "", size = text_size, x = c(1:2), y = c(rep((0.04),2))) +
  geom_signif(comparisons=list(c("Prefrontal","Limbic")), test = "t.test",
              test.args=list(alternative = "two.sided", var.equal = TRUE, paired=FALSE),
              map_signif_level=TRUE,
              y_position = c(0.005), tip_length = 0, vjust=-0.3, textsize=8) +
  theme(axis.title.y = element_blank()) + 
  scale_fill_manual(values=c(r_colors[2],r_colors[1],"green")) +
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 12,  face="bold")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  labs(fill = "Model") + theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("DATA/FLS_ples_AM2_betas_by_network.png"), width = 4, height = 3,dpi=300)


bar_plot = myplot + 
  annotate("text", label = "", size = 3, x = c(1:2), y = c(rep((0.03),2))) +
  geom_signif(comparisons=list(c("Prefrontal","Limbic")), test = "t.test",
              test.args=list(alternative = "greater", var.equal = TRUE, paired=FALSE),
              map_signif_level=TRUE,
              y_position = c(0.025), tip_length = 0, vjust=-0.3, textsize=3) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(hjust = 0.5, color="black", size=24, face="bold"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  scale_fill_manual(values=c(r_colors[2],r_colors[1],"green")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("%s/RSA_network_comparison_2PCS.png",out_dir), width = 4, height = 3,dpi=300)

## =================
## Plot Amplitude Modulated data by ROI
plotting_data = bar_plotting(psc_am_data,c('ROI','variable'),'value')
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

plot(myplot)
text_size = 8
bar_plot = myplot + 
  annotate("text", label = "", size = text_size, x = c(1:2), y = c(rep((0.2),2))) +
  geom_signif(comparisons=list(c("ples_AM","self_AM")), test = "t.test",
              test.args=list(alternative = "two.sided", var.equal = TRUE, paired=TRUE),
              map_signif_level=TRUE,
              y_position = c(0.10), tip_length = 0, vjust=-0.3, textsize=5) +
  theme(axis.title.y = element_blank()) + 
  scale_fill_manual(values=c("red","blue","green")) +
  theme(strip.text = element_text(hjust = 0.5, color="black", size=24, face="bold"))+
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 12,  face="bold")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  labs(fill = "Model") + theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("DATA/PSC_ples_vs_selfctrl_rating_modulated_by_roi.png"), width = 6, height = 4.5,dpi=300)

## =================
# Some t-tests

D = summary_stats(grouping_data,c('network','variable'),"value")

t.test(psc_data[psc_data$network=="Prefrontal","ples_AM"],psc_data[psc_data$network=="Prefrontal","self_AM"], paired = T)
# t = -4.8525, df = 279, p-value = 2.031e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.019648449 -0.008307587
# sample estimates:
#   mean of the differences 
# -0.01397802 

t.test(psc_data[psc_data$network=="Limbic","ples_AM"],psc_data[psc_data$network=="Limbic","self_AM"], paired = T)
# t = -0.4817, df = 359, p-value = 0.6303
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.008119480  0.004924453
# sample estimates:
#   mean of the differences 
# -0.001597514 
# 
t.test(psc_data[psc_data$network=="Visual","ples_AM"],psc_data[psc_data$network=="Visual","self_AM"], paired = T)
# t = -1.2458, df = 239, p-value = 0.2141
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.009267880  0.002086963
# sample estimates:
#   mean of the differences 
# -0.003590458 

summary(aov(value ~ variable*network+Error((subj)) , data = psc_am_data))
# Significant effect of network, condition, and condition x network interaction
# Error: subj
# Df Sum Sq  Mean Sq F value Pr(>F)
# Residuals 39 0.3391 0.008694               
# 
# Error: Within
# Df Sum Sq Mean Sq F value Pr(>F)    
#   variable            1 0.0157 0.01574   6.590 0.0104 *  
#   network             1 0.1766 0.17656  73.909 <2e-16 ***
#   variable:network    1 0.0121 0.01207   5.053 0.0248 *  
#   Residuals        1237 2.9551 0.00239                   
# ---
# 
# plotting_data = bar_plotting(grouping_data,c('network','variable'),'value')
# myplot = plotting_data$myplot
# D = plotting_data$df
# D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))
# 
# plot(myplot)
# text_size = 8
# bar_plot = myplot + 
#   annotate("text", label = "", size = text_size, x = c(1:2), y = c(rep((0.04),2))) +
#   geom_signif(comparisons=list(c("ples","selfctrl")), test = "t.test",
#               test.args=list(alternative = "two.sided", var.equal = TRUE, paired=TRUE),
#               map_signif_level=TRUE,
#               y_position = c(0.035), tip_length = 0, vjust=-0.3, textsize=5) +
#   theme(axis.title.y = element_blank()) + 
#   scale_fill_manual(values=c("red","blue","green")) +
#   theme(strip.text = element_text(hjust = 0.5, color="black", size=24, face="bold"))+
#   theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
#   theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank()) +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
#   theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
#   labs(fill = "Model") + theme(legend.position = "none")
# 
# plot(bar_plot)
# ggsave(sprintf("PSC_ples_vs_selfctrl_by_network.png"), width = 6, height = 4.5,dpi=300)
#    
## =================

psc_data$BMI = fls_demo_stats[as.character(psc_data$subj),"BMI"]
summary(aov(ples_UM ~ BMI*network + Error(subj), data = psc_data))
summary(aov(self_UM ~ BMI*network + Error(subj), data = psc_data))

## =================
## aggregate psc_data by network
network_data = aggregate(. ~ subj+network+variable,data=psc_am_data,mean)

plotting_data = bar_plotting(network_data,c('network','variable'),'value')
myplot = plotting_data$myplot
D = plotting_data$df
D$p.adj=p.adjust(D$p.val,method="fdr",n=length((D$p.val)/length(levels(D$variable))))

plot(myplot)
text_size = 8
bar_plot = myplot + 
  annotate("text", label = "", size = text_size, x = c(1:2), y = c(rep((0.06),2))) +
  geom_signif(comparisons=list(c("ples_AM","self_AM")), test = "t.test",
              test.args=list(alternative = "two.sided", var.equal = TRUE, paired=TRUE),
              map_signif_level=TRUE,
              y_position = c(0.04), tip_length = 0, vjust=-0.3, textsize=5) +
  theme(axis.title.y = element_blank()) + 
  scale_fill_manual(values=c("red","blue","green")) +
  theme(strip.text = element_text(hjust = 0.5, color="black", size=24, face="bold"))+
  theme(axis.text.y = element_text(hjust = 0.5, color="black", size=10, face="bold"))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 12,  face="bold")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  labs(fill = "Model") + theme(legend.position = "none")

plot(bar_plot)
ggsave(sprintf("DATA/PSC_ples_vs_selfctrl_rating_modulated_by_network.png"), width = 6, height = 4.5,dpi=300)


