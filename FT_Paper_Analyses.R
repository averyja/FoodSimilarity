## Analyses and Figure prep for Carrington et al., 2023 Food Similarity Paper
setwd("~/Documents/GitHub/FoodSimilarity/")
source("source_me.R")

task = "food_loc"; 
data_dir = "FoodTriplet_DATA";

## ========================================================================================== ##
## Read in Mturk data 
mturk_data1 = read.csv(sprintf("%s/food_triplet_dataset1.csv",data_dir), row.names = 1)
mturk_data2 = read.csv(sprintf("%s/food_triplet_dataset2.csv",data_dir), row.names = 1)
## ========================================================================================== ##
participant_stats = read.csv(sprintf("%s/participant_stats.csv",data_dir), row.names = 1)
food_frame = read.csv(sprintf("%s/foods.csv",data_dir), row.names = 1)
food_vars = read.csv(sprintf("%s/food_pics_data.csv",data_dir),row.names = 1, check.names = F)
## ========================================================================================== ##
# Create similarity matrix from mturk data
# combine set1 and set2 data
mturk_data = rbind(mturk_data1[,c("food1","food2","food3","answer")],mturk_data2[,c("food1","food2","food3","answer")])

# strip trailing numbers off foods in mturk data
food1 = sub("_[0-9]","",mturk_data$food1)
food2 = sub("_[0-9]","",mturk_data$food2)
food3 = sub("_[0-9]","",mturk_data$food3)
answer = sub("_[0-9]","",mturk_data$answer)
# create empty array for pvalues
pval_mat=array(0,c(length(foods),length(foods)),dimnames = list(foods,foods))

Y = pval_mat;

for (i in foods){
  for (j in foods){
    if (i == j){pval_mat[i,j] = 1; Y[i,j] = "NA"; next}
    # X[i,j] = 0; Y[i,j] = 0;
    # return indices of triplets containing both foods
    triplets = which((food1 == i | food2 == i | food3 == i) & (food1 == j | food2 == j | food3 == j))
    answers = answer[triplets] # subset answer column by those indices
    nval = length(answers) # how many triples are in this set
    num_targets = sum(answers == j | answers == i) #find the number of times either i or j is in the list of answers
    pval_mat[i,j] = 1 - (num_targets)/nval # proportion we're looking for
    Y[i,j] = nval # store this value in a separate matrix
  }
}

outfile = sprintf("%s/food_similarity_matrix.csv",data_dir) #set the output file

# write resulting p-value matrix to a csv file
write.csv(pval_mat, file = outfile, quote = TRUE, eol = "\n", na = "NA", row.names = TRUE)

## ========================================================================================== ##
# Run PCA on sim matrix
sim.pca = princomp(pval_mat)
scores=data.frame(sim.pca$scores) #assign pca scores to a df

#Make PC scatterplot using ggimage 
library(ggimage)

df=data.frame(scores[,c(1:2)])
names(df)=c("PC1","PC2")
df$image = sprintf("image_files/%s.png",row.names(df))
options(ggrepel.max.overlaps = Inf) #set max overlaps for ggrepel, for cluster plot

# Create clustered PC plot with images, with 5 clusters
final <- kmeans(df[,c(1:2)], 5, nstart = 25)
# final$cluster = clusters_vec[row.names(df)]

row.names(food_frame) = food_frame$foods
final$cluster = food_frame[row.names(df),"cluster"]

fviz_cluster(final, data = df, choose.vars = c("PC1","PC2"), labelsize = 0, pointsize = 0, 
             ellipse = T, ellipse.level = 0.75, repel = T) + theme_classic() +
  # geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  xlab("Principal Component 1") + theme(axis.title.x = element_blank()) +
  ylab("Principal Component 2") + theme(axis.title.y = element_blank()) +
  ggtitle("Behavioral Data: PC1 vs. PC2") + theme(plot.title = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_image(aes(image=df$image), size=.08)
ggsave(sprintf("%s/FLS_PC_Plot_5clusters_with_images_2.png",data_dir), width = 7, height = 5,dpi=300)

## ========================================================================================== ##
## New Scatterplots with PC corrs, using images!
food_vars$image = sprintf("image_files/%s.png",row.names(food_vars))

library(ggimage)

## ========================================================================================== ##
health_model = lm(Comp.1 ~ health, data = food_vars)
rsqr = signif(summary(health_model)$adj.r.squared, 2)

ggplot(food_vars, aes(x = health, y = Comp.1)) + 
  geom_image(aes(image=food_vars$image), size=.08) + 
  stat_smooth(method = "lm", col = "red") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 14, x = 5, y = 0.75) +
  ylab("Principal Component 1") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
  xlab("Average Healthiness Ratings") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold"))
ggsave(sprintf("%s/PC1_vs_Health_scatterplot_with_images.png",data_dir), width = 8, height = 6,dpi=300)
## ========================================================================================== ##
processed_model = lm(P_rating ~ Comp.1, data = food_vars)
rsqr = signif(summary(processed_model)$adj.r.squared, 2)

ggplot(food_vars, aes(x = P_rating, y = Comp.1)) + 
  geom_image(aes(image=food_vars$image), size=.08) + 
  stat_smooth(method = "lm", col = "red") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 14, x = 4, y = 0.75) +
  ylab("Principal Component 1") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
  xlab("Average Processed Ratings") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold"))
ggsave(sprintf("%s/PC1_vs_Processed_ratings_scatterplot_with_images.png",data_dir), width = 8, height = 6,dpi=300)

## ========================================================================================== ##
sugar_model = lm(Comp.2 ~ sugar_est, data = food_vars)
rsqr = signif(summary(sugar_model)$adj.r.squared, 2)

ggplot(food_vars, aes(x = sugar_est, y = Comp.2)) + 
  geom_image(aes(image=food_vars$image), size=.08) + 
  stat_smooth(method = "lm", col = "red") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 14, x = 5, y = 0.75) +
  ylab("Principal Component 2") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
  xlab("Estimated Sugar Ratings") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold"))
ggsave(sprintf("%s/PC2_vs_Sugar_scatterplot_with_images.png",data_dir), width = 8, height = 6,dpi=300)

## ========================================================================================== ##
fat_sugar_model = lm(Comp.2 ~ fat_est+sugar_est, data = food_vars)
rsqr = round(summary(fat_sugar_model)$adj.r.squared, digits = 2)

food_vars$fat_sugar = rescale(food_vars$fat_est - food_vars$sugar_est)

ggplot(food_vars, aes(x = fat_sugar, y = Comp.2)) + 
  geom_image(aes(image=food_vars$image), size=.08) + 
  stat_smooth(method = "lm", col = "red") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 14, x = 0.5, y = 1.5) +
  ylab("Principal Component 2") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
  xlab("(Fat - Sugar) Ratings, Scaled") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold"))
ggsave(sprintf("%s/Fat-Sugar_vs_PC2_scatterplot_with_images.png",out_dir), width = 8, height = 6,dpi=300)
## ========================================================================================== ##

## ========================================================================================== ##
## Test fit of objective and subjective models to PCs from similarity data
objective_frame = data.frame(food_vars[,c(1:7)])
subjective_frame = data.frame(food_vars[,c(1,2,8:11)])

PC1_obj_model=lm(PC1~.,data=objective_frame[,c(1,3:7)])
PC2_obj_model=lm(PC2~.,data=objective_frame[,c(2,3:7)])
summary(PC1_obj_model)
# Multiple R-squared:  0.8088,	Adjusted R-squared:  0.7769 
summary(PC2_obj_model)
# Multiple R-squared:  0.7161,	Adjusted R-squared:  0.6688 

PC1_subj_model=lm(PC1~.,data=subjective_frame[,c(1,3:6)])
PC2_subj_model=lm(PC2~.,data=subjective_frame[,c(2:6)])

summary(PC1_subj_model)
# Multiple R-squared:  0.9725,	Adjusted R-squared:  0.969 
summary(PC2_subj_model)
# Multiple R-squared:  0.913,	Adjusted R-squared:  0.9018 

AIC(PC1_obj_model)
# [1] 64.95284
AIC(PC1_subj_model)
# [1] -6.885613
AIC(PC2_obj_model)
# [1] 60.66431
AIC(PC2_subj_model)
# [1] 16.08172

# Now run Wilk's Likelihood-ratio test to compare goodness of fit
library(lmtest)

lrtest(PC1_subj_model,PC1_obj_model)
# Likelihood ratio test
# 
# Model 1: Comp.1 ~ P_rating + fat_est + sugar_est + health
# Model 2: Comp.1 ~ kcal + fat + carb + protein + sugar
# #Df   LogLik Df  Chisq Pr(>Chisq)    
# 1   6   9.4428                         
# 2   7 -25.4764  1 69.838  < 2.2e-16 ***

lrtest(PC2_subj_model,PC2_obj_model)
# Likelihood ratio test
# 
# Model 1: Comp.2 ~ P_rating + fat_est + sugar_est + health
# Model 2: Comp.2 ~ kcal + fat + carb + protein + sugar
# #Df   LogLik Df  Chisq Pr(>Chisq)    
# 1   6  -2.0409                         
# 2   7 -23.3322  1 42.583  6.776e-11 ***
## ========================================================================================== ##
## ========================================================================================== ##
## Compare linear and logarithmic models for fat and sugar ratings
fat_ratings = food_vars$fat_estimate
sugar_ratings = food_vars$sugar_estimate

fat_content = (food_vars$`Fat/100g`) + 1 
sugar_content = (food_vars$`Sugar/100g`) + 1 

fat_linear_model = lm(fat_ratings ~ fat_content)
fat_log_model = lm(fat_ratings ~ log(fat_content))
summary(fat_linear_model)
AIC(fat_linear_model)
summary(fat_log_model)
AIC(fat_log_model)

lrtest(fat_log_model,fat_linear_model)

sugar_linear_model = lm(sugar_ratings ~ sugar_content)
sugar_log_model = lm(sugar_ratings ~ log(sugar_content))

summary(sugar_linear_model)
AIC(sugar_linear_model)
summary(sugar_log_model)
AIC(sugar_log_model)
lrtest(sugar_log_model,sugar_linear_model)
## ========================================================================================== ##
## Make figures for fat/sugar vs fat/sugar_ratings

#Make PC scatterplot using ggimage 
library(ggimage)

df=data.frame(sugar_content, sugar_ratings)
names(df)=c("sugar_content","sugar_ratings")
row.names(df) = row.names(food_vars)
df$image = sprintf("image_files/%s.png",row.names(df))

options(ggrepel.max.overlaps = Inf) #set max overlaps for ggrepel, for cluster plot

sugar_model = lm(sugar_ratings ~ log(sugar_content), data = df)
rsqr = signif(summary(sugar_model)$adj.r.squared, 2)

## 
scatter_plot = ggplot(df, aes(x=sugar_content, y=sugar_ratings)) + 
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 14, x = 45, y = 5) +
  geom_image(aes(image=image), size=.08) +
  geom_smooth(method="lm",formula=y~log(x))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme_classic()+
  xlab("Sugar per 100g") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
  ylab("Average Sugar Ratings") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold"))
plot(scatter_plot)

ggsave(sprintf("%s/Sugar_vs_Sugar_ratings_Log_Model_scatterplot.png",out_dir), width = 8, height = 6,dpi=300)
## ========================================================================================== ##
## Now with fat ratings
df=data.frame(fat_content, fat_ratings)
names(df)=c("fat_content","fat_ratings")
row.names(df) = row.names(food_vars)
df$image = sprintf("image_files/%s.png",row.names(df))
options(ggrepel.max.overlaps = Inf) #set max overlaps for ggrepel, for cluster plot

fat_model = lm(fat_ratings ~ log(fat_content), data = df)
rsqr = signif(summary(fat_model)$adj.r.squared, 2)

## 
scatter_plot = ggplot(df, aes(x=fat_content, y=fat_ratings)) + 
  annotate("text", label = bquote(r^2~"="~.(rsqr)), size = 14, x = 20, y = 4.5) +
  geom_image(aes(image=image), size=.08) +
  geom_smooth(method="lm",formula=y~log(x))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme_classic()+
  xlab("Fat per 100g") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
  ylab("Average Fat Ratings") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold"))
plot(scatter_plot)

ggsave(sprintf("%s/Fat_vs_Fat_ratings_Log_Model_scatterplot.png",out_dir), width = 8, height = 6,dpi=300)
## ========================================================================================== ##


