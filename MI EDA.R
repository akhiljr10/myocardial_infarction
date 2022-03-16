#Setting the working directory

getwd()

setwd('D:/DOCUMENTS/process documents/IUPUI/Spring 2021/Applied Statistics/Project/MI')
getwd()

#installing the required packages

#install.packages('tidyverse')
#install.packages('ggplot2')
#install.packages("ggpubr")
#install.packages("missMDA")
#install.packages('ggfortify')
#install.packages("ISLR")
#install.packages("SmartEDA")
#install.packages("xfun")
#install.packages('nnet')


library("data.table")

library("nnet")

library("ISLR")

library("SmartEDA")

library('ggfortify')

library('tidyverse')

library('ggplot2')
library("ggpubr")

library('fastDummies')

library('mltools')

library('tidyverse')
library('ggplot2')
library("ggpubr")
library("corrplot")

library("plotROC")

library('pROC')

library('ResourceSelection')

library('caret')

library('MKmisc')

library('ROCR')

library('lmtest')

library('missMDA')

#reading the data

newdata <- data.frame(read.csv('Myocardial infarction complications Database.csv'))

#Subsetting dataframe with factors that can be measured at time of admission

my_data <- newdata[, c(2:91, 113:124)]

#Looking at the data

head(my_data)

tail(my_data)

descriptive <- ExpData(data=my_data,type=2, fun = c("mean", "median", "var"))
descriptive

#Looking for null values

sum(is.na(my_data))

mean(is.na(my_data))

summary(my_data)

#Handling the null values

names(which(colSums(is.na(my_data)) > 0))

#we have 89 columns with null values

my_data$AGE[is.na(my_data$AGE)] <- mean(my_data$AGE,na.rm=TRUE) 

## Remove columns with more than 30% NA

my_data <- my_data[, which(colMeans(!is.na(my_data)) > 0.30)]

#Handling the missing values

#fill the columns having numerical values with mean of the data

my_data$ROE[is.na(my_data$ROE)] <- mean(my_data$ROE,na.rm=TRUE) 

my_data$L_BLOOD[is.na(my_data$L_BLOOD)] <- mean(my_data$L_BLOOD,na.rm=TRUE) 

my_data$AST_BLOOD[is.na(my_data$AST_BLOOD)] <- mean(my_data$AST_BLOOD,na.rm=TRUE) 

my_data$NA_BLOOD[is.na(my_data$NA_BLOOD)] <- mean(my_data$NA_BLOOD,na.rm=TRUE) 

my_data$K_BLOOD[is.na(my_data$K_BLOOD)] <- mean(my_data$K_BLOOD,na.rm=TRUE) 

my_data$ALT_BLOOD[is.na(my_data$ALT_BLOOD)] <- mean(my_data$ALT_BLOOD,na.rm=TRUE) 

finaldf <- subset(my_data, select = -c(S_AD_KBRIG, D_AD_KBRIG))

finaldf <- sapply( finaldf, as.numeric)

finaldf[is.na(finaldf)] <- median(finaldf, na.rm = TRUE)

finaldf <- as.data.frame(finaldf)

summary(finaldf)

finaldf$LET_IS[finaldf$LET_IS != 0] <- 1
###############################################################################################

#EXPLORATORY DATA ANALYSIS
#Numerical Variables

numsum <- ExpNumStat(finaldf,by="A",gp=NULL,Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2,Nlim=10)
numsum

#Density plot (Univariate)

plot1 <- ExpNumViz(finaldf,target=NULL,nlim=10,Page=c(3,3),sample=9)
plot1[[1]]

#frequency for all categorical independent variables

catsum <- ExpCTable(finaldf,Target=NULL,margin=1,clim=10,nlim=3,round=2,bin=NULL,per=T)
catsum

#Bar plots for all categorical variables

plot4 <- ExpNumViz(finaldf,target="LET_IS",type=1,nlim=3,fname=NULL,col=c("darkgreen", 'blue', 'red'),Page=c(3,3),sample=9)
plot4[[1]]

#Correlation plot
# calculate correlation matrix

correlationMatrix <- cor(finaldf [,1:98])

# find attributes that are highly corrected (ideally >0.75)

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.25)

# print indexes of highly correlated attributes

newcorr <- cor(finaldf[, c(highlyCorrelated)])
corrplot(newcorr,order = 'original', hclust.method = "average", addrect = 12, tl.cex = 0.8)

#Statistical Test

et4 <- ExpCatStat(finaldf,Target="LET_IS",result = "Stat",clim=10,nlim=5,bins=10,Pclass="Yes",plot=FALSE,top=20,Round=2)

options(width = 150)
CData = finaldf
qqp <- ExpOutQQ(CData,nlim=10,fname=NULL,Page=c(3,2),sample=6)
qqp[[1]]

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.50)
exparcoord <- ExpParcoord(CData,Group=NULL,Stsize=NULL,Nvar=c(highlyCorrelated))
exparcoord

#################################################################################################

pca <- prcomp(finaldf[,1:98], center = TRUE,scale. = TRUE)
d <- print(pca)

pca.plot <- autoplot(pca, data = finaldf, colour = 'LET_IS', loadings = TRUE,
                     loadings.colour = 'red',
                     loadings.label = TRUE, loadings.label.size = 3)
pca.plot




Train <- createDataPartition(finaldf$LET_IS, p=0.8, list=FALSE)
training <- finaldf[ Train, ]
testing <- finaldf[ -Train, ]

#Modeling after PCA

components <- cbind(LET_IS = finaldf[, "LET_IS"], pca$x[, 1:15]) %>%
  as.data.frame()

fit_1 <- glm(LET_IS ~ AGE + ritm_ecg_p_02 + MP_TP_POST +  ritm_ecg_p_01
             + SEX + ZSN_A + ant_im + K_SH_POST + RAZRIV, data = finaldf)

fit_2 <- glm(finaldf$LET_IS ~ components$PC1 + components$PC2 + components$PC3 + 
               components$PC4 + components$PC5 + components$PC6 + components$PC7 
             + components$PC8 + components$PC9 + components$PC10)

summary(fit_1)

summary(fit_2)

# glance(fit_1)
# glance(fit_2)

predictions <- fit_1 %>% predict(training)

data.frame(
  R2 = R2(predictions, training$LET_IS),
  RMSE = RMSE(predictions, training$LET_IS),
  MAE = MAE(predictions, training$LET_IS)
)

predictions <- fit_2 %>% predict(finaldf)
data.frame(
  R2 = R2(predictions, finaldf$LET_IS),
  RMSE = RMSE(predictions, finaldf$LET_IS),
  MAE = MAE(predictions, finaldf$LET_IS)
)

#Hence some of the variables that we added in the model were not significant we dropped them and 
#replaced them using another variables

fit_1 <- glm(LET_IS ~ AGE +  ritm_ecg_p_01 + ZSN_A + ant_im + K_SH_POST +
               INF_ANAM + O_L_POST + L_BLOOD + STENOK_AN + RAZRIV, data = finaldf)

#Hence PC9 is not statistically significant, we dropped it from the model and added another PC's

fit_2 <- glm(finaldf$LET_IS ~ components$PC1 + components$PC2 + components$PC3 + 
               components$PC4 + components$PC5 + components$PC6 + components$PC7 
             + components$PC8 + components$PC10 + components$PC11 + components$PC12)


summary(fit_1)

summary(fit_2)

# glance(fit_1)
# glance(fit_2)

predictions <- fit_1 %>% predict(training)

data.frame(
  R2 = R2(predictions, training$LET_IS),
  RMSE = RMSE(predictions, training$LET_IS),
  MAE = MAE(predictions, training$LET_IS)
)

predictions <- fit_2 %>% predict(finaldf)
data.frame(
  R2 = R2(predictions, finaldf$LET_IS),
  RMSE = RMSE(predictions, finaldf$LET_IS),
  MAE = MAE(predictions, finaldf$LET_IS)
)