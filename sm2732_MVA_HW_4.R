# This data set dates from 1988 and consists of four databases: Cleveland, Hungary, 
# Switzerland, and Long Beach V. It contains 76 attributes, including the predicted 
# attribute, but all published experiments refer to using a subset of 14 of them. 
# The "target" field refers to the presence of heart disease in the patient. 
# It is integer valued 0 = no disease and 1 = disease.


# Attribute information:

# age
# sex
# chest pain type (4 values)
# resting blood pressure
# serum cholestoral in mg/dl
# fasting blood sugar > 120 mg/dl
# resting electrocardiographic results (values 0,1,2)
# maximum heart rate achieved
# exercise induced angina
# oldpeak = ST depression induced by exercise relative to rest
# the slope of the peak exercise ST segment
# number of major vessels (0-3) colored by flourosopy
# thal: 0 = normal; 1 = fixed defect; 2 = reversable defect


library(readr)
library(cluster)
library(factoextra)
library(magrittr)
library(NbClust)
library(data.table)
library(dplyr)


# Factor Analysis

library(psych)

heart_data <- read_csv("/Users/sarju/Desktop/MITA Sem 2/MVA/Homework/Week6/heart.csv")
heart_data

colnames(heart_data)
str(heart_data)


attach(heart_data)
#heart_data[1:6]
heart_data

fit.pc <- principal(heart_data, nfactors=4, rotate="varimax")
fit.pc

round(fit.pc$values, 3)
fit.pc$loadings

# Loadings with more digits
for (i in c(1,3,2,4)) { print(fit.pc$loadings[[1,i]])}

# Communalities
fit.pc$communality

# Rotated factor scores, Notice the columns ordering: RC1, RC3, RC2 and RC4
fit.pc$scores

# Play with FA utilities
fa.parallel(heart_data) # See factor recommendation
fa.plot(fit.pc) # See Correlations within Factors
fa.diagram(fit.pc) # Visualize the relationship
vss(heart_data) # See Factor recommendations for a simple structure


# Computing Correlation Matrix
corrm.emp <- cor(heart_data)
corrm.emp
plot(corrm.emp)

heart_data_pca <- prcomp(heart_data, scale=TRUE)
summary(heart_data_pca)
plot(heart_data_pca)

# A table containing eigenvalues and %'s accounted, follows. Eigenvalues are the sdev^2
(eigen_heart_data <- round(heart_data_pca$sdev^2,3))
round(fit.pc$values, 3)
names(eigen_heart_data) <- paste("PC",1:14,sep="")
eigen_heart_data

sumlambdas <- sum(eigen_heart_data)
sumlambdas

propvar <- round(eigen_heart_data/sumlambdas,2)
propvar

cumvar_heart_data <- cumsum(propvar)
cumvar_heart_data

matlambdas <- rbind(eigen_heart_data,propvar,cumvar_heart_data)
matlambdas

rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
rownames(matlambdas)

eigvec.emp <- heart_data_pca$rotation
print(heart_data_pca)

# Taking the first four PCs to generate linear combinations for all the variables with four factors
pcafactors.emp <- eigvec.emp[,1:4]
pcafactors.emp

# Multiplying each column of the eigenvector’s matrix by the square-root of the corresponding eigenvalue in order to get the factor loadings
unrot.fact.emp <- sweep(pcafactors.emp,MARGIN=2,heart_data_pca$sdev[1:4],`*`)
unrot.fact.emp

# Computing communalities
communalities.emp <- rowSums(unrot.fact.emp^2)
communalities.emp

# Performing the varimax rotation. The default in the varimax function is norm=TRUE thus, Kaiser normalization is carried out
rot.fact.emp <- varimax(unrot.fact.emp)

#View(unrot.fact.emp)
rot.fact.emp

# The print method of varimax omits loadings less than abs(0.1). In order to display all the loadings, it is necessary to ask explicitly the contents of the object $loadings
fact.load.emp <- rot.fact.emp$loadings[1:14,1:4]
fact.load.emp

# Computing the rotated factor scores for the 30 European Countries. Notice that signs are reversed for factors F2 (PC2), F3 (PC3) and F4 (PC4)
scale.emp <- scale(heart_data[1:14])
scale.emp
as.matrix(scale.emp)%*%fact.load.emp%*%solve(t(fact.load.emp)%*%fact.load.emp)







