library(randomForest)

#read in training dataset and the mito_score_matrix files
training <- read.table("all.training.csv", sep = ",", header = T, stringsAsFactors = F)
mitoscore <- read.table("mito_score_matrix", sep = ",", header = T, stringsAsFactors = F)

#train the model
rf <- randomForest(x = training[,-c(1,2,10,11)], y = as.factor(training$localization), importance = T, ntree = 500)

#predict 
predictedRF <- predict(rf, mitoscore[,-c(1,9,10)], type="response")  # predicted scores

#add column with prediction to matrix
predictedRF2 <- cbind(mitoscore, predictedRF)

#write matrix to file
write.csv(predictedRF2,file="rf.csv")





