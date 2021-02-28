# Example analysis of Golub et al 1999 gene expression dataset

#### Setup
# INSTALL THESE PACKAGES AND SET YOUR WORKING DIRECTORY WHERE THE DATA ARE

install.packages('e1071')

library(caret)
library(randomForest) # requires e1071
library(e1071)

# Read in data
outcomes <- read.csv('ML/actual.csv', header = TRUE, stringsAsFactors = FALSE)
train <- read.csv('ML/data_set_ALL_AML_train.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 'Gene.Accession.Number')
test <- read.csv('ML/data_set_ALL_AML_independent.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 'Gene.Accession.Number')

# Some data formatting 
# pull out only expression values only for genes that are in both train and test data
# format as subjects x genes for easy use in modeling
genesToUse <- intersect(rownames(train), rownames(test))
trainCleaned <- t(train[rownames(train) %in% genesToUse, 
                      grep('X', names(train))])
testCleaned <- t(test[rownames(test) %in% genesToUse, 
                    grep('X', names(test))])
outcomes$patient <- paste0('X', outcomes$patient)
geneDescriptions <- data.frame(Gene.Accession.Number = rownames(train)[rownames(train) %in% genesToUse],
                               Gene.Description = train[rownames(train) %in% genesToUse, 'Gene.Description'])

# CHALLENGE QUESTIONS
# Can you verify the expected numbers of subjects and genes in the cleaned train and test datasets?
# How many subjects have AML vs ALL in the training data? How about in the test data?
# How would you look up the gene description for a particular gene accession number of interest?


#### Identify associated biomarkers in training data
# Using univariate logistic regression
trainOutcome <- as.factor(outcomes$cancer[match(rownames(trainCleaned), outcomes$patient)])
univariateResults <- sapply(1:ncol(trainCleaned), function(i) {
  fit <- glm(trainOutcome ~ unlist(trainCleaned[,i]), family = 'binomial')
  results <- summary(fit)$coef['unlist(trainCleaned[, i])', c('Estimate', 'Pr(>|z|)')]
  c(colnames(trainCleaned)[i], results)
})

# Summarize and explore results - how many genes are significant? Which ones? Are they up or down regulated, and by how much?
histogram(as.numeric(univariateResults['Pr(>|z|)',]))
sum(as.numeric(univariateResults['Pr(>|z|)',]) < 0.05)
univariateResults[1,as.numeric(univariateResults['Pr(>|z|)',]) < 0.05]
histogram(as.numeric(univariateResults['Estimate',][as.numeric(univariateResults['Pr(>|z|)',]) < 0.05]))

# CHALLENGE QUESTIONS
# Instead of identifying significant genes, can you identify the genes most highly correlated with the outcome in this dataset?
# Can you display visually how significant or highly correlated genes differ in expression between AML and ALL?
# Can you apply unsupervised methods (e.g. a biclustered heatmap) to visually display sets genes that are up- or down- regulated together in groups of subjects, and whether these subjects group according to AML vs ALL diagnosis?





#### Account for multiple hypothesis testing
# What does this do to the number of significant genes?
bonferroni <- 0.05 / ncol(trainCleaned)
sum(as.numeric(univariateResults['Pr(>|z|)',]) < bonferroni)


permutedOutcome <- sapply(1:100, function(i) {
  sample(trainOutcome)
})

minPvalues <- apply(permutedOutcome, 2, function(outcome) {
  permutedResults <- sapply(1:ncol(trainCleaned), function(i) {
    fit <- glm(as.factor(outcome) ~ unlist(trainCleaned[,i]), family = 'binomial')
    results <- summary(fit)$coef['unlist(trainCleaned[, i])', c('Estimate', 'Pr(>|z|)')]
    c(colnames(trainCleaned)[i], results)
  })
  minP <- min(as.numeric(permutedResults['Pr(>|z|)',]))
  minP
})

histogram(minPvalues)
permCutoff <- quantile(minPvalues, 0.05)
sum(univariateResults['Pr(>|z|)',] < permCutoff)

# CHALLENGE QUESTIONS
# Which genes did you find were significant? What are their descriptions?
# What happens if you rerun the permutations here? Does anything change? Why?





#### Predictive model building with caret
# Set the seed
set.seed(1234)

trainPredict <- data.frame(outcome = trainOutcome,
                           trainCleaned)

# Train a model
control <- trainControl(method = 'cv',
                        number = 10)
rfFit <- train(outcome ~ .,
               data = trainPredict,
               method = 'rf',
               trControl = control)

# Look at the best model - what are its characteristics in the training data?
rfFit

# Apply the model to the validation set
testPredict <- data.frame(outcome = as.factor(outcomes$cancer[match(rownames(testCleaned), outcomes$patient)]),
                          testCleaned)
testPredictions <- data.frame(predicted = predict(rfFit, newdata = testPredict),
                              actual = testPredict$outcome)
table(testPredictions$predicted, testPredictions$actual)
confusionMatrix(data = testPredictions$predicted,
                reference = testPredictions$actual)

# Look at variable importance
varImp(rfFit)


# CHALLENGE QUESTIONS
# Which genes are most important in the predictive model? What do these genes do?
# What happens if you use logisitic regression instead of RF for the predictive model? Why?
# Can you use a smaller number of genes (perhaps the most important ones) in the logistic regression and build a highly performing (albeit overfit) model?
# Why might the number of significant genes found in the previous section be so small in spite of the high predictive accuracy?