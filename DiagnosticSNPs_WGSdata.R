# Load necessary libraries
library(SNPRelate)
library(dplyr)
library(randomForest)
library(caret)
library(rpart)
library(partykit)
library(doParallel)
library(parallel)

##############################################################################################
# Run Recursive Feature Elimination and rpart on full SNP dataset to extract diagnostic SNPs #
##############################################################################################


#### Set up dataset ####

# Set working directory
setwd("path/with/.vcf")

# Define paths and load data
vcf.fn <- paste0(getwd(), "/yourvcfname.vcf.gz")  # Path to VCF file
metadata <- read.table("Bd_metadata_trainingRF.txt", header = TRUE, sep = "\t")  # Load metadata

# Convert VCF to GDS format for efficient handling
snpgdsVCF2GDS(vcf.fn, "BdReferences.gds", method = "biallelic.only")
genofile <- snpgdsOpen("BdReferences.gds")  # Open GDS file

# Extract genotypes coded as 0, 1, 2
genotypes <- snpgdsGetGeno(genofile, snp.id = NULL, snpfirstdim = TRUE)

# Check distribution of genotypes (0, 1, 2, and NA)
table(c(genotypes), exclude = NULL)

# Extract sample names
sample_names <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Create a dataframe with genotypes and sample names
df <- data.frame(sample_names, t(genotypes))
colnames(df) <- c("ID", paste0(read.gdsn(index.gdsn(genofile, "snp.chromosome")), "_", 
                               read.gdsn(index.gdsn(genofile, "snp.position"))))
                  
# Sanity check: inspect the first few rows and columns
dim(df)
df[1:5, 1:5]

# Add group information from metadata
df <- df %>% mutate(Group = metadata$Group[match(df$ID, metadata$ID)])

# Convert Group and ID to factors
region.vec <- as.factor(df$Group)
id.vec <- df$ID

# Convert all columns to factors (may take time with many SNPs)
df[, 1:dim(df)[2]] <- lapply(df[, 1:dim(df)[2]], as.factor)


#### Impute missing values using random forest (rfImpute) ####

print("Start imputing missing values...")

# Function to set up parallel processing
amp_up_models <- function() {
  library(parallel)
  library(doParallel)
  no_cores <- parallel::detectCores() - 1  # Use all but one core
  cluster <- makePSOCKcluster(no_cores)
  registerDoParallel(cluster)
  cat("Model amped and ready to go with:", no_cores, "cores. \n")
}
amp_up_models()

# Prepare data for imputation
df2 <- subset(df, select = -c(Group, ID))  # Remove Group and ID columns
ncores <- min(parallel::detectCores() - 1)
n_parts <- ncores
n_cols <- ncol(df2)
parts <- split.default(df2, rep(1:n_parts, each = ceiling(n_cols / n_parts)))  # Split data into chunks

# Check number of levels per column
sapply(df2, function(x) if (is.factor(x)) length(levels(x)) else NA)

# Perform parallel imputation using rfImpute
results_list <- foreach(i = 1:length(parts), .packages = c("randomForest", "dplyr")) %dopar% {
  part <- parts[[i]]
  imp <- rfImpute(region.vec ~ ., data = part, iter = 5, ntree = 300)  # Impute missing values
  return(imp)
}

# Combine results from parallel processing
imputed_df <- do.call(cbind, results_list)

# Stop the parallel cluster
stopCluster(cluster)

# Inspect imputed data
imputed_df[1:5, 1:5]
dim(imputed_df)

# Clean up column names and add Group column
imputed_df2 <- imputed_df[, grep("NC", colnames(imputed_df))]
df_train <- cbind(imputed_df2, Group = region.vec)

# Save imputed genotypes to a file
print("Saving imputed genotypes...")
write.table(x = df_train, "bd_intercept_ingroup_miss95_minDP6_References_imputed.txt", quote = FALSE, row.names = FALSE)

# Load imputed genotypes if already saved
df_train <- read_delim("bd_intercept_ingroup_miss95_minDP6_References_imputed.txt", col_names = TRUE, num_threads = 6, delim = " ")

# Convert all columns to factors
df_train[, 1:dim(df_train)[2]] <- lapply(df_train[, 1:dim(df_train)[2]], as.factor)


#### Prepare random forest training and validation set ####

# Create training and validation sets
df_rf <- df_train

# Simplify Group labels
unique(df_rf$Group)
df_rf <- df_rf %>%
  mutate(Group = as.character(Group)) %>%
  mutate(Group = case_when(
    Group %in% c("East-Africa", "West-Africa") ~ "Africa",
    TRUE ~ Group
  ))

# Check number of levels in each column
level_check <- subset(df_rf, select = -c(Group)) %>%
  summarise(across(where(is.factor), ~ nlevels(.) == 3))
print(level_check)

# Split data into training and validation sets
train_index <- createDataPartition(df_rf$Group, p = 0.75, list = FALSE, times = 1)
df.train <- df_rf[train_index, ]
df.validate <- df_rf[-train_index, ]


#### Random Forest Analysis using parallel computing ####

# Set up parallel processing
amp_up_models()

# Feature selection using Recursive Feature Elimination (RFE)
subsets <- c(1, 2, 5, 10, 15, 20, 50, 80, 100, 150, 200)  # Number of SNPs to test
set.seed(123)
rfeCtrl <- rfeControl(functions = rfFuncs,
                      method = "cv",
                      verbose = FALSE,
                      number = 10,
                      repeats = 10,
                      p = 0.7,
                      allowParallel = TRUE)
rfProfile <- rfe(x = df.train[, !names(df.train) %in% "Group"],
                 y = as.factor(df.train$Group), 
                 sizes = subsets,
                 rfeControl = rfeCtrl)

# Plot RFE results
plot(rfProfile, type = c("g", "o"), xlim = subsets)
axis(1, subsets)

# Select best features based on RFE
best.features <- caret::predictors(rfProfile)  # Automatically pick the best features
MeanDecreaseGini <- rfProfile$fit$importance[, dim(rfProfile$fit$importance)[2]]

# Alternatively, manually select the number of features
n.features <- 2
best.features <- caret::predictors(rfProfile)[1:n.features]
MeanDecreaseGini <- rfProfile$fit$importance[1:n.features, dim(rfProfile$fit$importance)[2]]


# Important info:
# The best feature (highest rank) is determined by the Recursive Feature Elimination (RFE) process, 
# which selects features based on their contribution to the models predictive accuracy.
# On the other hand, Gini Importance measures how much a feature contributes to reducing impurity in the splits.
# These two metrics are related but not the same, and they can sometimes disagree.

# Plot Gini importance
gini <- cbind(best.features, as.numeric(MeanDecreaseGini)) %>% data.frame()
colnames(gini) <- c("snp", "MeanDecreaseGini")
sorted_gini <- arrange(gini, desc(MeanDecreaseGini))
norm_gini <- sorted_gini %>% mutate(norm = as.numeric(.[, 2]) / max(as.numeric(.[, 2])))
norm_gini$snp <- rownames(norm_gini)
plot(norm_gini$norm, ylab = "Importance")

# Save diagnostic SNPs
best.features_df <- data.frame(snp = best.features, stringsAsFactors = FALSE)
diagn_snps <- best.features_df %>%
  mutate(
    chrom = str_replace(snp, "_\\d+$", ""),  # Extract chromosome
    pos = str_replace(snp, ".*_(\\d+)$", "\\1")  # Extract position
  ) %>%
  dplyr::select(chrom, pos)  # Select only chromosome and position columns

write.table(diagn_snps, "rf_diagn.snps_2.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Create a .bed file
write.table(as.data.frame(cbind(diagn_snps$chrom, diagn_snps$pos, diagn_snps$pos)), 
            "rf_diagn.snps_2.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

## Run random forest using only extracted features and assess accuracy via permutation ##
df_purged <- data.frame(Group = df_rf$Group, df_rf[,best.features])
df_purged[,1:dim(df_purged)[2]] <- lapply(df_purged[,1:dim(df_purged)[2]], as.factor)
train.p_index <- createDataPartition(df_purged$Group, p = 0.75, list = FALSE, times = 1)
df.train.purged <- df_purged[train.p_index,] 
dim(df.train.purged)
df.validate.purged<-df_purged[-train.p_index,] 
dim(df.validate.purged)
# Do the permutations
set.seed(123)
# Define number of permutations
n_permutations <- 100
purged_performance <- list()

for (i in 1:n_permutations) {
  cat("Iteration:", i, "\n")
  
  # Shuffle data and reassign factors
  df_purged <- data.frame(Group = df_rf$Group, df_rf[, best.features])
  df_purged[, 1:dim(df_purged)[2]] <- lapply(df_purged[, 1:dim(df_purged)[2]], as.factor)
  
  # Partition the data
  train.p_index <- createDataPartition(df_purged$Group, p = 0.75, list = FALSE, times = 1)
  df.train.purged <- df_purged[train.p_index, ] 
  df.validate.purged <- df_purged[-train.p_index, ] 
  
  # Train the Random Forest model
  rf_model_purged <- randomForest(x = df.train.purged[, 2:ncol(df.train.purged)],
                                  y = df.train.purged$Group, ntree = 10000, 
                                  importance = TRUE, mtry = 2)
  
  # Store performance table
  purged.pred <- predict(rf_model_purged, df.validate.purged)
  performance_table <- as.data.frame(table(df.validate.purged[,1], purged.pred, dnn=c("Actual", "Predicted")))
  
  # Compute correct classification percentage
  total_counts <- aggregate(Freq ~ Actual, data = performance_table, sum)
  correct_counts <- performance_table %>% filter(Actual == Predicted)
  accuracy_df <- merge(correct_counts, total_counts, by = "Actual", suffixes = (c("_correct", "_total"))) %>%
    mutate(Accuracy = (Freq_correct / Freq_total) * 100, Iteration = i) %>%
    select(Actual, Accuracy, Iteration)
  
  purged_performance[[i]] <- accuracy_df
}

# Combine results into a single dataframe
performance_df <- bind_rows(purged_performance, .id = "Iteration")

# Compute mean and standard deviation of accuracy per class
summary_df <- performance_df %>%
  group_by(Actual) %>%
  summarise(Mean_Accuracy = mean(Accuracy), SD_Accuracy = sd(Accuracy), .groups = "drop")

# Output the summarized results
summary_df

#Save the model
save(rf_model_purged, file = "Bdorsalis_AfricaAsia_2SNPs_rfeModel.rfe")



#### Make a classification model using rpart: make decision tree ####
# Apply the rpart function and investigate results
tree <- rpart(Group~., data=cbind(Group = df_rf$Group,df_rf[,best.features]), method = "class")
plot(as.party(tree))
tree$cptable
rsq.rpart(tree)
print(tree)
summary(tree)
# Prune the tree in case of complex trees
pruned_tree <- prune(tree, cp=tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"])
# plot pruned tree
plot(as.party(pruned_tree))


#### Use diagnostic SNPs on query sample ####

# Path to query VCF file
vcf_file <- "bd_intercept_ingroup_miss95_minDP6_Belgium.vcf.gz"
bed_file <- "rf_diagn.snps_2.bed"  # Path to .bed file
gds_file <- "DiagnSnps.gds"  # Temporary GDS file

# Convert VCF to GDS format
snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only")
genofile <- snpgdsOpen(gds_file)  # Open GDS file

# Read .bed file
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_data) <- c("CHROM", "POS", "POS2")

# Extract SNP IDs from .bed file
bed_ids <- paste(bed_data$CHROM, bed_data$POS, sep = "_")

# Get genotype data for diagnostic SNPs
vcf_ids <- paste0(read.gdsn(index.gdsn(genofile, "snp.chromosome")), "_", 
                  read.gdsn(index.gdsn(genofile, "snp.position")))
geno <- snpgdsGetGeno(genofile, snpfirstdim = FALSE)
colnames(geno) <- vcf_ids
rownames(geno) <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Close GDS file
snpgdsClose(genofile)

# Filter genotype data to include only diagnostic SNPs
diagn.snps <- subset(geno, select = bed_ids) %>% data.frame()

# Convert all columns to factors
diagn.snps[, 1:dim(diagn.snps)[2]] <- lapply(diagn.snps[, 1:dim(diagn.snps)[2]], as.factor)

# Check dimensions and missing values
dim(diagn.snps)
table(is.na(diagn.snps))



########################################
# Use diagnostic SNPs on query sample #
#######################################

# Define file paths
vcf_file <- "/path/to/vcf/with/query/samples.vcf.gz"  # Path to query VCF file
bed_file <- "rf_diagn.snps_20.bed"  # Path to .bed file containing diagnostic SNPs
gds_file <- "DiagnSnps.gds"  # Temporary GDS file for storing genotype data

# Convert VCF to GDS format for efficient handling
snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only")

# Open the GDS file
genofile <- snpgdsOpen(gds_file)

# Read the .bed file containing diagnostic SNPs
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_data) <- c("CHROM", "POS", "POS2")  # Rename columns for clarity

# Extract SNP IDs from the .bed file
bed_ids <- paste(bed_data$CHROM, bed_data$POS, sep = "_")

# Get genotype data for all SNPs in the VCF file
vcf_ids <- paste0(read.gdsn(index.gdsn(genofile, "snp.chromosome")), "_", 
                  read.gdsn(index.gdsn(genofile, "snp.position")))
geno <- snpgdsGetGeno(genofile, snpfirstdim = FALSE)  # Extract genotypes
colnames(geno) <- vcf_ids  # Assign SNP IDs as column names
rownames(geno) <- read.gdsn(index.gdsn(genofile, "sample.id"))  # Assign sample IDs as row names

# Close the GDS file
snpgdsClose(genofile)

# Filter genotype data to include only diagnostic SNPs
diagn.snps <- subset(geno, select = bed_ids) %>% data.frame()

# Convert all columns to factors
diagn.snps[, 1:dim(diagn.snps)[2]] <- lapply(diagn.snps[, 1:dim(diagn.snps)[2]], as.factor)

# Check dimensions and missing values in the filtered data
dim(diagn.snps)
table(is.na(diagn.snps))


## Load the trained model ##

# Load the pre-trained random forest model
load("path/to/model/Bdorsalis_AfricaAsia_2SNPs_rfeModel.rfe", verbose = TRUE)
rf_model_purged  # Display the loaded model in the R environment

# Plot the random forest model
plot(rf_model_purged)

# Function to extract predictor SNPs from the model
predictors <- function(model) {
  if (class(model) == "randomForest") {
    importance_vars <- importance(model)  # Extract variable importance
    predictor_names <- rownames(importance_vars)  # Get predictor names
  } else if (class(model) == "ranger") {
    predictor_names <- model$forest$variables$names  # Get predictor names for ranger models
  } else {
    stop("Unsupported model type. Supported types: 'randomForest' and 'ranger'")
  }
  return(predictor_names)
}

# Extract predictor SNPs from the model
predictor_snps <- predictors(rf_model_purged)
print(predictor_snps)


#### Make predictions using the trained model ####

# Prepare the diagnostic SNPs data for prediction
rownames(diagn.snps) <- NULL  # Remove row names
diagn.snps[diagn.snps == "<NA>"] <- NA  # Replace "<NA>" with NA

# Make predictions for a single sample (e.g., the second sample)
predict(rf_model_purged, as.factor(diagn.snps[2, ]), type = "prob")

# Visualize prediction probabilities for all samples
par(mfrow = c(4, 5))  # Set up a 4x5 grid for plots
for (i in 1:nrow(diagn.snps)) {
  predict(rf_model_purged, diagn.snps[i, ], type = "prob") %>%
    barplot(main = rownames(diagn.snps)[i], cex.names = 0.2, col = "lightgrey")
}

# Print prediction probabilities for all samples
for (i in 1:nrow(diagn.snps)) {
  print(paste(rownames(diagn.snps)[i], ": ", 
              predict(rf_model_purged, diagn.snps[i, ], type = "prob")))
}