Recursive Feature Elimination and Random Forest Analysis for Diagnostic SNP Identification
This R script performs Recursive Feature Elimination (RFE) and Random Forest (RF) analysis on a Single Nucleotide Polymorphism (SNP) dataset to identify diagnostic SNPs for sample classification. It also includes data preprocessing, missing value imputation, and decision tree creation. The trained model can be applied to query samples for prediction.

Purpose
The script is designed to:

Identify diagnostic SNPs that are most informative for classifying samples into predefined groups (e.g., geographic regions).

Train a Random Forest model using the selected SNPs for classification.

Evaluate the model's performance using cross-validation and permutation tests.

Apply the trained model to a query sample for classification.

Prerequisites
Software and Libraries
R Environment: Ensure R is installed on your system.

R Libraries: Install the required R packages by running:
install.packages(c("SNPRelate", "dplyr", "randomForest", "caret", "rpart", "partykit", "doParallel", "parallel"))

Input Data
A VCF file containing SNP data.

A metadata file (Bd_metadata_trainingRF.txt) with sample IDs and corresponding group labels.

Computational Resources
The script uses parallel processing, so a multi-core system is recommended for efficient execution.

Use Case
This script is ideal for population genetics and genomic studies where:

You need to identify a small set of SNPs for accurate sample classification.

You want to build and evaluate a predictive model using Random Forest.

You need to classify new, unseen samples using the trained model.

Workflow
1. Data Preparation
Load the VCF file and convert it to GDS format for efficient handling.

Extract genotypes and merge them with metadata.

Convert data to factors for compatibility with machine learning algorithms.

2. Missing Value Imputation
Use Random Forest-based imputation (rfImpute) to handle missing genotype data.

Perform imputation in parallel to speed up the process.

3. Feature Selection
Use Recursive Feature Elimination (RFE) to identify the most informative SNPs.

Plot the results of RFE to visualize the importance of different numbers of SNPs.

4. Random Forest Model Training
Train a Random Forest model using the selected SNPs.

Evaluate the model's performance using cross-validation and permutation tests.

5. Decision Tree Creation
Build a decision tree using the rpart package to visualize classification rules.

Prune the tree to avoid overfitting.

6. Application to Query Sample
Apply the trained model to a new query sample to predict its classification.

Visualize the prediction probabilities.

Key Outputs
Diagnostic SNPs:

A list of SNPs identified as most informative for classification (rf_diagn.snps_2.txt and rf_diagn.snps_2.bed).

Trained Random Forest Model:

Saved as Bdorsalis_AfricaAsia_2SNPs_rfeModel.rfe.

Performance Metrics:

Mean accuracy and standard deviation of the model's predictions.

Decision Tree:

A pruned decision tree for visualizing classification rules.

Query Sample Predictions:

Classification probabilities for the query sample.

Example Usage
1. Set Up
Place your VCF file and metadata file in the working directory.

Update the file paths in the script (e.g., vcf.fn, metadata, vcf_file).

2. Run the Script
Execute the script in an R environment.

Monitor the output for progress and results.

3. Interpret Results
Use the diagnostic SNPs for further analysis.

Apply the trained model to classify new samples.

Repository Structure
Organize your repository as follows:

Copy
/Bdorsalis_RF_Analysis
├── /data
│   ├── yourvcfname.vcf.gz
│   ├── Bd_metadata_trainingRF.txt
├── /scripts
│   ├── RFE_RF_Analysis.R
├── /outputs
│   ├── rf_diagn.snps_2.txt
│   ├── Bdorsalis_AfricaAsia_2SNPs_rfeModel.rfe
├── README.md
Dependencies
List of required R packages:

SNPRelate

dplyr

randomForest

caret

rpart

partykit

doParallel

parallel

License
Include a license file (e.g., MIT, GPL) to specify usage terms.

