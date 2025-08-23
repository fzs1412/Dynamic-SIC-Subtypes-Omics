https://github.com/fzs1412/Dynamic-SIC-Subtypes-Omics.git
# ==============================================================================
# Project: Dynamic SIC Subtypes & Omics
# Description: This script performs a comprehensive analysis of Sepsis-Induced
#              Coagulopathy (SIC) patient trajectories. It includes patient
#              grouping, differential gene expression (DEG) analysis,
#              prognostic gene selection using LASSO regression, and dynamic
#              expression analysis of key pathways.
 

library(dplyr)
library(data.table)
library(DESeq2)
library(BiocManager)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(FactoMineR)
library(factoextra)
library(EnhancedVolcano)
library(gbmt)
library(magrittr)
library(tidyr)
library(stringr)
library(lcmm)
library(ggplot2)
library(survival)
library(survminer)
library(mice)
library(CBCgrps)
library(flextable)
library(DataExplorer)
library(glmnet)
library(pROC)
library(caret)
library(ReactomePA)
library(ggvenn)
library(ggpubr)
library(grid)
library(gridExtra)
library(pheatmap) # pheatmap is used but not listed in the user's header, so adding it here.

# ==============================================================================
# 1. Data Preprocessing and Patient Trajectory Analysis
# ==============================================================================

# 1.1 Load and prepare clinical data for SIC score calculation
data <- read.csv("CMAISE1.5v_x0.csv", row.names = 1)

# Calculate SIC score based on clinical criteria
data <- data %>%
  mutate(plt_score = case_when(
    plt < 100 ~ 2,
    plt >= 100 & plt <= 150 ~ 1,
    TRUE ~ 0 # Platelets > 150 get a score of 0
  )) %>%
  mutate(inr_score = case_when(
    inr > 1.4 ~ 2,
    inr >= 1.2 & inr <= 1.4 ~ 1,
    TRUE ~ 0 # INR < 1.2 get a score of 0
  )) %>%
  mutate(sofa_score_comp = case_when( # Rename to differentiate from total SOFA score
    SOFA >= 2 ~ 2,
    SOFA == 1 ~ 1,
    TRUE ~ 0 # SOFA < 1 get a score of 0
  )) %>%
  mutate(SIC = plt_score + inr_score + sofa_score_comp) # Calculate total SIC Score
data <- data %>% select(-plt_score, -inr_score, -sofa_score_comp)

# 1.2 Define patient trajectories based on Day 1 and Day 3 SIC scores
dat <- data %>% select(PtID, Days, SIC)

# Step 1: Extract Day 1 and Day 3 data for each patient
day1_data <- dat %>% filter(Days == 1) %>% select(PtID, SIC_Day1 = SIC)
day3_data <- dat %>% filter(Days == 3) %>% select(PtID, SIC_Day3 = SIC)

# Step 2: Merge data and define trajectory types
trajectory <- full_join(day1_data, day3_data, by = "PtID") %>%
  mutate(
    trajectory_type = case_when(
      SIC_Day1 >= 4 & SIC_Day3 >= 4 ~ "persistent",
      SIC_Day1 < 4 & SIC_Day3 < 4 ~ "nonTP",
      SIC_Day1 < 4 & SIC_Day3 >= 4 ~ "worsening",
      SIC_Day1 >= 4 & SIC_Day3 < 4 ~ "transient",
      TRUE ~ NA_character_  # Handle missing data
    )
  )
table(trajectory$trajectory_type)

# Step 3: Merge classification results back to the original data
dat_with_trajectory <- data %>%
  left_join(select(trajectory, PtID, trajectory_type), by = "PtID")
fwrite(dat_with_trajectory, "coldata.csv")

# 1.3 Prepare expression and clinical data for DEG analysis
coldata <- read.csv("coldata.csv", row.names = 1) %>%
  filter(Days == 1) %>% rename(class = trajectory_type)
coldata$class <- as.factor(coldata$class)

countdata <- read.csv("OMIX006457-02.csv", row.names = 1)
countdata_zero <- countdata
countdata_zero[is.na(countdata_zero)] <- 0
countdata_zero <- countdata_zero[, !grepl("d3$|d5$", colnames(countdata_zero))]

common_samples <- intersect(colnames(countdata_zero), rownames(coldata))
print(paste("Matching samples:", length(common_samples)))
# Retain 394 samples
countdata_zero <- countdata_zero[, common_samples]
coldata <- coldata[common_samples, ]

table(coldata$class)
# Extract row names, remove pipe and everything after it to keep only ENSEMBL ID
rownames(countdata_zero) <- gsub("\\|.*$", "", rownames(countdata_zero))
rownames(countdata_zero) <- sub("^[^|]+\\|", "", rownames(countdata_zero))
countdata_zero[1:20, 1:2]

write.csv(coldata, file = "coldata394.csv")
write.csv(countdata_zero, file = "countdata394.csv")

# 1.4 Longitudinal Trajectory Modeling (lcmm)
set.seed(123)
data <- fread("CMAISE1.5v_x0.csv")
dat <- read.csv("coldata394.csv") %>% select(PtID)
data <- data %>% inner_join(dat, by = "PtID")

# Convert PtID to a factor and then to a numeric ID for lcmm
data$ID <- as.numeric(as.factor(data$PtID))
data$Days <- as.numeric(data$Days)

# Recalculate SIC Score based on the provided criteria
data <- data %>%
  mutate(plt_score = case_when(
    plt < 100 ~ 2,
    plt >= 100 & plt <= 150 ~ 1,
    TRUE ~ 0 # Platelets > 150 get a score of 0
  )) %>%
  mutate(inr_score = case_when(
    inr > 1.4 ~ 2,
    inr >= 1.2 & inr <= 1.4 ~ 1,
    TRUE ~ 0 # INR < 1.2 get a score of 0
  )) %>%
  mutate(sofa_score_comp = case_when( # Rename to differentiate from total SOFA score
    SOFA >= 2 ~ 2,
    SOFA == 1 ~ 1,
    TRUE ~ 0 # SOFA < 1 get a score of 0
  )) %>%
  mutate(SIC = plt_score + inr_score + sofa_score_comp) # Calculate total SIC Score
data <- data %>% select(-plt_score, -inr_score, -sofa_score_comp)

# Run lcmm models with 1 to 6 groups
m1 <- hlme(SIC ~ poly(Days, degree = 2, raw = TRUE), subject = 'ID', ng = 1, data = data)
m2 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(SIC ~ poly(Days, degree = 2, raw = TRUE), mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 2, data = data))
m3 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(SIC ~ poly(Days, degree = 2, raw = TRUE), mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 3, data = data))
m4 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(SIC ~ poly(Days, degree = 2, raw = TRUE), mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 4, data = data))
m5 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(SIC ~ poly(Days, degree = 2, raw = TRUE), mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 5, data = data))
m6 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(SIC ~ poly(Days, degree = 2, raw = TRUE), mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 6, data = data))

summarytable(m1, m2, m3, m4, m5, m6,
             which = c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

dir.create("Table", recursive = TRUE)

# Plot the chosen model (m4)
plot(m4, which = "fit", var.time = "Days", break.times = 5,
     bty = "l", ylab = "SIC_Score",
     xlab = "days after ICU admission", lwd = 2,
     marg = TRUE, legend = NULL, shades = TRUE)

# Customize the plot for publication
plot(m4, which = "fit", var.time = "Days",
     break.times = 3,
     bty = "l",
     ylab = "SIC Score",
     xlab = "Days after ICU admission",
     lwd = 2.5, # Thicken lines for publication quality (0.5-1.5 pt)
     marg = TRUE,
     shades = TRUE,
     legend = NULL,
     col = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"), # Use new colors
     lty = rep(1, 4)) # Use solid lines for all trajectories

# Add a professional legend
legend("topright",
       legend = c("Group 1: Low-Risk",
                  "Group 2: Very-High-Risk",
                  "Group 3: Worsening",
                  "Group 4: Resolving"),
       col = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"),
       lty = rep(1, 4), # Use solid lines
       lwd = 2.5,
       cex = 0.85, # Set font size to 8-10pt range
       bty = "n", # No border
       inset = c(0.02, 0.14)) # Fine-tune position to avoid overlap

# Save patient classification and trajectory data
dtclass <- m4$pprob[, 1:2]
fwrite(dtclass, "dtclass.csv")
datcom <- merge(data, dtclass, by = 'ID')
fwrite(datcom, "datcom.csv")
datcomday1 <- datcom %>% filter(Days == 1)
fwrite(datcomday1, "datcomday1.csv")

# ==============================================================================
# 2. Inter-Group Differential Expression Analysis (Class 2 vs 3)
# ==============================================================================

# Load and prepare data for analysis
coldata <- read.csv("guiji_coldata394.csv", row.names = 1) %>%
  filter(class == "2" | class == "3")
coldata$class <- as.factor(coldata$class)

countdata <- read.csv("OMIX006457-02.csv", row.names = 1)
countdata_zero <- countdata
countdata_zero[is.na(countdata_zero)] <- 0
countdata_zero <- countdata_zero[, !grepl("d3$|d5$", colnames(countdata_zero))]

common_samples <- intersect(colnames(countdata_zero), rownames(coldata))
print(paste("Matching samples:", length(common_samples)))
countdata_zero <- countdata_zero[, common_samples]
coldata <- coldata[common_samples, ]

# Create a gene ID map from the original row names (format: ENSEMBL|SYMBOL)
original_rownames <- rownames(countdata_zero)
gene_map <- data.frame(
  ENSEMBL = gsub("\\|.*$", "", original_rownames),
  SYMBOL = gsub(".*\\|", "", original_rownames),
  stringsAsFactors = FALSE
)
# Filter out rows without a SYMBOL and ensure unique ENSEMBL IDs
gene_map <- gene_map %>% filter(SYMBOL != "" & !duplicated(ENSEMBL))
head(gene_map)

# Prepare the expression matrix for DESeq2, using pure ENSEMBL IDs as row names
countdata_for_deseq <- countdata_zero
rownames(countdata_for_deseq) <- gsub("\\|.*$", "", rownames(countdata_for_deseq))
# Filter for genes present in the gene map
countdata_for_deseq <- countdata_for_deseq[rownames(countdata_for_deseq) %in% gene_map$ENSEMBL, ]

# 2.1 DESeq2 Differential Expression Analysis
dds <- DESeqDataSetFromMatrix(countData = countdata_for_deseq,
                              colData = coldata,
                              design = ~ class)
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)
res_df_for_volcano <- merge(res_df, gene_map, by = "ENSEMBL", all.x = TRUE)

# Identify top genes to label
res_df_for_volcano <- res_df_for_volcano[order(res_df_for_volcano$padj), ]
top_genes <- bind_rows(
  res_df_for_volcano %>% filter(log2FoldChange > 1) %>% head(10),
  res_df_for_volcano %>% filter(log2FoldChange < -1) %>% head(10)
)
genes_to_label <- top_genes$SYMBOL

# 2.2 Plot the Volcano Plot
pdf("Volcano_14class.pdf", width = 8, height = 7)
p_volcano <- EnhancedVolcano(res_df_for_volcano,
                             lab = res_df_for_volcano$SYMBOL,
                             selectLab = genes_to_label,
                             x = 'log2FoldChange',
                             y = 'padj',
                             xlab = bquote(~Log[2]~ 'Fold Change'),
                             ylab = bquote(~-Log[10]~ 'Adjusted P-value'),
                             axisLabSize = 14, pCutoff = 0.05, FCcutoff = 1.0,
                             pointSize = 2.0, labSize = 4.0, labFace = 'bold',
                             title = NULL, subtitle = NULL,
                             col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                             colAlpha = 0.8, legendPosition = "right", legendLabSize = 12,
                             legendIconSize = 5.0,
                             legendLabels = c('Not significant', 'Log₂FC', 'p-value', 'p-value & Log₂FC'),
                             drawConnectors = TRUE, widthConnectors = 0.5, colConnectors = 'grey50',
                             gridlines.major = FALSE, gridlines.minor = FALSE,
                             border = 'full', borderWidth = 1.5, borderColour = 'black')

# Move legend to the top-right corner inside the plot
p_volcano <- p_volcano + theme(legend.position = c(0.95, 0.45),
                               legend.justification = c("right", "top"),
                               legend.background = element_rect(fill = alpha("white", 0.6)))
print(p_volcano)
dev.off()

# 2.3 Plot the Heatmap
# Filter for significant differentially expressed genes (DEGs)
res_df <- as.data.frame(subset(res, padj < 0.05))
DEG <- res_df
logFC_cutoff = 1
# Define UP/DOWN/NOT regulated genes
DEG$change <- ifelse(
  DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
  ifelse(DEG$log2FoldChange > logFC_cutoff, 'UP', 'DOWN'),
  'NOT'
)
DEG$change <- factor(DEG$change, levels = c("UP", "DOWN", "NOT"))
table(DEG$change)

# Select top genes for the heatmap
choose_gene <- subset(DEG, DEG$change != 'NOT')
choose_gene <- choose_gene[order(choose_gene$log2FoldChange, decreasing = TRUE), ]
choose_gene <- choose_gene[c(head(rownames(choose_gene), 40), tail(rownames(choose_gene), 40)), ]

# Get variance-stabilized transformed (VST) expression matrix
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)
choose_matrix <- vst_matrix[rownames(choose_gene), ]

# Change: Convert row names of the heatmap matrix from ENSEMBL ID to Gene Symbol
symbol_rownames <- gene_map$SYMBOL[match(rownames(choose_matrix), gene_map$ENSEMBL)]
# Handle potential NAs or duplicated Symbols to ensure unique row names
symbol_rownames[is.na(symbol_rownames)] <- rownames(choose_matrix)[is.na(symbol_rownames)]
rownames(choose_matrix) <- make.unique(symbol_rownames)

# Prepare annotation data
Group <- data.frame(Group = coldata$class, row.names = colnames(vst_matrix))
Up_Down <- data.frame(Up_Down = choose_gene$change, row.names = rownames(choose_matrix))

# Plot the heatmap
pdf("heatmap_14class.pdf", width = 10, height = 12)
# Replace group labels for better readability
Group$Group <- factor(ifelse(Group$Group == "2", "Low-Risk", "Severe"))
pheatmap_plot <- pheatmap(
  choose_matrix,
  scale = "row",
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_colnames = FALSE, show_rownames = TRUE, # Display row names (Gene Symbol)
  fontsize_row = 8,
  annotation_col = Group,
  annotation_colors = list(
    Group = c("Low-Risk" = "#e23d30", "Severe" = "#00dae0"),
    Up_Down = c("UP" = "#ff9289", "DOWN" = "#96ca00")
  ),
  annotation_row = Up_Down,
  border_color = 'white',
  main = paste0("Pheatmap ", length(rownames(choose_matrix)), " DEGs"),
  color = colorRampPalette(c(rep("blue", 2), "white", rep("red", 2)))(50)
)
print(pheatmap_plot)
dev.off()

# 2.4 Combine Volcano Plot and Heatmap
p_heatmap <- grid.grabExpr(draw(pheatmap_plot))
pdf("combined_volcano_heatmap.pdf", width = 16, height = 10)
grid.arrange(
  arrangeGrob(p_volcano, top = textGrob("A", x = unit(0, "npc"), y = unit(1, "npc"),
                                        just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))),
  arrangeGrob(p_heatmap, top = textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"),
                                        just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))),
  ncol = 2, # Arrange plots horizontally
  widths = c(1, 1.2)
)
dev.off()

# ==============================================================================
# 3. LASSO Regression for Prognostic Gene Selection
# ==============================================================================

# 3.1 Data preparation
# This section uses an older approach. Re-using the DESeq2 object for consistency.
dds <- DESeqDataSetFromMatrix(countData = countdata_for_deseq,
                              colData = coldata,
                              design = ~ class)
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(subset(res, padj < 0.05))
vsd <- vst(dds)

all_degs <- rownames(res_df)
X <- t(assay(vsd)[all_degs, ])
# Remove genes with near-zero variance
nzv <- nearZeroVar(X, saveMetrics = TRUE)
non_zero_var_genes <- rownames(nzv[nzv$nzv == FALSE, ])
X_filtered <- X[, non_zero_var_genes]
cat("Initial candidate gene count:", length(all_degs), "\n")
cat("Remaining genes after removing zero/near-zero variance:", ncol(X_filtered), "\n")

# Prepare the response variable y (patient group information)
y <- coldata$class
# Ensure sample order is consistent between X and y
stopifnot(identical(rownames(X_filtered), rownames(coldata)))

# Partition data into training and testing sets
set.seed(42) # Set random seed for reproducibility
train_indices <- createDataPartition(y, p = 0.7, list = FALSE, times = 1)
X_train <- X_filtered[train_indices, ]
y_train <- y[train_indices]
X_test <- X_filtered[-train_indices, ]
y_test <- y[-train_indices]
cat("Training sample count:", nrow(X_train), "\n")
cat("Testing sample count:", nrow(X_test), "\n")

# 3.2 Perform LASSO regression
X_train_matrix <- as.matrix(X_train)
y_train_numeric <- ifelse(y_train == "2", 1, 0)
y_test_numeric <- ifelse(y_test == "2", 1, 0)

# alpha = 1 for LASSO regression
cat("\n--- Running LASSO regression with cross-validation ---\n")
set.seed(42)
cv_lasso_model <- cv.glmnet(X_train_matrix, y_train_numeric,
                            family = "binomial",
                            alpha = 1,
                            type.measure = "auc")
# Plot the cross-validation curve
pdf("LASSO_CV_Curve.pdf", width = 8, height = 6)
plot(cv_lasso_model)
dev.off()
# Get the optimal lambda value
best_lambda <- cv_lasso_model$lambda.min * 1.5 # Or cv_lasso_model$lambda.1se
cat("Optimal Lambda (lambda.min):", best_lambda, "\n")

# 3.3 Extract feature genes and evaluate model performance
lasso_coeffs <- coef(cv_lasso_model, s = best_lambda)
lasso_selected_genes <- lasso_coeffs@Dimnames[[1]][which(lasso_coeffs != 0)]
# The first coefficient is the intercept and should be removed
lasso_selected_genes <- lasso_selected_genes[-1]
cat("\nNumber of optimal genes selected by LASSO:", length(lasso_selected_genes), "\n")

# Get gene symbols for the selected genes
selected_genes_df <- data.frame(ENSEMBL = lasso_selected_genes)
selected_genes_with_names <- merge(selected_genes_df, gene_map, by = "ENSEMBL")

# Evaluate model performance on the independent test set
cat("\n--- Evaluating model performance on the independent test set ---\n")
X_test_matrix <- as.matrix(X_test)
predictions_prob <- predict(cv_lasso_model, newx = X_test_matrix,
                            s = best_lambda, type = "response")
# Calculate and plot the ROC curve
roc_curve <- roc(response = y_test_numeric, predictor = as.vector(predictions_prob),
                 levels = c("0", "1"))
auc_value <- auc(roc_curve)
cat("\nAUC value of the model on the test set:", round(auc_value, 4), "\n")

# Plot the ROC curve
pdf("LASSO_ROC_Curve.pdf", width = 7, height = 7)
print(
  ggroc(roc_curve, colour = 'darkred', size = 1.5) +
    ggtitle(paste0("ROC Curve for LASSO Model (AUC = ", round(auc_value, 3), ")")) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "darkgrey") +
    theme_minimal() +
    annotate("text", x = 0.5, y = 0.25,
             label = paste(length(lasso_selected_genes), "genes selected"))
)
dev.off()

# ==============================================================================
# 4. SIC Status Transition Analysis (Sankey Plot)
# ==============================================================================

# Load data
data <- fread("datcom.csv")

# 4.1 Calculate the proportion of patients with SIC at each time point
sic_threshold <- 4
sic_proportions <- data %>%
  filter(Days %in% c(1, 3, 5)) %>% # Only consider Day 1, 3, 5
  group_by(Days) %>%
  summarise(
    total_patients = n_distinct(PtID),
    sic_patients = n_distinct(PtID[SIC >= sic_threshold], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(sic_proportion = sic_patients / total_patients)
cat("--- Proportion of SIC Patients on Day 1, 3, and 5 ---\n")
print(sic_proportions)

# 4.2 Visualize the proportion of SIC patients over time
sic_proportions$Days <- factor(sic_proportions$Days, levels = c(1, 3, 5))
plot_proportion <- ggplot(sic_proportions, aes(x = Days, y = sic_proportion, fill = Days)) +
  geom_bar(stat = "identity") +
  labs(x = "Days after ICU admission", y = "Proportion of patients with SIC (SIC Score >= 4)",
       title = "Proportion of SIC Patients on Day 1, 3, and 5") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(legend.position = "none")
print(plot_proportion)

# 4.3 Calculate patient SIC status transitions
data_sic_status <- data %>%
  filter(Days %in% c(1, 3, 5)) %>%
  mutate(sic_status = if_else(SIC >= sic_threshold, 1, 0)) %>%
  select(PtID, Days, sic_status) %>%
  pivot_wider(names_from = Days, values_from = sic_status, names_prefix = "Day")

status_labels <- c("No SIC", "SIC")

transitions_d1_d5 <- data_sic_status %>%
  filter(!is.na(Day1) & !is.na(Day3) & !is.na(Day5)) %>%
  group_by(Day1, Day3, Day5) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(
    status_d1 = status_labels[Day1 + 1],
    status_d3 = status_labels[Day3 + 1],
    status_d5 = status_labels[Day5 + 1],
    transition_path = paste(status_d1, "->", status_d3, "->", status_d5)
  ) %>%
  mutate(proportion = count / sum(count))
cat("\n--- SIC Status Transition Paths from Day 1 to Day 5 ---\n")
print(transitions_d1_d5)

# 4.4 Plot the Sankey diagram
data_for_alluvial <- data_sic_status %>%
  select(PtID, Day1, Day3, Day5) %>%
  pivot_longer(cols = starts_with("Day"), names_to = "Time", values_to = "Status") %>%
  filter(!is.na(Status)) %>%
  mutate(
    Time = factor(Time, levels = c("Day1", "Day3", "Day5")),
    Status = factor(status_labels[Status + 1], levels = status_labels)
  )

plot_alluvial <- ggplot(data_for_alluvial,
                        aes(x = Time, stratum = Status, alluvium = PtID,
                            fill = Status, label = Status)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .9) +
  geom_text(stat = "stratum", size = 3, color = "black") +
  labs(title = "SIC Status Transitions from Day 1 to Day 5",
       x = "Days after ICU admission",
       y = "Number of Patients") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(plot_alluvial)

# ==============================================================================
# 5. Patient Characteristics and Outcome Tables
# ==============================================================================

# Impute missing data using MICE
dat <- fread("datcomday1.csv")
aa <- mice(dat, seed = 123)
datmice <- complete(aa, action = 3)
datmice$infectionSite_SD <- ifelse(datmice$infectionSite_SD == 3 | datmice$infectionSite_SD == 7 | datmice$infectionSite_SD == 8, 6, datmice$infectionSite_SD)
fwrite(datmice, "datmice.csv")

# 5.1 Create Table 1 (Characteristics)
dat <- fread("datmice.csv") %>%
  select(class, age, sex, diabete, hyperten, myoinfarc, cardiofailure,
         SOFA, fluidin, fluidout, urine, hrmax, hrmin, mapmax, mapmin,
         sapmax, sapmin, rrmax, rrmin, tmax, tmin, wbc, hct, plt, pha,
         paco, pao, lac, pf, cr, crp)
skewvar2 <- c("lac", "paco", "pao", "pha", "cr", "crp", "urine", "wbc",
              "fluidin", "fluidout")
tab1 <- multigrps(dat, gvar = "class", norm.rd = 1, cat.rd = 1, sk.rd = 1, minfactorlevels = 10, skewvar = skewvar2)
colnames(tab1) <- c("Variable", "Group 1", "Group 2", "Group 3", "Group 4", "p-value", "Method")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1)), path = "tab1.docx")

# 5.2 Create Table for Outcomes
dat <- fread("datmice.csv") %>%
  select(class, infectionSite_SD, mort, MV_days, CRRT_days, VASO_days, Hospital_days, SIC)
tab_outcome <- multigrps(dat, gvar = "class", norm.rd = 1, cat.rd = 1, sk.rd = 1, minfactorlevels = 10,
                         sim = TRUE, workspace = 200000)
colnames(tab_outcome) <- c("Variable", "Group 1", "Group 2", "Group 3", "Group 4", "p-value", "Method")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab_outcome)), path = "tab_outcom2.docx")

# 5.3 Survival Analysis and Cox Regression
dd <- fread("datmice.csv") %>%
  mutate(Hospital_days = ifelse(Hospital_days > 30, 30, Hospital_days))
fit <- survfit(Surv(Hospital_days, mort) ~ class, dd)

ggsurvplot(fit,
           xlab = 'Days',
           pval = TRUE,
           pval.size = 3,
           pval.coord = c(0, 0.6),
           legend.labs = c('class1', 'class2', 'class3', 'class4'),
           surv.median.line = "hv", # Display median survival time
           ggtheme = theme_bw(), # Set ggplot2 theme
           palette = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"),
           risk.table = TRUE,
           risk.table.col = "strata",
           ylim = c(0.5, 1))

dd$class <- ifelse(dd$class == 1, "class1", ifelse(dd$class == 2, "class2", ifelse(dd$class == 3, "class3", "class4")))

# Cox regression models
cox_mode1 <- coxph(Surv(Hospital_days, mort) ~ class, dd)
tbl_regression(cox_mode1, exponentiate = TRUE)
cox_mode2 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex, dd)
tbl_regression(cox_mode2, exponentiate = TRUE)
cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex + hyperten + diabete +
                     SOFA + copd + lac + pha + pf + cr + crp + wbc + hct + plt + urine, dd)
tbl_regression(cox_mode3, exponentiate = TRUE)
save_as_docx(as_flex_table(tbl_regression(cox_mode1, exponentiate = TRUE)), path = "cox1.docx")
save_as_docx(as_flex_table(tbl_regression(cox_mode2, exponentiate = TRUE)), path = "cox2.docx")
save_as_docx(as_flex_table(tbl_regression(cox_mode3, exponentiate = TRUE)), path = "cox3.docx")

# ==============================================================================
# 6. Intra-Group Paired Differential Expression Analysis (Day 1 vs Day 5)
# ==============================================================================

# Helper function to prepare data for paired analysis
prepare_paired_data <- function(class_label) {
  coldata_raw <- read.csv("guiji_coldata394.csv", row.names = 1) %>%
    filter(class == class_label)
  
  countdata_raw <- read.csv("OMIX006457-02.csv", row.names = 1)
  countdata_zero <- countdata_raw
  countdata_zero[is.na(countdata_zero)] <- 0
  
  countdata_d1d5 <- countdata_zero[, grepl("d1$|d5$", colnames(countdata_zero))]
  count_samples <- colnames(countdata_d1d5)
  count_samples_base <- sub("_(d1|d5)$", "", count_samples)
  sample_info <- data.frame(full_name = count_samples, base_name = count_samples_base)
  
  both_timepoints <- sample_info %>%
    group_by(base_name) %>%
    filter(all(c("d1", "d5") %in% sub(".*_", "", full_name))) %>%
    pull(base_name) %>% unique()
  
  coldata_samples_base <- sub("_(d1|d5)$", "", rownames(coldata_raw))
  final_base_names <- intersect(both_timepoints, coldata_samples_base)
  
  final_samples <- sample_info %>%
    filter(base_name %in% final_base_names) %>%
    pull(full_name)
  
  countdata_final <- countdata_d1d5[, final_samples, drop = FALSE]
  
  coldata_matched <- data.frame()
  for (base in final_base_names) {
    coldata_row <- coldata_raw[grepl(paste0("^", base, "_"), rownames(coldata_raw)), , drop = FALSE]
    if (nrow(coldata_row) > 0) {
      coldata_d1 <- coldata_row[1, , drop = FALSE]
      rownames(coldata_d1) <- paste0(base, "_d1")
      coldata_d5 <- coldata_row[1, , drop = FALSE]
      rownames(coldata_d5) <- paste0(base, "_d5")
      coldata_matched <- rbind(coldata_matched, coldata_d1, coldata_d5)
    }
  }
  
  coldata_final <- coldata_matched[colnames(countdata_final), , drop = FALSE]
  
  # Prepare for DESeq2
  coldata_final$sample_id <- rownames(coldata_final)
  coldata_final$patient <- sub("_(d1|d5)$", "", coldata_final$sample_id)
  coldata_final$time <- ifelse(grepl("_d1$", coldata_final$sample_id), "d1", "d5")
  coldata_final$time <- factor(coldata_final$time, levels = c("d1", "d5"))
  
  countdata_for_deseq <- countdata_final
  rownames(countdata_for_deseq) <- gsub("\\|.*$", "", rownames(countdata_for_deseq))
  
  gene_map <- data.frame(
    ENSEMBL = gsub("\\|.*$", "", rownames(countdata_final)),
    SYMBOL = gsub(".*\\|", "", rownames(countdata_final)),
    stringsAsFactors = FALSE
  )
  gene_map <- gene_map %>% filter(SYMBOL != "" & !duplicated(ENSEMBL))
  countdata_for_deseq <- countdata_for_deseq[rownames(countdata_for_deseq) %in% gene_map$ENSEMBL, ]
  
  dds <- DESeqDataSetFromMatrix(countData = countdata_for_deseq,
                                colData = coldata_final,
                                design = ~ patient + time)
  dds <- DESeq(dds)
  res <- results(dds)
  res_df_for_volcano <- as.data.frame(res) %>%
    mutate(ENSEMBL = rownames(.), .before = 1) %>%
    merge(gene_map, by = "ENSEMBL", all.x = TRUE)
  
  return(list(res_df = res_df_for_volcano, dds = dds))
}

# --- Analyze "Persistent" (Class 2) group ---
cat("\n--- Running paired analysis for 'Persistent' (Class 2) group ---\n")
class2_results <- prepare_paired_data("2")
results_class2 <- class2_results$res_df
dds_class2 <- class2_results$dds

# --- Analyze "Resolving" (Class 4) group ---
cat("\n--- Running paired analysis for 'Resolving' (Class 4) group ---\n")
class4_results <- prepare_paired_data("4")
results_class4 <- class4_results$res_df
dds_class4 <- class4_results$dds

# ==============================================================================
# 7. Comparison and Causal Inference
# ==============================================================================

# 7.1 Extract "Resolving" and "Persistent" gene sets
P_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1.0

res_resolving <- results_class4 %>% filter(!is.na(padj))
resolving_up <- res_resolving %>% filter(padj < P_CUTOFF & log2FoldChange > LOGFC_CUTOFF) %>% pull(SYMBOL)
resolving_down <- res_resolving %>% filter(padj < P_CUTOFF & log2FoldChange < -LOGFC_CUTOFF) %>% pull(SYMBOL)
cat(sprintf("Resolving gene set: %d UP-regulated genes, %d DOWN-regulated genes\n", length(resolving_up), length(resolving_down)))

res_persistent <- results_class2 %>% filter(!is.na(padj))
persistent_up <- res_persistent %>% filter(padj < P_CUTOFF & log2FoldChange > LOGFC_CUTOFF) %>% pull(SYMBOL)
persistent_down <- res_persistent %>% filter(padj < P_CUTOFF & log2FoldChange < -LOGFC_CUTOFF) %>% pull(SYMBOL)
cat(sprintf("Persistent gene set: %d UP-regulated genes, %d DOWN-regulated genes\n", length(persistent_up), length(persistent_down)))

# 7.2 Visualize Gene Set Comparison: Venn Diagram
gene_lists <- list(
  `Resolving UP` = resolving_up,
  `Persistent UP` = persistent_up,
  `Resolving DOWN` = resolving_down,
  `Persistent DOWN` = persistent_down
)
pdf("Venn_Comparison_Signatures.pdf", width = 10, height = 8)
print(
  ggvenn(
    gene_lists,
    columns = c("Resolving UP", "Persistent UP", "Resolving DOWN", "Persistent DOWN"),
    fill_color = c("#ff9289", "#fbc1bb", "#96ca00", "#cce083"),
    stroke_size = 0.5, set_name_size = 4, text_size = 3
  ) + labs(title = "High-Risk Resolving vs. Persistently Severe: Differential Gene Expression")
)
dev.off()
cat("Venn diagram saved to Venn_Comparison_Signatures.pdf\n")

# 7.3 Identify Key Pathways: Pathway Enrichment Analysis
run_enrichment <- function(gene_symbol_list, title) {
  if (length(gene_symbol_list) == 0) {
    cat(paste("Warning: No genes available for enrichment in '", title, "', skipping.\n"))
    return(NULL)
  }
  entrez_ids <- bitr(gene_symbol_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  
  if (nrow(entrez_ids) == 0) {
    cat(paste("Warning: No genes in '", title, "' could be mapped to ENTREZID, skipping.\n"))
    return(NULL)
  }
  
  kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID,
                            organism = 'hsa',
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
  
  if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
    p_kegg <- barplot(kegg_enrich, showCategory = 20) + labs(title = paste("KEGG Enrichment:", title))
    ggsave(paste0("KEGG_Enrichment_", gsub("[ /]", "_", title), ".pdf"), p_kegg, width = 10, height = 8)
    cat(paste("KEGG enrichment plot saved:", paste0("KEGG_Enrichment_", gsub("[ /]", "_", title), ".pdf\n")))
  } else {
    cat(paste("No significant KEGG pathways found in '", title, "'.\n"))
  }
  return(kegg_enrich)
}

genes_resolving_only_down <- setdiff(resolving_down, persistent_down)
genes_resolving_only_up <- setdiff(resolving_up, persistent_up)
cat("\n--- Running pathway enrichment on key gene subsets ---\n")
enrich_resolving_down <- run_enrichment(genes_resolving_only_down, "Resolving Only DOWN")
enrich_resolving_up <- run_enrichment(genes_resolving_only_up, "Resolving Only UP")

# ==============================================================================
# 8. Final Figure Generation (Plots A and B)
# ==============================================================================

# 8.1 Preparation: Create a global DDS object with all Day 1 samples
cat("\n--- Preparing global DDS object with all Day 1 samples ---\n")
coldata_d1_full <- fread("guiji_coldata394.csv") %>%
  tibble::column_to_rownames("SampleName")
countdata_d1_full <- fread("countdata394.csv", header = TRUE) %>%
  rename(V1 = 1) %>% tibble::column_to_rownames("V1")

coldata_d1_full <- coldata_d1_full %>%
  filter(rownames(.) %in% colnames(countdata_d1_full))
countdata_d1_full <- countdata_d1_full[, rownames(coldata_d1_full)]

gene_map_full <- data.frame(
  ENSEMBL = gsub("\\|.*$", "", rownames(countdata_d1_full)),
  SYMBOL = gsub(".*\\|", "", rownames(countdata_d1_full)),
  stringsAsFactors = FALSE
) %>% filter(SYMBOL != "" & !duplicated(ENSEMBL))

countdata_for_deseq_full <- countdata_d1_full
rownames(countdata_for_deseq_full) <- gsub("\\|.*$", "", rownames(countdata_for_deseq_full))
countdata_for_deseq_full <- countdata_for_deseq_full[rownames(countdata_for_deseq_full) %in% gene_map_full$ENSEMBL, ]

dds_full_d1 <- DESeqDataSetFromMatrix(countData = countdata_for_deseq_full,
                                      colData = coldata_d1_full,
                                      design = ~ 1)
dds_full_d1 <- DESeq(dds_full_d1)
cat("--- Global DDS object created ---\n")

# 8.2 Plot A: Prognostic Value of a Key Gene
plot_prognostic_gene <- function(gene_symbol, dds_obj, coldata_obj, gene_map_obj) {
  cat(paste("\n--- Generating Plot A: Prognostic value of gene", gene_symbol, "---\n"))
  
  high_risk_coldata <- coldata_obj %>% filter(SIC >= 4)
  if (nrow(high_risk_coldata) == 0) {
    cat("Error: No samples found in the high-risk group. Skipping plot.\n")
    return(NULL)
  }
  
  ensembl_id <- gene_map_obj %>% filter(SYMBOL == gene_symbol) %>% pull(ENSEMBL)
  if (length(ensembl_id) == 0 || !ensembl_id %in% rownames(dds_obj)) {
    cat(paste("Error: Gene", gene_symbol, "not found in the DDS object.\n"))
    return(NULL)
  }
  
  counts_data <- plotCounts(dds_obj, gene = ensembl_id, intgroup = "mort", returnData = TRUE)
  plot_data <- counts_data %>%
    filter(rownames(.) %in% rownames(high_risk_coldata)) %>%
    mutate(Outcome = factor(ifelse(mort == 1, "Deceased", "Survived"), levels = c("Survived", "Deceased")))
  
  p <- ggplot(plot_data, aes(x = Outcome, y = count, fill = Outcome)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
    stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y.npc = 0.9) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_manual(values = c("Survived" = "#00BFC4", "Deceased" = "#F8766D")) +
    labs(
      title = paste("Prognostic Value of", gene_symbol, "at Admission"),
      subtitle = "in High-Risk SIC Patients (SIC Score >= 4)",
      x = "In-Hospital Outcome",
      y = "Normalized Expression (log10 scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none")
  
  ggsave(paste0("Prognostic_Plot_A_", gene_symbol, ".pdf"), p, width = 6, height = 7)
  cat(paste("Plot A saved as: Prognostic_Plot_A_", gene_symbol, ".pdf\n"))
  print(p)
}

# You can replace "DBNDD1" with a key gene from your LASSO model results
plot_prognostic_gene(gene_symbol = "DBNDD1", dds_obj = dds_full_d1, coldata_obj = coldata_d1_full, gene_map_obj = gene_map_full)

# 8.3 Plot B: Dynamic Expression of a Key Gene
plot_gene_expression <- function(gene_symbol, dds_class2, dds_class4, gene_map) {
  cat(paste("\n--- Generating Plot B: Dynamic expression of gene", gene_symbol, "---\n"))
  
  ensembl_id <- gene_map %>% filter(SYMBOL == gene_symbol) %>% pull(ENSEMBL)
  if (length(ensembl_id) == 0 || !ensembl_id %in% rownames(dds_class2) || !ensembl_id %in% rownames(dds_class4)) {
    cat(paste("Error: Gene", gene_symbol, "not found in one or both DDS objects. Skipping plot.\n"))
    return(NULL)
  }
  
  counts_class2 <- plotCounts(dds_class2, gene = ensembl_id, intgroup = c("time", "patient"), returnData = TRUE)
  counts_class2$group <- "Persistent"
  counts_class4 <- plotCounts(dds_class4, gene = ensembl_id, intgroup = c("time", "patient"), returnData = TRUE)
  counts_class4$group <- "Resolving"
  
  plot_data <- rbind(counts_class2, counts_class4)
  p <- ggplot(plot_data, aes(x = time, y = count, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 2) +
    facet_wrap(~group, scales = "free_x") +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_manual(values = c("Persistent" = "#e23d30", "Resolving" = "#7b61ff")) +
    labs(
      title = paste("Dynamic Expression of Gene", gene_symbol),
      x = "Time Point",
      y = "Normalized Expression (log10 scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold"))
  
  ggsave(paste0("Boxplot_", gene_symbol, "_Dynamics.pdf"), p, width = 8, height = 6)
  cat(paste("Boxplot for gene", gene_symbol, "dynamics saved.\n"))
  print(p)
}

# The user's code suggests two key genes for Plot B: VSIG4 and ITGB3.
# Let's run the function for both.
plot_gene_expression(gene_symbol = "VSIG4", dds_class2 = dds_class2, dds_class4 = dds_class4, gene_map = gene_map_full)
plot_gene_expression(gene_symbol = "ITGB3", dds_class2 = dds_class2, dds_class4 = dds_class4, gene_map = gene_map_full)

# End of script
# ==============================================================================
```