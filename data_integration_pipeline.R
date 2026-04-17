# ==========================================
# GDSC2 DATA INTEGRATION & CURATION PIPELINE
# Author: David Villafañe
# ==========================================

# ==========================================
# 1. LOAD LIBRARIES
# ==========================================
library(tidyverse)
library(readxl)
library(corrplot)

# ==========================================
# 2. LOAD DATA
# ==========================================
gdsc2 <- read_csv("data/GDSC2-dataset.csv")
compounds <- read_csv("data/Compounds-annotation.csv")
cell_lines <- read_excel("data/Cell_Lines_Details.xlsx")

# ==========================================
# 3. CLEAN COLUMN NAMES
# ==========================================
clean_names <- function(df) {
  names(df) <- names(df) %>%
    str_replace_all("[\r\n]", " ") %>%
    str_squish()
  return(df)
}

compounds <- clean_names(compounds)
cell_lines <- clean_names(cell_lines)

# ==========================================
# 4. PREPROCESS & MERGE DATASETS
# ==========================================
compounds_unique <- compounds %>%
  select(DRUG_NAME, TARGET, TARGET_PATHWAY) %>%
  distinct(DRUG_NAME, .keep_all = TRUE)

final_dataset <- gdsc2 %>%
  left_join(compounds_unique, by = "DRUG_NAME") %>%
  left_join(cell_lines, by = c("COSMIC_ID" = "COSMIC identifier")) %>%
  mutate(
    TARGET_FINAL  = coalesce(TARGET, PUTATIVE_TARGET),
    PATHWAY_FINAL = coalesce(TARGET_PATHWAY, PATHWAY_NAME)
  )

# ==========================================
# 5. INITIAL DATA QUALITY CHECK
# ==========================================
initial_rows <- nrow(final_dataset)
missing_ic50 <- sum(is.na(final_dataset$LN_IC50))

# ==========================================
# 6. CLEANING & DEDUPLICATION
# ==========================================
final_dataset <- final_dataset %>%
  filter(!is.na(LN_IC50)) %>%
  group_by(CELL_LINE_NAME, DRUG_NAME) %>%
  mutate(
    LN_IC50 = mean(LN_IC50, na.rm = TRUE),
    AUC = mean(AUC, na.rm = TRUE)
  ) %>%
  slice(1) %>%
  ungroup()

final_rows <- nrow(final_dataset)
duplicates_removed <- (initial_rows - missing_ic50) - final_rows

cat("\n--- DATA CLEANING SUMMARY ---\n")
cat("Initial rows:", initial_rows, "\n")
cat("Rows removed (missing LN_IC50):", missing_ic50, "\n")
cat("Duplicates removed:", duplicates_removed, "\n")
cat("Final dataset size:", final_rows, "\n")

# ==========================================
# 7. OUTLIER DETECTION
# ==========================================
final_dataset <- final_dataset %>%
  group_by(DRUG_NAME) %>%
  mutate(
    mean_drug = mean(LN_IC50, na.rm = TRUE),
    sd_drug   = sd(LN_IC50, na.rm = TRUE),
    Outlier_Flag = case_when(
      is.na(sd_drug) | sd_drug == 0 ~ "Normal",
      LN_IC50 > (mean_drug + 3 * sd_drug) ~ "Extreme Resistance",
      LN_IC50 < (mean_drug - 3 * sd_drug) ~ "Extreme Sensitivity",
      TRUE ~ "Normal"
    )
  ) %>%
  ungroup()

cat("\n--- OUTLIER SUMMARY ---\n")
print(table(final_dataset$Outlier_Flag))

# Top extreme sensitivities
cat("\n--- TOP 10 EXTREME SENSITIVITY CASES ---\n")
top_hits <- final_dataset %>%
  filter(Outlier_Flag == "Extreme Sensitivity") %>%
  select(
    CELL_LINE_NAME,
    Cancer_Type = `Cancer Type (matching TCGA label)`,
    DRUG_NAME,
    PATHWAY_FINAL,
    LN_IC50
  ) %>%
  arrange(LN_IC50) %>%
  head(10)

print(top_hits)

# ==========================================
# 8. MISSING DATA ANALYSIS
# ==========================================
missing_data <- as.data.frame(
  colSums(is.na(final_dataset)) / nrow(final_dataset) * 100
)
colnames(missing_data) <- "Percent_Missing"
missing_data$Variable <- rownames(missing_data)

ggplot(missing_data, aes(x = reorder(Variable, Percent_Missing), y = Percent_Missing)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Missing Data Overview", y = "% Missing", x = "Variables")

# ==========================================
# 9. CONTROL VISUALIZATIONS
# ==========================================

# A. Distribution
ggplot(final_dataset, aes(x = LN_IC50)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Distribution of LN_IC50")

# B. Cancer type distribution
final_dataset %>%
  count(`Cancer Type (matching TCGA label)`) %>%
  filter(!is.na(`Cancer Type (matching TCGA label)`)) %>%
  ggplot(aes(x = reorder(`Cancer Type (matching TCGA label)`, n), y = n)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Cell Line Distribution by Cancer Type")

# C. Correlation LN_IC50 vs AUC
ggplot(final_dataset, aes(x = LN_IC50, y = AUC)) +
  geom_point(alpha = 0.1, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(title = "LN_IC50 vs AUC Correlation")

# ==========================================
# 10. CORRELATION MATRIX
# ==========================================
num_vars <- final_dataset %>%
  select(LN_IC50, AUC, RMSE, Z_SCORE) %>%
  drop_na()

if (ncol(num_vars) > 1) {
  corrplot(cor(num_vars),
           method = "color",
           addCoef.col = "black",
           tl.col = "black",
           title = "\nCorrelation Matrix",
           mar = c(0,0,1,0))
}

# ==========================================
# 11. SAVE FINAL DATASET
# ==========================================
write_csv(final_dataset, "outputs/GDSC2_Final_Master_Clean.csv")

cat("\nDataset successfully saved!\n")
