
# ============================================================
# GLAUCIDIUM VOCAL DISCRIMINATION
# Linear Discriminant Analysis (LDA)
#
# Includes:
#   1. Data preparation
#   2. Stratified training/test split
#   3. LDA model fitting
#   4. Classification accuracy
#   5. Leave-one-out cross-validation
#   6. LD1-LD2 visualization
#   7. Export of results
# ============================================================

install.packages("caret")

library(MASS)
library(dplyr)
library(ggplot2)
library(caret)
library(psych)

# ============================================================
# 1. INPUT AND OUTPUT
# ============================================================

graphics.off()

input_file <- "nopo_only_lda.csv"

out_dir <- "LDA_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# ============================================================
# 2. READ DATA
# ============================================================

LDA_dat <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# Vocal variables
vocal_vars <- c(
  "avgfreq",
  "avgbw",
  "avgdur",
  "avginterval",
  "mininterval",
  "variance",
  "maxconsec",
  "pctsteady"
)

# Retain required columns
dat <- LDA_dat %>%
  select(taxonlump, all_of(vocal_vars)) %>%
  filter(complete.cases(.)) %>%
  mutate(taxonlump = factor(taxonlump)) %>%
  droplevels()

cat("\nSample sizes by taxon:\n")
print(table(dat$taxonlump))

# ============================================================
# 3. DEFINE TAXON COLORS
# ============================================================

custom_colors <- c(
  "brasilianum" = "#CC5500",
  "californicum" = "#D2B48C",
  "cobanense" = "#78866b",
  "gnoma" = "#964B00",
  "pinicola" = "#FFB6C1",
  "swarthi" = "#87CEEB"
)

# Verify color assignments
missing_colors <- setdiff(
  levels(dat$taxonlump),
  names(custom_colors)
)

if (length(missing_colors) > 0) {
  stop(
    "Missing colors for: ",
    paste(missing_colors, collapse = ", ")
  )
}

custom_colors <- custom_colors[levels(dat$taxonlump)]

# ============================================================
# 4. STRATIFIED TRAINING/TEST SPLIT
# ============================================================

# Approximately 60% training / 40% testing
# within each taxon.

train_index <- unlist(
  lapply(split(seq_len(nrow(dat)), dat$taxonlump),
         function(idx) {
           n_train <- max(
             1,
             min(length(idx) - 1,
                 round(length(idx) * 0.60))
           )
           sample(idx, n_train)
         }),
  use.names = FALSE
)

training_raw <- dat[train_index, ]
testing_raw  <- dat[-train_index, ]

training_raw <- droplevels(training_raw)
testing_raw$taxonlump <- factor(
  testing_raw$taxonlump,
  levels = levels(training_raw$taxonlump)
)

cat("\nTraining samples:\n")
print(table(training_raw$taxonlump))

cat("\nTesting samples:\n")
print(table(testing_raw$taxonlump))

# ============================================================
# 5. STANDARDIZE VOCAL VARIABLES
# ============================================================

# Estimate means and SDs using training data ONLY.
# Apply the same transformation to testing data.

train_means <- sapply(
  training_raw[vocal_vars], mean
)

train_sds <- sapply(
  training_raw[vocal_vars], sd
)

if (any(!is.finite(train_sds) | train_sds == 0)) {
  stop("At least one predictor has zero or invalid SD.")
}

training <- training_raw
testing <- testing_raw

training[vocal_vars] <- scale(
  training_raw[vocal_vars],
  center = train_means,
  scale = train_sds
)

testing[vocal_vars] <- scale(
  testing_raw[vocal_vars],
  center = train_means,
  scale = train_sds
)

# ============================================================
# 6. EXAMINE CORRELATIONS
# ============================================================

pdf(
  file.path(out_dir, "vocal_correlations.pdf"),
  width = 10,
  height = 10
)

pairs.panels(
  training[vocal_vars],
  gap = 0,
  pch = 21,
  bg = custom_colors[as.character(training$taxonlump)]
)

dev.off()

# ============================================================
# 7. FIT LDA MODEL
# ============================================================

linear <- lda(
  taxonlump ~ .,
  data = training,
  prior = rep(
    1 / nlevels(training$taxonlump),
    nlevels(training$taxonlump)
  )
)

print(linear)

# Proportion of between-group discriminant variance
ld_variance <- linear$svd^2 /
  sum(linear$svd^2)

cat("\nDiscriminant variance proportions:\n")
print(ld_variance)

# Export LDA coefficients
write.csv(
  linear$scaling,
  file.path(out_dir, "LDA_coefficients.csv")
)

# ============================================================
# 8. TRAINING CLASSIFICATION
# ============================================================

train_pred <- predict(linear, training)

train_confusion <- confusionMatrix(
  data = train_pred$class,
  reference = training$taxonlump
)

cat("\nTRAINING CLASSIFICATION:\n")
print(train_confusion)

# ============================================================
# 9. TEST CLASSIFICATION
# ============================================================

test_pred <- predict(linear, testing)

test_confusion <- confusionMatrix(
  data = test_pred$class,
  reference = testing$taxonlump
)

cat("\nTEST CLASSIFICATION:\n")
print(test_confusion)

# ============================================================
# 10. LEAVE-ONE-OUT CROSS-VALIDATION
# ============================================================

# Cross-validation performed on training data.
# Same class priors as main model.

cv_model <- lda(
  taxonlump ~ .,
  data = training,
  prior = rep(
    1 / nlevels(training$taxonlump),
    nlevels(training$taxonlump)
  ),
  CV = TRUE
)

cv_confusion <- confusionMatrix(
  data = cv_model$class,
  reference = training$taxonlump
)

cat("\nLOOCV CLASSIFICATION:\n")
print(cv_confusion)

# ============================================================
# 11. SUMMARIZE ACCURACY
# ============================================================

accuracy_results <- data.frame(
  Dataset = c(
    "Training",
    "Testing",
    "LOOCV"
  ),
  Accuracy = c(
    as.numeric(train_confusion$overall["Accuracy"]),
    as.numeric(test_confusion$overall["Accuracy"]),
    as.numeric(cv_confusion$overall["Accuracy"])
  )
)

accuracy_results$Percent <- 100 * accuracy_results$Accuracy

print(accuracy_results)

write.csv(
  accuracy_results,
  file.path(out_dir, "LDA_accuracy.csv"),
  row.names = FALSE
)

# ============================================================
# 12. EXTRACT DISCRIMINANT SCORES
# ============================================================

train_scores <- as.data.frame(train_pred$x)
train_scores$taxonlump <- training$taxonlump
train_scores$Dataset <- "Training"

test_scores <- as.data.frame(test_pred$x)
test_scores$taxonlump <- testing$taxonlump
test_scores$Dataset <- "Testing"

all_scores <- bind_rows(
  train_scores,
  test_scores
)

write.csv(
  all_scores,
  file.path(out_dir, "LDA_scores.csv"),
  row.names = FALSE
)

# ============================================================
# 13. LD1-LD2 PLOT
# ============================================================

ld1_pct <- round(ld_variance[1] * 100, 1)
ld2_pct <- round(ld_variance[2] * 100, 1)

p_lda <- ggplot(
  train_scores,
  aes(
    x = LD1,
    y = LD2,
    color = taxonlump
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  stat_ellipse(
    aes(group = taxonlump),
    level = 0.95,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = function(x) paste("G. g.", x)
  ) +
  labs(
    title = "Vocal differentiation among Glaucidium taxa",
    x = paste0("LD1 (", ld1_pct, "%)"),
    y = paste0("LD2 (", ld2_pct, "%)"),
    color = "Taxon"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

print(p_lda)

ggsave(
  file.path(out_dir, "LDA_LD1_LD2.pdf"),
  p_lda,
  width = 9,
  height = 7
)

ggsave(
  file.path(out_dir, "LDA_LD1_LD2.png"),
  p_lda,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 14. TRAINING VS TESTING PLOT
# ============================================================

p_test <- ggplot(
  all_scores,
  aes(
    x = LD1,
    y = LD2,
    color = taxonlump,
    shape = Dataset
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = function(x) paste("G. g.", x)
  ) +
  scale_shape_manual(
    values = c(
      "Training" = 16,
      "Testing" = 17
    )
  ) +
  labs(
    title = "LDA training and testing samples",
    x = paste0("LD1 (", ld1_pct, "%)"),
    y = paste0("LD2 (", ld2_pct, "%)"),
    color = "Taxon",
    shape = "Dataset"
  ) +
  theme_classic(base_size = 14)

print(p_test)

ggsave(
  file.path(out_dir, "LDA_training_testing.pdf"),
  p_test,
  width = 9,
  height = 7
)

# ============================================================
# 15. EXPORT INDIVIDUAL CLASSIFICATIONS
# ============================================================

test_results <- data.frame(
  Actual = testing$taxonlump,
  Predicted = test_pred$class,
  Correct = test_pred$class == testing$taxonlump,
  test_pred$posterior
)

write.csv(
  test_results,
  file.path(out_dir, "LDA_test_classifications.csv"),
  row.names = FALSE
)

# ============================================================
# 16. EXPORT CONFUSION MATRICES
# ============================================================

write.csv(
  as.data.frame(train_confusion$table),
  file.path(out_dir, "LDA_training_confusion.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(test_confusion$table),
  file.path(out_dir, "LDA_testing_confusion.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(cv_confusion$table),
  file.path(out_dir, "LDA_LOOCV_confusion.csv"),
  row.names = FALSE
)

cat("\nLDA analysis complete.\n")
cat("Results saved to:", out_dir, "\n")
