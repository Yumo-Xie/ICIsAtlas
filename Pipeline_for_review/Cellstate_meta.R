library(dplyr)
library(broom)
library(meta)

# Data Preparation and Normalization ---
# data0 contains the clinical response data per cohort
# data1 contains the cell abundance matrices
# All these data can be found in the supplemental table


analysis_output <- list()

for (i in 1:length(data0)) {
  current_cohort_name <- names(data0)[i]
  message(paste("Processing cohort:", current_cohort_name))
  
  # Extract clinical response (res) and ID
  clinical_df <- data0[[i]]
  if(ncol(clinical_df) > 27) {
    valid_samples <- clinical_df[clinical_df$P.value < 0.05, ]
  } else {
    valid_samples <- clinical_df
  }
  
  abundance_df <- as.data.frame(t(data1[[i]])) 
  colnames(abundance_df) <- paste(names(data1)[i], colnames(abundance_df), sep = "_")
  abundance_df$ID <- rownames(abundance_df)
  
  # Merge clinical and abundance data
  combined_data <- merge(abundance_df, valid_samples, by = "ID")
  
  # Isolate response variable and feature matrix
  response_vec <- as.factor(combined_data$res)
  feature_matrix <- combined_data %>% 
    select(-res, -ID) %>% 
    mutate(across(everything(), as.numeric))
  
  # Perform Z-score normalization
  scaled_features <- scale(feature_matrix)
  scaled_features[is.na(scaled_features)] <- 0
  
  # Logistic Regression (Intra-cohort) ---
  cohort_results <- data.frame()
  
  for (j in 1:ncol(scaled_features)) {
    # res ~ cell_state_abundance
    model <- glm(response_vec ~ scaled_features[, j], family = binomial(link = "logit"))
    tidy_model <- tidy(model)
    
    # Extract effect size (estimate) and uncertainty (std.error)
    # Removing intercept (row 1)
    stats_row <- tidy_model[-1, c("estimate", "std.error")]
    cohort_results <- rbind(cohort_results, stats_row)
    rownames(cohort_results)[j] <- colnames(scaled_features)[j]
  }
  
  analysis_output[[current_cohort_name]] <- cohort_results
}

# Random-Effects Meta-Analysis (Inter-cohort) ---
# Reshape data to group by cell state across all cohorts
cell_states <- rownames(analysis_output[[1]])
meta_summary_table <- data.frame()

for (state in cell_states) {
  state_data <- data.frame()
  for (cohort in names(analysis_output)) {
    row <- analysis_output[[cohort]][state, ]
    state_data <- rbind(state_data, cbind(Study = cohort, row))
  }
  
  # Run Random-Effects Meta-Analysis
  meta_model <- metagen(
    TE = estimate, 
    seTE = std.error, 
    data = state_data, 
    studlab = Study, 
    sm = "OR",
    random = TRUE
  )
  
  # Extract consolidated statistics
  pooled_stats <- data.frame(
    Cell_State = state,
    OR_random = exp(meta_model$TE.random),
    lower_CI = exp(meta_model$lower.random),
    upper_CI = exp(meta_model$upper.random),
    p_val_random = meta_model$pval.random,
    I2 = meta_model$I2
  )
  
  meta_summary_table <- rbind(meta_summary_table, pooled_stats)
}
