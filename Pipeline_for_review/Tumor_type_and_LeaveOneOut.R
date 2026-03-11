library(meta)
library(tidyverse)

anno = read_csv("~/Rproject/ICI/ICIs/R-data/R-csv/Anno.CSV")
Index = read_csv("~/Rproject/ICI/Index.CSV")


table(anno$Type)

# output_new is the list by study with cellstate meta-analysis results 

# Match the cellstate names
output_new = lapply(output_new, function(x){
  x$Names = Index$Names[match(rownames(x),Index$CellState)]

  x = x[!(is.na(x$Names)),]
  rownames(x) = NULL
  x = x %>% column_to_rownames("Names")
  return(x)
})



# cohort metadata  ----
cohort_info = data.frame(
  cohort = names(output_new),
  tumor_type = anno$Type[match(names(output_new),anno$File)]
)

cohort_info = cohort_info %>%
  mutate(
    tumor_group = ifelse(tolower(tumor_type) == "melanoma", "melanoma", "carcinoma"),
    tumor_group = factor(tumor_group, levels = c("carcinoma", "melanoma"))
  )


coef_long = bind_rows(lapply(names(output_new), function(cc) {
  df = output_new[[cc]]
  df$feature = rownames(df)
  df$cohort  = cc
  df
})) %>%
  select(cohort, feature, estimate, std.error, statistic, p.value) %>%
  left_join(cohort_info, by = "cohort") %>%
  filter(!is.na(tumor_group))  


# Extract summary from metagen 
extract_meta = function(m, which = c("fixed","random")) {
  which = match.arg(which)
  if (is.null(m)) return(
    data.frame(
      TE = NA_real_, seTE = NA_real_,
      lower = NA_real_, upper = NA_real_, p = NA_real_,
      I2 = NA_real_, pQ = NA_real_, H = NA_real_,
      OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_
    )
  )
  if (which == "fixed") {
    TE    = m$TE.fixed
    seTE  = m$seTE.fixed
    lower = m$lower.fixed
    upper = m$upper.fixed
    p     = m$pval.fixed
  } else {
    TE    = m$TE.random
    seTE  = m$seTE.random
    lower = m$lower.random
    upper = m$upper.random
    p     = m$pval.random
  }
  data.frame(
    TE = TE, seTE = seTE, lower = lower, upper = upper, p = p,
    I2 = m$I2, pQ = m$pval.Q, H = m$H,
    OR = exp(TE), OR_lower = exp(lower), OR_upper = exp(upper)
  )
}

get_metareg_p = function(m) {
  out = tryCatch(metareg(m, ~ tumor_group), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)

  if (!is.null(out$pval) && length(out$pval) >= 2) return(as.numeric(out$pval[2]))
  
  ss = tryCatch(summary(out), error = function(e) NULL)
  if (!is.null(ss)) {
    if (!is.null(ss$pval) && length(ss$pval) >= 2) return(as.numeric(ss$pval[2]))
    if (!is.null(ss$coef) && nrow(ss$coef) >= 2 && ncol(ss$coef) >= 4) return(as.numeric(ss$coef[2,4]))
  }
  NA_real_
}


features = sort(unique(coef_long$feature))

meta_table = bind_rows(lapply(features, function(f) {
  df = coef_long %>% filter(feature == f)
  
  # Need at least 2 cohorts overall for meta
  if (nrow(df) < 2) return(NULL)
  
  # overall meta
  m_all = metagen(TE = df$estimate, seTE = df$std.error,
                   studlab = df$cohort, sm = "OR", data = df,
                   comb.fixed = TRUE, comb.random = TRUE)
  
  # subgroup metas
  df_mel = df %>% filter(tumor_group == "melanoma")
  df_car = df %>% filter(tumor_group == "carcinoma")
  
  m_mel = if (nrow(df_mel) >= 2) metagen(TE=df_mel$estimate, seTE=df_mel$std.error,
                                          studlab=df_mel$cohort, sm="OR",
                                          comb.fixed=TRUE, comb.random=TRUE) else NULL
  
  m_car = if (nrow(df_car) >= 2) metagen(TE=df_car$estimate, seTE=df_car$std.error,
                                          studlab=df_car$cohort, sm="OR",
                                          comb.fixed=TRUE, comb.random=TRUE) else NULL
  
  # interaction test
  p_int = if (length(unique(df$tumor_group)) >= 2) get_metareg_p(m_all) else NA_real_
  
  out = data.frame(feature = f,
                    k_all = nrow(df),
                    k_melanoma = nrow(df_mel),
                    k_carcinoma = nrow(df_car),
                    p_interaction = p_int)
  
  # random-effects
  all_re = extract_meta(m_all, "random")
  mel_re = extract_meta(m_mel, "random")
  car_re = extract_meta(m_car, "random")
  
  colnames(all_re) = paste0("all_", colnames(all_re))
  colnames(mel_re) = paste0("mel_", colnames(mel_re))
  colnames(car_re) = paste0("car_", colnames(car_re))
  
  cbind(out, all_re, mel_re, car_re)
}))

# FDR correction on interaction tests
meta_table = meta_table %>%
  mutate(p_interaction_fdr = p.adjust(p_interaction, method = "fdr"))

sig_ove = meta_table %>% filter(all_p < 0.05) %>% pull(feature)

# Save 

meta_table =  cbind(sig_state = ifelse(meta_table$feature %in% sig_ove,"Yes","No"),meta_table)

save(meta_table, file = "~/Rproject/ICIs/New_data/meta_melanoma_vs_carcinoma.Rdata")
write.csv(meta_table, file = "~/Rproject/ICIs/New_tables/meta_melanoma_vs_carcinoma.csv", row.names = FALSE)



# Direction discordance ---------
dir_df = meta_table %>%
  filter(!is.na(mel_TE) & !is.na(car_TE)) %>%
  mutate(direction_same = sign(mel_TE) == sign(car_TE))

table(dir_df$direction_same)


library(ggplot2)
library(ggpubr)
plot_df = meta_table %>%
  filter(!is.na(mel_TE) & !is.na(car_TE))


pdf("~/Rproject/ICIs/New_plots/Tumor-group heterogeneity melanoma vs carcinoma.pdf",width = 4.5,height = 4.5)
ggscatter(plot_df,x = "car_TE", y = "mel_TE",add = "reg.line",cor.coef = T)+
  labs(x = "Carcinoma pooled log(OR) (random-effects)",
       y = "Melanoma pooled log(OR) (random-effects)",
       title = "Melanoma vs carcinoma heterogeneity")
dev.off()

# leave one out ==========

# tumor_types = sort(unique(na.omit(coef_long$tumor_type)))
# 
# loto_results = bind_rows(lapply(features, function(f) {
#   df0 = coef_long %>% filter(feature == f)
#   if (nrow(df0) < 3) return(NULL)
# 
#   base = metagen(TE=df0$estimate, seTE=df0$std.error, studlab=df0$cohort,
#                   sm="OR", comb.random=TRUE)
#   base_TE = base$TE.random
# 
#   bind_rows(lapply(tumor_types, function(tt) {
#     df = df0 %>% filter(tumor_type != tt)
#     if (nrow(df) < 2) return(NULL)
#     m = metagen(TE=df$estimate, seTE=df$std.error, studlab=df$cohort,
#                  sm="OR", comb.random=TRUE)
#     data.frame(feature=f, left_out=tt,
#                TE_random=m$TE.random, OR_random=exp(m$TE.random),
#                sign_flip = sign(m$TE.random) != sign(base_TE))
#   }))
# }))


# safe leave one out ===========

tumor_types = sort(unique(na.omit(coef_long$tumor_type)))
features    = sort(unique(coef_long$feature))

safe_metagen = function(df, random = TRUE) {
  df = df %>%
    filter(is.finite(estimate),
           is.finite(std.error),
           std.error > 0)
  
  if (nrow(df) < 2) return(NULL)
  
  if (random) {
   
    out = tryCatch(
      metagen(TE = df$estimate, seTE = df$std.error, studlab = df$cohort,
              sm = "OR",
              comb.fixed = TRUE, comb.random = TRUE,
              method.tau = "DL", hakn = FALSE),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
  
}

# extract TE preferentially from random
get_TE = function(m) {
  if (is.null(m)) return(NA_real_)
  if (!is.null(m$TE.random) && is.finite(m$TE.random)) return(m$TE.random)
  if (!is.null(m$TE.fixed)  && is.finite(m$TE.fixed))  return(m$TE.fixed)
  NA_real_
}

loto_results = dplyr::bind_rows(lapply(features, function(f) {
  df0 = coef_long %>% filter(feature == f)
  
  df0c = df0 %>% filter(is.finite(estimate), is.finite(std.error), std.error > 0)
  if (nrow(df0c) < 3) return(NULL)
  
  base_m  = safe_metagen(df0c, random = TRUE)
  base_TE = get_TE(base_m)
  if (!is.finite(base_TE)) return(NULL)
  
  bind_rows(lapply(tumor_types, function(tt) {
    df = df0c %>% filter(tumor_type != tt)
    if (nrow(df) < 2) return(NULL)
    
    m  = safe_metagen(df, random = TRUE)
    TE = get_TE(m)
    if (!is.finite(TE)) return(NULL)
    
    data.frame(
      feature   = f,
      left_out  = tt,
      k_left    = nrow(df),
      TE        = TE,
      OR        = exp(TE),
      sign_flip = (sign(TE) != sign(base_TE))
    )
  }))
}))

loto_results$sig_state = ifelse(loto_results$feature %in% sig_ove,"Yes","No")

write.csv(loto_results, "~/Rproject/ICIs/New_tables/leave_one_tumor_type_out.csv", row.names = FALSE)


# Pick the baseline TE from the full meta-analysis ----

meta_base = meta_table %>%
  transmute(
    feature,
    base_TE = all_TE,     # pooled log(OR) from overall meta-analysis
    base_OR = all_OR,     # pooled OR
    base_p  = all_p,      # overall p-value
    base_I2 = all_I2      # overall I2
  )

 # Join LOTO runs to baseline and compute deviations
loto_aug = loto_results %>%
  left_join(meta_base, by = "feature") %>%
  mutate(
    delta_TE = TE - base_TE,
    abs_delta_TE = abs(delta_TE),
    # sign(0) = 0; treat exact zeros as "no flip" unless you prefer otherwise
    sign_flip = (sign(TE) != sign(base_TE))
  )

# Compressed  summary 
loto_summary = loto_aug %>%
  group_by(feature) %>%
  summarise(
    n_runs = n(),
    min_k_left = min(k_left, na.rm = TRUE),
    max_k_left = max(k_left, na.rm = TRUE),
    
    OR_min = min(OR, na.rm = TRUE),
    OR_max = max(OR, na.rm = TRUE),
    OR_ratio = OR_max / OR_min,
    
    max_abs_delta_TE = max(abs_delta_TE, na.rm = TRUE),
    any_sign_flip = any(sign_flip, na.rm = TRUE),
    
    # which left-out tumor type caused the largest deviation
    worst_left_out = left_out[which.max(abs_delta_TE)],
    worst_k_left   = k_left[which.max(abs_delta_TE)],
    worst_TE       = TE[which.max(abs_delta_TE)],
    worst_OR       = OR[which.max(abs_delta_TE)],
    
    .groups = "drop"
  ) %>%
  left_join(meta_base, by = "feature") %>%
  arrange(base_p)

loto_summary = loto_summary %>%
  mutate(
    robust_no_flip = !any_sign_flip,
  )

loto_summary = loto_summary %>% mutate(sig_state = ifelse(loto_summary$feature %in% sig_ove,"Yes","No"),.after = 1)

# Write
write.csv(loto_summary, "~/Rproject/ICIs/New_tables/leave_one_tumor_type_out_summary.csv", row.names = FALSE)
