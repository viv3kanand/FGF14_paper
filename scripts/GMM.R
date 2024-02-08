process_file <- function(file_path) {
  require(tidyverse)
  require(mclust) # For fitting GMM
  
  df <- read.table(file_path, header = TRUE)
  
  df_filtered <- df %>%
    filter(count > 3, score_prefix >= 4, score_suffix >= 4) %>% 
    mutate(id = str_remove(basename(file_path), ".hg19.strique.tsv"))
  
  set.seed(666)
  
  fit <- Mclust(df_filtered$count)
  prob <- predict(fit)
  
  df_filtered$allele <- apply(prob$z, 1, which.max)
  
  mode_df <- df_filtered %>%
    group_by(id, allele) %>%
    summarise(
      allele_freq = n(),
      mode_repeat = as.numeric(names(sort(table(count), decreasing = TRUE)[1])),
      repeat_freq = max(table(count)),
      prop = repeat_freq/allele_freq,
      .groups = 'drop'
    ) %>% 
    arrange(desc(allele_freq))
  
  return(list(mode = mode_df, df = df_filtered))
}