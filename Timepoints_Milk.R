# Load libraries
library(tidyverse)

# Load your two CSV files
final72 <- read.csv("C:/Users/singh/Documents/R studio Genome Analysis/Milk analysis 2hr vs 72hrs/72merged.bracken_S.counts+percent.csv",
                    check.names = FALSE)

alltimes <- read.csv("C:/Users/singh/Documents/R studio Genome Analysis/Milk analysis 2hr vs 72hrs/all_times_merged.bracken_S.counts+percent.csv",
                     check.names = FALSE)

# Peek at the first few rows
head(final72)
head(alltimes)






# Keep only species + percent columns
final72_perc <- final72 %>%
  select(species, contains("percent"))

alltimes_perc <- alltimes %>%
  select(species, contains("percent"))

# Peek at the first few rows
head(final72_perc)
head(alltimes_perc)





alltimes_long <- alltimes_perc %>%
  pivot_longer(-species, names_to = "sample", values_to = "abundance")

head(alltimes_long)






library(stringr)

alltimes_long <- alltimes_long %>%
  mutate(
    time_sec = as.numeric(str_extract(sample, "\\d+(?=s_start)")),
    time_hr  = time_sec / 3600
  )

head(alltimes_long)






alltimes_long <- alltimes_long %>%
  mutate(
    replicate = str_extract(sample, "barcode\\d+")
  )

head(alltimes_long)




alltimes_avg <- alltimes_long %>%
  group_by(species, time_hr) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  ungroup()

head(alltimes_avg)




final72_avg <- final72_perc %>%
  pivot_longer(-species, names_to = "replicate", values_to = "abundance") %>%
  group_by(species) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  mutate(time_hr = 72) %>%
  ungroup()

head(final72_avg)


combined <- bind_rows(alltimes_avg, final72_avg)

head(combined)


unique(combined$time_hr)




combined_filtered <- combined %>%
  filter(mean_abundance >= 0.01)

head(combined_filtered)






library(tidyr)

# reshape into wide: each row = timepoint, each column = species
community_matrix <- combined_filtered %>%
  pivot_wider(
    names_from = species,
    values_from = mean_abundance,
    values_fill = 0
  ) %>%
  arrange(time_hr)

head(community_matrix)



library(vegan)
# Abundance matrix with time_hr as row names
abundances <- community_matrix %>%
  column_to_rownames("time_hr")

# Bray–Curtis dissimilarity matrix
bray_dist <- vegdist(abundances, method = "bray")
bray_mat <- as.matrix(bray_dist)

# Extract distances to 72h
dist_to_72 <- bray_mat[ , "72"]

# Save results
bray_results <- tibble(
  time_hr = as.numeric(names(dist_to_72)),
  bray_vs_72 = dist_to_72
)

print(bray_results, n = Inf)
write.csv(bray_results, "bray_curtis_vs72.csv", row.names = FALSE)





# Convert abundances to presence/absence (0/1)
abund_pa <- abundances
abund_pa[abund_pa > 0] <- 1

# Compute Jaccard dissimilarity
jaccard_dist <- vegdist(abund_pa, method = "jaccard", binary = TRUE)
jaccard_mat <- as.matrix(jaccard_dist)

# Extract distances to 72h
jaccard_to_72 <- jaccard_mat[ , "72"]

# Combine with Bray-Curtis results
beta_results <- tibble(
  time_hr = as.numeric(names(jaccard_to_72)),
  bray_vs_72 = bray_results$bray_vs_72,
  jaccard_vs_72 = jaccard_to_72
)

print(beta_results, n = Inf)
write.csv(beta_results, "bray_jaccard_vs72.csv", row.names = FALSE)





#SPEARMAN CORRELATION
# Extract 72h abundances as a numeric vector
abund_72 <- as.numeric(abundances["72", ])

# Spearman correlation per timepoint vs 72h
spearman_results <- apply(abundances, 1, function(x) {
  cor(as.numeric(x), abund_72, method = "spearman")
})

spearman_table <- tibble(
  time_hr = as.numeric(rownames(abundances)),
  spearman_vs_72 = spearman_results
)

print(spearman_table, n = Inf)
write.csv(spearman_table, "spearman_vs72.csv", row.names = FALSE)




#PERMANOVA
sum(is.na(abundances_repl))



abundances_repl <- abundances_repl %>% replace(is.na(.), 0)

sum(is.na(abundances_repl))   # should now return 0


library(tidyverse)
library(vegan)

### 1. Prepare replicate-level data (with 1% cutoff applied)

# From your long-format alltimes data
replicate_matrix <- alltimes_long %>%
  filter(abundance >= 0.01) %>%
  select(species, replicate, time_hr, abundance) %>%
  pivot_wider(
    names_from = species,
    values_from = abundance,
    values_fill = list(abundance = 0)   # fill missing with 0
  )

# Prepare 72h replicate-level data
replicate_72 <- final72_perc %>%
  pivot_longer(-species, names_to = "replicate", values_to = "abundance") %>%
  mutate(time_hr = 72) %>%
  filter(abundance >= 0.01) %>%
  pivot_wider(
    names_from = species,
    values_from = abundance,
    values_fill = list(abundance = 0)
  )

# Combine both
replicate_matrix <- bind_rows(replicate_matrix, replicate_72)

### 2. Split metadata and abundance matrix
metadata <- replicate_matrix %>% select(replicate, time_hr)
abundances_repl <- replicate_matrix %>% select(-replicate, -time_hr)

# Ensure no NAs remain
abundances_repl <- abundances_repl %>% replace(is.na(.), 0)

### 3. Run PERMANOVA
perm_results <- adonis2(
  abundances_repl ~ factor(metadata$time_hr),
  method = "bray",
  permutations = 999
)

print(perm_results)

### 4. Save results as CSV
perm_df <- as.data.frame(perm_results)
write.csv(perm_df, "permanova_results.csv")




install.packages("patchwork")
library(patchwork)















library(tidyverse)
library(vegan)

### 1. Prepare replicate-level data (with 1% cutoff applied)

# From your long-format alltimes data
replicate_matrix <- alltimes_long %>%
  filter(abundance >= 0.01) %>%
  select(species, replicate, time_hr, abundance) %>%
  pivot_wider(
    names_from = species,
    values_from = abundance,
    values_fill = list(abundance = 0)   # fill missing with 0
  )

# Prepare 72h replicate-level data
replicate_72 <- final72_perc %>%
  pivot_longer(-species, names_to = "replicate", values_to = "abundance") %>%
  mutate(time_hr = 72) %>%
  filter(abundance >= 0.01) %>%
  pivot_wider(
    names_from = species,
    values_from = abundance,
    values_fill = list(abundance = 0)
  )

# Combine both
replicate_matrix <- bind_rows(replicate_matrix, replicate_72)

### 2. Split metadata and abundance matrix
metadata <- replicate_matrix %>% select(replicate, time_hr)
abundances_repl <- replicate_matrix %>% select(-replicate, -time_hr)

# Ensure no NAs remain
abundances_repl <- abundances_repl %>% replace(is.na(.), 0)

### 3. Run PERMANOVA for Bray–Curtis
perm_bray <- adonis2(
  abundances_repl ~ factor(metadata$time_hr),
  method = "bray",
  permutations = 999
)

### 4. Run PERMANOVA for Jaccard
perm_jaccard <- adonis2(
  abundances_repl ~ factor(metadata$time_hr),
  method = "jaccard",
  permutations = 999
)

### 5. Combine results into one table
perm_df_bray <- as.data.frame(perm_bray) %>% 
  rownames_to_column("Term") %>%
  mutate(Method = "Bray-Curtis")

perm_df_jaccard <- as.data.frame(perm_jaccard) %>% 
  rownames_to_column("Term") %>%
  mutate(Method = "Jaccard")

perm_combined <- bind_rows(perm_df_bray, perm_df_jaccard)

### 6. Save results
write.csv(perm_combined, "permanova_bray_jaccard_results.csv", row.names = FALSE)

### 7. Print to console
print(perm_combined)









ggplot(bray_results, aes(x = factor(time_hr), y = bray_vs_72)) +
  geom_boxplot(fill = "firebrick", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(x = "Time (hours)", y = "Bray–Curtis dissimilarity")
ggplot(bray_results, aes(x = factor(time_hr), y = bray_vs_72)) +
  geom_boxplot(fill = "firebrick", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(x = "Time (hours)", y = "Bray–Curtis dissimilarity") +
  theme_minimal()





library(dplyr)
library(tidyr)

# Use your filtered combined dataset (the one with ≥1% applied)
species_abund_diff <- combined_filtered %>%
  filter(time_hr %in% c(2, 72)) %>%
  group_by(species, time_hr) %>%
  summarise(mean_abundance = mean(mean_abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = time_hr,
    values_from = mean_abundance,
    names_prefix = "hr_"
  ) %>%
  # Only keep species that differ between 2h and 72h
  filter(hr_2 != hr_72)

# Print all differences
print(species_abund_diff, n = Inf)

# Save to CSV
write.csv(species_abund_diff, "species_diff_2h_vs_72h_filtered.csv", row.names = FALSE)




#TRIAL
library(dplyr)
library(vegan)
library(tidyr)

# ---- Prepare community matrix (species as columns, timepoints as rows) ----
community_matrix <- combined_filtered %>%
  select(species, time_hr, mean_abundance) %>%
  pivot_wider(names_from = species, values_from = mean_abundance, values_fill = 0)

# Extract the abundance-only matrix
abundances <- community_matrix %>% select(-time_hr)

# ---- Calculate alpha diversity indices ----

# Species richness (number of non-zero species per timepoint)
richness <- rowSums(abundances > 0)

# Shannon index
shannon <- diversity(abundances, index = "shannon")

# Simpson index
simpson <- diversity(abundances, index = "simpson")

# ---- Combine results into one tibble ----
alpha_div <- tibble(
  time_hr = community_matrix$time_hr,
  richness = richness,
  shannon = shannon,
  simpson = simpson
)

# Print results
print(alpha_div, n = Inf)

# Save to CSV
write.csv(alpha_div, "alpha_diversity_filtered.csv", row.names = FALSE)








library(dplyr)

# Get species present at 2h
species_2h <- combined_filtered %>%
  filter(time_hr == 2, mean_abundance > 0) %>%
  pull(species)

# Get species present at 72h
species_72h <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  pull(species)

# Species present at 72h but not at 2h
new_species <- setdiff(species_72h, species_2h)

# Print them
print(new_species)

# Save to CSV
write.csv(data.frame(new_species), "species_present_at_72h_not_2h.csv", row.names = FALSE)





# Species present at 2h but not at 72h (lost species)
lost_species <- setdiff(species_2h, species_72h)

print(lost_species)
write.csv(data.frame(lost_species), "species_present_at_2h_not_72h.csv", row.names = FALSE)





library(dplyr)

# Species at 2h
species_2h_abund <- combined_filtered %>%
  filter(time_hr == 2, mean_abundance > 0) %>%
  select(species, abundance_2h = mean_abundance)

# Species at 72h
species_72h_abund <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  select(species, abundance_72h = mean_abundance)

# Species gained at 72h (present at 72h but not at 2h)
gained <- anti_join(species_72h_abund, species_2h_abund, by = "species") %>%
  mutate(status = "Gained at 72h") %>%
  rename(relative_abundance = abundance_72h)

# Species lost by 72h (present at 2h but not at 72h)
lost <- anti_join(species_2h_abund, species_72h_abund, by = "species") %>%
  mutate(status = "Lost by 72h") %>%
  rename(relative_abundance = abundance_2h)

# Combine into one summary table
summary_table <- bind_rows(gained, lost) %>%
  arrange(status, desc(relative_abundance))

# Print full table
print(summary_table, n = Inf)

# Save to CSV
write.csv(summary_table, "species_gain_loss_2h_vs_72h.csv", row.names = FALSE)




# Species present at 2h
species_2h <- combined_filtered %>%
  filter(time_hr == 2, mean_abundance > 0) %>%
  pull(species)

# Species present at 72h
species_72h <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  pull(species)

# Net new species (in 72h but not in 2h, minus those that disappear)
species_net_72 <- setdiff(species_72h, species_2h)

# Create table with their relative abundance at 72h
species_net_table <- combined_filtered %>%
  filter(time_hr == 72, species %in% species_net_72) %>%
  select(species, mean_abundance) %>%
  arrange(desc(mean_abundance))

print(species_net_table, n = Inf)

write.csv(species_net_table, "species_net_gain_72h_vs_2h.csv", row.names = FALSE)








# Species present at 2h
species_2h <- combined_filtered %>%
  filter(time_hr == 2, mean_abundance > 0) %>%
  pull(species)

# Species present at 72h
species_72h <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  pull(species)

# Species lost (present at 2h, absent at 72h)
lost_species <- setdiff(species_2h, species_72h)

# Species gained (present at 72h, absent at 2h)
gained_species <- setdiff(species_72h, species_2h)

# Net = gained minus lost = 19 species
net_species <- setdiff(gained_species, lost_species)

# Table of the 19 species with abundances at 72h
species_net_table <- combined_filtered %>%
  filter(time_hr == 72, species %in% net_species) %>%
  select(species, mean_abundance) %>%
  arrange(desc(mean_abundance))

print(species_net_table, n = Inf)

write.csv(species_net_table, "species_net_19_gain_72h_vs_2h.csv", row.names = FALSE)







library(dplyr)

# Species at 2h
species_2h <- combined_filtered %>%
  filter(time_hr == 2, mean_abundance > 0) %>%
  pull(species)

# Species at 72h
species_72h <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  pull(species)

# Gained (72h not in 2h)
gained_species <- setdiff(species_72h, species_2h)

# Lost (2h not in 72h)
lost_species <- setdiff(species_2h, species_72h)

# Counts
summary_counts <- tibble(
  Category = c("Species at 2h", "Species at 72h", 
               "Gained at 72h", "Lost by 72h", "Net difference"),
  Count = c(length(species_2h), length(species_72h), 
            length(gained_species), length(lost_species),
            length(species_72h) - length(species_2h))
)

print(summary_counts)

# Save summary
write.csv(summary_counts, "species_richness_summary_2h_vs_72h.csv", row.names = FALSE)




library(dplyr)

# ---- Species sets ----
species_2h <- combined_filtered %>%
  filter(time_hr == 2, mean_abundance > 0) %>%
  select(species, abundance_2h = mean_abundance)

species_72h <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  select(species, abundance_72h = mean_abundance)

# ---- Gained species (72h not in 2h) ----
gained_species <- anti_join(species_72h, species_2h, by = "species") %>%
  mutate(status = "Gained at 72h")

# ---- Lost species (2h not in 72h) ----
lost_species <- anti_join(species_2h, species_72h, by = "species") %>%
  mutate(status = "Lost by 72h")

# ---- Summary counts ----
summary_counts <- tibble(
  Category = c("Species at 2h", "Species at 72h", 
               "Gained at 72h", "Lost by 72h", "Net difference"),
  Count = c(nrow(species_2h), nrow(species_72h),
            nrow(gained_species), nrow(lost_species),
            nrow(species_72h) - nrow(species_2h))
)

# ---- Print results ----
print(summary_counts)
print(gained_species, n = Inf)
print(lost_species, n = Inf)

# ---- Save to CSV ----
write.csv(summary_counts, "summary_species_counts.csv", row.names = FALSE)
write.csv(gained_species, "species_gained_at_72h.csv", row.names = FALSE)
write.csv(lost_species, "species_lost_by_72h.csv", row.names = FALSE)








library(dplyr)
library(purrr)

# Get all unique timepoints except 72
timepoints <- combined_filtered %>% 
  pull(time_hr) %>% 
  unique() %>% 
  setdiff(72) %>% 
  sort()

# Species at 72h
species_72h <- combined_filtered %>%
  filter(time_hr == 72, mean_abundance > 0) %>%
  select(species, abundance_72h = mean_abundance)

# Function to compare one timepoint vs 72h
compare_to_72 <- function(tp) {
  
  # Species at current timepoint
  species_tp <- combined_filtered %>%
    filter(time_hr == tp, mean_abundance > 0) %>%
    select(species, abundance_tp = mean_abundance)
  
  # Gained = in 72h but not in tp
  gained <- anti_join(species_72h, species_tp, by = "species") %>%
    mutate(status = "Gained at 72h",
           time_hr = tp)
  
  # Lost = in tp but not in 72h
  lost <- anti_join(species_tp, species_72h, by = "species") %>%
    mutate(status = "Lost by 72h",
           time_hr = tp)
  
  bind_rows(gained, lost)
}

# Apply to all timepoints
summary_all <- map_dfr(timepoints, compare_to_72)

# Reorder columns for clarity
summary_all <- summary_all %>%
  select(time_hr, species, status, everything())

# Print a sample
print(summary_all, n = 50)

# Save full table
write.csv(summary_all, "species_gain_loss_all_vs_72h.csv", row.names = FALSE)












library(dplyr)

# ---- Top 20 taxa per timepoint ----
top20_by_time <- combined_filtered %>%
  group_by(time_hr) %>%
  arrange(time_hr, desc(mean_abundance)) %>%
  slice_head(n = 20) %>%   # take top 20 for each timepoint
  ungroup()

# View the top 20 lists
print(top20_by_time, n = Inf)

# ---- Compare overlaps across timepoints ----
# Count how many unique species are in the top 20 at each timepoint
top20_summary <- top20_by_time %>%
  group_by(time_hr) %>%
  summarise(unique_species = n_distinct(species),
            top_species = paste(species, collapse = "; "))

print(top20_summary, n = Inf)

# ---- Optionally: find "core" top taxa across all timepoints ----
core_top_species <- top20_by_time %>%
  group_by(species) %>%
  summarise(n_timepoints = n_distinct(time_hr)) %>%
  filter(n_timepoints == length(unique(combined_filtered$time_hr)))

print(core_top_species, n = Inf)

# ---- Save to CSV ----
write.csv(top20_by_time, "top20_taxa_per_timepoint.csv", row.names = FALSE)
write.csv(top20_summary, "top20_summary_overlap.csv", row.names = FALSE)
write.csv(core_top_species, "core_top20_species_all_timepoints.csv", row.names = FALSE)
