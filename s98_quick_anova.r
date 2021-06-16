df = read_csv('./data/final_csv/predictions_resolution.csv') %>% 
  filter(resolution == 0.01) %>% 
  select(method, location_error, orientation_error, x_source, y_source, orientation_source)
  mutate('source_id' = group_indices(., x_source, y_source, orientation_source))

# Summary stats
df %>%
  group_by(method) %>%
  get_summary_stats(location_error, orientation_error, type = "mean_sd")

# Create boxplots
bxp_loc <- ggboxplot(df, x = "method", y = "location_error", add = "point")
bxp_or <- ggboxplot(df, x = "method", y = "orientation_error", add = "point")


# There are extreme outliers
df %>%
  group_by(method) %>%
  identify_outliers(location_error) %>%
  filter(is.extreme)
df %>%
  group_by(method) %>%
  identify_outliers(orientation_error) %>%
  filter(is.extreme)

# Check normality (fails because to large sample size)
df %>%
  group_by(method) %>%
  shapiro_test(location_error)
df %>%
  group_by(method) %>%
  shapiro_test(orientation_error)

# Data is not normally distributed...
ggqqplot(df, "location_error", facet.by = "method")

# Still try it
res.aov <- anova_test(data = df, dv = location_error, wid = source_id, within = method)
get_anova_table(res.aov)

# Pairwise tests
pwc <- df %>%
  pairwise_t_test(
    location_error ~ method, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
