library(dplyr)
library(ggplot2)


###=============================================================================
# All Populations
###=============================================================================

all_likes <- read.delim('ngsadmix_all_likelihoods.tsv')

all_delta_k <- all_likes %>%
  group_by(K) %>%
  summarise(
    mean_like = mean(logLik),
    sd_like = sd(logLik),
    .groups = "drop"
  ) %>%
  arrange(K) %>%
  mutate(
    delta1 = mean_like - lag(mean_like),
    delta2 = lead(delta1) - delta1,
    deltaK = abs(delta2) / sd_like
  )

print(all_delta_k)

ggplot(all_delta_k, aes(x = K, y = deltaK)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 3, na.rm = TRUE) +
  scale_x_continuous(breaks = 1:7) +
  labs(
    title = "Evanno-style Delta K all populations",
    x = "Number of ancestral populations (K)",
    y = expression(Delta*K)
  ) +
  theme_classic()

###=============================================================================
# Nopo Only
###=============================================================================

nopo_likes <- read.delim('ngsadmix_nopo_only_likelihoods.tsv')

nopo_delta_k <- nopo_likes %>%
  group_by(K) %>%
  summarise(
    mean_like = mean(logLik),
    sd_like = sd(logLik),
    .groups = "drop"
  ) %>%
  arrange(K) %>%
  mutate(
    delta1 = mean_like - lag(mean_like),
    delta2 = lead(delta1) - delta1,
    deltaK = abs(delta2) / sd_like
  )

print(nopo_delta_k)

ggplot(nopo_delta_k, aes(x = K, y = deltaK)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 3, na.rm = TRUE) +
  scale_x_continuous(breaks = 1:7) +
  labs(
    title = "Evanno-style Delta K only G. gnoma populations",
    x = "Number of ancestral populations (K)",
    y = expression(Delta*K)
  ) +
  theme_classic()
