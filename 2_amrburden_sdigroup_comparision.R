library(dplyr)
library(ggplot2)
library(lmerTest)
library(lme4)
library(emmeans)
library(performance)

df <- read.csv("dataset/df_final.csv") %>%
   mutate(priority_group = factor(
      priority_group,
      levels = c("Critical", "High", "Medium")
   )) %>%
   mutate(sdi_cat = factor(
      sdi_cat,
      levels = c("low", "middle_low", "middle", "middle_high", "high"))) 


df_attdrate <- df %>%
   filter(metric == "Attributable_Deaths")

df_assdrate <- df %>%
   filter(metric == "Associated_Deaths")

df_attdalys <- df %>%
   filter(metric == "Attributable_DALYs (Disability-Adjusted Life Years)")

df_assdalys <- df %>%
   filter(metric == "Associated_DALYs (Disability-Adjusted Life Years)")



# test
df_attdrate$log_rate <- log(df_attdrate$ASR)

model <- lmer(log_rate ~ sdi_cat * priority_group + (1|Location), data = df_attdrate)
summary(model)

anova(model, type = 3)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
performance::r2(model)

check_model(model)


# emmeans
emm_options(pbkrtest.limit = 20000, lmerTest.limit = 20000)

attdeath_emm_obj <- emmeans(model, ~ sdi_cat | priority_group)
attdeath_emm_df <- as.data.frame(attdeath_emm_obj)

attdeath_emm_df$emmean_back <- exp(attdeath_emm_df$emmean)
attdeath_emm_df$upper_back <- exp(attdeath_emm_df$upper.CL)
attdeath_emm_df$lower_back <- exp(attdeath_emm_df$lower.CL)


# getting emm pairwise comparision
attdeath_pairs_obj <- pairs(attdeath_emm_obj, adjust = "tukey", infer = c(TRUE, TRUE))
attdeath_pairs_df <- as.data.frame(attdeath_pairs_obj)

attdeath_pairs_df$estimate_back <- exp(attdeath_pairs_df$estimate)
attdeath_pairs_df$upper_back <- exp(attdeath_pairs_df$upper.CL)
attdeath_pairs_df$lower_back <- exp(attdeath_pairs_df$lower.CL)

# visualisation
ggplot(attdeath_emm_df, aes(x = sdi_cat, y = emmean_back, 
                            color = priority_group, 
                            group = priority_group)) +
   geom_point(size = 3) +
   geom_line(size = 1) +
   geom_errorbar(aes(ymin = lower_back, ymax = upper_back), width = 0.2) +
   theme_minimal(base_size = 14) + 
   theme(legend.position = "top") +
   ggtitle("Interaction Between SDI Category and Priority Group") +
   ylab("Attributable Rate (Back-transformed)") +
   xlab("SDI Category") +
   labs(color = "Priority Group")
