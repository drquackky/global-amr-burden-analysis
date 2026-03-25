# ===============================================
# LOAD REQUIRED LIBRARIES
# ===============================================
library(dplyr)        # data manipulation
library(ggplot2)      # plotting
library(lme4)         # mixed model engine
library(lmerTest)     # the add-on to mixed models providing p-values, interaction tools
library(emmeans)      # estimated marginal means
library(performance)  # model diagnostics (R²)



# ===============================================
# LOAD AND PREPARE DATA
# ===============================================
df <- read.csv("dataset/asr_lmm_4run.csv") %>%
   
   # Convert variables to ordered factors
   mutate(priority_group = factor(
      priority_group,
      levels = c("Critical", "High", "Medium")
   )) %>%
   
   mutate(sdi_cat = factor(
      sdi_cat,
      levels = c("low", "middle_low", "middle", "middle_high", "high")
   ))


# ===============================================
# INITIAL MODEL CHECK (RAW DATA)
# ===============================================

# → Residuals vs fitted: should be randomly scattered around 0
# → QQ plot: points should follow straight line
# which violates the assumption of the mixed linear model

model <- lmer(
   ASR ~ sdi_cat * priority_group + Year + (1 | Location),
   data = df[df$metric == "Associated_Deaths", ]
)

# Residual diagnostics
plot(model) # → Residuals not randomly scattered around 0, patterns suggest non-linearity or unequal variance

qqnorm(resid(model)) # → Deviations indicate non-normal residuals
qqline(resid(model))

# Conclusion: Residuals not well behaved → transformation needed




# ===============================================
# LOG TRANSFORMATION
# ===============================================
# → Log transformation reduces skewness and stabilizes variance
# → Helps meet normality assumption for linear models

df_log_rate <- df %>%
   mutate(log_rate = log(ASR))


# ===============================================
# WRITE A FUNCTION FOR FULL ANALYSIS PIPELINE
# ===============================================

run_analysis <- function(data, metric_name){

   # 1. filter data - select one of the 4 burden indicators (ex. Attributable_Deaths)
   df_filtered <- data %>%
      filter(metric == metric_name)
   
   # 2. fit mixed linear model on log-transformed rate
   # Fixed effects: SDI, priority group, year
   # Random effect: variation across locations
   model <- lmer(
      log_rate ~ sdi_cat * priority_group + Year + (1 | Location),
      data = df_filtered)

   # 3.model output
   # test if the interaction effect (sdi_cat x priority_group) is significant
   model_anova <- anova(model, type = 3)

   # R² represents the proportion of variance in the outcome explained by the model
   model_r2 <- performance::r2(model)
   
   # → Marginal R²: only variance explained by fixed effects (SDI, priority group, year)
   # → Conditional R²: variance explained by full model (fixed + random(location differences))
   
   
   # 4. estimated marginal means - model-adjusted means for each SDI category within each priority group
   emm_obj <- emmeans(model, ~ sdi_cat | priority_group)
   emm_df  <- as.data.frame(emm_obj)
   # Back-transform from log scale by exponentiation for interpretation
   emm_df <- emm_df %>%
      mutate(
         emmean_back = exp(emmean),
         lower_back  = exp(lower.CL),
         upper_back  = exp(upper.CL))
   

   # 5. pairwise comparison - compare SDI groups within each priority group

   pairs_obj <- pairs(
      emm_obj,
      adjust = "tukey",
      infer = c(TRUE, TRUE))
   pairs_df <- as.data.frame(pairs_obj)
   
   # Back-transform contrasts
   pairs_df <- pairs_df %>%
      mutate(
         estimate_back = exp(estimate),
         lower_back    = exp(lower.CL),
         upper_back    = exp(upper.CL))
   

   # 6. return results
   return(list(
      model    = model,
      anova    = model_anova,
      r2       = model_r2,
      emm      = emm_df,
      contrast = pairs_df
   ))
}

# ===============================================
# RUN ANALYSIS
# ===============================================

# Increase limit for degrees-of-freedom calculations in emmeans
emm_options(pbkrtest.limit = 20000, lmerTest.limit = 20000)

# check the metric name to put for the function
unique(df_log_rate$metric)

# run the function
assdeath_result <- run_analysis(df_log_rate, "Associated_Deaths")

assdeath_result$anova
# → Significant interaction: interpret sdi_cat × priority_group together, not separately

assdeath_result$r2
# → Fixed effects (sdi_cat, priority_group, Year) explain ~8% of variation in AMR burden
# → while including Location as a random effect increases explained variation to ~77%
# → This means that SDI and priority group alone do not fully explain AMR burden
# → and there are substantial differences in burden across countries

assdeath_result$emm
# → show model-adjusted means for each SDI category within each priority group

assdeath_result$contrast
# compare burden rate between SDI groups within each priority group
# estimate_back is interpreted as ratio
# ex. low - high (Critical) - 2.41 meaning that low sdi has 141 percent higher burden than high sdi


# ===============================================
# DATA VISUALISATION
# ===============================================
# extract estimated marginal means for plotting
assdeath_emm_df <- assdeath_result$emm

# plot interaction using ggplot
library(ggplot2)
library(RColorBrewer)

p1 <- ggplot(assdeath_emm_df, aes(
   x = sdi_cat, 
   y = emmean_back,          
   color = priority_group, 
   group = priority_group)) +
   
   geom_point(size = 4, shape = 19) +
   geom_line(size = 1.2) +
   geom_errorbar(aes(ymin = lower_back, ymax = upper_back), width = 0.2, size = 0.8, alpha = 0.7) +
   
   scale_color_brewer(palette = "Set1") +
   
   theme_minimal(base_size = 14) +
   theme(
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
   ) +
   
   ggtitle("Interaction Between SDI Category and Priority Group influencing the difference in Associated Death rate") +
   ylab("Age-standardized rate per 100.000 (Back-transformed)") +
   xlab("SDI Category") +
   labs(color = "Priority Group")

ggsave("figures/AMR_interaction_plot.png",     
       width = 10, height = 6, dpi = 300)   






# box plot of associated death rate across different SDI levels
# superimpose the estimated mean from the model (emmeans) to the box plot of rate distribution


ass_death_data <- df_log_rate %>%
   filter(metric == "Associated_Deaths")

p2 <- ggplot(ass_death_data, aes(x = sdi_cat, y = log_rate)) +
   geom_boxplot(outlier.colour = "red",
                outlier.size = 2,
                outlier.shape = 8,
                #notch = TRUE,
                varwidth = TRUE,
                size = 0.8) +
   geom_violin(alpha = 0.01) +
   ggforce::geom_sina(
      data = ass_death_data %>% sample_n(700), 
      alpha = 0.15, 
      size = 3, 
      aes(color = sdi_cat)
      ) +
   stat_summary(fun.data = mean_cl_boot, geom = "point",
                color = "blue", size = 3) +
   stat_summary(fun.data = mean_cl_boot, geom = "errorbar",
                color = "blue", size = 1, width = 0.4) + 
   theme_minimal() +
   theme(legend.position = "none",
         plot.title = element_text(size = 16, face = "bold"),
         axis.text.x = element_text(angle = 45, hjust = 1),
         axis.title.x = element_text(face = "bold"),
         axis.title.y = element_text(face = "bold")) +
   labs(
      title = "Distribution of Associated Death Rate across SDI categories",
      subtitle = "Log Transformation on data grouped by WHO pathogen priority group",
      caption = "",
      x = "Social Demographic Index Category",
      y = "Log transformed of Death Rate per 100,000"
   ) +
   facet_grid(~priority_group) +
   geom_hline (yintercept = mean(ass_death_data$log_rate),
               linetype = "dashed", color = "green",
               size = 1)


# superimpose a model predicted mean to a raw rate
p3 <- p2 +
   geom_line(
      data = assdeath_emm_df,
      aes(x = sdi_cat, y = emmean, group = priority_group),
      inherit.aes = FALSE,
      color = "red",
      size = 1.3,
      position = position_dodge(width = 0.6)
   ) +
   geom_point(
      data = assdeath_emm_df,
      aes(x = sdi_cat, y = emmean, color = priority_group),
      inherit.aes = FALSE,
      size = 3,
      color = "purple",
      position = position_dodge(width = 0.6)
   ) +
   geom_errorbar(
      data = assdeath_emm_df,
      aes(x = sdi_cat, ymin = lower.CL, ymax = upper.CL, color = priority_group),
      inherit.aes = FALSE,
      width = 0.2,
      color = "purple",
      size = 1,
      position = position_dodge(width = 0.6)
   ) 

ggsave("figures/AMR_distribution.png",     
       width = 17, height = 11, dpi = 300)   
