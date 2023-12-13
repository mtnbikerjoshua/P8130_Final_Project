library(tidyverse)
library(jcreg)
set.seed(1111)

#Import data
cancer <- read_csv("data/Project_2_data.csv", 
                   col_types = "iffffffffiffiiif") %>%
  janitor::clean_names() %>%
  mutate(
    status = case_match(status, "Alive" ~ 0, "Dead" ~ 1) %>%
      as.factor(),
    id = 1:nrow(.)
  )

#Separate into training and test data
cancer_train <- cancer %>%
  sample_frac(0.80)
cancer_test <- anti_join(cancer, cancer_train, by = 'id')

#Visualize correlation between predictors and response
comp_bar <- function(data, x) {
  ggplot(data = data, mapping = aes(x = .data[[x]], fill = status)) +
    geom_bar(position = "fill") +
    scale_fill_viridis_d(labels = c("Alive", "Dead")) +
    scale_x_discrete(guide = guide_axis(angle = -45)) +
    labs(y = "Proportion", fill = "Status") +
    theme_classic()
}

cancer_cat_varnames <- c("race", "marital_status", "t_stage", "n_stage", 
                         "x6th_stage", "differentiate", "grade", "a_stage",
                         "estrogen_status", "progesterone_status")
comp_bar_plots <- map(cancer_cat_varnames, \(var) comp_bar(cancer_train, var))
do.call(gridExtra::grid.arrange, comp_bar_plots)




# cancer %>%
#   select(where(is.numeric) | where(is.factor)) %>%
#   cor_graphic()
# 
# model <- glm(status ~ age + race+ t_stage + grade, 
#              family=binomial(link='logit'), data=cancer)
