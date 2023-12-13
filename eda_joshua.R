library(tidyverse)
library(car)
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

cancer_cat_varnames <- cancer_train %>%
  select(where(is.factor)) %>%
  colnames()

comp_bar_plots <- map(cancer_cat_varnames, \(var) comp_bar(cancer_train, var))
do.call(gridExtra::grid.arrange, comp_bar_plots)


#Check for multicolinearity
cancer_train %>% 
  select(where(is.numeric)) %>%
  cor()

cancer %>%
  select(where(is.numeric)) %>%
  cor_graphic()


cancer_cramerV <- matrix(nrow = length(cancer_cat_varnames), 
                         ncol = length(cancer_cat_varnames))
colnames(cancer_cramerV) <- cancer_cat_varnames
rownames(cancer_cramerV) <- cancer_cat_varnames
for(rowvar in cancer_cat_varnames) {
  for(colvar in cancer_cat_varnames) {
    cancer_cramerV[rowvar, colvar] <- 
      table(cancer_train[[rowvar]], cancer_train[[colvar]]) %>%
      rcompanion::cramerV()
  }
}
cancer_cramerV

cancer_con_varnames <- cancer_train %>%
  select(-id) %>%
  select(where(is.numeric)) %>%
  colnames()
cancer_eta2 <- matrix(nrow = length(cancer_cat_varnames), 
                         ncol = length(cancer_con_varnames))
colnames(cancer_eta2) <- cancer_con_varnames
rownames(cancer_eta2) <- cancer_cat_varnames
for(rowvar in cancer_cat_varnames) {
  for(colvar in cancer_con_varnames) {
    cancer_eta2[rowvar, colvar] <-
      lm(cancer_train[[colvar]] ~ cancer_train[[rowvar]]) %>%
      effectsize::eta_squared(partial = FALSE) %>%
      .$Eta2
  }
}
cancer_eta2

cancer_glm <- glm(status ~ age + race + marital_status + x6th_stage +
                    differentiate + a_stage + tumor_size +
                    estrogen_status + progesterone_status +
                    regional_node_examined,
                  family=binomial(link='logit'), data=cancer_train)
vif(cancer_glm)
