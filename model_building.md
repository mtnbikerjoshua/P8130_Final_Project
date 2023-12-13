Model Building
================
2023-12-12

``` r
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(dplyr)
library(MASS)
```

    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(ggplot2)
library(corrplot)
```

    ## corrplot 0.92 loaded

``` r
library(leaps)
library(glmnet)
```

    ## Loading required package: Matrix
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack
    ## 
    ## Loaded glmnet 4.1-8

``` r
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     %--%, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union
    ## 
    ## The following objects are masked from 'package:purrr':
    ## 
    ##     compose, simplify
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing
    ## 
    ## The following object is masked from 'package:tibble':
    ## 
    ##     as_data_frame
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(arules)
```

    ## 
    ## Attaching package: 'arules'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     abbreviate, write

# Data Cleaning

``` r
# import data
breastcancer = read_csv("./data/Project_2_data.csv")

#Data Cleaning
breastcancer_1 = breastcancer|>
  janitor::clean_names()|>
   mutate(
    race = as_factor(race),
    marital_status = factor(marital_status, levels = c("Single", "Married", "Divorced", "Separated", "Widowed")),
    t_stage = factor(t_stage, levels = c("T1", "T2", "T3", "T4")),
    n_stage = factor(n_stage, levels = c("N1", "N2", "N3")),
    x6th_stage = factor(x6th_stage, levels = c("IIA", "IIB", "IIIA", "IIIB", "IIIC")),
    differentiate = factor(differentiate, levels = c("Moderately differentiated", "Poorly differentiated", "Undifferentiated", "Well differentiated")),
    grade = factor(grade, levels = c("1", "2", "3", "anaplastic; Grade IV")),
    a_stage = factor(a_stage, levels = c("Distant", "Regional")),
    estrogen_status = as_factor(estrogen_status),
    progesterone_status = as_factor(progesterone_status),
    status = ifelse(status == "Dead", 1, 0),
    status = as_factor(status))

breastcancer_clean = breastcancer_1|>
  mutate(node_positive_prop = reginol_node_positive/regional_node_examined,
         node_positive_prop = round(node_positive_prop, 4))
  
breastcancer_clean
```

    ## # A tibble: 4,024 × 17
    ##      age race  marital_status t_stage n_stage x6th_stage differentiate     grade
    ##    <dbl> <fct> <fct>          <fct>   <fct>   <fct>      <fct>             <fct>
    ##  1    68 White Married        T1      N1      IIA        Poorly different… 3    
    ##  2    50 White Married        T2      N2      IIIA       Moderately diffe… 2    
    ##  3    58 White Divorced       T3      N3      IIIC       Moderately diffe… 2    
    ##  4    58 White Married        T1      N1      IIA        Poorly different… 3    
    ##  5    47 White Married        T2      N1      IIB        Poorly different… 3    
    ##  6    51 White Single         T1      N1      IIA        Moderately diffe… 2    
    ##  7    51 White Married        T1      N1      IIA        Well differentia… 1    
    ##  8    40 White Married        T2      N1      IIB        Moderately diffe… 2    
    ##  9    40 White Divorced       T4      N3      IIIC       Poorly different… 3    
    ## 10    69 White Married        T4      N3      IIIC       Well differentia… 1    
    ## # ℹ 4,014 more rows
    ## # ℹ 9 more variables: a_stage <fct>, tumor_size <dbl>, estrogen_status <fct>,
    ## #   progesterone_status <fct>, regional_node_examined <dbl>,
    ## #   reginol_node_positive <dbl>, survival_months <dbl>, status <fct>,
    ## #   node_positive_prop <dbl>

# Exploratory Analysis

## Checking Association Between Numerical Variables

``` r
# exploratory
pairs(breastcancer_clean)
```

![](model_building_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# correlation plot
breastcancer_num = breastcancer_clean|>
  dplyr::select(age, tumor_size, regional_node_examined, reginol_node_positive, node_positive_prop)

corrplot(cor(breastcancer_num), type = "upper", diag = FALSE)
```

![](model_building_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
# Boxplots for each variable
par(mfrow=c(2,3))
boxplot(breastcancer_clean$age, main='Age')
boxplot(breastcancer_clean$tumor_size, main='Tumor Size')
boxplot(breastcancer_clean$regional_node_examined,main='Node Examined' )
boxplot(breastcancer_clean$reginol_node_positive, main='Positive Node')
boxplot(breastcancer_clean$node_positive_prop, main='Proportion of Positive Nodes')
```

![](model_building_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

- `tumor_size` has substantial amounts of outliers. <br>
- `reginol_node_positive` and `node_positive_prop` are highly
  correlated.

## Checking Association Between Categorical Variables

``` r
breastcancer_cag = breastcancer_clean|>
  dplyr::select(-age, -tumor_size, -regional_node_examined, -reginol_node_positive, -node_positive_prop, -survival_months)

rules <- apriori(breastcancer_cag, parameter = list(supp = 0.001, conf = 0.8))
```

    ## Apriori
    ## 
    ## Parameter specification:
    ##  confidence minval smax arem  aval originalSupport maxtime support minlen
    ##         0.8    0.1    1 none FALSE            TRUE       5   0.001      1
    ##  maxlen target  ext
    ##      10  rules TRUE
    ## 
    ## Algorithmic control:
    ##  filter tree heap memopt load sort verbose
    ##     0.1 TRUE TRUE  FALSE TRUE    2    TRUE
    ## 
    ## Absolute minimum support count: 4 
    ## 
    ## set item appearances ...[0 item(s)] done [0.00s].
    ## set transactions ...[36 item(s), 4024 transaction(s)] done [0.00s].
    ## sorting and recoding items ... [36 item(s)] done [0.00s].
    ## creating transaction tree ... done [0.00s].
    ## checking subsets of size 1 2 3 4 5 6 7 8 9 10 done [0.08s].
    ## writing ... [424628 rule(s)] done [0.10s].
    ## creating S4 object  ... done [0.19s].

``` r
# Inspect the top 5 rules
inspect(head(sort(rules, by = "confidence"), 5))
```

    ##     lhs                                 rhs                                 support confidence   coverage       lift count
    ## [1] {grade=anaplastic; Grade IV}     => {differentiate=Undifferentiated} 0.00472167          1 0.00472167 211.789474    19
    ## [2] {differentiate=Undifferentiated} => {grade=anaplastic; Grade IV}     0.00472167          1 0.00472167 211.789474    19
    ## [3] {grade=anaplastic; Grade IV}     => {a_stage=Regional}               0.00472167          1 0.00472167   1.023398    19
    ## [4] {differentiate=Undifferentiated} => {a_stage=Regional}               0.00472167          1 0.00472167   1.023398    19
    ## [5] {x6th_stage=IIIB}                => {t_stage=T4}                     0.01665010          1 0.01665010  39.450980    67

# Checking Logistic Regression Assumption

# Fitting Model

## Backward Selection

``` r
breastcancer_clean1 <- breastcancer_clean|>
  drop_na()

full_model <- glm(status ~ ., data = breastcancer_clean, family = binomial())

summary(full_model)
```

    ## 
    ## Call:
    ## glm(formula = status ~ ., family = binomial(), data = breastcancer_clean)
    ## 
    ## Coefficients: (4 not defined because of singularities)
    ##                                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                        -0.440384   0.548792  -0.802 0.422287    
    ## age                                 0.028147   0.006515   4.320 1.56e-05 ***
    ## raceBlack                           0.454071   0.190407   2.385 0.017092 *  
    ## raceOther                          -0.422436   0.236111  -1.789 0.073592 .  
    ## marital_statusMarried              -0.060428   0.157866  -0.383 0.701883    
    ## marital_statusDivorced              0.104212   0.206690   0.504 0.614122    
    ## marital_statusSeparated             0.585305   0.480555   1.218 0.223232    
    ## marital_statusWidowed               0.147030   0.256999   0.572 0.567251    
    ## t_stageT2                           0.236842   0.229656   1.031 0.302405    
    ## t_stageT3                           0.744156   0.372432   1.998 0.045706 *  
    ## t_stageT4                           1.351247   0.571092   2.366 0.017978 *  
    ## n_stageN2                           0.649642   0.281493   2.308 0.021008 *  
    ## n_stageN3                           0.547902   0.356740   1.536 0.124573    
    ## x6th_stageIIB                       0.266978   0.268017   0.996 0.319191    
    ## x6th_stageIIIA                     -0.206690   0.341201  -0.606 0.544666    
    ## x6th_stageIIIB                     -0.017167   0.663988  -0.026 0.979373    
    ## x6th_stageIIIC                            NA         NA      NA       NA    
    ## differentiatePoorly differentiated  0.433317   0.122956   3.524 0.000425 ***
    ## differentiateUndifferentiated       1.717588   0.795600   2.159 0.030861 *  
    ## differentiateWell differentiated   -0.599953   0.207922  -2.885 0.003908 ** 
    ## grade2                                    NA         NA      NA       NA    
    ## grade3                                    NA         NA      NA       NA    
    ## gradeanaplastic; Grade IV                 NA         NA      NA       NA    
    ## a_stageRegional                     0.161619   0.325207   0.497 0.619208    
    ## tumor_size                         -0.003307   0.004731  -0.699 0.484495    
    ## estrogen_statusNegative             0.372422   0.227943   1.634 0.102294    
    ## progesterone_statusNegative         0.519439   0.152536   3.405 0.000661 ***
    ## regional_node_examined             -0.019948   0.011963  -1.667 0.095442 .  
    ## reginol_node_positive               0.058731   0.023440   2.506 0.012224 *  
    ## survival_months                    -0.061493   0.002764 -22.250  < 2e-16 ***
    ## node_positive_prop                  0.468375   0.375561   1.247 0.212349    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3444.7  on 4023  degrees of freedom
    ## Residual deviance: 2230.9  on 3997  degrees of freedom
    ## AIC: 2284.9
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
backward_model <- step(full_model, direction = "backward")
```

    ## Start:  AIC=2284.88
    ## status ~ age + race + marital_status + t_stage + n_stage + x6th_stage + 
    ##     differentiate + grade + a_stage + tumor_size + estrogen_status + 
    ##     progesterone_status + regional_node_examined + reginol_node_positive + 
    ##     survival_months + node_positive_prop
    ## 
    ## 
    ## Step:  AIC=2284.88
    ## status ~ age + race + marital_status + t_stage + n_stage + x6th_stage + 
    ##     differentiate + a_stage + tumor_size + estrogen_status + 
    ##     progesterone_status + regional_node_examined + reginol_node_positive + 
    ##     survival_months + node_positive_prop
    ## 
    ##                          Df Deviance    AIC
    ## - marital_status          4   2234.1 2280.1
    ## - x6th_stage              3   2233.6 2281.6
    ## - a_stage                 1   2231.1 2283.1
    ## - tumor_size              1   2231.4 2283.4
    ## - node_positive_prop      1   2232.4 2284.4
    ## <none>                        2230.9 2284.9
    ## - estrogen_status         1   2233.5 2285.5
    ## - regional_node_examined  1   2233.8 2285.8
    ## - t_stage                 3   2238.0 2286.0
    ## - n_stage                 1   2236.3 2288.3
    ## - reginol_node_positive   1   2237.3 2289.3
    ## - race                    2   2240.6 2290.6
    ## - progesterone_status     1   2242.1 2294.1
    ## - age                     1   2249.9 2301.9
    ## - differentiate           3   2261.3 2309.3
    ## - survival_months         1   2948.6 3000.6
    ## 
    ## Step:  AIC=2280.08
    ## status ~ age + race + t_stage + n_stage + x6th_stage + differentiate + 
    ##     a_stage + tumor_size + estrogen_status + progesterone_status + 
    ##     regional_node_examined + reginol_node_positive + survival_months + 
    ##     node_positive_prop
    ## 
    ##                          Df Deviance    AIC
    ## - x6th_stage              3   2236.7 2276.7
    ## - a_stage                 1   2234.3 2278.3
    ## - tumor_size              1   2234.5 2278.5
    ## - node_positive_prop      1   2235.7 2279.7
    ## <none>                        2234.1 2280.1
    ## - estrogen_status         1   2236.7 2280.7
    ## - regional_node_examined  1   2236.8 2280.8
    ## - t_stage                 3   2241.2 2281.2
    ## - n_stage                 1   2239.3 2283.3
    ## - reginol_node_positive   1   2240.5 2284.5
    ## - race                    2   2245.0 2287.0
    ## - progesterone_status     1   2245.8 2289.8
    ## - age                     1   2255.5 2299.5
    ## - differentiate           3   2264.2 2304.2
    ## - survival_months         1   2956.3 3000.3
    ## 
    ## Step:  AIC=2276.73
    ## status ~ age + race + t_stage + n_stage + differentiate + a_stage + 
    ##     tumor_size + estrogen_status + progesterone_status + regional_node_examined + 
    ##     reginol_node_positive + survival_months + node_positive_prop
    ## 
    ##                          Df Deviance    AIC
    ## - a_stage                 1   2237.0 2275.0
    ## - tumor_size              1   2237.1 2275.1
    ## - node_positive_prop      1   2238.3 2276.3
    ## <none>                        2236.7 2276.7
    ## - estrogen_status         1   2239.5 2277.5
    ## - regional_node_examined  1   2239.5 2277.5
    ## - n_stage                 2   2241.6 2277.6
    ## - reginol_node_positive   1   2243.3 2281.3
    ## - race                    2   2248.3 2284.3
    ## - t_stage                 3   2251.9 2285.9
    ## - progesterone_status     1   2248.3 2286.3
    ## - age                     1   2257.3 2295.3
    ## - differentiate           3   2266.7 2300.7
    ## - survival_months         1   2958.1 2996.1
    ## 
    ## Step:  AIC=2274.97
    ## status ~ age + race + t_stage + n_stage + differentiate + tumor_size + 
    ##     estrogen_status + progesterone_status + regional_node_examined + 
    ##     reginol_node_positive + survival_months + node_positive_prop
    ## 
    ##                          Df Deviance    AIC
    ## - tumor_size              1   2237.3 2273.3
    ## - node_positive_prop      1   2238.6 2274.6
    ## <none>                        2237.0 2275.0
    ## - n_stage                 2   2241.7 2275.7
    ## - regional_node_examined  1   2239.7 2275.7
    ## - estrogen_status         1   2239.7 2275.7
    ## - reginol_node_positive   1   2243.6 2279.6
    ## - race                    2   2248.6 2282.6
    ## - t_stage                 3   2252.1 2284.1
    ## - progesterone_status     1   2248.6 2284.6
    ## - age                     1   2257.7 2293.7
    ## - differentiate           3   2267.3 2299.3
    ## - survival_months         1   2958.2 2994.2
    ## 
    ## Step:  AIC=2273.28
    ## status ~ age + race + t_stage + n_stage + differentiate + estrogen_status + 
    ##     progesterone_status + regional_node_examined + reginol_node_positive + 
    ##     survival_months + node_positive_prop
    ## 
    ##                          Df Deviance    AIC
    ## - node_positive_prop      1   2238.9 2272.9
    ## <none>                        2237.3 2273.3
    ## - n_stage                 2   2241.8 2273.8
    ## - estrogen_status         1   2239.9 2273.9
    ## - regional_node_examined  1   2240.0 2274.0
    ## - reginol_node_positive   1   2243.9 2277.9
    ## - race                    2   2248.8 2280.8
    ## - progesterone_status     1   2248.9 2282.9
    ## - t_stage                 3   2259.0 2289.0
    ## - age                     1   2258.3 2292.3
    ## - differentiate           3   2267.6 2297.6
    ## - survival_months         1   2958.2 2992.2
    ## 
    ## Step:  AIC=2272.88
    ## status ~ age + race + t_stage + n_stage + differentiate + estrogen_status + 
    ##     progesterone_status + regional_node_examined + reginol_node_positive + 
    ##     survival_months
    ## 
    ##                          Df Deviance    AIC
    ## <none>                        2238.9 2272.9
    ## - estrogen_status         1   2241.6 2273.6
    ## - n_stage                 2   2245.8 2275.8
    ## - race                    2   2250.4 2280.4
    ## - progesterone_status     1   2250.7 2282.7
    ## - regional_node_examined  1   2254.4 2286.4
    ## - t_stage                 3   2261.7 2289.7
    ## - reginol_node_positive   1   2258.5 2290.5
    ## - age                     1   2260.1 2292.1
    ## - differentiate           3   2269.2 2297.2
    ## - survival_months         1   2961.8 2993.8

``` r
summary(backward_model)
```

    ## 
    ## Call:
    ## glm(formula = status ~ age + race + t_stage + n_stage + differentiate + 
    ##     estrogen_status + progesterone_status + regional_node_examined + 
    ##     reginol_node_positive + survival_months, family = binomial(), 
    ##     data = breastcancer_clean)
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                        -0.143815   0.402348  -0.357 0.720763    
    ## age                                 0.028709   0.006300   4.557 5.20e-06 ***
    ## raceBlack                           0.485245   0.185249   2.619 0.008808 ** 
    ## raceOther                          -0.452583   0.236372  -1.915 0.055530 .  
    ## t_stageT2                           0.351189   0.129003   2.722 0.006482 ** 
    ## t_stageT3                           0.525114   0.170942   3.072 0.002127 ** 
    ## t_stageT4                           1.279758   0.301169   4.249 2.14e-05 ***
    ## n_stageN2                           0.387559   0.148726   2.606 0.009164 ** 
    ## n_stageN3                           0.486693   0.274316   1.774 0.076030 .  
    ## differentiatePoorly differentiated  0.436269   0.122503   3.561 0.000369 ***
    ## differentiateUndifferentiated       1.684447   0.781285   2.156 0.031084 *  
    ## differentiateWell differentiated   -0.580959   0.206311  -2.816 0.004863 ** 
    ## estrogen_statusNegative             0.375306   0.226598   1.656 0.097667 .  
    ## progesterone_statusNegative         0.530750   0.151917   3.494 0.000476 ***
    ## regional_node_examined             -0.030898   0.008047  -3.839 0.000123 ***
    ## reginol_node_positive               0.078963   0.017838   4.427 9.58e-06 ***
    ## survival_months                    -0.061417   0.002750 -22.329  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3444.7  on 4023  degrees of freedom
    ## Residual deviance: 2238.9  on 4007  degrees of freedom
    ## AIC: 2272.9
    ## 
    ## Number of Fisher Scoring iterations: 6

## Forward Selection

``` r
null_model <- glm(status ~ 1, data = breastcancer_clean, family = binomial())
forward_model <- step(null_model, scope = list(lower = null_model, upper = full_model), direction = "forward")
```

    ## Start:  AIC=3446.68
    ## status ~ 1
    ## 
    ##                          Df Deviance    AIC
    ## + survival_months         1   2533.3 2537.3
    ## + x6th_stage              4   3197.4 3207.4
    ## + n_stage                 2   3214.0 3220.0
    ## + reginol_node_positive   1   3236.8 3240.8
    ## + node_positive_prop      1   3237.8 3241.8
    ## + progesterone_status     1   3335.1 3339.1
    ## + estrogen_status         1   3338.7 3342.7
    ## + differentiate           3   3337.3 3345.3
    ## + grade                   3   3337.3 3345.3
    ## + t_stage                 3   3349.2 3357.2
    ## + tumor_size              1   3380.1 3384.1
    ## + a_stage                 1   3415.7 3419.7
    ## + race                    2   3418.9 3424.9
    ## + marital_status          4   3419.2 3429.2
    ## + age                     1   3432.0 3436.0
    ## + regional_node_examined  1   3439.9 3443.9
    ## <none>                        3444.7 3446.7
    ## 
    ## Step:  AIC=2537.32
    ## status ~ survival_months
    ## 
    ##                          Df Deviance    AIC
    ## + x6th_stage              4   2382.1 2394.1
    ## + n_stage                 2   2395.4 2403.4
    ## + reginol_node_positive   1   2404.6 2410.6
    ## + node_positive_prop      1   2406.8 2412.8
    ## + differentiate           3   2461.0 2471.0
    ## + grade                   3   2461.0 2471.0
    ## + t_stage                 3   2471.5 2481.5
    ## + progesterone_status     1   2483.6 2489.6
    ## + estrogen_status         1   2498.6 2504.6
    ## + tumor_size              1   2498.9 2504.9
    ## + age                     1   2517.2 2523.2
    ## + race                    2   2518.5 2526.5
    ## + a_stage                 1   2523.4 2529.4
    ## + marital_status          4   2520.2 2532.2
    ## + regional_node_examined  1   2530.1 2536.1
    ## <none>                        2533.3 2537.3
    ## 
    ## Step:  AIC=2394.05
    ## status ~ survival_months + x6th_stage
    ## 
    ##                          Df Deviance    AIC
    ## + differentiate           3   2339.7 2357.7
    ## + grade                   3   2339.7 2357.7
    ## + progesterone_status     1   2345.6 2359.6
    ## + node_positive_prop      1   2355.2 2369.2
    ## + estrogen_status         1   2359.5 2373.5
    ## + age                     1   2363.0 2377.0
    ## + race                    2   2366.6 2382.6
    ## + reginol_node_positive   1   2370.2 2384.2
    ## + n_stage                 1   2374.5 2388.5
    ## + regional_node_examined  1   2374.9 2388.9
    ## + marital_status          4   2373.1 2393.1
    ## <none>                        2382.1 2394.1
    ## + t_stage                 3   2376.5 2394.5
    ## + a_stage                 1   2381.3 2395.3
    ## + tumor_size              1   2381.8 2395.8
    ## 
    ## Step:  AIC=2357.69
    ## status ~ survival_months + x6th_stage + differentiate
    ## 
    ##                          Df Deviance    AIC
    ## + node_positive_prop      1   2310.9 2330.9
    ## + progesterone_status     1   2314.8 2334.8
    ## + age                     1   2315.3 2335.3
    ## + reginol_node_positive   1   2326.6 2346.6
    ## + race                    2   2326.1 2348.1
    ## + estrogen_status         1   2328.1 2348.1
    ## + regional_node_examined  1   2331.8 2351.8
    ## + n_stage                 1   2332.8 2352.8
    ## + marital_status          4   2331.4 2357.4
    ## <none>                        2339.7 2357.7
    ## + a_stage                 1   2339.2 2359.2
    ## + tumor_size              1   2339.5 2359.5
    ## + t_stage                 3   2336.0 2360.0
    ## 
    ## Step:  AIC=2330.89
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop
    ## 
    ##                          Df Deviance    AIC
    ## + progesterone_status     1   2285.8 2307.8
    ## + age                     1   2289.0 2311.0
    ## + estrogen_status         1   2298.5 2320.5
    ## + race                    2   2298.2 2322.2
    ## + reginol_node_positive   1   2307.3 2329.3
    ## + n_stage                 1   2308.4 2330.4
    ## <none>                        2310.9 2330.9
    ## + marital_status          4   2303.1 2331.1
    ## + tumor_size              1   2310.5 2332.5
    ## + a_stage                 1   2310.5 2332.5
    ## + regional_node_examined  1   2310.9 2332.9
    ## + t_stage                 3   2307.1 2333.1
    ## 
    ## Step:  AIC=2307.82
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop + 
    ##     progesterone_status
    ## 
    ##                          Df Deviance    AIC
    ## + age                     1   2265.8 2289.8
    ## + race                    2   2273.7 2299.7
    ## + reginol_node_positive   1   2281.6 2305.6
    ## + n_stage                 1   2283.2 2307.2
    ## <none>                        2285.8 2307.8
    ## + estrogen_status         1   2284.4 2308.4
    ## + marital_status          4   2279.2 2309.2
    ## + tumor_size              1   2285.4 2309.4
    ## + a_stage                 1   2285.5 2309.5
    ## + regional_node_examined  1   2285.8 2309.8
    ## + t_stage                 3   2281.8 2309.8
    ## 
    ## Step:  AIC=2289.83
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop + 
    ##     progesterone_status + age
    ## 
    ##                          Df Deviance    AIC
    ## + race                    2   2254.2 2282.2
    ## + reginol_node_positive   1   2261.9 2287.9
    ## + n_stage                 1   2263.1 2289.1
    ## + estrogen_status         1   2263.3 2289.3
    ## <none>                        2265.8 2289.8
    ## + tumor_size              1   2264.9 2290.9
    ## + t_stage                 3   2261.2 2291.2
    ## + a_stage                 1   2265.7 2291.7
    ## + regional_node_examined  1   2265.8 2291.8
    ## + marital_status          4   2261.2 2293.2
    ## 
    ## Step:  AIC=2282.25
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop + 
    ##     progesterone_status + age + race
    ## 
    ##                          Df Deviance    AIC
    ## + reginol_node_positive   1   2249.9 2279.9
    ## + n_stage                 1   2251.8 2281.8
    ## + estrogen_status         1   2251.9 2281.9
    ## <none>                        2254.2 2282.2
    ## + tumor_size              1   2253.5 2283.5
    ## + t_stage                 3   2249.6 2283.6
    ## + regional_node_examined  1   2254.2 2284.2
    ## + a_stage                 1   2254.2 2284.2
    ## + marital_status          4   2251.1 2287.1
    ## 
    ## Step:  AIC=2279.88
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop + 
    ##     progesterone_status + age + race + reginol_node_positive
    ## 
    ##                          Df Deviance    AIC
    ## + estrogen_status         1   2247.6 2279.6
    ## + regional_node_examined  1   2247.8 2279.8
    ## <none>                        2249.9 2279.9
    ## + t_stage                 3   2244.2 2280.2
    ## + n_stage                 1   2248.5 2280.5
    ## + tumor_size              1   2248.6 2280.6
    ## + a_stage                 1   2249.8 2281.8
    ## + marital_status          4   2246.8 2284.8
    ## 
    ## Step:  AIC=2279.62
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop + 
    ##     progesterone_status + age + race + reginol_node_positive + 
    ##     estrogen_status
    ## 
    ##                          Df Deviance    AIC
    ## + regional_node_examined  1   2245.4 2279.4
    ## <none>                        2247.6 2279.6
    ## + t_stage                 3   2241.7 2279.7
    ## + n_stage                 1   2246.3 2280.3
    ## + tumor_size              1   2246.4 2280.4
    ## + a_stage                 1   2247.6 2281.6
    ## + marital_status          4   2244.5 2284.5
    ## 
    ## Step:  AIC=2279.39
    ## status ~ survival_months + x6th_stage + differentiate + node_positive_prop + 
    ##     progesterone_status + age + race + reginol_node_positive + 
    ##     estrogen_status + regional_node_examined
    ## 
    ##                  Df Deviance    AIC
    ## <none>                2245.4 2279.4
    ## + t_stage         3   2239.7 2279.7
    ## + n_stage         1   2243.7 2279.7
    ## + tumor_size      1   2244.2 2280.2
    ## + a_stage         1   2245.3 2281.3
    ## + marital_status  4   2242.2 2284.2

``` r
summary(forward_model)
```

    ## 
    ## Call:
    ## glm(formula = status ~ survival_months + x6th_stage + differentiate + 
    ##     node_positive_prop + progesterone_status + age + race + reginol_node_positive + 
    ##     estrogen_status + regional_node_examined, family = binomial(), 
    ##     data = breastcancer_clean)
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                        -0.389517   0.434745  -0.896 0.370270    
    ## survival_months                    -0.061291   0.002741 -22.358  < 2e-16 ***
    ## x6th_stageIIB                       0.453360   0.162581   2.789 0.005295 ** 
    ## x6th_stageIIIA                      0.566815   0.170656   3.321 0.000896 ***
    ## x6th_stageIIIB                      1.388888   0.373839   3.715 0.000203 ***
    ## x6th_stageIIIC                      0.762973   0.292623   2.607 0.009124 ** 
    ## differentiatePoorly differentiated  0.444345   0.122099   3.639 0.000273 ***
    ## differentiateUndifferentiated       1.824929   0.752241   2.426 0.015267 *  
    ## differentiateWell differentiated   -0.601322   0.207578  -2.897 0.003769 ** 
    ## node_positive_prop                  0.593517   0.365126   1.626 0.104053    
    ## progesterone_statusNegative         0.530136   0.151884   3.490 0.000482 ***
    ## age                                 0.027713   0.006268   4.421 9.82e-06 ***
    ## raceBlack                           0.475487   0.185035   2.570 0.010178 *  
    ## raceOther                          -0.471807   0.235491  -2.004 0.045123 *  
    ## reginol_node_positive               0.058102   0.023163   2.508 0.012128 *  
    ## estrogen_statusNegative             0.353398   0.226096   1.563 0.118042    
    ## regional_node_examined             -0.017351   0.011789  -1.472 0.141091    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3444.7  on 4023  degrees of freedom
    ## Residual deviance: 2245.4  on 4007  degrees of freedom
    ## AIC: 2279.4
    ## 
    ## Number of Fisher Scoring iterations: 6

## Discussion

AIC of backward selection is 2992.2 compared to AIC = 2996.4 for forward
selection. The difference in AIC is approximately 4 for the backward and
forward selction models. This difference suggests less support for the
model with higher AIC. Also,
