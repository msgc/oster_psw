#Checking for packages
list.of.packages <- c("ggplot2", "broom", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(broom)
library(dplyr)
library(ggplot2)

# The dataset for this example contains treated cases and control cases: the treated cases are those who actually received job training under the NSW experiment, and the control cases include a large number of people from the PSID, a nationally representative survey of income and employment. Since the PSID is a representative sample of the whole population while the NSW experiment was a program targeted at poor workers, you should expect that the treatment group will look very different from the control group. Remember that this is exactly the problem that matching estimators are designed to correct.

a<-read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", "10X6BkKJaEvD8ntEMHfyqFqk4goW7PNII"))

# Data
# Income in 1978 (one or two years after the job training) is our outcome variable here, and the
# treatment is whether the individual received job training. The data contain a number of demographic
# and work-related variables designed to control for other factors that might affect income. The
# variables included are:
# Variable Description
# treatment 1 if in NSW Job Training Program; 0 if PSID control case
# age Age
# educ Education (# of years)
# black 1 if black
# hisp 1 if Hispanic
# married 1 if married
# nodegree 1 if no high school degree
# re74 1974 earnings Baseline, before training
# re75 1975 earnings assume during training
# re78 1978 earnings (outcome)
# age_sq Age squared
# educ_sq Education squared
# re74_sq 1974 earnings squared
# re75_sq 1975 earnings squared
# u74 Unemployed in 1974
# u74_black Black and Unemployed in 1974
# u75 Unemployed in 1975
# u74_hisp Hispanic and Unemployed in 1974
# w PSW obtained with TWANG

o_delta_psw <- 
function (y, x, con, w, m = "none", id = "none", time = "none", 
    beta = 0, R2max, type, data) 
{
    if (type == "lm") {
        data <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))
    }
    if (type == "plm") {
        data <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w), 
            id_var = all_of(id), time_var = all_of(time))
        data_plm <- pdata.frame(data, index = c("id_var", 
            "time_var"), drop.index = TRUE, row.names = TRUE)
    }
    if (m == "none") {
        model0_formula <- as.formula("y_var ~ x_var")
        model1_formula <- as.formula(paste("y_var ~ x_var +", 
            con))
        aux_model_formula <- as.formula(paste("x_var ~", 
            con))
        if (type == "plm") {
            sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))
        }
    }
    else {
        model0_formula <- as.formula(paste("y_var ~ x_var +", 
            m))
        model1_formula <- as.formula(paste("y_var ~ x_var +", 
            con, "+", m))
        aux_model_formula <- as.formula(paste("x_var ~", 
            con, "+", m))
        if (type == "lm") {
            sigma_xx_model_formula <- as.formula(paste("x_var ~ ", 
                m))
        }
        if (type == "plm") {
            sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var) +", 
                m))
        }
    }
    if (type == "lm") {
        model0 <- lm(model0_formula, data = data, na.action = na.exclude, weights=w_var)
        model1 <- lm(model1_formula, data = data, na.action = na.exclude, weights=w_var)
        aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude, weights=w_var)
        if (m != "none") {
            model_xx <- lm(sigma_xx_model_formula, data = data, 
                na.action = na.exclude, weights=w_var)
        }
    }
    if (type == "plm") {
        model0 <- plm(model0_formula, data = data_plm, model = "within", 
            na.action = na.exclude, weights=w_var)
        model1 <- plm(model1_formula, data = data_plm, model = "within", 
            na.action = na.exclude, weights=w_var)
        aux_model <- plm(aux_model_formula, data = data_plm, 
            model = "within", na.action = na.exclude, weights=w_var)
        model_xx <- lm(sigma_xx_model_formula, data = data, na.action = na.exclude, weights=w_var)
    }
    if (type == "lm") {
        b0 = as.numeric(tidy(model0)[2, 2])
    }
    else if (type == "plm") {
        b0 = as.numeric(tidy(model0)[1, 2])
    }
    if (type == "lm") {
        b1 = as.numeric(tidy(model1)[2, 2])
    }
    else if (type == "plm") {
        b1 = as.numeric(tidy(model1)[1, 2])
    }
    if (type == "lm") {
        R20 = summary(model0)$r.squared
    }
    else if (type == "plm") {
        R20 = as.numeric(summary(model0)$r.squared[1])
    }
    if (type == "lm") {
        R21 = summary(model1)$r.squared
    }
    else if (type == "plm") {
        R21 = as.numeric(summary(model1)$r.squared[1])
    }
    sigma_yy = var(data$y_var, na.rm = T)
    if (m == "none") {
        if (type == "lm") {
            sigma_xx = var(data$x_var, na.rm = T)
        }
        else if (type == "plm") {
            sigma_xx = var(model_xx$residuals)
        }
    }
    else {
        if (type == "lm") {
            sigma_xx = var(model_xx$residuals)
        }
        else if (type == "plm") {
            sigma_xx = var(model_xx$residuals)
        }
    }
    t_x = var(aux_model$residuals)
    bt_m_b = b1 - beta
    rt_m_ro_t_syy = (R21 - R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    num = num1 + num2 + num3 + num4
    den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
    den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    den = den1 + den2 + den3 + den4
    delta_star = num/den
    result_delta <- tribble(~Name, ~Value, "delta*", round(delta_star, 
        6), "Uncontrolled Coefficient", b0, "Controlled Coefficient", 
        b1, "Uncontrolled R-square", R20, "Controlled R-square", 
        R21, "Max R-square", R2max, "beta hat", beta)
    if (R21 > R2max) 
        warning("The max R-square value is smaller than the R-square of the controlled model")
    return(result_delta)
}

oster_psw_dist<-function (y, x, con, w, m = "none", id = "none", time = "none", 
    beta = 0, type, data) 
{
    if (type == "lm") {
        data_aux <- data %>% dplyr::rename(y_var = all_of(y), 
            x_var = all_of(x), w_var = all_of(w))
    }
    if (type == "plm") {
        data_aux <- data %>% dplyr::rename(y_var = all_of(y), 
            x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
        data_plm_aux <- pdata.frame(data_aux, index = c("id_var", 
            "time_var"), drop.index = TRUE, row.names = TRUE)
    }
    if (m == "none") {
        model1_formula <- as.formula(paste("y_var ~ x_var +", 
            con))
    }
    else {
        model1_formula <- as.formula(paste("y_var ~ x_var +", 
            con, "+", m))
    }
    if (type == "lm") {
        model1 <- lm(model1_formula, data = data_aux, na.action = na.exclude, weights=w_var)
    }
    if (type == "plm") {
        model1 <- plm(model1_formula, data = data_plm_aux, model = "within", 
            na.action = na.exclude, weights=w_var)
    }
    if (type == "lm") {
        R21_aux = summary(model1)$r.squared
    }
    else if (type == "plm") {
        R21_aux = as.numeric(summary(model1)$r.squared[1])
    }
    r_start <- ceiling(R21_aux/0.01) * 0.01
    n_runs <- length(seq(r_start:1, by = 0.01))
    delta_over_rsq_results <- tibble(x = 1:n_runs, `Max R-square` = 0, 
        `delta*` = 0)
    for (i in 1:n_runs) {
        R2max = seq(r_start:1, by = 0.01)[i]
        delta_results_aux <- o_delta_psw(y = y, x = x, con = con, w,
            m = m, id = id, time = time, beta = beta, R2max = R2max, 
            type = type, data = data)
        delta_over_rsq_results[i, 2] <- R2max
        delta_over_rsq_results[i, 3] <- round(delta_results_aux[1, 
            2], 6)
    }
    return(delta_over_rsq_results)
}

o_delta_rsq_viz_psw <- function (y, x, w, con, m = "none", id = "none", time = "none", 
    beta = 0, type, data) 
{
    if (type == "lm") {
        data_aux <- data %>% dplyr::rename(y_var = all_of(y), 
            x_var = all_of(x), w_var = all_of(w))
    }
    if (type == "plm") {
        data_aux <- data %>% dplyr::rename(y_var = all_of(y), 
            x_var = all_of(x), w_var = all_of(w), id_var = all_of(id), time_var = all_of(time))
        data_plm_aux <- pdata.frame(data_aux, index = c("id_var", 
            "time_var"), drop.index = TRUE, row.names = TRUE)
    }
    if (m == "none") {
        model1_formula <- as.formula(paste("y_var ~ x_var +", 
            con))
    }
    else {
        model1_formula <- as.formula(paste("y_var ~ x_var +", 
            con, "+", m))
    }
    if (type == "lm") {
        model1 <- lm(model1_formula, data = data_aux, na.action = na.exclude, weights=w_var)
    }
    if (type == "plm") {
        model1 <- plm(model1_formula, data = data_plm_aux, model = "within", 
            na.action = na.exclude, weights=w_var)
    }
    if (type == "lm") {
        R21_aux = summary(model1)$r.squared
    }
    else if (type == "plm") {
        R21_aux = as.numeric(summary(model1)$r.squared[1])
    }
    r_start <- ceiling(R21_aux/0.01) * 0.01
    n_runs <- length(seq(r_start:1, by = 0.01))
    delta_over_rsq_results <- tibble(x = 1:n_runs, rmax = 0, 
        delta = 0)
    for (i in 1:n_runs) {
        R2max = seq(r_start:1, by = 0.01)[i]
        delta_results_aux <- o_delta_psw(y = y, x = x, w = w, con = con, 
            m = m, id = id, time = time, beta = beta, R2max = R2max, 
            type = type, data = data)
        delta_over_rsq_results[i, 2] <- R2max
        delta_over_rsq_results[i, 3] <- round(delta_results_aux[1, 
            2], 6)
    }
    theme_set(theme_bw())
    result_plot <- ggplot(data = delta_over_rsq_results, aes(x = rmax, 
        y = delta)) + geom_line(size = 1.3) + scale_y_continuous(name = expression(delta^"*")) + 
        scale_x_continuous(name = expression("maximum R"^2)) + 
        theme(axis.title = element_text(size = 15), axis.text = element_text(size = 13))
    return(result_plot)
}

# Examples
# No PSW in the calculation of delta
# lm1<-lm(re78 ~ treatment + age + age_sq + educ + educ_sq + married + nodegree + black + hisp + re74 + re75 + re74_sq + re75_sq + u74 + u75 + u74_hisp + u74_black, data=lalonde, weights=NULL)
# a$w1<-1 #assigns a weight of 1
# o_delta_psw(y = "re78",           # dependent variable
# x = "treatment",            # independent treatment variable
# con = "age + age_sq + educ + educ_sq + married + nodegree + black + hisp + re74 + re75 + re74_sq + re75_sq + u74 + u75 + u74_hisp + u74_black",   # other control variables
# w = "w1",            # PSW
# beta = 0,            # beta
# type = "lm",         # model type
# R2max = summary(lm1)$r.squared*1.3,         # maximum R-square
# data = a)   # dataset

# # PSW in the calculation of delta
# lm1<-lm(re78 ~ treatment + age + age_sq + educ + educ_sq + married + nodegree + black + hisp + re74 + re75 + re74_sq + re75_sq + u74 + u75 + u74_hisp + u74_black, data=lalonde, weights=w)
# a$w1<-1 #assigns a weight of 1
# o_delta_psw(y = "re78",           # dependent variable
# x = "treatment",            # independent treatment variable
# con = "age + age_sq + educ + educ_sq + married + nodegree + black + hisp + re74 + re75 + re74_sq + re75_sq + u74 + u75 + u74_hisp + u74_black",   # other control variables
# w = "w",            # PSW
# beta = 0,            # beta
# type = "lm",         # model type
# R2max = summary(lm1)$r.squared*1.3,         # maximum R-square
# data = a)   # dataset

# # Distribution of delta
# oster_psw_dist(y = "re78",           # dependent variable
# x = "treatment",            # independent treatment variable
# con = "age + age_sq + educ + educ_sq + married + nodegree + black + hisp + re74 + re75 + re74_sq + re75_sq + u74 + u75 + u74_hisp + u74_black",   # other control variables
# w = "w",            # PSW
# beta = 0,            # beta
# type = "lm",         # model type
# data = a)   # dataset

# # Visualization of distribution of delta
# o_delta_rsq_viz_psw(y = "re78",           # dependent variable
# x = "treatment",            # independent treatment variable
# con = "age + age_sq + educ + educ_sq + married + nodegree + black + hisp + re74 + re75 + re74_sq + re75_sq + u74 + u75 + u74_hisp + u74_black",   # other control variables
# w = "w",            # PSW
# beta = 0,            # beta
# type = "lm",         # model type
# data = a)   
