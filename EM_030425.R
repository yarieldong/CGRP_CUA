#3/4/25: Base-case model comparing eptinezumab and erenumab with divalproex ER

library(diagram)
library(ggplot2)
library(rgho)
library(dplyr)
library(plyr)
library(tidyverse)
library(bayestestR)
library(darthtools)
library(heemod)
library(gridExtra)

#Parameters Calculation##########################

###Life table & cohort demographics####
age_init = 39.8
age_max = 50
sex <- 0.84 #84% of the patient are females, male = 0, female = 1
life_table = read.csv("life_table.csv", header = TRUE)
mr_by_age <- life_table %>% #return a data frame with age and corresponding mortality of female, male and total mortality
  dplyr::filter(Age >= 39 & Age < age_max) %>%
  dplyr::select(Age, Female, Male) %>%
  as.data.frame()

#convert probability of dying for 1 year to annual rate
annual_mr <- data.frame(
  age = mr_by_age$Age,
  female = prob_to_rate(mr_by_age$Female, t = 12),
  male = prob_to_rate(mr_by_age$Male, t = 12)
)

#convert yearly mortality rate  to monthly probability
p_mr_monthly <- data.frame(
  age = annual_mr$age,
  female = rate_to_prob(annual_mr$female, t = 1),
  male = rate_to_prob(annual_mr$male, t= 1)
)


MMD_baseline = 6.3*(122/(122+223+286)) + 8.7*(223/(122+223+286)) + 8.1*(286/(122+223+286))
MMD_dival = MMD_baseline - 1.7 #MMD on divalproex (SOC)
MMD_eptin = MMD_baseline - 3.9 #MMD on eptinezumab
MMD_erenu = MMD_baseline - 2.9 #MMD on erenumab

###Treatment intervention costs####
c_dival = 21.66/100 * 28 * 2 #cost of divalproex ER per cycle
c_admin_dival = 0 #administration cost of divalproex ER per cycle

c_eptin = 1585.68 / 3 #cost of eptinezumab per cycle
c_admin_eptin = (62.58 + 68.37) / 3 #administration cost of eptinezumab per cycle

c_erenu = 670.53 #cost of eptinezumab per cycle
c_admin_erenu = 0 #administration cost of eptinezumab per cycle

###Healthcare utilization costs####
c_rescue_med = 9 #cost of rescue medication per migraine day

c_hospitalization = 7364 #cost of migraine-related hospitalization
p_hospitalization_female = 0.002 #probability of migraine-related hospitalization per migraine attack for female
p_hospitalization_male = 0.003#probability of migraine-related hospitalization per migraine attack for male

c_outpt = 253 #cost of migraine-related outpatient provider visit
p_outpt_female = 0.064 #probability of outpatient provider visit per migraine attack for female
p_outpt_male = 0.031 ##probability of outpatient provider visit per migraine attack for male

c_ed = 1115 #cost of migraine-related ED visit
p_ed_female = 0.013 #probability of ED visit per migraine attack for female
p_ed_male = 0.007 #probability of ED visit per migraine attack for male

c_other_cost_on_dival = MMD_dival * (c_rescue_med + sex * (p_ed_female * c_ed + p_outpt_female * c_outpt + p_hospitalization_female * c_hospitalization) + (1-sex) * (p_ed_male * c_ed + p_outpt_male* c_outpt + p_hospitalization_male * c_hospitalization)) #other costs than medication and administration costs while on divalproex per cycle

c_other_cost_on_eptin = MMD_eptin * (c_rescue_med + sex * (p_ed_female * c_ed + p_outpt_female * c_outpt + p_hospitalization_female * c_hospitalization) + (1-sex) * (p_ed_male * c_ed + p_outpt_male* c_outpt + p_hospitalization_male * c_hospitalization)) #other costs than medication and administration costs while on eptinezumab per cycle

c_other_cost_on_erenu = MMD_erenu * (c_rescue_med + sex * (p_ed_female * c_ed + p_outpt_female * c_outpt + p_hospitalization_female * c_hospitalization) + (1-sex) * (p_ed_male * c_ed + p_outpt_male* c_outpt + p_hospitalization_male * c_hospitalization)) #other costs than medication and administration costs while on erenumab per cycle

c_other_cost_npt = MMD_baseline * (c_rescue_med + sex * (p_ed_female * c_ed + p_outpt_female * c_outpt + p_hospitalization_female * c_hospitalization) + (1-sex) * (p_ed_male * c_ed + p_outpt_male* c_outpt + p_hospitalization_male * c_hospitalization)) #other costs than medication and administration costs while not on preventive treatments
###Decision Tree parameters####
dr_dt = 0.015 #discount rate of the decision tree portion of the model - 3% annual discount rate -> 24 weeks

##Utility####
#baseline utility score
u_baseline = (0.7938 - 0.0189 * MMD_baseline) / 12

#utility score on divalproex
u_on_dival = (0.7938 - 0.0189 * MMD_dival) / 12

#utility scores of eptin
u_on_eptin = (0.7938 - 0.0189 * MMD_eptin) / 12

#utility scores of erenu
u_on_erenu = (0.7938 - 0.0189 * MMD_erenu) / 12

##transition probabilities####
p_dc_dival = 0.0028
p_dc_erenu = 0.0055
p_dc_eptin = 0.0091

p_clinical_response_dival = 0.51 #Probability of achieving 50% reduction in MMD during first 6 months of tx
p_clinical_response_eptin = 0.75

#Base Case Model##########################

###Defining Model Parameters####
param <- define_parameters(

  ##demographics
  age_init = 39.8, #age of patients at cycle 0 (end of decision tree model)
  sex = 0.84, ##male = 0, female = 1, 84% of patients are female

  age = age_init + model_time/12, #current age of patient, model_time default to 1 in heemod

  #mortality
  p_mr = sex*look_up(p_mr_monthly, age = age, value = "female", bin = TRUE) + (1-sex)*look_up(p_mr_monthly, age = age, value = "male", bin = TRUE),

  ##transition probabilities of divalproex ER
  p_dc_dival = 0.0028,
  p_stay_on_tx_dival = 1 - p_dc_dival - p_mr,
  p_stay_off_tx_dival = 1 - p_mr,

  ##transition probabilities of eptin
  p_dc_eptin = 0.0091,
  p_stay_on_tx_eptin = 1 - p_dc_eptin - p_mr,
  p_stay_off_tx_eptin = 1 - p_mr,

  ##transition probability of erenu
  p_dc_erenu = 0.0055,
  p_stay_on_tx_erenu = 1 - p_dc_erenu - p_mr,
  p_stay_off_tx_erenu = 1 - p_mr,

  ##baseline utility score
  u_baseline = (0.7938 - 0.0189 * MMD_baseline) / 12,

  ##utility score on divalproex
  u_on_dival = (0.7938 - 0.0189 * MMD_dival) / 12,

  ##utility scores of eptin
  u_on_eptin = (0.7938 - 0.0189 * MMD_eptin) / 12,

  ##utility scores of erenu
  u_on_erenu = (0.7938 - 0.0189 * MMD_erenu) / 12,

  ##costs
  c_dival = c_dival,
  c_admin_dival = c_admin_dival,
  c_eptin = c_eptin,
  c_admin_eptin = c_admin_eptin,
  c_erenu = c_erenu,
  c_other_cost_on_erenu = c_other_cost_on_erenu,
  c_other_cost_npt = c_other_cost_npt

)

#defining state names
state_names <- c("state_on_preventive_tx", "state_off_preventive_tx", "state_death")

###Divalproex model####

#Define transition matrix
p_dival <- define_transition(
  p_stay_on_tx_dival, p_dc_dival, p_mr,
  0, p_stay_off_tx_dival, p_mr,
  0, 0, 1,
  state_names = state_names,
)

strat_dival <- define_strategy(
  state_on_preventive_tx = define_state(
    c_drug = discount(c_dival, 0.03/12),
    c_dt = (5.5*(c_dival+ c_admin_dival) + 6*(c_other_cost_on_dival))/(1+dr_dt),
    #cost of medications during the decision tree phase, patient only received half dose the first week
    cost_total = discount(c_dival + c_admin_dival + c_other_cost_on_dival, 0.03/12) + c_dt/120,
    utility = discount(u_on_dival, 0.03/12)
  ),

  state_off_preventive_tx = define_state(
    c_drug = 0,
    c_dt = (5.5*(c_dival+ c_admin_dival) + 6*(c_other_cost_on_dival))/(1+dr_dt),
    cost_total = discount(c_other_cost_npt, 0.03/12) + c_dt/120,
    utility = discount(u_baseline, 0.03/12)
  ),

  state_death = define_state(
    c_drug = 0,
    c_dt = (5.5*(c_dival+ c_admin_dival) + 6*(c_other_cost_on_dival))/(1+dr_dt),
    cost_total = c_dt/120,
    utility = 0
  ),

  transition = p_dival
)

#run model
dival_mod <- run_model(
  dival = strat_dival,
  init = c(0.51, 0.49, 0), #number of patients entering each state at the beginning of Markov model)
  cycles = 120,
  cost = cost_total,
  parameters = param,
  effect = utility
)

###Eptinezumab model####
#Eptinezumab
#define transition matrix
p_eptin <- define_transition(
  p_stay_on_tx_eptin, p_dc_eptin, p_mr,
  0, p_stay_off_tx_eptin, p_mr,
  0, 0, 1,
  state_names = state_names,
)


#define the cost and utility associated with each state
strat_eptin <- define_strategy(
  state_on_preventive_tx = define_state(
    c_drug = c_eptin,
    c_dt = 6*(c_eptin + c_admin_eptin + c_other_cost_on_eptin)/(1+dr_dt), #cost of medications during the decision tree phase
    cost_total = discount(c_eptin + c_admin_eptin + c_other_cost_on_eptin, 0.03/12) + c_dt/120,
    utility = discount(u_on_eptin, 0.03/12)
  ),

  state_off_preventive_tx = define_state(
    c_drug = 0,
    c_dt = 6*(c_eptin + c_admin_eptin + c_other_cost_on_eptin)/(1+dr_dt),
    cost_total = discount(c_other_cost_npt, 0.03/12) + c_dt/120,
    utility = discount(u_baseline, 0.03/12)
  ),

  state_death = define_state(
    c_drug = 0,
    c_dt = 6*(c_eptin + c_admin_eptin + c_other_cost_on_eptin)/(1+dr_dt),
    cost_total = c_dt/120,
    utility = 0
  ),

  transition = p_eptin
)

#run model
eptin_mod <- run_model(
  eptin = strat_eptin,
  init = c(0.75, 0.25, 0), #number of patients entering each state at the beginning of Markov model)
  cycles = 120,
  cost = cost_total,
  parameters = param,
  effect = utility
)

###Erenumab model####
#erenumab
#Transition Matrix
p_erenu <- define_transition(
  p_stay_on_tx_erenu, p_dc_erenu, p_mr,
  0, p_stay_off_tx_erenu, p_mr,
  0, 0, 1,
  state_names = state_names,
)

#define the cost and utility associated with each state
strat_erenu <- define_strategy(
  state_on_preventive_tx = define_state(
    c_drug = c_erenu,
    c_dt = 6*(c_erenu + c_admin_erenu + c_other_cost_on_erenu)/(1+dr_dt),
    cost_total = discount(c_erenu + c_other_cost_on_erenu, 0.03/12) + c_dt/120,
    utility = discount(u_on_erenu, 0.03/12)
  ),

  state_off_preventive_tx = define_state(
    c_drug = 0,
    c_dt = 6*(c_erenu + c_admin_erenu + c_other_cost_on_erenu)/(1+dr_dt),
    cost_total = discount(c_other_cost_npt, 0.03/12) + c_dt/120,
    utility = discount(u_baseline, 0.03/12)
  ),

  state_death = define_state(
    c_drug = 0,
    c_dt = 6*(c_erenu + c_admin_erenu + c_other_cost_on_erenu)/(1+dr_dt),
    cost_total = 0 + c_dt/120,
    utility = 0
  ),

  transition = p_erenu
)

#run model
erenu_mod <- run_model(
  erenu = strat_erenu,
  init = c(0.64, 0.36, 0),
  cycles = 120,
  cost = cost_total,
  parameters = param,
  effect = utility
)

#Base Case Result##################################
###CE Frontier####
ce_frontier <- data.frame(
  strategy = c("Divalproex","Eptinezumab", "Erenumab"),
  cost = c(dival_mod$run_model$cost_total, eptin_mod$run_model$cost_total, erenu_mod$run_model$cost_total),
  utility = c(dival_mod$run_model$utility, eptin_mod$run_model$utility, erenu_mod$run_model$utility)
)
ce_frontier

#plotting the cost-effectiveness frontier
ggplot(ce_frontier) +
  geom_point(aes(x=utility, y=cost, shape = strategy), size=2.5) +
  geom_line(aes(x=utility, y=cost, group=1),color ="black") +
  theme_bw()

###Determistic Results####
deterministic_result <- data.frame(
  strategy = c("Eptinezumab", "Erenumab", "Divalproex"),
  cost = c(eptin_mod$run_model$cost_total, erenu_mod$run_model$cost_total, dival_mod$run_model$cost_total),
  drug_cost = c(eptin_mod$run_model$c_drug + (6*(c_eptin + c_admin_eptin))/(1 + dr_dt), erenu_mod$run_model$c_drug + (6*(c_erenu + c_admin_eptin))/(1 + dr_dt), dival_mod$run_model$c_drug + (5.5*(c_dival + c_admin_dival))/(1 + dr_dt)),
  QALY = c(eptin_mod$run_model$utility, erenu_mod$run_model$utility, dival_mod$run_model$utility),
  incremental_cost = c(eptin_mod$run_model$cost_total - dival_mod$run_model$cost_total, erenu_mod$run_model$cost_total - dival_mod$run_model$cost_total, NA),
  incremental_qaly = c(eptin_mod$run_model$utility - dival_mod$run_model$utility, erenu_mod$run_model$utility - dival_mod$run_model$utility, NA)
)

deterministic_result$ICER <- c(deterministic_result[1,5]/deterministic_result[1,6], deterministic_result[2,5]/deterministic_result[2,6],  NA)

deterministic_result$iNMB <- c(deterministic_result[1,6]*150000 - deterministic_result[1,5], deterministic_result[2,6]*150000 - deterministic_result[2,5], NA)

deterministic_result

#Probablistic Sensitivity Analysis#####################################
set.seed(123)
##Defining PSA parameters########
psa_param <- define_psa(

  c_dival ~ gamma(mean = c_dival, sd = 0.25*c_dival),
  c_eptin ~ gamma(mean = c_eptin, sd = 0.25*c_eptin),
  c_erenu ~ gamma(mean = c_erenu, sd = 0.25*c_erenu),

  p_dc_dival ~ beta(15.95, 5681.33),
  p_dc_eptin ~ beta(15.8453, 1725.396),
  p_dc_erenu ~ beta(15.9065, 2876.184),

  u_on_dival ~ beta(15.04, 252.55),
  u_on_eptin ~ beta(15.00, 238.78),
  u_on_erenu ~ beta(15.02, 246.14)
)

pm_dival <- run_psa(
  model = dival_mod,
  psa = psa_param,
  N = 10000 #10000 simulations
)

pm_eptin <- run_psa(
  model = eptin_mod,
  psa = psa_param,
  N = 10000
)

pm_erenu <- run_psa(
  model = erenu_mod,
  psa = psa_param,
  N = 10000
)

#PSA Result###############################
psa_result <- data.frame(
  eptin_cost = c(pm_eptin$psa$cost_total),
  eptin_qaly = c(pm_eptin$psa$utility),
  dival_cost = c(pm_dival$psa$cost_total),
  dival_qaly = c(pm_dival$psa$utility),
  incr_cost = c(pm_eptin$psa$cost_total-pm_dival$psa$cost_total),
  incr_qaly = c(pm_eptin$psa$utility-pm_dival$psa$utility),
  inmb = c((pm_eptin$psa$utility - pm_dival$psa$utility)*150000 - (pm_eptin$psa$cost_total-pm_dival$psa$cost_total))
)

psa_result$icer <- c(psa_result$incr_cost/psa_result$incr_qaly)

##Scatterplot of PSA##########################
ylim <- max(psa_result$incr_cost) * 1.2
xlim <- max(psa_result$incr_qaly) * 1.2
ggplot(psa_result,
       aes(x = incr_qaly, y = incr_cost)) +
  geom_jitter(size = .5)  +
  xlab("Incremental QALYs") +
  ylab("Incremental cost") +
  scale_y_continuous(limits = c(-ylim, ylim),
                     labels = scales::dollar_format()) +
  scale_x_continuous(limits = c(-xlim, xlim)) +
  geom_abline(slope = 150000, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  annotate("text", x = 0.23, y = 75000, label = "WTP Threshold of $150,000", angle = 0, size = 3) +
  theme_classic()

##cost-effectiveness acceptability curve####
wtp_threshold <- 0:300 * 1000

compute_p_ceac <- function(psa_result, threshold){
  nmb <- psa_result$incr_qaly * threshold - psa_result$incr_cost
  p = sum(nmb >= 0)/10000
  return(p)
}

probability_ceac <- matrix(NA, nrow = 301, ncol=1)
for (i in 1:301) {
  probability_ceac[i] <- compute_p_ceac(psa_result, wtp_threshold[i])
}

ceac <- cbind(data.frame(wtp_threshold), data.frame(probability_ceac))
ceac

##plotting CEAC
ceac_plot <- ggplot(ceac, mapping = aes(x = wtp_threshold)) +
  labs(color = NULL, caption = "(a) Cost-effectiveness acceptability curve") +
  xlab("Willingness to pay ($/QALY)") +
  ylab("Probability of Cost-effectiveness")+
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE)) +
  geom_line(mapping = aes(x = wtp_threshold, y = probability_ceac, color = "Eptinezumab"), linewidth = 0.8) +
  geom_line(mapping = aes(x = wtp_threshold, y = 1 - probability_ceac, color = "Divalproex ER"), linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 150000, aes(color = "WTP Threshold"), linetype = "dashed") +
  annotate("text", x = 230000, y = 1.1, label = "WTP Threshold of $150,000", angle = 0, size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size=8),legend.text = element_text(size=8), plot.caption = element_text(hjust = 0), plot.caption.position = "plot" ) +
  scale_color_manual(values = c("lightblue","navy")) +
  theme(legend.key.width = unit(1.5, 'cm'))
ceac_plot

#One-way sensitivity analysis#####################
##Defining OSA Parameters####
osa_param <- define_dsa(
  c_dival, 5.2, 7.0,
  c_eptin, 449.28, 607.84,

  p_dc_dival, 0.00238, 0.00322,
  p_dc_eptin, 0.00774, 0.0105,

  u_on_dival, 0.04779, 0.06465,
  u_on_eptin, 0.05073, 0.06863,
) #all parameter is modified +/- 15%

##run OSA####
osa_eptin <- run_dsa(
  model = eptin_mod,
  dsa = osa_param)


osa_erenu <- run_dsa(
  model = erenu_mod,
  dsa = osa_param)


osa_dival <- run_dsa(
  model = dival_mod,
  dsa = osa_param)

#create a table with parameters and low/high range
osa_param_list <- data.frame(
  parameters = c("Cost of divalproex", "Cost of eptinezumab", "Discontinuation rate of divalproex", "Discontinuation rate of eptinezumab", "Utility on divalproex", "Utility on eptinezumab", "Probability of achieving 50% MMD reduction on first 6 months of divaproex", "Probability of achieving 50% MMD reduction on first 6 months of eptinezumab"),
  value = c(c_dival, c_eptin, p_dc_dival, p_dc_eptin, u_on_dival, u_on_eptin, p_clinical_response_dival, p_clinical_response_eptin)
)

osa_param_list$low <- osa_param_list$value * 0.85
osa_param_list$high <- osa_param_list$value * 1.15

##Calculating incremental net monetary benefit####
WTP=150000
osa_nmb <- function(i){
  nmb <- (150000 * (osa_eptin$dsa$utility[i]-osa_dival$dsa$utility[i])) - (osa_eptin$dsa$cost_total[i]-osa_dival$dsa$cost_total[i])
  return(nmb)
}

low_nmb <- matrix(NA, nrow = 6, ncol = 1)
high_nmb <- matrix(NA, nrow = 6, ncol = 1)

for(i in seq(1,11,by = 2)){
  low_nmb[i] <- osa_nmb(i)
}
low_nmb <- as.data.frame(low_nmb)
low_nmb <- low_nmb[seq(1, nrow(low_nmb), by = 2), ]
low_nmb <- c(low_nmb, 5850, -1386) #

for(i in seq(2,12,by = 2)){
  high_nmb[i] <- osa_nmb(i)
}
high_nmb <- as.data.frame(high_nmb)
high_nmb <- high_nmb[seq(2, nrow(high_nmb), by = 2), ]
high_nmb <- c(high_nmb, -599, 6637)

osa_result <- cbind.data.frame(osa_param_list, low_nmb, high_nmb)
osa_result <- as.data.frame(osa_result)

osa_result

base_case_nmb <- (WTP*deterministic_result$incremental_qaly[1]) - deterministic_result$incremental_cost[1]

#calculate spread
osa_result$spread <- abs(osa_result$high_nmb - osa_result$low_nmb)
osa_result <- osa_result[order(osa_result$spread, decreasing = TRUE),]
rownames(osa_result) <- NULL
osa_result$parameters <- factor(osa_result$parameters, levels = rev(osa_result$parameters))

width <- 0.5

osa_plot_low <- data.frame(
  parameters = osa_result$parameters,
  ymin = pmin(base_case_nmb, osa_result$low_nmb),
  ymax = pmax(base_case_nmb, osa_result$low_nmb),
  xmin = as.numeric(osa_result$parameters)-width/2,
  xmax = as.numeric(osa_result$parameters)+width/2
)

order_parameters <- osa_result %>% arrange(spread) %>%
  mutate(parameters=factor(x=parameters, levels=parameters)) %>%
  select(parameters) %>% unlist() %>% levels()

osa_plot_high <- data.frame(
  parameters = osa_result$parameters,
  ymin = pmin(base_case_nmb, osa_result$high_nmb),
  ymax = pmax(base_case_nmb, osa_result$high_nmb),
  xmin = as.numeric(osa_result$parameters)-width/2,
  xmax = as.numeric(osa_result$parameters)+width/2
)

labels <- c("Low Input", "High Input")
legend_colors <- setNames(c("navy","lightblue"),labels)
osa_plot <- ggplot() +
  labs(caption = "(b) Tornado diagram of the net monetary benefit changes between eptinezumab and divaproex") +
  geom_rect(data = osa_plot_low, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "Low Input")) +
  geom_rect(data = osa_plot_high, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "High Input")) +
  geom_hline(yintercept = base_case_nmb) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_x_continuous(breaks = c(1:length(order_parameters)),
                     labels = order_parameters) +
  ylab("Net Monetary Benefit ($)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right", legend.text =element_text(size=8, color = "black"),axis.title = element_text(size=8, color = "black"), plot.caption = element_text(hjust = 0), plot.caption.position = "plot") +
  scale_fill_manual("", values = legend_colors) +
  coord_flip()
osa_plot

#Figure 1##########################
grid.arrange(ceac_plot,osa_plot,nrow=2)

