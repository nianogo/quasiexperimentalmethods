#------------------------------------------------------------------------------#
#-----Quasi-Experimental Methods Based on the Timing of Natural Experiments-----
#-------------------------A Guide for Epidemiologists"--------------------------
#------------------------------------------------------------------------------#
#Authors:
#Roch Nianogo (niaroch@ucla.edu) 
#Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu)
#O'Neill Stephen (stephen.oneill@lshtm.ac.uk)

#This script provides sample codes for replicating some of the results in 
#manuscripts and for applying the methods to different settings

#Outline
#I:   DATA GENERATING PROCESS
#II:  DATA EXPLORATION
#III: ANALYSIS
# 1) Pre-post analysis
# 2) Controlled pre-post analysis
# 3) Interrupted Time series
# 4a) 2x2 Difference-in-difference
# 4b) Difference-in-difference/Controlled ITS (CITS)
# 5) Traditional Synthetic Control Method (TSCM) (Abadie et al)
# 6) Generalized Synthetic Control Method (GSCM) (Xu)

###############################################################################-
#################################0: SET-UP######################################
###############################################################################-


if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
} # a nice package to load several packages simultaneously


p_load("tidyverse","magrittr","broom",        #Manipulate data
       "Synth",                               #Traditional SCM
       "gsynth",                              #Generalized SCM
       "panelView",                           #Exploring the data
       "lme4", "estimatr")                    #multi-level model


###############################################################################-
###########################I: DATA GENERATING PROCESS###########################
###############################################################################-

# This is a simulated data where a fictitious policy ban was implemented
# The policy was enacted in one state:California (or state number 5)
# The policy was enacted in year 40
# The unit of analysis is the state
# y  = is the outcome
# xi = unit-varying variable (varies across states but constant across years)
# xt = time-varying variable (variess across years but constant across states)
# xit = unit-time-variable (i.e. varies by years and statess)


  set.seed(123)
  
      # PTrend=T : Whether the parallel trend assumption is met
      # Cshock=T : Whether the common shock assumption is met
      # Valid_controls=T, Whether invalid controls have been introduced
      # LinearTrends=T, Whether the outcome trends are linear

      #-----creating IDs-----------
      n_state     =50
      n_treated   =n_treated
      
      year_start  =1
      year_end    =50
      year_policy <<-40
      n_year      =length(seq(year_start, year_end, 1))
      states      = state.name[1:n_state]
      years       = 1:n_year
      state_nums  = 1:n_state
      
      state       = rep(states, each = n_year)
      state_num   = as.numeric(as.factor(state))
      
      year        = rep(years, times = n_state)
      
      violated_crtl <-  ifelse(state_num %in% c(10, 27:50), 1, 0)

      #-----Parameters-----------
      b0=1    #intercept
      b1=1    #ui1
      b2=20   #lambda1t
      b3=1     #linear trend
      b4=0.15  #Parallel trend
      b5=100   #Common shock
      b6=0.001 #invalid controls
      b7=1    #xi:unit-varying
      b8=1    #xt:time-varying
      b9=1    #xit:unit-time-varying
      tau=100 #ATT
      
      treated <- ifelse(state %in% treated_states, 1, 0)
      
      post <-  ifelse(year >= year_policy, 1, 0)
      treatedpost <- treated*post
      
      convexhullupper = as.numeric(state_num==3)
      convexhulllower = as.numeric(state_num==4)
  
      #unit-varying covariates
      
      xi  = 2*treated + rep(rnorm(n_state, mean=5, sd=1), each=n_year)  +
        -5*convexhulllower  + 5*convexhullupper
      ui1 = 2*treated + rep(rnorm(n_state, mean=5, sd=1), each=n_year)  +
        -5*convexhulllower  + 5*convexhullupper
      ui2 = 2*treated + rep(rnorm(n_state, mean=5, sd=1), each=n_year)  +
        -5*convexhulllower  + 5*convexhullupper
      ui3 = 1         + rep(rnorm(n_state, mean=5, sd=1), each=n_year)  +
        -5*convexhulllower  + 5*convexhullupper 
  

       #time-varying covariates
      xt = rep(rnorm(n_year, mean = 1 + 0.5*years, sd=1), times=n_state)     
      
      lambda1t = year
      lambda1nonlinear =        (0*I(year %in%  0:9) +
                                 5*I(year %in% 10:29)  +
                                 -10*I(year %in%  30:39)  +
                                 50*I(year %in% 40:44) +
                                 100*I(year %in% 45:50))

      lambda2t = (1-treated)*(3 + 10*year*I(year<=40))
      lambda3t = year^3
      
      #unit-time varying covariates
      xit = rnorm(n_year*n_state, mean=0, sd=1) 
      
      uit = rnorm(n_year*n_state, mean=0, sd=1)
      
      #epsilon
      epsilon = rnorm(n_state*n_year, mean = 0, sd=1)
      
      #Potential outcomes
      y_cf0 = b0 + 
        b1*ui1 + 
        b2*lambda1t + 
        b3*lambda1nonlinear*(1-LinearTrends) + 
        b4*ui2*lambda2t*(1-PTrend) + 
        b5*treatedpost*(1-Cshock)  + 
        b6*ui3*lambda3t*violated_crtl*(1-Valid_controls) + 
        b7*xi + b8*xt + b9*xit + epsilon
      
      y_cf1 = y_cf0 + tau
      
      y = treatedpost*y_cf1 + (1-treatedpost)*y_cf0
      
      
      mydata <- tibble(state, violated_crtl, state_num, year, treated, post, treatedpost, y, y_cf0,y_cf1,
                   xit,xt, xi) 


###############################################################################-
#############################II: DATA EXPLORATION###############################
###############################################################################-

#-------------------------------------------------#
#--------------Panel view of the data--------------
#--------------------------------------------------#


mydata <- mydata %>% 
  as.data.frame() # need to convert to a data.frame for some function to work

p_load("panelView")
panelView(y ~ treatedpost, data = mydata, index = c("state","year"), pre.post = TRUE) 


###############################################################################-
###################################III: ANALYSIS###############################
###############################################################################-

#---------------------------#
#----1)Pre-post analysis----
#---------------------------#

##Step1:Create the data---
dt <- mydata %>%
  select(state,year, post, xi, xt, xit, y) %>% 
  filter(state=="California",
         year %in% c(39, 45)) 


##Step2:Preview the data---
head(dt,4)


##Step3:Plot the data---
dt %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Pre-post"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  scale_x_continuous(breaks= scales::pretty_breaks()) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 


##Step4:Analyze the data---
fit <- lm(y ~ post + xit, data=dt)
dta <- tidy(fit, conf.int = T)

pre_post <- round(data.frame(ATT  = dta[dta$term=="post", "estimate"],
                             se       = dta[dta$term=="post", "std.error"],
                             low_ci   = dta[dta$term=="post", "conf.low"],
                             high_ci  = dta[dta$term=="post", "conf.high"]),2)

pre_post
#Cannot estimate SE as expected
#-----------------------------------#
#----2)Controlled Pre-post analysis----
#-----------------------------------#
##Step1:Create the data---
dt <- mydata %>%
  select(state,year, post, treated, xi, xt, xit, y) %>% 
  filter(state=="California" | state=="Georgia",
         year %in% c(39, 45)) 


##Step2:Preview the data---
head(dt,4)


##Step3:Plot the data---
dt %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Controlled pre-post"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  scale_x_continuous(breaks= scales::pretty_breaks()) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 


##Step4:Analyze the data---
fit <- lm(y ~ treated*post + xit, data=dt)
dta <- tidy(fit, conf.int = T)
dta
crtl_pre_post <- round(data.frame(ATT = dta[dta$term=="treated:post", "estimate"],
                                  se       = dta[dta$term=="treated:post", "std.error"],
                                  low_ci   = dta[dta$term=="treated:post", "conf.low"],
                                  high_ci  = dta[dta$term=="treated:post", "conf.high"]),2)

crtl_pre_post
#Cannot estimate SE and Point estimate if adjusting for another unit-time varying
#variable

#-----------------------------------------#
#---3)Interrupted Time Series analysis----
#-----------------------------------------#
##Step1:Create the data---
dt <- mydata %>%
  select(state,year, post, treated, xi, xt, xit, y) %>% 
  filter(state=="California") 


##Step2:Preview the data---
head(dt)

##Step3:Plot the data---
dt %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Interrupted Time Series"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  scale_x_continuous(breaks= scales::pretty_breaks()) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 


##Step4:Analyze the data---
#Best to recenter year to the time of the policy for better interpretability
fit <- lm(y ~ post*I(year-40) + xit, data=dt)
dta <- tidy(fit, conf.int = T)
dta
its <- round(data.frame(ATT      = dta[dta$term=="post", "estimate"],
                        se       = dta[dta$term=="post", "std.error"],
                        low_ci   = dta[dta$term=="post", "conf.low"],
                        high_ci  = dta[dta$term=="post", "conf.high"]),2)

its

#---------------------------------------#
#-----4a)2x2 Difference in difference-----
#---------------------------------------#
##Step1:Create the data---
dt <- mydata %>%
  select(state,year, post, treated, xi, xt, xit, y) %>% 
  filter(state=="California" | state=="Georgia") 


##Step2:Preview the data---
head(dt)

##Step3:Plot the data---
dt %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("2x2 DID"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  scale_x_continuous(breaks= scales::pretty_breaks()) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 


##Step4:Analyze the data---
p_load("estimatr")
dta <- lm_robust(y ~ treated*post + xit, 
                 data = dt, 
                 se_type = "stata")
dta
twobytwo_did <- round(data.frame(ATT     = dta$coefficients["treated:post"], 
                                 se      = dta$std.error["treated:post"],
                                 low_ci  = dta$conf.low["treated:post"],
                                 high_ci = dta$conf.hig["treated:post"]),2)
twobytwo_did

#-----------------------------------------#
#------4b)CITS/Difference in difference-----
#-----------------------------------------#
##Step1:Create the data---
#Step 1,2 and 3 will be the same for the CITS/DID, SCM and GSCM
dt <- mydata %>%
  select(state,state_num, year, post, treated, treatedpost, xi, xt, xit, y) %>% 
  as.data.frame() #this is necessary for the SCM using the Synth package


##Step2:Preview a sample of the the data---
head(dt)

##Step3:Plot the data---
#For all the data
dt %>% 
  ggplot(aes(x=year, y=y, group=state, 
             size=factor(treated), 
             color=factor(treated))) + 
  geom_line() +
  scale_size_discrete(range = c(0.1, 1),labels=c("Control", "Treated")) +
  labs(title = "Outcome by states",
       x = "Year", 
       y = "Outcome",
       color = "Treatment") +
  scale_color_discrete(labels=c("Control", "Treated")) +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(size = "none")


#On average
dt %>% 
  group_by(year, treated) %>% 
  summarise(y=mean(y),.groups="keep") %>% 
  ggplot(aes(x=year, y=y, group=treated, color = factor(treated))) + 
  labs(title = "Average outcome by group",
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  scale_color_discrete(labels=c("Controls", "Treated")) +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##Step4:Analyze the data---

#Checking for the parallel trend asssumption

#Check the parallel trends assumptions
#To do so
#restrict the data to the period before the policy
#Run a linear regression of the outcome on time-varying
#covariates and on an interaction term between treated indicator and year

p_load("estimatr") #for the lm_robust() function
pretrend_data <- dt %>% 
  filter(post == 0)

#the lm_robust procedure is best because SE are correctly estimated
res_pretrend <- lm_robust(y ~ treated*year + xit, 
                          data = pretrend_data,
                          fixed_effects=state,
                          clusters = state, se_type = "stata")
summary(res_pretrend)
#Pay particular attention to the value and p-value value of
# treated:year. If the coefficient of the interaction is close to 0 or 
#the p-value large, then
# can say that the unit-specific confounders do not have time-varying effects
#in other words, the parallel trend assumption might be okay
#recall that this is only true if we assume that the trend is linear

# If okay, now we can implement the DID analysis method
# Fit a linear model with the treated, post and their interaction
# Don't forget to add any other time-varying covariates
#Method 1: Lm_robust()
p_load("estimatr")
dta <- lm_robust(y ~ treated*post + xit, 
                 data = dt, 
                 fixed_effects=state,
                 clusters = state, 
                 se_type = "stata")
dta
did <- round(data.frame(ATT     = dta$coefficients["treated:post"], 
                        se      = dta$std.error["treated:post"],
                        low_ci  = dta$conf.low["treated:post"],
                        high_ci = dta$conf.hig["treated:post"]),2)
did

#Method 2: lmer(): Multilevel
p_load("lmerTest")
fit <- lmerTest::lmer(y ~ treated*post + xit + 
                        (1| state) + 
                        (1 | year), 
                      data = dt, 
                      REML = T)
summary(fit)

#-----------------------------------------------#
#------5) Synthetic Control Method (Abadie)-----
#-----------------------------------------------#
##Step4:Analyze the data---
# create a list of treated states
list_int <- dt %>% 
  filter(treated==1) %>% 
  select(state_num) %>% 
  distinct() %>% 
  as.matrix() %>% 
  as.vector()

# create a list of control states
list_control <- dt %>% 
  filter(treated==0) %>% 
  select(state_num) %>% 
  distinct() %>% 
  as.matrix() %>% 
  as.vector()

# create some handy object to be used in the Synth package
year_start  = 1
year_end    = 50
year_policy = 40
year_policy_prior = year_policy - 1

# Prep the data
dataprep.out <-
  dataprep(foo = dt,
           predictors    = "xit", #include unit-time-varying covariates only
           special.predictors = list(list("y", year_start:year_policy_prior, "mean")),
           predictors.op = "mean",
           dependent     = "y",
           unit.variable = "state_num",
           time.variable = "year",
           unit.names.variable = "state",
           treatment.identifier  = 5, #cannot be several states
           controls.identifier   = list_control,
           time.predictors.prior = c(year_start:year_policy_prior),
           time.optimize.ssr     = c(year_start:year_policy_prior),
           time.plot             = c(year_start:year_end))


# Run the Synth command
synth.out <- synth(dataprep.out)

# Get result tables
print(synth.tables   <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res    = synth.out)
)

# see the weights
round(synth.out$solution.w,2)


# Obtain the estimate using a DID
##--Estimate the average (weighted) outcome in the synthetic control prior to the policy

y_ctrl_pre <- dataprep.out$Y0plot%*%synth.out$solution.w %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  rename(control=w.weight) %>% 
  filter(year %in% c(year_start:year_policy_prior)) %>% 
  group_by() %>% 
  summarise(mean=mean(control)) %>% 
  as.matrix() %>% 
  as.vector()

##--Estimate the average (weighted) outcome in the synthetic control after the policy
y_ctrl_post <- dataprep.out$Y0plot%*%synth.out$solution.w %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  rename(control=w.weight) %>% 
  filter(year %in% c(year_policy:year_end)) %>% 
  group_by() %>% 
  summarise(mean=mean(control)) %>% 
  as.matrix() %>% 
  as.vector()

##--Estimate the average outcome in the treated state prior to the policy
y_treat_pre <- dataprep.out$Y1plot %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  filter(year %in% c(year_start:year_policy_prior)) %>% 
  select(-c(year)) %>% 
  as.matrix() %>% 
  as.vector() %>% mean()

##--Estimate the average outcome in the treated state after the policy
y_treat_post <- dataprep.out$Y1plot %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  filter(year %in% c(year_policy:year_end)) %>% 
  select(-c(year)) %>% 
  as.matrix() %>% 
  as.vector() %>% mean()

##--Estimate a simple difference-in-difference
diff_post = y_treat_post - y_ctrl_post
diff_pre = y_treat_pre - y_ctrl_pre

DID = diff_post - diff_pre
DID

# Plot the results
##--Plot treated and synthetic state
##--Plot treated and synthetic state
path.plot(synth.res    = synth.out,
          dataprep.res = dataprep.out,
          Ylab         = c("Y"),
          Xlab         = c("Year"),
          Legend.position = c("topleft"))

abline(v   = 40,
       lty = 2)


##--Plot gap between treated and synthetic state
gaps.plot(synth.res    = synth.out,
          dataprep.res = dataprep.out,
          Ylab         = c("Gap"),
          Xlab         = c("Year"),
          Main         = ""
)
abline(v   = 40,
       lty = 2)



#Generate placebo/falsification test
p_load("SCtools")
# run the synth command to create the synthetic control
synth.out <- synth(dataprep.out, Sigf.ipop=2)

## run the generate.placebos command to reassign treatment status
## to each unit listed as control, one at a time, and generate their
## synthetic versions. Sigf.ipop = 2 for faster computing time. 
## Increase to the default of 5 for better estimates. 
tdf <- generate.placebos(dataprep.out,synth.out, Sigf.ipop = 2, strategy='multiprocess')

## Plot the gaps in outcome values over time of each unit --
## treated and placebos -- to their synthetic controls
p <- plot_placebos(tdf,discard.extreme=T, mspe.limit=10, xlab='Year')
p
#Treated/California should be the bolded line. 
#There must be a bug in the package SCtools


#-------------------------------------------------------#
#------6) Generalized Synthetic Control Method (Xu)-----
#-------------------------------------------------------#

##Step4:Analyze the data---
y <- gsynth(y ~ treatedpost + xit, data = dt,  EM = F, index = c("state","year"), 
            inference = "parametric", se = TRUE,
            nboots = 200,  r = c(0, 5), CV = TRUE, force = "two-way", parallel = FALSE)
y1 <- round(data.frame(y$est.avg),2)
y1
p0 <- plot(y)
p0

# END
