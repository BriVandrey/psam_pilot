# function to plot and analyse PSAM data
# developed by Brianna Vandrey, last updated 05-2022

# import modules
library(tidyverse)
library(rstatix)
library(gridExtra)
library(ggpubr)

# read file
path <- "~/Desktop/psam_pilot/plots" # specify location to save plots 
df <- read.csv(file="~/Desktop/psam_pilot/psam_df.csv", head=TRUE, sep=",") #csv to read in 

#add absolute values for ipsp 
df$ipsp_abs <- abs(df$ipsp)
df$ipsp_max_abs <- abs(df$ipsp_max)

#subset based on patterns PSAM expression
psam_df <- subset(df, psam_expression_cell=="Y") #PSAM+ neurons
control_df <- subset(df, psam_expression_cell=="N") # PSAM- neurons

psam_pc_df <- subset(psam_df, cell_type != "l2_in"); psam_pc_df <- subset(psam_pc_df, cell_type != "l1_in") # PSAM+ PCs
psam_in_df <- subset(psam_df, cell_type == "l2_in" | cell_type == "l1_in") # PSAM+ interneurons

control_pc_df <- subset(control_df, cell_type != "l2_in"); control_pc_df <- subset(control_pc_df, cell_type != "l1_in") # PSAM- PCs
control_in_df <- subset(control_df, cell_type != "l2_pc"); control_in_df <- subset(control_in_df, cell_type != "l2_sc") # PSAM- INs

chr2_psam <- subset(control_df, psam_expression_slice == "Y") #data with ChR2 excitation & PSAM expression in the slice

# --------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------Plot subthreshold/spiking data ----------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------

# plotting functions - adapted from lineplot code originally written by Sau Tsoi, 2021

# --------------------------------------------------------------------------------------------------------------------------------------------------------
#function for producing lineplot with individual observations for x 3 conditions
# --------------------------------------------------------------------------------------------------------------------------------------------------------
lp_3_conds <- function(df, col, drug_dose, xlabel1, xlabel2, ylabel, filename){
  
  df <- subset(df, condition==drug_dose) 
  
  # handle col variable
  col <- enquo(col)
  
  #convert condition into a factor with the levels in the correct order
  df$drug <- as.character(df$drug)
  df$drug <- factor(df$drug, levels=unique(df$drug))
  
  #summarise data 
  sum <-group_by(df, drug, id) %>% 
    summarise(mean = mean(!!col)) 
  sum_mean<-group_by(sum, drug) %>% 
    get_summary_stats(type = "mean_se")
  
  #lineplot 
  lp <- ggplot(sum_mean, aes(x=drug, y=mean))+
    theme_classic()+
    geom_line(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean, group = id), size =0.6, colour="grey74")+ #individual traces
    geom_point(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean, group = id), size = 2, shape=111, colour="grey74") + #individual points
    geom_point(size = 1) + #mean points
    geom_line(aes(group=1), size = 0.6) + #mean trace
    geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width = 0.2, size = 0.6) + #error bars 
    scale_x_discrete(labels=c("Baseline", xlabel1, xlabel2)) +
    xlab("") + ylab(ylabel)
  
  lp + theme(axis.text.x = element_text(size = 8, angle=-50, color='black'), 
                  axis.text.y = element_text(size = 8, color='black'),
                  axis.title.y = element_text(size=12, margin= margin(t=0,r=5,b=0,l=0))) 
  
  ggsave(path=path, file=paste(filename, "lp.eps", sep="_"), width =50, height =60, units="mm") #save plot. 
}

# plot data for PSAM-expressing principal neurons with 10 nM uPSEM-792 application
lp_3_conds(psam_pc_df, vm, 10, "uPSEM-792 (10 nM)", "Washout", "RMP (mV)", "psam_pc_rmp_10nM")
lp_3_conds(psam_pc_df, ir_neg, 10, "uPSEM-792 (10 nM)", "Washout", "IR-", "psam_pc_irneg_10nM")
lp_3_conds(psam_pc_df, sag, 10, "uPSEM-792 (10 nM)", "Washout", "Sag", "psam_pc_sag_10nM")
lp_3_conds(psam_pc_df, rheobase, 10, "uPSEM-792 (10 nM)", "Washout", "Rheobase (pA)", "psam_pc_rheobase_10nM")
lp_3_conds(psam_pc_df, ap_thresh, 10, "uPSEM-792 (10 nM)", "Washout", "AP Thresh (mV)", "psam_pc_apthresh_10nM")

# plot data for PSAM-expressing interneurons with 10 nM uPSEM-792 application
lp_3_conds(psam_in_df, vm, 10, "uPSEM-792 (10 nM)", "Washout", "RMP (mV)", "psam_in_rmp_10nM")
lp_3_conds(psam_in_df, ir_neg, 10, "uPSEM-792 (10 nM)", "Washout", "IR-", "psam_in_irneg_10nM")
lp_3_conds(psam_in_df, sag, 10, "uPSEM-792 (10 nM)", "Washout", "Sag", "psam_in_sag_10nM")
lp_3_conds(psam_in_df, rheobase, 10, "uPSEM-792 (10 nM)", "Washout", "Rheobase (pA)", "psam_in_rheobase_10nM")
lp_3_conds(psam_in_df, ap_thresh, 10, "uPSEM-792 (10 nM)", "Washout", "AP Thresh (mV)", "psam_in_apthresh_10nM")

# plot control data - 10 nM
lp_3_conds(control_pc_df, vm, 10, "uPSEM-792 (10 nM)", "Washout", "RMP (mV)", "control_pc_rmp_10nM")
lp_3_conds(control_pc_df, ir_neg, 10, "uPSEM-792 (10 nM)", "Washout", "IR-", "control_pc_irneg_10nM")
lp_3_conds(control_pc_df, sag, 10, "uPSEM-792 (10 nM)", "Washout", "Sag", "control_pc_sag_10nM")
lp_3_conds(control_pc_df, rheobase, 10, "uPSEM-792 (10 nM)", "Washout", "Rheobase (pA)", "control_pc_rheobase_10nM")
lp_3_conds(control_pc_df, ap_thresh, 10, "uPSEM-792 (10 nM)", "Washout", "AP Thresh (mV)", "control_pc_apthresh_10nM")
lp_3_conds(control_in_df, vm, 10, "uPSEM-792 (10 nM)", "Washout", "RMP (mV)", "control_in_rmp_10nM")
lp_3_conds(control_in_df, ir_neg, 10, "uPSEM-792 (10 nM)", "Washout", "IR-", "control_in_irneg_10nM")
lp_3_conds(control_in_df, sag, 10, "uPSEM-792 (10 nM)", "Washout", "Sag", "control_in_sag_10nM")
lp_3_conds(control_in_df, rheobase, 10, "uPSEM-792 (10 nM)", "Washout", "Rheobase (pA)", "control_in_rheobase_10nM")
lp_3_conds(control_in_df, ap_thresh, 10, "uPSEM-792 (10 nM)", "Washout", "AP Thresh (mV)", "control_in_apthresh_10nM")

# plot data for PSAM-expressing principal neurons with 20 nM uPSEM-792 application
lp_3_conds(psam_pc_df, vm, 20, "uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "psam_pc_rmp_20nM")
lp_3_conds(psam_pc_df, ir_neg, 20, "uPSEM-792 (20 nM)", "Washout", "IR-", "psam_pc_irneg_20nM")
lp_3_conds(psam_pc_df, sag, 20, "uPSEM-792 (20 nM)", "Washout", "Sag", "psam_pc_sag_20nM")

# plot data for PSAM-expressing interneurons with 20 nM uPSEM-792 application
lp_3_conds(psam_in_df, vm, 20, "uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "psam_in_rmp_20nM")
lp_3_conds(psam_in_df, ir_neg, 20, "uPSEM-792 (20 nM)", "Washout", "IR-", "psam_in_irneg_20nM")
lp_3_conds(psam_in_df, sag, 20, "uPSEM-792 (20 nM)", "Washout", "Sag", "psam_in_sag_20nM")
lp_3_conds(psam_in_df, rheobase, 20, "uPSEM-792 (20 nM)", "Washout", "Rheobase (pA)", "psam_in_rheobase_20nM")
lp_3_conds(psam_in_df, ap_thresh, 20, "uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)", "psam_in_apthresh_20nM")

# plot control data - 20 nM
lp_3_conds(control_pc_df, vm, 20, "uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "control_pc_rmp_20nM")
lp_3_conds(control_pc_df, ir_neg, 20, "uPSEM-792 (20 nM)", "Washout", "IR-", "control_pc_irneg_20nM")
lp_3_conds(control_pc_df, sag, 20, "uPSEM-792 (20 nM)", "Washout", "Sag", "control_pc_sag_20nM")
lp_3_conds(control_pc_df, rheobase, 20, "uPSEM-792 (20 nM)", "Washout", "Rheobase (pA)", "control_pc_rheobase_20nM")
lp_3_conds(control_pc_df, ap_thresh, 20, "uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)", "control_pc_apthresh_20nM")


# --------------------------------------------------------------------------------------------------------------------------------------------------------
#function for producing lineplot with individual observations for x 4 conditions
# --------------------------------------------------------------------------------------------------------------------------------------------------------
lp_4_conds <- function(df, col, drug_dose, xlabel1, xlabel2, xlabel3, ylabel, filename){
  
  df <- subset(df, condition==drug_dose) 
  
  # handle col variable
  col <- enquo(col)
  
  #convert condition into a factor with the levels in the correct order
  df$drug <- as.character(df$drug)
  df$drug <- factor(df$drug, levels=unique(df$drug))
  
  #summarise data 
  sum <-group_by(df, drug, id) %>% 
    summarise(mean = mean(!!col)) 
  sum_mean<-group_by(sum, drug) %>% 
    get_summary_stats(type = "mean_se")
  
  #lineplot 
  lp <- ggplot(sum_mean, aes(x=drug, y=mean))+
    theme_classic()+
    geom_line(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean, group = id), size =0.6, colour="grey74")+ #individual traces
    geom_point(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean, group = id), size = 2, shape=111, colour="grey74") + #individual points
    geom_point(size = 1) + #mean points
    geom_line(aes(group=1), size = 0.6) + #mean trace
    geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width = 0.2, size = 0.6) + #error bars 
    scale_x_discrete(labels=c("Baseline", xlabel1, xlabel2, xlabel3)) +
    xlab("") + ylab(ylabel)
  
  lp + theme(axis.text.x = element_text(size = 8, angle=-50, color='black'), 
             axis.text.y = element_text(size = 8, color='black'),
             axis.title.y = element_text(size=12, margin= margin(t=0,r=5,b=0,l=0))) 
  
  ggsave(path=path, file=paste(filename, "lp.eps", sep="_"), width =60, height =60, units="mm") #save plot
}

# plot data for PSAM-expressing principal neurons with 10-20 nM uPSEM-792 application
lp_4_conds(psam_pc_df, vm, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "psam_pc_rmp_10-20nM")
lp_4_conds(psam_pc_df, ir_neg, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "IR-", "psam_pc_irneg_10-20nM")
lp_4_conds(psam_pc_df, sag, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "Sag", "psam_pc_sag_10-20nM")
lp_4_conds(psam_pc_df, rheobase, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "Rheobase", "psam_pc_rheobase_10-20nM")
lp_4_conds(psam_pc_df, ap_thresh, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)", "psam_pc_apthresh_10-20nM")

# plot data for PSAM-expressing interneurons with 10-20 nM uPSEM-792 application
lp_4_conds(psam_in_df, vm, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "psam_in_rmp_10-20nM")
lp_4_conds(psam_in_df, ir_neg, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "IR-", "psam_in_irneg_10-20nM")
lp_4_conds(psam_in_df, sag, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "Sag", "psam_in_sag_10-20nM")
lp_4_conds(psam_in_df, rheobase, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "Rheobase", "psam_in_rheobase_10-20nM")
lp_4_conds(psam_in_df, ap_thresh, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)", "psam_in_apthresh_10-20nM")

# plot control data - 10 to 20 nM
lp_4_conds(control_pc_df, vm, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "control_pc_rmp_10-20nM")
lp_4_conds(control_pc_df, ir_neg, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "IR-", "control_pc_irneg_10-20nM")
lp_4_conds(control_pc_df, sag, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "Sag", "control_pc_sag_10-20nM")
lp_4_conds(control_pc_df, rheobase, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "Rheobase", "control_pc_rheobase_10-20nM")
lp_4_conds(control_pc_df, ap_thresh, 1020, "uPSEM-792 (10 nM)","uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)", "control_pc_apthresh_10-20nM")


# --------------------------------------------------------------------------------------------------------------------------------------------------------
#function for producing lineplot with individual observations for x 5 conditions
# --------------------------------------------------------------------------------------------------------------------------------------------------------
lp_5_conds <- function(df, col, drug_dose, xlabel1, xlabel2, xlabel3, xlabel4, ylabel, filename){
  
  df <- subset(df, condition==drug_dose) 
  
  # handle col variable
  col <- enquo(col)
  
  #convert condition into a factor with the levels in the correct order
  df$drug <- as.character(df$drug)
  df$drug <- factor(df$drug, levels=unique(df$drug))
  
  #summarise data 
  sum <-group_by(df, drug, id) %>% 
    summarise(mean = mean(!!col)) 
  sum_mean<-group_by(sum, drug) %>% 
    get_summary_stats(type = "mean_se")
  
  #lineplot 
  lp <- ggplot(sum_mean, aes(x=drug, y=mean))+
    theme_classic()+
    geom_line(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean, group = id), size =0.6, colour="grey74")+ #individual traces
    geom_point(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean, group = id), size = 2, shape=111, colour="grey74") + #individual points
    geom_point(size = 1) + #mean points
    geom_line(aes(group=1), size = 0.6) + #mean trace
    geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width = 0.2, size = 0.6) + #error bars 
    scale_x_discrete(labels=c("Baseline", xlabel1, xlabel2, xlabel3, xlabel4)) +
    xlab("") + ylab(ylabel)
  
  lp + theme(axis.text.x = element_text(size = 8, angle=-50, color='black'), 
             axis.text.y = element_text(size = 8, color='black'),
             axis.title.y = element_text(size=12, margin= margin(t=0,r=5,b=0,l=0))) 
  
  ggsave(path=path, file=paste(filename, "lp.eps", sep="_"), width=70, height =60, units="mm") #save plot
}

# plot control data - 10 to 20 nM  
lp_5_conds(control_pc_df, vm, 1, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "RMP (mV)", "control_pc_rmp_1-20nM")
lp_5_conds(control_pc_df, ir_neg, 1, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "IR-", "control_pc_irneg_1-20nM")
lp_5_conds(control_pc_df, sag, 1, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "Sag", "control_pc_sag_1-20nM")
lp_5_conds(control_pc_df, rheobase, 1, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "Rheobase", "control_pc_rheobase_1-20nM")
lp_5_conds(control_pc_df, ap_thresh, 1, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)", "control_pc_apthresh_1-20nM")



# --------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------ Stats subthreshold/spiking data -----------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------
#function to perform t-test
run_t_test <- function (df, col, filter){
  
  df_filt <- subset(df, drug != 'washout') #remove washout data 
  df_filt <- subset(df_filt, condition==filter) #subset drug concentration
  res <- t.test(eval(substitute(col~drug)), data=df_filt) #t-test
  
  return(res)
}


#function to perform ANOVA analysis w pairwise comparisons
run_aov <- function (df, col, filter){
  
  #filter based on data
  df_filt <- subset(df, condition==filter)
  
  #handle col data
  col1 <- enquo(col)
  
  #summary statistics
  sum_stats <- df_filt %>%
    group_by(drug) %>%
    summarise(count = n(),
            mean = mean(!!col1, na.rm=TRUE), 
            sd = sd(!!col1, na.rm=TRUE))
  
  #anova
  aov <- aov(eval(substitute(col~drug)), data=df_filt) #need to use eval/substitute to read col into Base R
  aov_res <-summary(aov)
  pw <- TukeyHSD(aov) #pairwise comparisons
  
  #output results
  output <- c(sum_stats, aov_res, pw) 
  return(output)
  
}

# stats for PSAM-expressing neurons
psam_pc_vm_10 <- run_t_test(psam_pc_df, vm, 10) #rmp
psam_pc_vm_20 <- run_t_test(psam_pc_df, vm, 20)
psam_in_vm_10 <- run_t_test(psam_in_df, vm, 10)

psam_pc_ir_10 <- run_t_test(psam_pc_df, ir_neg, 10) #ir
psam_pc_ir_20 <- run_t_test(psam_pc_df, ir_neg, 20)
psam_in_ir_10 <- run_t_test(psam_in_df, ir_neg, 10) #** p = 0.036

psam_pc_sag_10 <- run_t_test(psam_pc_df, sag, 10) #sag
psam_pc_sag_20 <- run_t_test(psam_pc_df, sag, 20)
psam_in_sag_10 <- run_t_test(psam_in_df, sag, 10) #** p = 0.019

psam_pc_rheobase_10 <- run_t_test(psam_pc_df, rheobase, 10) #rheobase
psam_in_rheobase_10 <- run_t_test(psam_in_df, rheobase, 10)
psam_in_rheobase_20 <- run_t_test(psam_in_df, rheobase, 20)

psam_pc_ap_10 <- run_t_test(psam_pc_df, ap_thresh, 10) #ap thresh
psam_in_ap_10 <- run_t_test(psam_in_df, ap_thresh, 10)
psam_in_ap_20 <- run_t_test(psam_in_df, ap_thresh, 20)

#compile psam+ outputs
psam_stats <- list(psam_pc_vm_10, psam_pc_vm_20, psam_in_vm_10,
                   psam_pc_ir_10, psam_pc_ir_20, psam_in_ir_10,
                   psam_pc_sag_10, psam_pc_sag_20, psam_in_sag_10,
                   psam_pc_rheobase_10, psam_in_rheobase_10, psam_in_rheobase_20,
                   psam_pc_ap_10, psam_in_ap_10, psam_in_ap_20)
                   
                   
# stats for control neurons
cont_vm_1020 <- run_aov(control_pc_df, vm, 1020) #rmp
cont_pc_vm_10 <- run_t_test(control_pc_df, vm, 10) 
cont_pc_vm_20 <- run_t_test(control_pc_df, vm, 20)
cont_in_vm_10 <- run_t_test(control_in_df, vm, 10)

cont_ir_1020 <- run_aov(control_pc_df, ir_neg, 1020) #ir
cont_pc_ir_10<- run_t_test(control_pc_df, ir_neg, 10) 
cont_pc_ir_20 <- run_t_test(control_pc_df, ir_neg, 20)
cont_in_ir_10 <- run_t_test(control_in_df, ir_neg, 10)

cont_sag_1020 <- run_aov(control_pc_df, sag, 1020) #sag
cont_pc_sag_10 <- run_t_test(control_pc_df, sag, 10) 
cont_pc_sag_20 <- run_t_test(control_pc_df, sag, 20)
cont_in_sag_10 <- run_t_test(control_in_df, sag, 10)

cont_rheobase_1020 <- run_aov(control_pc_df, rheobase, 1020) #rheobase
cont_pc_rheobase_10 <- run_t_test(control_pc_df, rheobase, 10) 
cont_pc_rheobase_20 <- run_t_test(control_pc_df, rheobase, 20)
cont_in_rheobase_10 <- run_t_test(control_in_df, rheobase, 10)

cont_apthresh_1020 <- run_aov(control_pc_df, ap_thresh, 1020) #rheobase
cont_pc_apthresh_10 <- run_t_test(control_pc_df, ap_thresh, 10) 
cont_pc_apthresh_20 <- run_t_test(control_pc_df, ap_thresh, 20)
cont_in_apthresh_10 <- run_t_test(control_in_df, ap_thresh, 10)


#compile psam- outputs
control_psam_stats <- list(cont_vm_1020, cont_pc_vm_10, cont_pc_vm_20, cont_in_vm_10,
                      cont_ir_1020, cont_pc_ir_10, cont_pc_ir_20, cont_in_ir_10,
                      cont_sag_1020, cont_pc_sag_10, cont_pc_sag_20, cont_in_sag_10,
                      cont_rheobase_1020, cont_pc_rheobase_10, cont_pc_rheobase_20, cont_in_rheobase_10,
                      cont_apthresh_1020, cont_pc_apthresh_10, cont_pc_apthresh_20, cont_in_apthresh_10)


#save output to .txt file(s)
capture.output(psam_stats, file='psam_stats.txt')
capture.output(control_psam_stats, file='control_psam_stats.txt')

# --------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------PSAM-ChR2 analysis-------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------
