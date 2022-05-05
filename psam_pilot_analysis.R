# function to plot and analyse PSAM pilot data
# reads in dataframe containing subthreshold, spike data, and PSP measurements (where relevant) for EC L1 & L2 neurons in PSAM+ and control animals
# columns are included specifying whether neuron, slice, or animal had expression of PSAM-GFP
# if a neuron stopped firing after application of uPSEM-792, this is recorded as 'no_spikes'

# import modules
library(tidyverse)
library(rstatix)
library(gridExtra)
library(ggpubr)

# read file
path <- "~/Desktop/psam_pilot/in_plots"
df <- read.csv(file="~/Desktop/psam_pilot/psam_df.csv", head=TRUE, sep=",")

#add absolute values for ipsp 
df$ipsp_abs <- abs(df$ipsp)
df$ipsp_max_abs <- abs(df$ipsp_max)

#subset based on patterns PSAM expression
psam_df <- subset(df, psam_expression_cell=="Y") #PSAM+ neurons
control_df <- subset(df, psam_expression_cell=="N") # PSAM- neurons
psam_pc_df <- subset(psam_df, cell_type != "l2_in"); psam_pc_df <- subset(psam_pc_df, cell_type != "l1_in") # PSAM+ principal neurons 
psam_in_df <- subset(psam_df, cell_type == "l2_in" | cell_type == "l1_in") # PSAM+ interneurons
control_pc_df <- subset(control_df, cell_type != "l2_in"); control_pc_df <- subset(control_pc_df, cell_type != "l1_in") # data from slices with PSAM-GFP expression
control_in_df <- subset(control_df, cell_type != "l2_pc"); control_in_df <- subset(control_in_df, cell_type != "l2_sc")


# --------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------- Plot comparative data ----------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------

# plotting functions - adapted from lineplot code originally written by Sau Tsoi, 2021
# TODO make col input flexible - can't figure out why this doesn't work? Atm specify at lines **

#function for producing lineplot with individual observations for x 3 conditions
lp_3_conds <- function(df, xlabel1, xlabel2, ylabel){
  
  df <- subset(df, condition==20) # extract 3 cond mice

  #convert condition into a factor with the levels in the correct order
  df$drug <- as.character(df$drug)
  df$drug <- factor(df$drug, levels=unique(df$drug))
  
  #summarise data 
  sum <-group_by(df, drug, id) %>% 
    summarise(mean = mean(vm)) #** specify here
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
  
  ggsave(path=path, file=paste("rmp","expressing", "20nM_rmp.eps", sep="_"), width =50, height =60, units="mm") #save plot. **specify here
}

#lp_3_conds(control_in_df, "uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)")
lp_3_conds(psam_in_df, "uPSEM-792 (20 nM)", "Washout", "RMP (mV)")


#function for producing lineplot with individual observations for x 4 conditions
lp_4_conds <- function(df, xlabel1, xlabel2, xlabel3, ylabel){
  
  df <- subset(df, condition==1020) # extract 3 cond mice
  
  #convert condition into a factor with the levels in the correct order
  df$drug <- as.character(df$drug)
  df$drug <- factor(df$drug, levels=unique(df$drug))
  
  #summarise data 
  sum <-group_by(df, drug, id) %>% 
    summarise(mean = mean(vm)) #** specify here
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
  
  ggsave(path=path, file=paste("rmp","expressing", "10-20nM.eps", sep="_"), width =60, height =60, units="mm") #save plot. **specify here
}

#lp_4_conds(control_in_df, "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)")
#lp_4_conds(psam_in_df, "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "RMP (mV)")


#function for producing lineplot with individual observations for x 5 conditions
lp_5_conds <- function(df, xlabel1, xlabel2, xlabel3, xlabel4, ylabel){
  
  df <- subset(df, condition==1) 
  
  #convert condition into a factor with the levels in the correct order
  df$drug <- as.character(df$drug)
  df$drug <- factor(df$drug, levels=unique(df$drug))
  
  #summarise data 
  sum <-group_by(df, drug, id) %>% 
    summarise(mean = mean(ap_thresh)) #** specify here
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
  
  ggsave(path=path, file=paste("apthresh","expressing", "1-20nM.eps", sep="_"), width=70, height =60, units="mm") #save plot. **specify here
}
#lp_5_conds(control_pc_df, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "AP Thresh (mV)")
#lp_5_conds(psam_pc_df, "uPSEM-792 (1 nM)", "uPSEM-792 (10 nM)", "uPSEM-792 (20 nM)", "Washout", "RMP (mV)")
