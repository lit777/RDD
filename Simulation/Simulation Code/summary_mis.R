library(tidyverse)
library(gridExtra)

type_ZIP_trial100 = read.csv('100_type_ZIP_trial_500_mis.csv')
type_ZIP_trial300 = read.csv('300_type_ZIP_trial_500_mis.csv')
type_ZIP_trial500 = read.csv('500_type_ZIP_trial_500_mis.csv')
type_ZIP_trial1000 = read.csv('1000_type_ZIP_trial_500_mis.csv')
type_ZIP_trial3000 = read.csv('3000_type_ZIP_trial_500_mis.csv')



make_result = function(dataframe, truevalue, sample){
  
  dataframe = dataframe[complete.cases(dataframe),]
  ind = which(dataframe$zipest!= Inf) 

  mse_local = round(mean((dataframe$zipest[ind]-truevalue)^2), 3)
  mse_lin_IK = round(mean((dataframe$linest_IK_r-truevalue)^2), 3)
  mse_lin = round(mean((dataframe$linest-truevalue)^2), 3)
  
  abs_local = round(mean(abs(dataframe$zipest[ind]-truevalue)/truevalue), 3)
  abs_lin_IK = round(mean(abs(dataframe$linest_IK_r-truevalue)/truevalue), 3)
  abs_lin = round(mean(abs(dataframe$linest - truevalue)/truevalue), 3)
  
  ci_local = round(mean(dataframe$zip_ci_left[ind] <= truevalue & dataframe$zip_ci_right[ind] >= truevalue, na.rm=T), 3)
  ci_lin_IK = round(mean(dataframe$lin_ci_left_IK_r <= truevalue & dataframe$lin_ci_right_IK_r >= truevalue, na.rm=T), 3)
  ci_lin = round(mean(dataframe$lin_ci_left <= truevalue & dataframe$lin_ci_right >= truevalue, na.rm=T), 3)
  
  result = data.frame('num_sam' = sample, 'mselocal' = mse_local, 'abslocal' = c(abs_local), 'cilocal' = c(ci_local), 'mselin_IK' = mse_lin_IK, 'abslin_IK' = c(abs_lin_IK), 'cilin_IK' = c(ci_lin_IK), 'mselin' = mse_lin, 'abslin' = abs_lin, 'cilin' = ci_lin)
  
  return(result)
}



truediff_zip = exp(1.1)/(1+exp(-2.3)) - exp(-0.1)/(1+exp(-0.5))
zip100 = make_result(type_ZIP_trial100 , truediff_zip, 100)
zip300 = make_result(type_ZIP_trial300 , truediff_zip, 300)
zip500 = make_result(type_ZIP_trial500, truediff_zip, 500)
zip1000 = make_result(type_ZIP_trial1000, truediff_zip, 1000)
zip3000 = make_result(type_ZIP_trial3000, truediff_zip, 3000)

zip = rbind(zip100, zip300, zip500, zip1000, zip3000)
colnames(zip) = c('N', 'local likelihood mse', 'local likelihood arb', 'local likelihood coverage', 'local linear IK mse', 'local linear IK arb', 'local linear IK coverage', 'local linear mse', 'local linear arb', 'local linear coverage')
zip


library(ggpubr)
make_graph = function(data){
  
  mse = data %>% gather(model, mse, contains("mse")) %>% mutate(.data = ., model = str_sub(string = model, start = -30, end = -4)) %>% 
    ggplot() + geom_line(mapping = aes(x=factor(N), y= mse, color=model, group=model), size=1) +
    geom_point(mapping = aes(x = factor(N), y=mse, color=model, group=model), size=3) +theme_bw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
          #axis.ticks.x=element_blank(), 
          #axis.ticks.y=element_blank(), 
          axis.text.y=element_text(size=8), 
          axis.text.x=element_text(size=8),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          legend.title = element_blank(), legend.text = element_text(size=10), legend.position='bottom') + scale_color_brewer(palette= "Set1")+ xlab('Sample size') + ylab('MSE')
  
  arb = data %>% gather(model, arb, contains("arb")) %>%  mutate(.data = ., model = str_sub(string = model, start = -30, end = -4)) %>% 
    ggplot() + geom_line(mapping = aes(x=factor(N), y= arb, color=model, group=model), size=1) +
    geom_point(mapping = aes(x = factor(N), y=arb, color=model, group=model), size=3) +theme_bw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
          #axis.ticks.x=element_blank(), 
          #axis.ticks.y=element_blank(), 
          axis.text.y=element_text(size=8), 
          axis.text.x=element_text(size=8),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          legend.title = element_blank(), legend.text = element_text(size=10), legend.position='bottom') + scale_color_brewer(palette= "Set1")+ xlab('Sample size') + ylab('ARB')
  
  coverage = data %>% gather(model, coverage, contains("coverage")) %>%  mutate(.data = ., model = str_sub(string = model, start = -35, end = -10)) %>% 
    ggplot() + geom_line(mapping = aes(x=factor(N), y= coverage, color=model, group=model), size=1) +
    geom_point(mapping = aes(x = factor(N), y=coverage, color=model, group=model), size=3) +theme_bw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
          #axis.ticks.x=element_blank(), 
          #axis.ticks.y=element_blank(), 
          axis.text.y=element_text(size=8), 
          axis.text.x=element_text(size=8),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          legend.title = element_blank(), legend.text = element_text(size=10), legend.position='bottom') + ylim(0, 1) +  scale_color_brewer(palette= "Set1")+ geom_hline(yintercept=0.95, linetype=2) + xlab('Sample size') + ylab('Coverage')
  
  return(ggarrange(mse, arb, coverage, nrow=1, common.legend = TRUE, legend="bottom"))
  
}


make_graph(zip)


