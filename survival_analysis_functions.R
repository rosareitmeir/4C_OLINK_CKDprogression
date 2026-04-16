### functions for survival analysis:

#### plot distribution of egfr/proteinuria for the different i-score groups: 
plot_boxplot_iscore <- function(data, proteinset,proteinnames,prefix="above_residual_",cutoff,grouped=T,xlabel=1.3,psize=3.5){
  ## title for plot
  title <- paste0("set: ", paste0(proteinnames,collapse=","))
  title <- paste0("cut-off:",cutoff, "\n", title)
  
  ## calculate i-score
  pnames <- paste0(prefix, proteinset)
  data<- data  %>% group_by(patid) %>% mutate(i_score= as.character(sum(c_across(pnames), na.rm = TRUE))) %>% ungroup() 
  ### grouping i-scores in three groups: 0, 1-2, 3-4
  if(grouped){
    data <- data %>% mutate(i_score=case_when(i_score=="1"|i_score=="2"~"1-2", 
                                              i_score=="3"|i_score=="4"~"3-4",
                                              T~"0"
    ))
  }
  ## plot: for each i-score grp a boxplot for eGFR and log(alburie)
  
  boxplot <- ggplot(data %>% mutate(log_alburie = log(alburie + 0.01)) %>%
                      gather(egfr, log_alburie, key = "key", value = "value") 
                    ,
                    aes(x = i_score, y = value, fill = i_score)) +
    geom_boxplot() +
    xlab("inflammation score") + theme_classic() +
    stat_compare_means(aes(label = paste(..method.., ", p =", ..p.format..)),
                       label.x=xlabel, size = psize) +
    ylab("") +
    theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
    facet_wrap(~ key, scales = "free_y") +
    ggtitle(title)
  ### if grouped , then overwrite colors for the filling
  if(grouped){
    boxplot <- boxplot + scale_fill_manual(values=c("0"="#FFD353FF","1-2"="#FFB242FF","3-4"="#BB292CFF"))
  } 
  else{
    boxplot <- boxplot + scale_fill_viridis_d()
  }
  
  return(boxplot)
}


####function KM curves plot

plot_kaplan_meier_iscore <- function(data, proteinset,proteinnames,cutoff,category=F,prefix="above_residual_",survline="none",timelimit=NULL,pvalsize=3,conf.int=F,title="",color,savepath=NULL,fontsize=5){
  ## calculate i-score 
  pnames <- paste0(prefix, proteinset)
  iscore_max<- as.character(length(proteinset))
  data<- data  %>% group_by(patid) %>% mutate(i_score= as.character(sum(c_across(pnames), na.rm = TRUE)))
  if (category){
    data <- data %>% mutate(i_score= factor(case_when(i_score==0 ~"0", 
                                                      i_score>=1 & i_score <3 ~"1-2",
                                                      T ~ paste0("3-",iscore_max)), levels=c("0","1-2",paste0("3-",iscore_max))))
    
  }
  # km model 
  km_model <- surv_fit( as.formula("Surv( survival_time,reachep)~ i_score"), data= data)
  
  ## 50% survival time of all strata, saving in table
  print(summary(km_model)$table[,"median"])
  
  ## for log-rank test for i-score grp 0 and max i-score grp (5)
  if(category){
    iscore_max <- paste0("3-",iscore_max)
  }
  formula<- as.formula("Surv( survival_time,reachep)~ i_score")
  logrank_test <- survdiff(formula, data =  data %>% filter(i_score=="0" | i_score==iscore_max))
  p_value <- round(1 - pchisq(logrank_test$chisq, df = 1),3)
  if (p_value< 0.01){
    p_value= "< 0.01"}
  
  if(title==""){
    title <- paste0("set: ", paste0(proteinnames,collapse=","))
    title <- paste0("cut-off:",cutoff, "\n", title)
  }
  ## plot curves
  if(is.null(timelimit)){
    km_plot <-  ggsurvplot(km_model, pval=paste0("grp_0 vs grp_max\nlog-rank p ",p_value),surv.median.line=survline, pval.size=pvalsize,conf.int=conf.int,risk.table=T, title=title, fontsize=fontsize, ggtheme=theme_classic(),palette=color,print=F)
  }
  else{
    km_plot <-  ggsurvplot(km_model, pval=paste0("grp_0 vs grp_max\nlog-rank p ",p_value),surv.median.line=survline,xlim=timelimit ,pval.size=pvalsize,conf.int=conf.int,risk.table=T, title=title, fontsize=fontsize, ggtheme=theme_classic(),palette=color,print=F)
    
  }
  
  if(! is.null(savepath)){
    fname <- paste0(savepath, "KM_cutoff",cutoff*100, ".pdf")
    ggsave(fname,arrange_ggsurvplots(list(km_plot)),width=7,height = 8)
  }
  
  return(km_plot)
}




### iterate cut off 

## plots kaplan-meier for each i-score cutoff and set of proteins 
# min/maxcutoff and step, which cutoffs should be used to define the i-scores
# proteinset: on which proteins shoudl be the iscore be calculated, # proteinnames: plot names for proteins
# which color for the curves and boxplotsm with confidene intervals as shadows or not, and save path 
iterate_cutoffs <- function(survdata, proteinset,proteinnames, mincutoff=0.5, maxcutoff=0.9,category=F, step=0.1,conf.int=F,color,savepath="~/Projects/4cstudy/Plots/figure5_survival/inflammation_score/kaplanmeier_cutoffs/"){
  
  km_list <- list()
  
  for (cutoff in seq(mincutoff, maxcutoff, step)){
    km_plot <- plot_kaplan_meier_iscore(data=surv_data,
                                        prefix = paste0("high_",cutoff*100,"_residual_"),
                                        proteinset=proteinset,proteinnames=proteinnames, cutoff = cutoff,
                                        conf.int = conf.int, color=color,savepath=savepath)
    
    km_list[[as.character(cutoff)]] <- km_plot
    km_list[[paste0(as.character(cutoff),"_bp")]] <- plot_boxplot_iscore(data=surv_data,
                                                                         prefix = paste0("high_",cutoff*100,"_residual_"),
                                                                         proteinset=proteinset,cutoff=cutoff)
    
  }
  
  return(km_list)
}

### plot residual
plot_residual <- function(protein,protein_name,survdata=surv_data,olinkdata=olink_data,savepath=NULL,
                          prefix="high_80_residual_",cutoff="80"){
  nativeprotein <- paste0("native_",protein)
  ylab <- paste0(protein_name, " NPX values")
  plot <-  survdata %>% left_join(olinkdata,by="patid")  %>%
    ggplot(aes(x=egfr,y=!!sym(nativeprotein))) +
    geom_point(aes(color=as.character(!!sym(paste0(prefix,protein)))))+ 
    geom_smooth(method="lm",se=F,color="gray33") + 
    theme_classic() + theme(legend.position = "top") + ylab(ylab) + xlab("eGFR")+
    scale_color_manual(values=c("lightblue", "darkblue"), labels=c("no","yes"), name=paste0("above ", cutoff,"% percentile") ) 
  
  if(!is.null(savepath)){
    ggsave(paste0(savepath,"residual_egfr_",protein_name,"_i80.pdf"),width=6.45,height=3.95 )}
  
  return(plot)
}

plot_npx_survtime <- function(protein,protein_name,survdata=surv_data,olinkdata=olink_data,savepath=NULL){
  nativeprotein <- paste0("native_",protein)
  ylab <- paste0(protein_name, " NPX values")
  
  ### plot npx values vs egfr and coloring survival time
  plot <- survdata %>% left_join(olinkdata,by="patid") %>% mutate(reachep=if_else(reachep==1,"reached","censored"))%>%  #filter(reachep==1) %>%
    ggplot(aes(x=egfr,y=!!sym(nativeprotein),color=survival_time,shape=reachep)) + geom_point() 
  
  ## aesthetics 
  plot <- plot+ ylab(ylab) + xlab("eGFR")+ scale_shape_manual(values=c(1,17),name="endpoint") +   theme_classic() + scale_colour_viridis_c(trans="log10", name="time to CKD endpoint") 
  
  if(!is.null(savepath)){
    ggsave(paste0(savepath,"npx_survtime_egfr_",protein_name,".pdf"),width=6.45,height=3.95 )}
  
  return(plot)
}

##  mutate(survtime_category=case_when(survival_time<=1~"1", survival_time>1 | survival_time <= 2 ~ "2",survival_time>2 | survival_time <= 3 ~ "3",T ~ "4")) 
