# 04_simulation_figures.R
# Publication/thesis figures from completed simulation results. No refitting.

dir.create("thesis_analysis/figures_v9006", recursive=TRUE, showWarnings=FALSE)
stopifnot(requireNamespace("ggplot2",quietly=TRUE), requireNamespace("scales",quietly=TRUE))
library(ggplot2)
df <- read.csv("thesis_analysis/tables_v9006/simulation_full_factorial_summary.csv",check.names=FALSE)
df$criterion <- factor(df$criterion,levels=c("AIC","BIC","CAIC","ICOMP_IFIM","CICOMP"))
df$cov_method <- toupper(df$cov_method); df$true_d<-factor(df$true_d); df$snr<-factor(df$snr); df$rho_x<-factor(df$rho_x)

save_plot <- function(p,nm,w=11,h=7){
  ggsave(file.path("thesis_analysis/figures",paste0(nm,".pdf")),p,width=w,height=h,units="in")
  ggsave(file.path("thesis_analysis/figures",paste0(nm,".png")),p,width=w,height=h,units="in",dpi=300)
}

p1 <- ggplot(df,aes(factor(p),dimension_recovery_rate,group=criterion,linetype=criterion))+
  geom_line()+geom_point()+
  facet_grid(paste0("d0=",true_d,", rho=",rho_x) ~ paste0("n=",n,", SNR=",snr,", ",cov_method))+
  scale_y_continuous(limits=c(0,1),labels=scales::label_percent())+
  labs(x="Number of predictors (p)",y="P(hat(d)=d0)",linetype="Criterion",
       title="Structural-dimension recovery across simulation conditions")+theme_bw(base_size=9)+theme(legend.position="bottom")
save_plot(p1,"fig_dimension_recovery",15,10)

heat_dim <- aggregate(dimension_recovery_rate~criterion+cov_method+n+p,df,mean)
p2 <- ggplot(heat_dim,aes(criterion,factor(p),fill=dimension_recovery_rate))+
  geom_tile()+geom_text(aes(label=sprintf("%.2f",dimension_recovery_rate)),size=3)+
  facet_grid(cov_method~n)+scale_fill_viridis_c(limits=c(0,1))+
  labs(x="Criterion",y="p",fill="Recovery",title="Mean structural-dimension recovery probability")+
  theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot(p2,"fig_dimension_recovery_heatmap",10,6)

p3 <- ggplot(df,aes(factor(p),mean_rmse,group=criterion,linetype=criterion))+
  geom_line()+geom_point()+
  facet_grid(paste0("d0=",true_d,", rho=",rho_x) ~ paste0("n=",n,", SNR=",snr,", ",cov_method))+
  labs(x="Number of predictors (p)",y="Mean test RMSE",linetype="Criterion",
       title="Predictive RMSE across simulation conditions")+theme_bw(base_size=9)+theme(legend.position="bottom")
save_plot(p3,"fig_prediction_rmse",15,10)

heat_sub <- aggregate(mean_subspace_distance~criterion+cov_method+n+p,df,
                      function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE))
p4 <- ggplot(heat_sub,aes(criterion,factor(p),fill=mean_subspace_distance))+
  geom_tile()+geom_text(aes(label=ifelse(is.na(mean_subspace_distance),"NA",sprintf("%.2f",mean_subspace_distance))),size=3)+
  facet_grid(cov_method~n)+scale_fill_viridis_c(direction=-1,na.value="grey90")+
  labs(x="Criterion",y="p",fill="Distance",title="Mean subspace distance conditional on correct dimension recovery")+
  theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot(p4,"fig_subspace_distance_heatmap",10,6)

p5 <- ggplot(df,aes(factor(p),mean_f1_variable,group=criterion,linetype=criterion))+
  geom_line()+geom_point()+
  facet_grid(paste0("d0=",true_d,", rho=",rho_x) ~ paste0("n=",n,", SNR=",snr,", ",cov_method))+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Number of predictors (p)",y="Mean variable-selection F1",linetype="Criterion",
       title="Feature-selection performance across simulation conditions")+theme_bw(base_size=9)+theme(legend.position="bottom")
save_plot(p5,"fig_variable_selection_f1",15,10)

p6 <- ggplot(df,aes(factor(p),mean_n_selected,group=criterion,linetype=criterion))+
  geom_hline(yintercept=5,linewidth=.4)+geom_line()+geom_point()+facet_grid(cov_method~n)+
  labs(x="Number of predictors (p)",y="Mean number selected",linetype="Criterion",
       title="Selected-model cardinality; horizontal line denotes five true signals")+
  theme_bw()+theme(legend.position="bottom")
save_plot(p6,"fig_selected_model_size",10,6)

p7 <- ggplot(df,aes(mean_recall,mean_precision,shape=criterion))+geom_point()+facet_grid(cov_method~p)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(x="Mean recall",y="Mean precision",shape="Criterion",title="Feature-selection precision-recall trade-off")+
  theme_bw()+theme(legend.position="bottom")
save_plot(p7,"fig_precision_recall",10,6)

overall <- read.csv("thesis_analysis/tables_v9006/table_overall_method_summary.csv",check.names=FALSE)
overall$criterion <- factor(overall$criterion,levels=c("AIC","BIC","CAIC","ICOMP_IFIM","CICOMP")); overall$cov_method<-toupper(overall$cov_method)
p8 <- ggplot(overall,aes(criterion,dimension_recovery_rate,shape=cov_method,group=cov_method))+
  geom_point(size=3)+geom_line()+scale_y_continuous(limits=c(0,1),labels=scales::label_percent())+
  labs(x="Criterion",y="Overall dimension-recovery rate",shape="Covariance",
       title="Overall criterion comparison across all simulation scenarios")+
  theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot(p8,"fig_overall_dimension_recovery",8,5)
cat("Simulation figures generated in thesis_analysis/figures.\n")
