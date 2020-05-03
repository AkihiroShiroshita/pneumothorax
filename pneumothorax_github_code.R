#Multiple imputation using mice package
#Generating 100 imputed datasets
m <- 100 　　　　
imp <- mice(data_analysis, m = m ,seed=12345 ,maxit=50 ,printFlag=FALSE)
#Calculating propensity score
imp_stack <- complete(imp, action="long") %>% 
  as_data_frame
imp_stack<-imp_stack %>%
  group_by(.imp) %>%
  nest() %>%
  mutate(ps=map(data,function(df){
    ps_model<-glm(midclavicular_approach ~ place+department+age+gender+bmi+copd+ild+cancer+bronchiectasis,family=binomial,data=df)
    ps<-predict(ps_model,type="response")
    return(ps)
  }))%>%
  unnest()
#Propensity score matching using inverse probability of treatment weighting and matching weight
imp_stack<-imp_stack %>%
  group_by(.imp)%>%
  mutate(iptw=if_else(midclavicular_approach ==1,
                      1/ps,
                      1/(1-ps)),
         mw=pmin(ps,1-ps)*iptw,
         min_ps_treated=min(ps[midclavicular_approach==1]),
         max_ps_treated=max(ps[midclavicular_approach==1]),
         within_restriction=between(ps,unique(min_ps_treated),unique(max_ps_treated)),
         within_restriction=TRUE)
#Sample size calculation
fit.alt1 <-glm(malposition ~ 1,data=subset(imp_stack,midclavicular_approach==1),family=gaussian)
summary(fit.alt1)
fit.alt2 <-glm(malposition ~ 1,data=subset(imp_stack,midclavicular_approach==0),family=gaussian)
summary(fit.alt2)
#Balance check with stacked complete dataset
imp_data_wt_range_restricted<-imp_stack%>%
  filter(within_restriction)
myVars<-c("age", "gender", "bmi","place","department", "copd", "ild", "cancer", "bronchiectasis", "malposition")
svy_imp_stack_iptw <- svydesign(ids=~1,data=imp_data_wt_range_restricted,weights=~iptw)
svy_imp_stack_mw <-svydesign(ids=~1,data=imp_data_wt_range_restricted,weights =~mw)
tab_iptw_imp_stack <- svyCreateTableOne(vars=myVars,strata="midclavicular_approach",data=svy_imp_stack_iptw,test=FALSE)
print(tab_iptw_imp_stack,smd=TRUE)
tab_mw_imp_stack<-svyCreateTableOne(vars=myVars,strata="midclavicular_approach",data=svy_imp_stack_mw,test=FALSE)
print(tab_mw_imp_stack,smd=TRUE)
#Creating weighted datasets with range restriction
imp_stack<-imp_stack%>%
  group_by(.imp) %>%
  nest()%>%
  mutate(svy_data_iptw=map(data,function(df){
    svydesign(ids=~1,data=filter(df,within_restriction),weights=~iptw)
  }),
  svy_data_mw=map(data,function(df){
    svydesign(ids=~1,data=filter(df,within_restriction),weights=~mw)
  }))
#Balance check within each imputed dataset
imp_stack %>%
  mutate(smd_iptw=map(svy_data_iptw,function(svy_data) {
    smd<-ExtractSmd(svyCreateTableOne(vars=c("place","department","age", "gender","copd","ild","cancer","bronchiectasis","malposition"),strata="midclavicular_approach",data=svy_data, test=FALSE))
    df<-as_data_frame(smd) %>%
      mutate(name=rownames(smd))
    names(df)<-c("smd","name")
    df
  }),
  smd_mw=map(svy_data_mw,function(svy_data){
    smd<-ExtractSmd(svyCreateTableOne(vars = c("place","department","age", "gender","copd","ild","cancer","bronchiectasis","malposition"),strata="midclavicular_approach",data=svy_data,test=FALSE))
    df<-as_data_frame(smd) %>%
      mutate(name=rownames(smd))
    names(df)<-c("smd","name")
    df
  })) %>%
  select(.imp, smd_iptw, smd_mw) %>%
  unnest() %>%
  gather(key=key,value=smd,smd,smd1) %>%
  mutate(key=if_else(key=="smd",
                     "iptw",
                     "mw")) %>%
  ggplot(mapping=aes(x=name,y=smd,color=key,group=.imp))+
  geom_line()+
  geom_hline(yintercept=0.1)+
  facet_grid(.~key)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5),
        legend.key=element_blank(),
        plot.title=element_text(hjust=0.5),
        strip.background=element_blank())
#Weighted outcome analysis
imp_stack_results<- imp_stack %>%
  mutate(outcome_model_iptw=map(svy_data_iptw,function(svy_data){
    svyglm(malposition ~ midclavicular_approach, family = binomial(link = "logit"),design=svy_data)
  }),
  outcome_model_mw=map(svy_data_mw,function(svy_data){
    svyglm(malposition ~ midclavicular_approach, family = binomial(link = "logit"),design=svy_data)
  })) %>%
  summarize(iptw=list(MIcombine(outcome_model_iptw)),
            mw=list(MIcombine(outcome_model_mw))) %>%
  gather()
names(imp_stack_results$value) <- imp_stack_results$key
imp_stack_results$value
#Weighted outcome analysis with double-adjustment
imp_stack_results <- imp_stack %>%
  mutate(outcome_model_iptw=map(svy_data_iptw,function(svy_data){
    svyglm(malposition ~ midclavicular_approach+department+age+gender+bmi+bronchiectasis, family = binomial(link = "logit"),design=svy_data)
  }),
  outcome_model_mw=map(svy_data_mw,function(svy_data){
    svyglm(malposition ~ midclavicular_approach, family = binomial(link = "logit"),design=svy_data)
  }))%>%
  summarize(iptw=list(MIcombine(outcome_model_iptw)),
            mw=list(MIcombine(outcome_model_mw))) %>%
  gather()
names(imp_stack_results$value) <- imp_stack_results$key
imp_stack_results$value


