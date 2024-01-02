#functions I stole from Sam
library(tidyverse)
library(lubridate)

zoop_predict<-function(model, data, confidence=95){
  prob_lower<-(100-confidence)/200
  probs<-c(prob_lower, 1-prob_lower)
  
  quantiles<-paste0("Q", probs*100)
  
  newdata<-expand_grid(SalSurfsc=quantile(data$SalSurfsc, probs=seq(0.05, 0.95, by=0.05), na.rm = T),
                       Secchisc=quantile(data$Secchisc, probs=seq(0.05, 0.95, by=0.05), na.rm = T),
                       logChlsc=quantile(data$logChlsc, probs=seq(0.05, 0.95, by=0.05), na.rm = T),
                       Month =1:12,
                       Year=2000)
  
  pred<-fitted(model, newdata=newdata, re_formula=NA, scale="response", probs=probs)
  
  newdata_pred<-newdata%>%
    mutate(Pred=pred[,"Estimate"],
           lowerCI=pred[,quantiles[1]],
           upperCI=pred[,quantiles[2]])%>%
    mutate(Month2=month(Month, label=T))
  
  return(newdata_pred)
}

zoop_plot<-function(data, type){
  
  require(ggplot2)
  require(dplyr)
  require(stringr)
  
  if(!type%in%c("Chla", "Secchi", "Month", "Salinity")){
    stop('Valid types are ""Chla", "Secchi", "Month", or "Salinity"')
  }

  
  if(type=="Chla"){
    data<-filter(data, SalSurfsc%in%unique(data$SalSurfsc)[seq(1,19, by=6)] & Year%in%plot_years1)
  } else{
    if(type=="year"){
      data<-filter(data, Salinity%in%unique(data$Salinity)[seq(1,19, by=6)] & Day==16)
    }else{
      data<-filter(data, Year%in%plot_years2 & Day==16)
    }
  }
  
  
  xvar<-case_when(type=="season" ~ "Julian_day", 
                  type=="year" ~ "Year", 
                  type=="salinity" ~ "Salinity")
  
  fillvar<-case_when(type=="season" ~ "Salinity", 
                     type=="year" ~ "Salinity", 
                     type=="salinity" ~ "Year")
  
  facetvar<-case_when(type=="season" ~ "Year", 
                      type=="year" ~ "Month2", 
                      type=="salinity" ~ "Month2")
  
  xlabel<-str_replace(xvar, "_", " ")
  
  data_orphans<-data%>%
    group_by(.data[[fillvar]], .data[[facetvar]])%>%
    mutate(Pred_lag=lag(Pred, order_by = .data[[xvar]]), Pred_lead=lead(Pred, order_by = .data[[xvar]]))%>%
    ungroup()%>%
    filter(is.na(Pred_lag) & is.na(Pred_lead) & !is.na(Pred))%>%
    select(-Pred_lag, -Pred_lead)
  
  if(type=="season"){
    scales<-list(scale_color_viridis_c(aesthetics = c("color", "fill"), trans="log", 
                                       breaks=unique(data$Salinity), 
                                       labels=round(unique(data$Salinity), 3),
                                       limits=range(data$Salinity)),
                 scale_x_continuous(breaks=c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335,
                                             15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349), 
                                    labels=c(rep("", 12), as.character(month(1:12, label = T))), limits=c(0,366),
                                    expand=expansion(0,0)),
                 theme(axis.ticks.x = element_line(color = c(rep("black", 12), rep(NA, 12))), axis.text.x=element_text(angle=45, hjust=1),
                       panel.grid.minor=element_blank(), panel.grid.major.x=element_line(color = c(rep("grey92", 12), rep(NA, 12)))))
  }else{
    if(type=="year"){
      scales<-list(scale_color_viridis_c(aesthetics = c("color", "fill"), trans="log", 
                                         breaks=unique(data$Salinity), 
                                         labels=round(unique(data$Salinity), 3),
                                         limits=range(data$Salinity)),
                   theme(axis.text.x=element_text(angle=45, hjust=1)))
      
    }else{
      scales<-list(scale_color_viridis_c(aesthetics = c("color", "fill")),
                   scale_x_continuous(trans="log", breaks=round(exp(seq(log(min(data$Salinity)), log(max(data$Salinity)), length.out=5)), 3),  minor_breaks = NULL),
                   theme(axis.text.x=element_text(angle=45, hjust=1)))
    }
    
  }
  
  p<-ggplot(data, aes(x=.data[[xvar]], y=Pred, ymin=lowerCI, ymax=upperCI, fill=.data[[fillvar]], group=.data[[fillvar]]))+
    geom_ribbon(alpha=0.4)+
    geom_line(aes(color=.data[[fillvar]]))+
    geom_pointrange(data=data_orphans, aes(color=.data[[fillvar]]), shape=21, position=position_dodge(width=2), size=0.2)+
    facet_wrap(~.data[[facetvar]], scales = "free_y")+
    scale_y_continuous(expand=c(0,0), limits = c(0, NA))+
    ylab("CPUE")+
    xlab(xlabel)+
    theme_bw()+
    scales
  
  
  return(p)
}

zoop_vario<-function(model, data, yvar, resid_type="standardized", cores=4){
  require(sp)
  require(gstat)
  require(spacetime)
  require(brms)
  require(dplyr)
  require(sf)
  
  if(resid_type=="standardized"){
    resids<-residuals(model, method="posterior_predict")/sd(data[[yvar]])
  }else{
    if(resid_type=="deviance"){
      resids<-residuals(model, method="posterior_predict")^2 # This seems to be the method used to calculate deviance residuals by mgcv
    }else
      stop("Only 'standardized' or 'deviance' resid_type values are accepted")
  }
  
  Data_vario<-data%>%
    mutate(Resid=resids[,"Estimate"])
  
  Data_coords<-Data_vario%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
    st_transform(crs=26910)%>%
    st_coordinates()%>%
    as_tibble()%>%
    mutate(across(c(X,Y), ~(.x-mean(.x))/1000))
  
  Data_vario<-bind_cols(Data_vario%>%
                          select(Date, Resid), Data_coords)
  sp<-SpatialPoints(coords=data.frame(X=Data_vario$X, Y=Data_vario$Y))
  sp2<-STIDF(sp, time=Data_vario$Date, 
             data=data.frame(Residuals=Data_vario$Resid))
  mb2M_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=cores, tlags=seq(0,30, by=2))
  
  return(mb2M_vario)
}