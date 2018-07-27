# Load Librarys
library("xtable")
library("loo")
library("rstan")
library("ggplot2")
library("fludata")
library("dplyr")

dat <- data.frame(y=OG$dat$y,
                  week=OG$dat$x,
                  num=OG$dat$num,
                  group=OG$dat$og_group,
                  seas=OG$dat$og_seas)

load("~/Downloads/ForecastRdata/Common_OG.RData")
load("~/Downloads/ForecastRdata/GroupCS_OG.RData")
load("~/Downloads/ForecastRdata/Indep_OG.RData")

res_common   <- extract(Common_OG)
res_vague    <- extract(Indep_OG)
res_groupCS  <- extract(GroupCS_OG)

###### Comparison via WAIC and LOO #####
ll_Common      <- extract_log_lik(Common_OG)
ll_Vague       <- extract_log_lik(Indep_OG)
ll_LKJGroupCS  <- extract_log_lik(GroupCS_OG)

waicCommon <- waic(ll_Common)
waicVague  <- waic(ll_Vague)
waicGroupCS <- waic(ll_LKJGroupCS)


get_res <- function(y, week, num, group, seas, res, model){
  resdat = data.frame(week = week, num = num, y = y, group = group, seas = seas, model = model)
  resdat$est <- resdat$upper <- resdat$lower <- resdat$pred <- resdat$pred_upper <- resdat$pred_lower <- MSE <- NA
  binoms <- matrix(data=NA, nrow=nrow(res$lpsi), ncol=nrow(resdat))
  for(i in 1:nrow(resdat)){
    # Est of mean
    resdat$est[i]   <- median(boot::inv.logit(res$lpsi[,i]))
    resdat$lower[i] <- quantile(boot::inv.logit(res$lpsi[,i]), probs=c(0.025))
    resdat$upper[i] <- quantile(boot::inv.logit(res$lpsi[,i]), probs=c(0.975))
    # Predictios
    binoms[,i] <- rbinom(length(res$lpsi[,i]), resdat$num[i], boot::inv.logit(res$lpsi[,i]))
    resdat$pred[i]       <- median(binoms[,i])
    resdat$pred_lower[i] <- quantile(binoms[,i], probs=c(0.025))
    resdat$pred_upper[i] <- quantile(binoms[,i], probs=c(0.975))
  }
  
  for(m in 1:nrow(res$lpsi)){
    MSE[m] <- mean((binoms[m,] - resdat$y)^2)
  }
  
  return(list(resdat = resdat, MSE = MSE))
}

resCommon   <- get_res(dat$y, dat$week, dat$num, dat$group, dat$seas, res_common, "Common")
resVague    <- get_res(dat$y, dat$week, dat$num, dat$group, dat$seas, res_vague, "Indep")
resGroupCS  <- get_res(dat$y, dat$week, dat$num, dat$group, dat$seas, res_groupCS, "RegMnSd")

cdat <- rbind(resCommon$resdat,resVague$resdat,resGroupCS$resdat)
cdat <- subset(cdat, seas=="15-16" & (group=="Region 9" | group=="Region 7"))


ggplot(cdat) + 
  geom_point(aes(x=week, y=y/num, size=num), shape=1) + 
  geom_line(aes(x=week, y=est, color=model)) + 
  geom_ribbon(aes(x=week, ymin=lower, ymax=upper, fill=model), alpha=0.1) + 
  facet_grid(~group)+
  theme_bw() + labs(x="Week", y="ILI Percentage")

##### MSE and WAIC #####

# Check Prediction inteval lengths
cdat %>% group_by(model, group) %>% summarise(PIlength = mean(pred_upper-pred_lower),
                                       CIlength = mean(upper-lower),
                                       MSE = sum(y-pred)^2)


resdf <- data.frame(cdat %>% group_by(model, seas) %>% summarise(PIlength = mean(pred_upper-pred_lower),
                                                                 CIlength = mean(upper-lower),
                                                                 MSE = mean(pred-y)^2))

ggplot(subset(resdf, model!="Common")) + geom_point(aes(x=seas, y=MSE, color=model))

ggplot(cdat, aes(x=week, y=num, color=seas)) + geom_point() + geom_smooth(method="lm") + facet_grid(~group)

resdf <- data.frame(MSE=c(mean((subset(cdat, model=="Common")$est - subset(cdat, model=="Common")$y/subset(cdat, model=="Common")$num)^2),
                          mean((subset(cdat, model=="Indep")$est - subset(cdat, model=="Indep")$y/subset(cdat, model=="Indep")$num)^2),
                          mean((subset(cdat, model=="RegMnSd")$est - subset(cdat, model=="RegMnSd")$y/subset(cdat, model=="RegMnSd")$num)^2))*1000000,
                    MSPE=c(mean((subset(cdat, model=="Common")$est - subset(cdat, model=="Common")$y)^2),
                           mean((subset(cdat, model=="Indep")$est - subset(cdat, model=="Indep")$y)^2),
                           mean((subset(cdat, model=="RegMnSd")$est - subset(cdat, model=="RegMnSd")$y)^2))
)


resdf$WAIC <- c(waicCommon$waic, waicVague$waic, waicGroupCS$waic)
resdf <- t(resdf)
colnames(resdf) <- c("Common", "Indep", "RegMnSd")

mresdf <- reshape2::melt(resdf)
ggplot(subset(mresdf, Var2!="Common")) + 
  geom_point(aes(x="none", y=value, color=Var2)) + 
  facet_wrap(~Var1, scales="free_y") + labs(y="", color="Model", x="") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

minres <- apply(resdf, 1, FUN=min)

resdf[1,] <- resdf[1,] - minres[1]
resdf[2,] <- resdf[2,] - minres[2]

print(xtable(resdf, caption="This table compares the different hierarchy structres accross MSE and WAIC;
             the minimum value of each statistic has been subtracted from each row to make comparisons easier. 
             The hierarchy structure using a region means outpreforms the rest as seen via MSE though when taking 
             into consideration model complexity via WAIC, it just slightly outpreforms the vague model.",
             digits=c(0,2,2,2,2,2,2,2,2)))
