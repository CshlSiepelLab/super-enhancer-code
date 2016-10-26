rm(list=ls())
dir.create("plots")
plot.prefix="wap"
source("enhancerLib.R")

### WARNING: MUST SET WORKING DIRECTORY TO DIRECTORY OF THIS SCRIPT
## Read in Data
enh=fread("WAP_full_data.csv")
enh.dat=melt(enh)
setnames(enh.dat,c("variable","value"),c("condition","expression"))
enh.dat[,condition:=gsub("\\?","",enh.dat$condition)]
enh.dat=enh.dat[!is.na(expression)]
## enh.dat=enh.dat[!grepl("m|4",condition)]
enh.dat[,E1:=as.numeric(!grepl("1",condition))]
enh.dat[,E2:=as.numeric(!grepl("2",condition))]
enh.dat[,E3:=as.numeric(!grepl("3",condition))]

## Test linear model with log-normal error
linear.ln=deOptim.mod(response = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
                      ll.fun=linear.mod, parms = runif(5),lower = c(rep(-200,4),10^-2),
                      upper = 200,error.mod="log-normal",refine=TRUE)
pred.linear=linear.predict(design.mat = enh.dat[,3:5,with=FALSE],x = linear.ln$par)
pdf(file.path("plots",paste0(plot.prefix,"_linear_ln_pvr.pdf")))
pvr.plot(pred.linear$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_linear_ln_model.pdf")))
    plot.model(model = linear.ln,observed = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
           val.fun = linear.val,error.type="log-normal")
dev.off()

## Test multiplicative model with log-normal error
multi.ln=deOptim.mod(response = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
                      ll.fun=multiplicative.model, parms = runif(5),lower = c(rep(-200,4),10^-2),
                      upper = 200,error.mod="log-normal",refine=TRUE)
pred.multi=multiplicative.predict(design.mat = enh.dat[,3:5,with=FALSE],
                                  x = multi.ln$par)
pdf(file.path("plots",paste0(plot.prefix,"_multi_ln_pvr.pdf")))
pvr.plot(pred.multi$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_multi_ln_model.pdf")))
plot.model(model = multi.ln,observed = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
           val.fun = multiplicative.val,error.type="log-normal")
dev.off()

## Test logistic model with log-normal error
log.ln=deOptim.mod(response = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
                   ll.fun=logistic.model, parms = runif(6),lower = c(rep(-200,5),10^-2),
                   upper = 200,error.mod="log-normal",refine=TRUE)
pred.log=logistic.predict(design.mat = enh.dat[,3:5,with=FALSE],x = log.ln$par)
pdf(file.path("plots",paste0(plot.prefix,"_logistic_ln_pvr.pdf")))
pvr.plot(pred.log$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_logistic_ln_model.pdf")))
plot.model(model = log.ln,observed = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
           val.fun = logistic.val,scale=log.ln$par[5],error.type="log-normal")
dev.off()

## Plot of all residuals
all.mod=rbind(data.table(enh.dat$expression,pred.linear$expression,"linear"),
              data.table(enh.dat$expression,pred.multi$expression,"independant"),
              data.table(enh.dat$expression,pred.log$expression,"logistic"))
setnames(all.mod,colnames(all.mod),c("obs","pred","mod"))
pdf(file.path("plots",paste0(plot.prefix,"_ln_residuals.pdf")))
plot.residual(observed = all.mod$obs, all.mod$pred,all.mod$mod)
dev.off()

## Test linear model with normal error
linear.norm=deOptim.mod(response = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
                        ll.fun=linear.mod, parms = runif(5),lower = c(rep(-200,4),10^-2),
                        upper = 200,error.mod="gaussian",refine=TRUE)
pred.linear=linear.predict(design.mat = enh.dat[,3:5,with=FALSE],x = linear.norm$par)
pdf(file.path("plots",paste0(plot.prefix,"_linear_norm_pvr.pdf")))
pvr.plot(pred.linear$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_linear_norm_model.pdf")))
plot.model(model = linear.norm,observed = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
           val.fun = linear.val,error.type="gaussian")
dev.off()


## Test multiplicative model with normal error
multi.norm=deOptim.mod(response = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
                       ll.fun=multiplicative.model, parms = runif(5),lower = c(rep(-200,4),10^-2),
                       upper = 200,error.mod="gaussian",refine=TRUE)
pred.multi=multiplicative.predict(design.mat = enh.dat[,3:5,with=FALSE],
                                  x = multi.norm$par)
pdf(file.path("plots",paste0(plot.prefix,"_multi_norm_pvr.pdf")))
pvr.plot(pred.multi$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_multi_norm_model.pdf")))
plot.model(model = multi.norm,observed = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
           val.fun = multiplicative.val,error.type="gaussian")
dev.off()

## Test logistic model with normal error
log.norm=deOptim.mod(response = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
                           ll.fun=logistic.model, parms = runif(6),lower = c(rep(-200,5),10^-2),
                           upper = 200,error.mod="gaussian",refine=TRUE)
pred.log=logistic.predict(design.mat = enh.dat[,3:5,with=FALSE],x = log.norm$par)
pdf(file.path("plots",paste0(plot.prefix,"_logistic_norm_pvr.pdf")))
pvr.plot(pred.log$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_logistic_norm_model.pdf")))
plot.model(model = log.norm,observed = enh.dat$expression,design.mat = enh.dat[,3:5,with=FALSE],
           val.fun = logistic.val,scale=log.norm$par[5],error.type="gaussian")
dev.off()

## Plot of all residuals
all.mod=rbind(data.table(enh.dat$expression,pred.linear$expression,"linear"),
              data.table(enh.dat$expression,pred.multi$expression,"independant"),
              data.table(enh.dat$expression,pred.log$expression,"logistic"))
setnames(all.mod,colnames(all.mod),c("obs","pred","mod"))
pdf(file.path("plots",paste0(plot.prefix,"_norm_residuals.pdf")))
plot.residual(observed = all.mod$obs, all.mod$pred,all.mod$mod)
dev.off()

############################
## Hierarchical model
hier2.des=cbind(enh.dat[,3:5,with=F],apply(enh.dat[,c(3,4),with=FALSE],1,prod))
log.hier2.ln=deOptim.mod(response = enh.dat$expression,design.mat = hier2.des,
                        ll.fun=logistic.model, parms = runif(7),lower = c(rep(-200,6),10^-2),
                        upper = 200,error.mod="log-normal",refine=TRUE)
pred.log=logistic.predict(design.mat = hier2.des,x = log.hier2.ln$par)
pdf(file.path("plots",paste0(plot.prefix,"_logistic_hier12_ln_pvr.pdf")))
pvr.plot(pred.log$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_logistic_hier12_ln_model.pdf")))
plot.model(model = log.hier2.ln,observed = enh.dat$expression,design.mat = hier2.des,
           val.fun = logistic.val,scale=log.hier2.ln$par[6],error.type="log-normal")
dev.off()
bic(log.hier2.ln,nrow(enh.dat))

hier4.des=cbind(enh.dat[,3:5,with=F],apply(enh.dat[,c(4,5),with=FALSE],1,prod))
log.hier4.ln=deOptim.mod(response = enh.dat$expression,design.mat = hier4.des,
                         ll.fun=logistic.model, parms = runif(7),lower = c(rep(-200,6),10^-2),
                         upper = 200,error.mod="log-normal",refine=TRUE)
pred.log=logistic.predict(design.mat = hier4.des,x = log.hier4.ln$par)
pdf(file.path("plots",paste0(plot.prefix,"_logistic_hier23_ln_pvr.pdf")))
pvr.plot(pred.log$expression,enh.dat$expression,enh.dat$condition)
dev.off()
pdf(file.path("plots",paste0(plot.prefix,"_logistic_hier23_ln_model.pdf")))
plot.model(model = log.hier4.ln,observed = enh.dat$expression,design.mat = hier4.des,
           val.fun = logistic.val,scale=log.hier4.ln$par[6],error.type="log-normal")
dev.off()
bic(log.hier4.ln,nrow(enh.dat))

#############################


## Calculate relative BIC of independant-logistic models
## Evidence against higher BIC is as follows:
## 0 to 2   Not worth more than a bare mention
## 2 to 6   Positive
## 6 to 10   Strong
## >10   Very Strong

bic(multi.ln,nrow(enh.dat))-bic(log.ln,nrow(enh.dat))
bic(linear.ln,nrow(enh.dat))-bic(log.ln,nrow(enh.dat))

bic(multi.norm,nrow(enh.dat))-bic(log.norm,nrow(enh.dat))
bic(linear.norm,nrow(enh.dat))-bic(log.norm,nrow(enh.dat))

min(bic(multi.norm,nrow(enh.dat)),bic(multi.ln,nrow(enh.dat)))-
  min(bic(log.norm,nrow(enh.dat)),bic(log.ln,nrow(enh.dat)))

## Some output writing
linear.ln$bic=bic(linear.ln,nrow(enh.dat))
multi.ln$bic=bic(multi.ln,nrow(enh.dat))
log.ln$bic=bic(log.ln,nrow(enh.dat))
linear.norm$bic=bic(linear.norm,nrow(enh.dat))
multi.norm$bic=bic(multi.norm,nrow(enh.dat))
log.norm$bic=bic(log.norm,nrow(enh.dat))

mod.list=list(linear.ln,multi.ln,log.ln,linear.norm,multi.norm,log.norm)
unlist(lapply(mod.list,function(x) return(x$bic)))

lapply(mod.list,function(x) return(x$par))

save.image("wap.RData")
