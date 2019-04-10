## ----setup, include=FALSE------------------------------------------------
# Package loacking
library(readr)
library(stringr)
library(knitr)
library(ggplot2)
library(dplyr)
library(lmtest)
library(sandwich)
library(broom)
library(car)
library(tidyr)
library(splines)
library(survival)
library(ROCR)
library(caret)

# Paths and other defined variables
data_dir <- "../Data/"

logistic <- function(x) 1/(1 + exp(-x))

## ----include=F-----------------------------------------------------------
get_ests_table <- function(model_obj){
  model_ests <- data.frame(names(coef(model_obj)),
  round(coef(model_obj), digits=2),
  round(sqrt(diag(sandwich(model_obj))), digits=2),
      round(exp(coef(model_obj)), 2),
      exp(confint(model_obj)))
model_ests$ci <- paste0("(",
                             round(model_ests[, 5], 2),
                             ",",
                             round(model_ests[, 6], 2),
                             ")")
model_ests <- model_ests[, -5:-6]
colnames(model_ests) <- c("Model Term", "Estimate", "SE",
                              "Exp(Est.)","95% CI")
return(model_ests)
}

## ----include=FALSE-------------------------------------------------------
pcs <- read_csv(paste0(data_dir, "pcs.csv"))

## ----echo=TRUE-----------------------------------------------------------
psa_m1 <- glm(tumor~psa, data=pcs, family=binomial)

## ----message=F, include=F------------------------------------------------
psa_m1_ests <- get_ests_table(psa_m1)
kable(psa_m1_ests, row.names=F)

## ----out.width="85%"-----------------------------------------------------
psa_m1_preds <- predict(psa_m1, type="response")

ggplot() + theme_bw() + 
   geom_point(aes(x=pcs$psa,
                          y=psa_m1_preds)) +
   ylab("Predicted Probability of Cancer Spread") +
   xlab("PSA Concentration")

## ----output.lines=1:5, size="footnotesize"-------------------------------
ggplot() + theme_bw() + 
   geom_point(aes(x=pcs$psa,
                          y=predict(psa_m1, type="response"))) +
   geom_point(aes(x=psa,
                  y=as.numeric(predict(psa_m1, type="response")>0.5)),
              col="red",
              data=pcs) + 
   ylab("Predicted Probability of Cancer Spread") +
   xlab("PSA Concentration")

## ----output.lines=1:5, size="footnotesize"-------------------------------
ggplot() + theme_bw() + 
   geom_point(aes(x=pcs$psa,
                          y=predict(psa_m1, type="response"))) +
   geom_point(aes(x=psa,
                  y=as.numeric(predict(psa_m1, type="response")>0.75)),
              col="red",
              data=pcs) + 
   ylab("Predicted Probability of Cancer Spread") +
   xlab("PSA Concentration")

## ----eval=F, echo=TRUE---------------------------------------------------
## cutoff <- 0.5
## yhat_c <- ifelse(predict(psa_m1,type="response") > cutoff, 1, 0)
## confusionMatrix(data=factor(yhat_c),
##                 reference=factor(pcs$tumor),
##                 positive="1")

## ----echo=F, output.lines=1:12-------------------------------------------
cutoff <- 0.5
yhat_c <- ifelse(predict(psa_m1,type="response") > cutoff, 1, 0)
confusionMatrix(data=factor(yhat_c),
                reference=factor(pcs$tumor),
                positive="1")

## ----echo=F, output.lines=-1:-12-----------------------------------------
confusionMatrix(data=factor(yhat_c),
                reference=factor(pcs$tumor),
                positive="1")

## ------------------------------------------------------------------------
calc_specsens <- function(yhat, y){
   agree <- yhat==y
   sens <- mean(agree[y==1])
   spec <- mean(agree[y==0])
   c(spec=spec, sens=sens)
}
calc_yhat <- function(phat, cut){
   yhat <- as.numeric(phat>=cut)
   yhat
}
calc_roc <- function(phat, y, cuts){
   
   rocs <- sapply(cuts, function(w) {yhats <- calc_yhat(phat=phat, cut=w); calc_specsens(yhat=yhats, y=y)} )
   out <- as.data.frame(t(rocs))
   out$cut <- cuts
   out
}

psa_m1_roc <- calc_roc(phat=psa_m1_preds, y=pcs$tumor, cuts=seq(0, 1, length=500))


ggplot() + theme_bw() +
   geom_point(aes(x=cut,
                  y=spec,
              col="spec"),
              data=psa_m1_roc) +
   geom_point(aes(x=cut,
                  y=sens,
                  col="sens"),
              data=psa_m1_roc) + ylab("") + xlab("Cut Point") + 
   scale_color_manual(values=c("red", "blue"), breaks=c("sens", "spec"), labels=c("Sensitivity", "Specificity"), name="")

## ----include=FALSE-------------------------------------------------------

table(true=pcs$tumor, psa_m1_preds>0.5)
table(pcs$tumor)
table(true=pcs$tumor, psa_m1_preds>0.75)

## ----out.width="90%"-----------------------------------------------------
ggroc <- ggplot() + theme_bw() + geom_abline(aes(intercept=0, slope=1), lty=2, col="grey50") + 
   scale_x_continuous(limits=0:1,
                     breaks=seq(0, 1, by=0.2),
                     expand=c(0, 0)) + 
  scale_y_continuous(limits=0:1,
                     breaks=seq(0, 1, by=0.2),
                     expand=c(0, 0)) + 
   ylab("Sensitivity") + xlab("1 - Specificity")
ggroc

## ------------------------------------------------------------------------
ggroc + geom_path(aes(x=1-spec, y=sens), data=psa_m1_roc)

## ----echo=TRUE, out.width="80%"------------------------------------------
predobj <- prediction(predictions = predict(psa_m1,
                                              type="response"),
                        labels=pcs$tumor)
perfobj <- performance(predobj, 'tpr', "fpr")
plot(perfobj)
abline(0, 1, col="grey")

## ----echo=TRUE, size="scriptsize"----------------------------------------
psa_m1 <- glm(tumor~psa, data=pcs, family=binomial)
psa_m2 <- glm(tumor~gleason, data=pcs, family=binomial)
psa_m3 <- glm(tumor~psa + gleason, data=pcs, family=binomial)
psa_m4 <- glm(tumor~psa + dpros + dcaps, data=pcs, family=binomial)
psa_m5 <- glm(tumor~race + age, data=pcs, family=binomial)
psa_m6 <- glm(tumor~psa + gleason +  dpros + dcaps, data=pcs, family=binomial)

## ----echo=TRUE-----------------------------------------------------------

pred_list <- prediction(list(m1=predict(psa_m1, type="response"),
                             m2=predict(psa_m2, type="response"),
                             m3=predict(psa_m3, type="response"),
                             m4=predict(psa_m4, type="response"),
                             m5=predict(psa_m5, type="response"),
                             m6=predict(psa_m6, type="response")),
           labels=list(pcs$tumor,
                       pcs$tumor,
                       pcs$tumor,
                       pcs$tumor,
                       pcs$tumor[!is.na(pcs$race)],
                       pcs$tumor))
perf_list <- performance(pred_list, 'tpr', "fpr")

## ----echo=TRUE-----------------------------------------------------------
plot(perf_list, colorize=F)
abline(0,1, col="grey")

## ----echo=TRUE-----------------------------------------------------------
aucs_obj <- performance(pred_list, 'auc')
aucs <- unlist(aucs_obj@y.values)
kable(data.frame(Model=1:6, AUC=aucs), digits=3)

## ----echo=TRUE-----------------------------------------------------------
perf_list <- performance(predobj, 'sens', 'spec')
plot(perf_list, colorize=T)

## ----echo=TRUE, out.width="80%"------------------------------------------
perf_list <- performance(predobj, 'prec', "rec")
plot(perf_list)

## ------------------------------------------------------------------------
set.seed(11)
x <- rep(1:10, each=10)
mufun <- function(x) {
   mu <- -2 + 0.5*x - 0.03 *x^2
   1/(1 + exp(-mu))
}
# y0 <- rep(rep(c(0, 1), times=5),
#           times=c(9, 1, 8, 2, 7, 3, 7, 3, 5, 5))
y0 <- rbinom(n=length(x), size=1, prob=mufun(x))

# table(x, y0)
g1 <- glm(y0~x + I(x^2), family=binomial)
g2 <- glm(y0~factor(x), family=binomial)

## ------------------------------------------------------------------------
ggplot() + theme_bw() +
   geom_point(aes(x=x, y=mufun(x), col="g1", shape="g1")) +
   geom_line(aes(x=x, y=mufun(x), col="g1")) +
   geom_point(aes(x=unique(x), y=tapply(y0, x, mean), col="g2", shape="g2")) +
geom_point(aes(x=x, y=predict(g1, type="response"), col="g3", shape="g3")) +
   geom_line(aes(x=x, y=predict(g1, type="response"), col="g3")) +
   scale_color_manual(breaks=c("g1", "g2", "g3"), labels=c("True Mean", "Observed Mean", "Modeled Mean"), values=c("red", "blue", "green"), name="Values") +
   scale_shape_manual(breaks=c("g1", "g2", "g3"), labels=c("True Mean", "Observed Mean", "Modeled Mean"), values=c(16, 17, 15), name="Values") + ylab("")


## ----eval=F, include=F---------------------------------------------------
## 
## cutoff <- 0.5
## yhat1_c <- ifelse(predict(g1,type="response") > cutoff, 1, 0)
## confusionMatrix(data=factor(yhat1_c),
##                 reference=factor(y0),
##                 positive="1")
## # Acc - 0.69
## 
## yhat2_c <- ifelse(predict(g2,type="response") > cutoff, 1, 0)
## confusionMatrix(data=factor(yhat2_c),
##                 reference=factor(y0),
##                 positive="1")
## # acc - 0.73
## 
## set.seed(12)
## 
## B <- 100
## acc1 <- numeric(B)
## acc2 <- acc1
## sens1 <- acc1
## sens2 <- acc1
## for (i in 1:B){
## ynew <- rbinom(n=length(x), size=1, prob=mufun(x))
## cm1 <- confusionMatrix(data=factor(yhat1_c),
##                 reference=factor(ynew),
##                 positive="1")
## acc1[i] <- cm1$overall["Accuracy"]
## sens1[i] <- cm1$byClass["Sensitivity"]
## cm2 <- confusionMatrix(data=factor(yhat2_c),
##                 reference=factor(ynew),
##                 positive="1")
## acc2[i] <- cm2$overall["Accuracy"]
## sens2[i] <- cm2$byClass["Sensitivity"]
## }
## mean(acc1) # 0.592
## mean(acc2) # 0.597
## mean(sens1) # 0.372
## mean(sens2) # 0.377

## ----echo=TRUE, size="footnotesize"--------------------------------------
set.seed(5)
k <- 10
flds <- createFolds(pcs$tumor, k = k, list = TRUE, returnTrain = FALSE)
str(flds)

## ----echo=TRUE-----------------------------------------------------------
# Create lists for storing predictions
# One element for each model
pred_list <- list(numeric(nrow(pcs)),
                  numeric(nrow(pcs)),
                  numeric(nrow(pcs)))

## ----echo=TRUE, size="footnotesize"--------------------------------------
for (i in 1:k){
   # Subset training data
   train <- pcs[-flds[[i]],]
   # Subset test data
   test <- pcs[flds[[i]],]
   # Fit each model
   psa_cvm1 <- glm(tumor~psa, data=train, family=binomial)
   psa_cvm2 <- glm(tumor~gleason, data=train, family=binomial)
   psa_cvm3 <- glm(tumor~psa + gleason, data=train, family=binomial)
   # Store predictions
   pred_list[[1]][flds[[i]]]=predict(psa_cvm1, type="response",
                                     newdata=test)
   pred_list[[2]][flds[[i]]]=predict(psa_cvm2, type="response",
                                     newdata=test)
   pred_list[[3]][flds[[i]]]=predict(psa_cvm3, type="response",
                                     newdata=test)
}

## ----echo=T--------------------------------------------------------------
predobj <- prediction(pred_list,
           labels=list(pcs$tumor,
                       pcs$tumor,
                       pcs$tumor))
perfobj <- performance(predobj, 'tpr', 'fpr')
plot(perfobj)

## ----echo=T--------------------------------------------------------------
aucs_cv <- performance(predobj, 'auc')
aucs <- unlist(aucs_cv@y.values)
aucs

## ----eval=F, include=F---------------------------------------------------
## # Create empty lists
## pred_list <- list(numeric(nrow(pcs)),
##                   numeric(nrow(pcs)),
##                   numeric(nrow(pcs)),
##                   numeric(nrow(pcs)),
##                   numeric(nrow(pcs)),
##                   numeric(nrow(pcs)))
## for (i in 1:10){
##    train <- pcs[-flds[[i]],]
##    test <- pcs[flds[[i]],]
##    psa_cvm1 <- glm(tumor~psa, data=train, family=binomial)
##    psa_cvm2 <- glm(tumor~gleason, data=train, family=binomial)
##    psa_cvm3 <- glm(tumor~psa + gleason, data=train, family=binomial)
##    psa_cvm4 <- glm(tumor~psa + dpros + dcaps, data=train, family=binomial)
##    psa_cvm5 <- glm(tumor~race + age, data=train, family=binomial)
##    psa_cvm6 <- glm(tumor~psa + gleason +  dpros + dcaps, data=train, family=binomial)
## 
## 
##    pred_list[[1]][flds[[i]]]=predict(psa_cvm1, type="response", newdata=test)
##    pred_list[[2]][flds[[i]]]=predict(psa_cvm2, type="response", newdata=test)
##    pred_list[[3]][flds[[i]]]=predict(psa_cvm3, type="response", newdata=test)
##    pred_list[[4]][flds[[i]]]=predict(psa_cvm4, type="response", newdata=test)
##    pred_list[[5]][flds[[i]]]=predict(psa_cvm5, type="response", newdata=test)
##    pred_list[[6]][flds[[i]]]=predict(psa_cvm6, type="response", newdata=test)
## }
## 
## predobj <- prediction(pred_list,
##            labels=list(pcs$tumor,
## pcs$tumor,
## pcs$tumor,
## pcs$tumor,
## pcs$tumor,
## pcs$tumor))
##    perfobj <- performance(predobj, 'tpr', 'fpr')
## plot(perfobj)
## 
## aucs_cv <- performance(predobj, 'auc')
## aucs <- unlist(aucs_cv@y.values)
## aucs

## ----include=FALSE-------------------------------------------------------
glow <- read_csv(paste0(data_dir, "glow.csv"))

## ----echo=TRUE-----------------------------------------------------------
set.seed(10)
trainfull_inds <- createDataPartition(glow$fracture, p=0.75)
trainfull <- glow[trainfull_inds[[1]], ]
test_df <- glow[-trainfull_inds[[1]], ]
train_inds <- createDataPartition(trainfull$fracture, p=0.66)
train_df <- trainfull[train_inds[[1]], ]
validate_df <- trainfull[-train_inds[[1]], ]
c(nrow(train_df), nrow(validate_df), nrow(test_df))

## ----echo=TRUE, size="footnotesize"--------------------------------------
glow_m1 <- glm(fracture~age + smoke + pre_meno +
                  prior_frac + mom_frac,
               data=train_df, family=binomial)
glow_m2 <- glm(fracture~age + prior_frac*mom_frac,
               data=train_df, family=binomial)
glow_m3 <- glm(fracture~age + bmi + smoke +
                  prior_frac + mom_frac + arm_assist,
               data=train_df, family=binomial)
glow_m4 <- glm(fracture~ns(age, 4) +  bmi + pre_meno +
                  prior_frac*mom_frac + arm_assist +
                  height + weight + smoke,
               data=train_df, family=binomial)

## ----echo=TRUE, size="footnotesize"--------------------------------------
predobj <- prediction(list(m1=predict(glow_m1,
                                      newdata=validate_df,
                                      type="response"),
                           m2=predict(glow_m2,
                                      newdata=validate_df,
                                      type="response"),
                           m3=predict(glow_m3,
                                      newdata=validate_df,
                                      type="response"),
                           m4=predict(glow_m4,
                                      newdata=validate_df,
                                      type="response")),
                      labels=list(validate_df$fracture,
                                  validate_df$fracture,
                                  validate_df$fracture,
                                  validate_df$fracture))

## ----echo=TRUE-----------------------------------------------------------
aucs_obj <- performance(predobj, 'auc')
aucs <- unlist(aucs_obj@y.values)
kable(data.frame(Model=1:4, AUC=aucs), digits=3)

## ----echo=TRUE-----------------------------------------------------------
predobj2 <- prediction(predict(glow_m1, 
                               newdata=test_df, 
                               type="response"),
                       test_df$fracture)
aucs_obj2 <- performance(predobj2, 'auc')
aucs2 <- unlist(aucs_obj2@y.values)
aucs2

## ----echo=TRUE-----------------------------------------------------------
roc_glow <- performance(predobj2, 'tpr', 'fpr')
plot(roc_glow)

