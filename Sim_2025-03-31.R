##Load package
#install.packages('hamlet')
#install.packages('sanon')
library(hamlet)
library(stddiff)
library(lme4)
library(lmerTest)
library(sanon)
library(gtools)
source('/Users/jrigdon/OneDrive - Wake Forest Baptist Health/Useful_Functions/Tables_v2.R')
source('/Users/jrigdon/OneDrive - Wake Forest Baptist Health/Useful_Functions/Figures.R')
source('/Users/jrigdon/OneDrive - Wake Forest Baptist Health/Useful_Functions/Functions.R')


##Function for calculating D
calcD = function(data, var, rand) {
data$var = data[, names(data)==var]
data$rand = data[, names(data)==rand]
s = sd(data$var) #overall SD of weight
m = as.numeric(by(data$var, data$rand, mean)) #record group means
z = (m-mean(m)) / s #record z-scores
D = sqrt(sum(z^2)) #sqrt of sum of z^2s
D
}

##Function to get ANOVA p-value
getP = function(data, var, rand) {
data$var = data[, names(data)==var]
data$rand = data[, names(data)==rand]
m1 = lm(var ~ factor(rand), data=data)
an = anova(m1)
pval = an["factor(rand)", "Pr(>F)"]
pval
}

##Set seed for reproducibility
set.seed(22)

##############################
##Set up a simple simulation##
##############################
##Only 3 batches for now; can extend
##Try different distributions for baseline weight
##1. rnorm(nb*b, 200, 50)
##2. runif(nb*b, 50, 350)
##3. rexp(nb*b, 1/200)


##Try more per group (eliminate hamlet)
##Loop over B batches now
simF2 = function(nsim, nb, B=3, alpha, beta, sd_e=0.05, icc=0, nperm=10000) {
##Solve for sd_b for data generating process
v_e = sd_e^2
v_b = (icc*v_e)/(1-icc)
sd_b = sqrt(v_b)
    
##Matrices in which to record values
##Simple randomization
d0 = matrix(NA, nrow=1, ncol=nsim)
d1 = matrix(NA, nrow=1, ncol=nsim)
p1 = matrix(NA, nrow=1, ncol=nsim)

s_unLM_est = matrix(NA, nrow=1, ncol=nsim)
s_unLM_se = matrix(NA, nrow=1, ncol=nsim)
s_unLM_l = matrix(NA, nrow=1, ncol=nsim)
s_unLM_u = matrix(NA, nrow=1, ncol=nsim)
s_unLM_p = matrix(NA, nrow=1, ncol=nsim)
s_unLM_icc = matrix(NA, nrow=1, ncol=nsim)
    
s_aLM_est = matrix(NA, nrow=1, ncol=nsim)
s_aLM_se = matrix(NA, nrow=1, ncol=nsim)
s_aLM_l = matrix(NA, nrow=1, ncol=nsim)
s_aLM_u = matrix(NA, nrow=1, ncol=nsim)
s_aLM_p = matrix(NA, nrow=1, ncol=nsim)
s_aLM_icc = matrix(NA, nrow=1, ncol=nsim)

##Hamlet
d0h = matrix(NA, nrow=1, ncol=nsim)
d1h = matrix(NA, nrow=1, ncol=nsim)
p1h = matrix(NA, nrow=1, ncol=nsim)

sh_unLM_est = matrix(NA, nrow=1, ncol=nsim)
sh_unLM_se = matrix(NA, nrow=1, ncol=nsim)
sh_unLM_l = matrix(NA, nrow=1, ncol=nsim)
sh_unLM_u = matrix(NA, nrow=1, ncol=nsim)
sh_unLM_p = matrix(NA, nrow=1, ncol=nsim)
sh_unLM_icc = matrix(NA, nrow=1, ncol=nsim)

sh_aLM_est = matrix(NA, nrow=1, ncol=nsim)
sh_aLM_se = matrix(NA, nrow=1, ncol=nsim)
sh_aLM_l = matrix(NA, nrow=1, ncol=nsim)
sh_aLM_u = matrix(NA, nrow=1, ncol=nsim)
sh_aLM_p = matrix(NA, nrow=1, ncol=nsim)
sh_aLM_icc = matrix(NA, nrow=1, ncol=nsim)
    
##Constrained ANOVA P
d0p = matrix(NA, nrow=1, ncol=nsim)
d1p = matrix(NA, nrow=1, ncol=nsim)
p1p = matrix(NA, nrow=1, ncol=nsim)

sp_unLM_est = matrix(NA, nrow=1, ncol=nsim)
sp_unLM_se = matrix(NA, nrow=1, ncol=nsim)
sp_unLM_l = matrix(NA, nrow=1, ncol=nsim)
sp_unLM_u = matrix(NA, nrow=1, ncol=nsim)
sp_unLM_p = matrix(NA, nrow=1, ncol=nsim)
sp_unLM_icc = matrix(NA, nrow=1, ncol=nsim)
    
sp_aLM_est = matrix(NA, nrow=1, ncol=nsim)
sp_aLM_se = matrix(NA, nrow=1, ncol=nsim)
sp_aLM_l = matrix(NA, nrow=1, ncol=nsim)
sp_aLM_u = matrix(NA, nrow=1, ncol=nsim)
sp_aLM_p = matrix(NA, nrow=1, ncol=nsim)
sp_aLM_icc = matrix(NA, nrow=1, ncol=nsim)
    
##Constrained distance
d0d = matrix(NA, nrow=1, ncol=nsim)
d1d = matrix(NA, nrow=1, ncol=nsim)
p1d = matrix(NA, nrow=1, ncol=nsim)

sd_unLM_est = matrix(NA, nrow=1, ncol=nsim)
sd_unLM_se = matrix(NA, nrow=1, ncol=nsim)
sd_unLM_l = matrix(NA, nrow=1, ncol=nsim)
sd_unLM_u = matrix(NA, nrow=1, ncol=nsim)
sd_unLM_p = matrix(NA, nrow=1, ncol=nsim)
sd_unLM_icc = matrix(NA, nrow=1, ncol=nsim)
    
sd_aLM_est = matrix(NA, nrow=1, ncol=nsim)
sd_aLM_se = matrix(NA, nrow=1, ncol=nsim)
sd_aLM_l = matrix(NA, nrow=1, ncol=nsim)
sd_aLM_u = matrix(NA, nrow=1, ncol=nsim)
sd_aLM_p = matrix(NA, nrow=1, ncol=nsim)
sd_aLM_icc = matrix(NA, nrow=1, ncol=nsim)
    
##Make matrix of assignments to sample from
assign = combinations(nb, nb/2)
if (dim(assign)[1]>nperm) {assign = assign[sample(1:dim(assign)[1], nperm), ]}

for (i in 1:nsim) {
##Set up simple simulation with confounding
dta = data.frame(id=1:(nb*B), batch=rep(seq(1:B), each=nb), wt=rnorm(nb*B, 200, 50))

##Set up potential outcomes (with ICC)
dta$y0 = rnorm(nb*B, mean=0.1+alpha*dta$wt, sd=sd_e) + rep(rnorm(B, 0, sd=sd_b), each=nb)
dta$y1 = dta$y0+beta #this is the issue

##Create batch 1 with simple randomization (for all approaches)
##Batch 1
a1 = sample(c(1:dim(assign)[1]), 1) #assignment for batch 1
b1 = dta[dta$batch==1, ] #subset to batch 1
b1$Z = 0 #set simple rand assignments
b1$Z[assign[a1, ]] = 1 #append assignment to batch 1 data frame

b1$Zp = b1$Z
b1$Zd = b1$Z

##Hamlet
b1D = as.matrix(dist(b1$wt, diag=TRUE, upper=TRUE))
sol = match.bb(b1D, g=2)
submatches = paste("Submatch_", LETTERS[1:nb][sol$solution], sep="")
names(submatches) = names(sol$solution)
b1$Submatch = submatches
b1$AllocatedGroups =  match.allocate(b1$Submatch)
b1$Zh = 0
b1$Zh[b1$AllocatedGroups=="Group_B"] = 1
    
all = b1
    
##Then, loop thru batches, randomize, estimate
for (b in 2:B) {
   
##Simple randomization
ab = sample(c(1:dim(assign)[1]), 1) #assignment for batch b
bb = dta[dta$batch==b, ] #subset to batch b
bb$Z = 0 #set simple rand assignments
bb$Z[assign[ab, ]] = 1 #append assignment to batch 1 data frame

bb$Zp = bb$Z
bb$Zd = bb$Z

##Hamlet
bD = as.matrix(dist(bb$wt, diag=TRUE, upper=TRUE))
sol = match.bb(bD, g=2)
submatches = paste("Submatch_", LETTERS[1:nb][sol$solution], sep="")
names(submatches) = names(sol$solution)
bb$Submatch = submatches
bb$AllocatedGroups =  match.allocate(bb$Submatch)
bb$Zh = 0
bb$Zh[bb$AllocatedGroups=="Group_B"] = 1

##Sequential constrained: pick an assignment randomly for batch 1, then constrain for subsequent ones (using ANOVA or D)

#Vector to record P or D for each assignment in current batch
Pb = matrix(NA, nrow=1, ncol=dim(assign)[1]) 
Db = matrix(NA, nrow=1, ncol=dim(assign)[1])
    
for (j in 1:dim(assign)[1]) {
temp = rep(0, nb)
temp[assign[j, ]] = 1

cb = dta[dta$batch==b, ]
cb$Z = temp #placeholder
cb$Zp = temp
cb$Zd = temp

allS = rbind(all[, !names(all) %in% c("Submatch", "AllocatedGroups", "Zh")], cb)

Pb[, j] = getP(allS, "wt", "Zp") #Zp to incorporate previous assignments
Db[, j] = calcD(allS, "wt", "Zd") #Zd to incorporate previous assignments
}

##Assign proportional to ANOVA p-value
asP = sample(c(1:dim(assign)[1]), 1, prob=Pb)
tempP = rep(0, nb)
tempP[assign[asP, ]] = 1

##Assign randomly out of top 25% of assignments in terms of balance metric D
subD = which(Db<=quantile(Db, probs=c(0.25)))
asD = sample(subD, 1)
tempD = rep(0, nb)
tempD[assign[asD, ]] = 1

##Update matrix
bb$Zp = tempP
bb$Zd = tempD

all = rbind(all, bb)
}    
    
##Analytic approaches
##Outcome reveals
##Simple
all$Y = all$y0
all$Y[all$Z==1] = all$y1[all$Z==1]

##Hamlet
all$Yh = all$y0
all$Yh[all$Zh==1] = all$y1[all$Zh==1]
    
##Sequential constrained: ANOVA P
all$Yp = all$y0
all$Yp[all$Zp==1] = all$y1[all$Zp==1]

##Sequential constrained: Distance matrix
all$Yd = all$y0
all$Yd[all$Zd==1] = all$y1[all$Zd==1]

##Balance calculations
##Simple
s1 = stddiff.numeric(data=all, gcol=which(names(all)=="Z"), vcol=which(names(all)=="wt"))
d0[1, i] = s1[colnames(s1)=="stddiff"]
d1[1, i] = calcD(all, "wt", "Z")
p1[1, i] = getP(all, "wt", "Z")

##Hamlet
s1h = stddiff.numeric(data=all, gcol=which(names(all)=="Zh"), vcol=which(names(all)=="wt"))
d0h[1, i] = s1h[colnames(s1h)=="stddiff"]
d1h[1, i] = calcD(all, "wt", "Zh")
p1h[1, i] = getP(all, "wt", "Zh")
    
##Constrained ANOVA p
s1p = stddiff.numeric(data=all, gcol=which(names(all)=="Zp"), vcol=which(names(all)=="wt"))
d0p[1, i] = s1p[colnames(s1p)=="stddiff"]
d1p[1, i] = calcD(all, "wt", "Zp")
p1p[1, i] = getP(all, "wt", "Zp")

##Constrained distance
s1d = stddiff.numeric(data=all, gcol=which(names(all)=="Zd"), vcol=which(names(all)=="wt"))
d0d[1, i] = s1d[colnames(s1d)=="stddiff"]
d1d[1, i] = calcD(all, "wt", "Zd")
p1d[1, i] = getP(all, "wt", "Zd")

##Unadj LMM
##Simple
m1 = lmer(Y ~ Z + (1|batch), data=all)
m1sum = summary(m1, ddf="Kenward-Roger")$coeff #same stderr as above
m1ci = confint(m1, parm="Z") #about the same
m1v = as.data.frame(summary(m1)$varcor)
s_unLM_est[1, i] = m1sum[rownames(m1sum)=="Z", "Estimate"]
s_unLM_se[1, i] = m1sum[rownames(m1sum)=="Z", "Std. Error"]
s_unLM_l[1, i] = m1ci[1]
s_unLM_u[1, i] = m1ci[2]
s_unLM_p[1, i] = m1sum[rownames(m1sum)=="Z", "Pr(>|t|)"]
s_unLM_icc[1, i] = m1v$vcov[m1v$grp=="batch"]/sum(m1v$vcov)

##Hamlet
all$match = paste(all$Submatch, all$batch, sep="_")
m1H = lmer(Yh ~ Zh + (1|match) + (1|batch), data=all)
m1Hsum = summary(m1H, ddf="Kenward-Roger")$coeff #same stderr as above
m1Hci = confint(m1H, parm="Zh") #about the same
m1Hv = as.data.frame(summary(m1H)$varcor)    
sh_unLM_est[1, i] = m1Hsum[rownames(m1Hsum)=="Zh", "Estimate"]
sh_unLM_se[1, i] = m1Hsum[rownames(m1Hsum)=="Zh", "Std. Error"]
sh_unLM_l[1, i] = m1Hci[1]
sh_unLM_u[1, i] = m1Hci[2]
sh_unLM_p[1, i] = m1Hsum[rownames(m1Hsum)=="Zh", "Pr(>|t|)"]
sh_unLM_icc[1, i] = m1Hv$vcov[m1Hv$grp=="batch"]/sum(m1Hv$vcov)
    
##Constrained ANOVA p 
m1P = lmer(Yp ~ Zp + (1|batch), data=all)
m1Psum = summary(m1P, ddf="Kenward-Roger")$coeff #same stderr as above
m1Pci = confint(m1P, parm="Zp") #about the same
m1Pv = as.data.frame(summary(m1P)$varcor)
sp_unLM_est[1, i] = m1Psum[rownames(m1Psum)=="Zp", "Estimate"]
sp_unLM_se[1, i] = m1Psum[rownames(m1Psum)=="Zp", "Std. Error"]
sp_unLM_l[1, i] = m1Pci[1]
sp_unLM_u[1, i] = m1Pci[2]
sp_unLM_p[1, i] = m1Psum[rownames(m1Psum)=="Zp", "Pr(>|t|)"]
sp_unLM_icc[1, i] = m1Pv$vcov[m1Pv$grp=="batch"]/sum(m1Pv$vcov)
    
##Constrained distance
m1D = lmer(Yd ~ Zd + (1|batch), data=all)
m1Dsum = summary(m1D, ddf="Kenward-Roger")$coeff #same stderr as above
m1Dci = confint(m1D, parm="Zd") #about the same
m1Dv = as.data.frame(summary(m1D)$varcor)
sd_unLM_est[1, i] = m1Dsum[rownames(m1Dsum)=="Zd", "Estimate"]
sd_unLM_se[1, i] = m1Dsum[rownames(m1Dsum)=="Zd", "Std. Error"]
sd_unLM_l[1, i] = m1Dci[1]
sd_unLM_u[1, i] = m1Dci[2]
sd_unLM_p[1, i] = m1Dsum[rownames(m1Dsum)=="Zd", "Pr(>|t|)"]
sd_unLM_icc[1, i] = m1Dv$vcov[m1Dv$grp=="batch"]/sum(m1Dv$vcov)
    
##Adj LMM (for weight)
##Simple
m2 = lmer(Y ~ Z + wt + (1|batch), data=all)
m2sum = summary(m2, ddf="Kenward-Roger")$coeff
m2ci = confint(m2, parm="Z")
m2v = as.data.frame(summary(m2)$varcor)
s_aLM_est[1, i] = m2sum[rownames(m2sum)=="Z", "Estimate"]
s_aLM_se[1, i] = m2sum[rownames(m2sum)=="Z", "Std. Error"]
s_aLM_l[1, i] = m2ci[1]
s_aLM_u[1, i] = m2ci[2]
s_aLM_p[1, i] = m2sum[rownames(m2sum)=="Z", "Pr(>|t|)"]
s_aLM_icc[1, i] = m2v$vcov[m2v$grp=="batch"]/sum(m2v$vcov)

##Hamlet
m2H = lmer(Yh ~ Zh + wt + (1|match) + (1|batch), data=all)
m2Hsum = summary(m2H, ddf="Kenward-Roger")$coeff #same stderr as above
m2Hci = confint(m2H, parm="Zh") #about the same
m2Hv = as.data.frame(summary(m2H)$varcor) 
sh_aLM_est[1, i] = m2Hsum[rownames(m2Hsum)=="Zh", "Estimate"]
sh_aLM_se[1, i] = m2Hsum[rownames(m2Hsum)=="Zh", "Std. Error"]
sh_aLM_l[1, i] = m2Hci[1]
sh_aLM_u[1, i] = m2Hci[2]
sh_aLM_p[1, i] = m2Hsum[rownames(m2Hsum)=="Zh", "Pr(>|t|)"]
sh_aLM_icc[1, i] = m2Hv$vcov[m2Hv$grp=="batch"]/sum(m2Hv$vcov)    

##Constrained ANOVA p
m2P = lmer(Yp ~ Zp + wt + (1|batch), data=all)
m2Psum = summary(m2P, ddf="Kenward-Roger")$coeff #same stderr as above
m2Pci = confint(m2P, parm="Zp") #about the same
m2Pv = as.data.frame(summary(m2P)$varcor)
sp_aLM_est[1, i] = m2Psum[rownames(m2Psum)=="Zp", "Estimate"]
sp_aLM_se[1, i] = m2Psum[rownames(m2Psum)=="Zp", "Std. Error"]
sp_aLM_l[1, i] = m2Pci[1]
sp_aLM_u[1, i] = m2Pci[2]
sp_aLM_p[1, i] = m2Psum[rownames(m2Psum)=="Zp", "Pr(>|t|)"]
sp_aLM_icc[1, i] = m2Pv$vcov[m2Pv$grp=="batch"]/sum(m2Pv$vcov)
    
##Constrained distance
m2D = lmer(Yd ~ Zd + wt + (1|batch), data=all)
m2Dsum = summary(m2D, ddf="Kenward-Roger")$coeff #same stderr as above
m2Dci = confint(m2D, parm="Zd") #about the same
m2Dv = as.data.frame(summary(m2D)$varcor)
sd_aLM_est[1, i] = m2Dsum[rownames(m2Dsum)=="Zd", "Estimate"]
sd_aLM_se[1, i] = m2Dsum[rownames(m2Dsum)=="Zd", "Std. Error"]
sd_aLM_l[1, i] = m2Dci[1]
sd_aLM_u[1, i] = m2Dci[2]
sd_aLM_p[1, i] = m2Dsum[rownames(m2Dsum)=="Zd", "Pr(>|t|)"]
sd_aLM_icc[1, i] = m2Dv$vcov[m2Dv$grp=="batch"]/sum(m2Dv$vcov)
}

tab = rbind(
c(round(c(mean(d0), mean(d1), mean(p1<0.05)), 4), rep(NA, 10)), #Simple rand
c(rep(NA, 3), round(c(mean(s_unLM_est), abs(beta-mean(s_unLM_est)), mean(abs(beta-s_unLM_est)), mean(s_unLM_se), sd(s_unLM_est), mean(s_unLM_l<=beta & beta<=s_unLM_u), mean(s_unLM_p<0.05), mean(s_unLM_u<beta), mean(s_unLM_l>beta), mean(s_unLM_icc)), 4)), #Unadj LM
c(rep(NA, 3), round(c(mean(s_aLM_est), abs(beta-mean(s_aLM_est)), mean(abs(beta-s_aLM_est)), mean(s_aLM_se), sd(s_aLM_est), mean(s_aLM_l<=beta & beta<=s_aLM_u), mean(s_aLM_p<0.05), mean(s_aLM_u<beta), mean(s_aLM_l>beta), mean(s_aLM_icc)), 4)), #Adj LM
c(round(c(mean(d0h), mean(d1h), mean(p1h<0.05)), 4), rep(NA, 10)), #Hamlet
c(rep(NA, 3), round(c(mean(sh_unLM_est), abs(beta-mean(sh_unLM_est)), mean(abs(beta-sh_unLM_est)), mean(sh_unLM_se), sd(sh_unLM_est), mean(sh_unLM_l<=beta & beta<=sh_unLM_u), mean(sh_unLM_p<0.05), mean(sh_unLM_u<beta), mean(sh_unLM_l>beta), mean(s_unLM_icc)), 4)), #Unadj LM
c(rep(NA, 3), round(c(mean(sh_aLM_est), abs(beta-mean(sh_aLM_est)), mean(abs(beta-sh_aLM_est)), mean(sh_aLM_se), sd(sh_aLM_est), mean(sh_aLM_l<=beta & beta<=sh_aLM_u), mean(sh_aLM_p<0.05), mean(sh_aLM_u<beta), mean(sh_aLM_l>beta), mean(sh_aLM_icc)), 4)), #Adj LM
c(round(c(mean(d0p), mean(d1p), mean(p1p<0.05)), 4), rep(NA, 10)), #Seq constr, ANOVA p
c(rep(NA, 3), round(c(mean(sp_unLM_est), abs(beta-mean(sp_unLM_est)), mean(abs(beta-sp_unLM_est)), mean(sp_unLM_se), sd(sp_unLM_est), mean(sp_unLM_l<=beta & beta<=sp_unLM_u), mean(sp_unLM_p<0.05), mean(sp_unLM_u<beta), mean(sp_unLM_l>beta), mean(sp_unLM_icc)), 4)), #Unadj LM
c(rep(NA, 3), round(c(mean(sp_aLM_est), abs(beta-mean(sp_aLM_est)), mean(abs(beta-sp_aLM_est)), mean(sp_aLM_se), sd(sp_aLM_est), mean(sp_aLM_l<=beta & beta<=sp_aLM_u), mean(sp_aLM_p<0.05), mean(sp_aLM_u<beta), mean(sp_aLM_l>beta), mean(sp_aLM_icc)), 4)), #Adj LM
c(round(c(mean(d0d), mean(d1d), mean(p1d<0.05)), 4), rep(NA, 10)), #Seq constr, D metric
c(rep(NA, 3), round(c(mean(sd_unLM_est), abs(beta-mean(sd_unLM_est)), mean(abs(beta-sd_unLM_est)), mean(sd_unLM_se), sd(sd_unLM_est), mean(sd_unLM_l<=beta & beta<=sd_unLM_u), mean(sd_unLM_p<0.05), mean(sd_unLM_u<beta), mean(sd_unLM_l>beta), mean(sd_unLM_icc)), 4)), #Unadj LM
c(rep(NA, 3), round(c(mean(sd_aLM_est), abs(beta-mean(sd_aLM_est)), mean(abs(beta-sd_aLM_est)), mean(sd_aLM_se), sd(sd_aLM_est), mean(sd_aLM_l<=beta & beta<=sd_aLM_u), mean(sd_aLM_p<0.05), mean(sd_aLM_u<beta), mean(sd_aLM_l>beta), mean(sd_aLM_icc)), 4)) #Adj LM
)

rownames(tab) = c("Simple rand", "Unadj LM", "Adj LM", "Hamlet", "Unadj LM", "Adj LM", "Seq, ANOVA", "Unadj LM", "Adj LM", "Seq, D", "Unadj LM", "Adj LM")
colnames(tab) = c("Sdiff", "D", "A<.05", "Est", "Bias", "MAD", "ASE", "ESE", "Cov", "Power", "Left", "Right", "ICC")
output.all = list(tab=tab, simple=s_unLM_est, simpleA=s_aLM_est, hamlet=sh_unLM_est, hamletA=sh_aLM_est, anova=sp_unLM_est, anovaA=sp_aLM_est, d=sd_unLM_est, dA=sd_aLM_est, s_unLM_p=s_unLM_p, s_aLM_p=s_aLM_p, sh_unLM_p=sh_unLM_p, sh_aLM_p=sh_aLM_p, sp_unLM_p=sp_unLM_p, sp_aLM_p=sp_aLM_p, sd_unLM_p=sd_unLM_p, sd_aLM_p=sd_aLM_p, d0=d0, d0h=d0h, d0p=d0p, d0d=d0d, s_icc=s_unLM_icc, s_iccA=s_aLM_icc, h_icc=sh_unLM_icc, h_iccA=sh_aLM_icc, p_icc=sp_unLM_icc, p_iccA=sp_aLM_icc, d_icc=sd_unLM_icc, d_iccA=sd_aLM_icc)
return(output.all)
}



##What happens for the "bad" randomizations (imbalance)
##Also look at SEs for each approach? 

##No ICC, equal confounder effect size, increasing sample size
sim1 = simF2(nsim=500, nb=6, B=3, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim1b = simF2(nsim=500, nb=6, B=6, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim1c = simF2(nsim=500, nb=6, B=12, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim1d = simF2(nsim=500, nb=6, B=18, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim1e = simF2(nsim=500, nb=6, B=24, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)

##Increase number per group
sim2 = simF2(nsim=500, nb=12, B=3, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim2b = simF2(nsim=500, nb=12, B=6, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim2c = simF2(nsim=500, nb=12, B=12, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim2d = simF2(nsim=500, nb=12, B=18, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)
sim2e = simF2(nsim=500, nb=12, B=24, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0)


##Save results (3/17)
mkDF = function(obj) {
 obj2 = obj[!names(obj) %in% c("tab", "dP")]
 obj3 = data.frame(t(matrix(unlist(obj2), nrow=length(obj2), byrow=TRUE)))
 names(obj3) = names(obj2)
 obj3
}

summary(mkDF(sim1))
summary(mkDF(sim1b))

write.csv(mkDF(sim1), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim1b), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1b_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim1c), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1c_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim1d), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1d_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim1e), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1e_2025-03-17.csv", row.names=FALSE)

write.csv(mkDF(sim2), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim2b), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2b_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim2c), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2c_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim2d), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2d_2025-03-17.csv", row.names=FALSE)
write.csv(mkDF(sim2e), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2e_2025-03-17.csv", row.names=FALSE)



##Increase ICC to 0.05
sim3 = simF2(nsim=500, nb=6, B=3, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim3b = simF2(nsim=500, nb=6, B=6, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim3c = simF2(nsim=500, nb=6, B=12, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim3d = simF2(nsim=500, nb=6, B=18, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim3e = simF2(nsim=500, nb=6, B=24, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)

##Increase number per group
sim4 = simF2(nsim=500, nb=12, B=3, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim4b = simF2(nsim=500, nb=12, B=6, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim4c = simF2(nsim=500, nb=12, B=12, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim4d = simF2(nsim=500, nb=12, B=18, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)
sim4e = simF2(nsim=500, nb=12, B=24, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.05)


##Save
write.csv(mkDF(sim3), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim3b), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3b_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim3c), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3c_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim3d), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3d_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim3e), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3e_2025-03-18.csv", row.names=FALSE)

write.csv(mkDF(sim4), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim4b), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4b_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim4c), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4c_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim4d), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4d_2025-03-18.csv", row.names=FALSE)
write.csv(mkDF(sim4e), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4e_2025-03-18.csv", row.names=FALSE)


##Increase ICC to 0.1
sim5 = simF2(nsim=500, nb=6, B=3, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim5b = simF2(nsim=500, nb=6, B=6, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim5c = simF2(nsim=500, nb=6, B=12, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim5d = simF2(nsim=500, nb=6, B=18, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim5e = simF2(nsim=500, nb=6, B=24, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)

##Increase number per group
sim6 = simF2(nsim=500, nb=12, B=3, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim6b = simF2(nsim=500, nb=12, B=6, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim6c = simF2(nsim=500, nb=12, B=12, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim6d = simF2(nsim=500, nb=12, B=18, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)
sim6e = simF2(nsim=500, nb=12, B=24, alpha=0.001, beta=-0.05, sd_e=0.05, icc=0.1)


##Save
write.csv(mkDF(sim5), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim5b), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5b_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim5c), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5c_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim5d), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5d_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim5e), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5e_2025-03-19.csv", row.names=FALSE)

write.csv(mkDF(sim6), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim6b), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6b_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim6c), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6c_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim6d), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6d_2025-03-19.csv", row.names=FALSE)
write.csv(mkDF(sim6e), "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6e_2025-03-19.csv", row.names=FALSE)



##Read in simulated data and make smaller figures
s0_n6_B3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1_2025-03-17.csv", header=TRUE)
s0_n6_B6 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1b_2025-03-17.csv", header=TRUE)
s0_n6_B12 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1c_2025-03-17.csv", header=TRUE)

s0_n12_B3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2_2025-03-17.csv", header=TRUE)
s0_n12_B6 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2b_2025-03-17.csv", header=TRUE)
s0_n12_B12 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2c_2025-03-17.csv", header=TRUE)

s05_n6_B3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3_2025-03-18.csv", header=TRUE)
s05_n6_B6 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3b_2025-03-18.csv", header=TRUE)
s05_n6_B12 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3c_2025-03-18.csv", header=TRUE)

s05_n12_B3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4_2025-03-18.csv", header=TRUE)
s05_n12_B6 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4b_2025-03-18.csv", header=TRUE)
s05_n12_B12 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4c_2025-03-18.csv", header=TRUE)


##Function to extract necessary info from saved results
getInfo = function(d0, icc=0, nb=6, B=3) {
 df = data.frame(   
 icc=icc,
 nb=nb, 
 B=B,
 
 meanDiff_s = mean(abs(d0$simple-d0$simpleA)),
 meanDiff_h = mean(abs(d0$hamlet-d0$hamletA)),
 meanDiff_a = mean(abs(d0$anova-d0$anovaA)),
 meanDiff_d = mean(abs(d0$d-d0$dA)),
 
 maxDiff_s = max(abs(d0$simple-d0$simpleA)),
 maxDiff_h = max(abs(d0$hamlet-d0$hamletA)),
 maxDiff_a = max(abs(d0$anova-d0$anovaA)),
 maxDiff_d = max(abs(d0$d-d0$dA)),

 bias_s = mean(d0$simple+0.05),
 bias_sA = mean(d0$simpleA+0.05), 
 bias_h = mean(d0$hamlet+0.05),
 bias_hA = mean(d0$hamletA+0.05),
 bias_a = mean(d0$anova+0.05),
 bias_aA = mean(d0$anovaA+0.05),
 bias_d = mean(d0$d+0.05),
 bias_dA = mean(d0$dA+0.05), 
 
 mad_s = mean(abs(d0$simple+0.05)),
 mad_sA = mean(abs(d0$simpleA+0.05)), 
 mad_h = mean(abs(d0$hamlet+0.05)),
 mad_hA = mean(abs(d0$hamletA+0.05)),
 mad_a = mean(abs(d0$anova+0.05)),
 mad_aA = mean(abs(d0$anovaA+0.05)),
 mad_d = mean(abs(d0$d+0.05)),
 mad_dA = mean(abs(d0$dA+0.05)),

 ese_s = sd(d0$simple),
 ese_sA = sd(d0$simpleA), 
 ese_h = sd(d0$hamlet),
 ese_hA = sd(d0$hamletA),
 ese_a = sd(d0$anova),
 ese_aA = sd(d0$anovaA),
 ese_d = sd(d0$d),
 ese_dA = sd(d0$dA), 

 pow_s = mean(d0$s_unLM_p<0.05),
 pow_sA = mean(d0$s_aLM_p<0.05),
 pow_h = mean(d0$sh_unLM_p<0.05),
 pow_hA = mean(d0$sh_aLM_p<0.05),
 pow_a = mean(d0$sp_unLM_p<0.05),
 pow_aA = mean(d0$sp_aLM_p<0.05),
 pow_d = mean(d0$sd_unLM_p<0.05),
 pow_dA = mean(d0$sd_aLM_p<0.05) 
 )

df 
}

##All data
d1 = rbind(
 getInfo(s0_n6_B3),
 getInfo(s0_n6_B6, B=6),
 getInfo(s0_n6_B12, B=12),
 getInfo(s0_n12_B3, nb=12),
 getInfo(s0_n12_B6, nb=12, B=6),
 getInfo(s0_n12_B12, nb=12, B=12),
 getInfo(s05_n6_B3, icc=0.05),
 getInfo(s05_n6_B6, icc=0.05, B=6),
 getInfo(s05_n6_B12, icc=0.05, B=12),
 getInfo(s05_n12_B3, icc=0.05, nb=12),
 getInfo(s05_n12_B6, icc=0.05, nb=12, B=6),
 getInfo(s05_n12_B12, icc=0.05, nb=12, B=12)
)


##Figure 1: differences (mean, max)
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/unadj_adj_diffs_2025-04-01.pdf")
par(mfrow=c(2, 2))
plot(d1$B[1:3], d1$meanDiff_s[1:3], type="n", ylim=c(0, 0.025), xaxt="n", xlab="Number of batches", ylab="Mean abs diff unadj vs. adj", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
lines(d1$B[1:3], d1$meanDiff_s[1:3], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$meanDiff_h[1:3], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$meanDiff_a[1:3], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$meanDiff_d[1:3], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(d1$B[1:3], d1$meanDiff_s[7:9], type="n", ylim=c(0, 0.025), xaxt="n", xlab="Number of batches", ylab="Mean abs diff unadj vs. adj", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
lines(d1$B[1:3], d1$meanDiff_s[7:9], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$meanDiff_h[7:9], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$meanDiff_a[7:9], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$meanDiff_d[7:9], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(d1$B[1:3], d1$maxDiff_s[1:3], type="n", ylim=c(0, 0.12), xaxt="n", xlab="Number of batches", ylab="Max abs diff unadj vs. adj", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
lines(d1$B[1:3], d1$maxDiff_s[1:3], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$maxDiff_h[1:3], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$maxDiff_a[1:3], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$maxDiff_d[1:3], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(d1$B[1:3], d1$maxDiff_s[7:9], type="n", ylim=c(0, 0.12), xaxt="n", xlab="Number of batches", ylab="Max abs diff unadj vs. adj", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
lines(d1$B[1:3], d1$maxDiff_s[7:9], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$maxDiff_h[7:9], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$maxDiff_a[7:9], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$maxDiff_d[7:9], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)
#dev.off()


##Figure 2: Bias and MAD
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/bias_mad_2025-04-01.pdf")
par(mfrow=c(2, 2))
plot(d1$B[1:3], d1$bias_s[1:3], type="n", ylim=c(-0.01, 0.015), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
abline(h=0, col="lightgrey")
lines(d1$B[1:3], d1$bias_s[1:3], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$bias_sA[1:3], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$bias_h[1:3], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$bias_hA[1:3], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$bias_a[1:3], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$bias_aA[1:3], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$bias_d[1:3], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$bias_dA[1:3], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(d1$B[1:3], d1$bias_s[7:9], type="n", ylim=c(-0.01, 0.015), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
abline(h=0, col="lightgrey")
lines(d1$B[1:3], d1$bias_s[7:9], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$bias_sA[7:9], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$bias_h[7:9], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$bias_hA[7:9], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$bias_a[7:9], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$bias_aA[7:9], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$bias_d[7:9], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$bias_dA[7:9], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(d1$B[1:3], d1$mad_s[1:3], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Mean abs diff", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
abline(h=0, col="lightgrey")
lines(d1$B[1:3], d1$mad_s[1:3], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$mad_sA[1:3], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$mad_h[1:3], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$mad_hA[1:3], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$mad_a[1:3], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$mad_aA[1:3], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$mad_d[1:3], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$mad_dA[1:3], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(d1$B[1:3], d1$mad_s[7:9], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Mean abs diff", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
abline(h=0, col="lightgrey")
lines(d1$B[1:3], d1$mad_s[7:9], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$mad_sA[7:9], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$mad_h[7:9], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$mad_hA[7:9], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$mad_a[7:9], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$mad_aA[7:9], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$mad_d[7:9], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$mad_dA[7:9], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()


##Figure 3: ESE and power
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/ese_power_2025-04-01.pdf")
par(mfrow=c(2, 2))
plot(d1$B[1:3], d1$ese_s[1:3], type="n", ylim=c(0.01, 0.045), xaxt="n", xlab="Number of batches", ylab="Emp std err", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
lines(d1$B[1:3], d1$ese_s[1:3], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$ese_sA[1:3], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$ese_h[1:3], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$ese_hA[1:3], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$ese_a[1:3], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$ese_aA[1:3], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$ese_d[1:3], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$ese_dA[1:3], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(d1$B[1:3], d1$ese_s[7:9], type="n", ylim=c(0.01, 0.045), xaxt="n", xlab="Number of batches", ylab="Emp std err", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
lines(d1$B[1:3], d1$ese_s[7:9], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$ese_sA[7:9], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$ese_h[7:9], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$ese_hA[7:9], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$ese_a[7:9], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$ese_aA[7:9], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$ese_d[7:9], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$ese_dA[7:9], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(d1$B[1:3], d1$pow_s[1:3], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
abline(h=0.8, col="lightgrey")
lines(d1$B[1:3], d1$pow_s[1:3], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$pow_sA[1:3], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$pow_h[1:3], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$pow_hA[1:3], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$pow_a[1:3], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$pow_aA[1:3], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$pow_d[1:3], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$pow_dA[1:3], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(d1$B[1:3], d1$pow_s[7:9], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12), c(3, 6, 12))
abline(h=0.8, col="lightgrey")
lines(d1$B[1:3], d1$pow_s[7:9], type="o", pch=16, lty=1)
lines(d1$B[1:3], d1$pow_sA[7:9], type="o", pch=16, lty=3)
lines(d1$B[1:3], d1$pow_h[7:9], type="o", pch=16, lty=1, col=2)
lines(d1$B[1:3], d1$pow_hA[7:9], type="o", pch=16, lty=3, col=2)
lines(d1$B[1:3], d1$pow_a[7:9], type="o", pch=16, lty=1, col=3)
lines(d1$B[1:3], d1$pow_aA[7:9], type="o", pch=16, lty=3, col=3)
lines(d1$B[1:3], d1$pow_d[7:9], type="o", pch=16, lty=1, col=4)
lines(d1$B[1:3], d1$pow_dA[7:9], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()



##Supplementary material
##Read in extra results
s0_n6_B18 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1d_2025-03-17.csv", header=TRUE)
s0_n6_B24 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim1e_2025-03-17.csv", header=TRUE)

s0_n12_B18 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2d_2025-03-17.csv", header=TRUE)
s0_n12_B24 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim2e_2025-03-17.csv", header=TRUE)

s05_n6_B18 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3d_2025-03-18.csv", header=TRUE)
s05_n6_B24 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim3e_2025-03-18.csv", header=TRUE)

s05_n12_B18 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4d_2025-03-18.csv", header=TRUE)
s05_n12_B24 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim4e_2025-03-18.csv", header=TRUE)

s10_n6_B3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5_2025-03-19.csv", header=TRUE)
s10_n6_B6 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5b_2025-03-19.csv", header=TRUE)
s10_n6_B12 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5c_2025-03-19.csv", header=TRUE)
s10_n6_B18 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5d_2025-03-19.csv", header=TRUE)
s10_n6_B24 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim5e_2025-03-19.csv", header=TRUE)

s10_n12_B3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6_2025-03-19.csv", header=TRUE)
s10_n12_B6 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6b_2025-03-19.csv", header=TRUE)
s10_n12_B12 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6c_2025-03-19.csv", header=TRUE)
s10_n12_B18 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6d_2025-03-19.csv", header=TRUE)
s10_n12_B24 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/sim6e_2025-03-19.csv", header=TRUE)


##Tables of means and max, with bootstrap CI for mean
bootCI = function(x, nsamp=1000) {
 estMean = double()
 estMax = double()
 for (i in 1:nsamp) {   
  x2 = sample(x, length(x), replace=TRUE)
  estMean = c(estMean, mean(x2))
  estMax = c(estMax, max(x2))
 }
 lMean = round(quantile(estMean, 0.025), 5)
 uMean = round(quantile(estMean, 0.975), 5)
 lMax = round(quantile(estMax, 0.025), 5)
 uMax = round(quantile(estMax, 0.975), 5)
 presMean = paste(paste(paste(paste(round(mean(x), 5), " (", sep=""), lMean, sep=""), uMean, sep=", "), ")", sep="")
 presMax = paste(paste(paste(paste(round(max(x), 5), " (", sep=""), lMax, sep=""), uMax, sep=", "), ")", sep="") 
 output.all = list(presMean=presMean, presMax=presMax)
 return(output.all)
}

#bootCI(rnorm(500))

##Do this for each simulation result across the four methods
getEst = function(dta) {
 s = bootCI(abs(dta$simple-dta$simpleA))
 h = bootCI(abs(dta$hamlet-dta$hamletA))
 a = bootCI(abs(dta$anova-dta$anovaA))
 d = bootCI(abs(dta$d-dta$dA))
 me = data.frame(s$presMean, h$presMean, a$presMean, d$presMean)
 names(me) = c("s", "h", "a", "d")
 mx = data.frame(s$presMax, h$presMax, a$presMax, d$presMax)
 names(mx) = c("s", "h", "a", "d") 
 output.all = list(me=me, mx=mx)
 return(output.all)
}
 

##Table of means
tabMean = rbind(getEst(s0_n6_B3)$me,
                getEst(s0_n6_B6)$me,
                getEst(s0_n6_B12)$me,                
                getEst(s0_n6_B18)$me,
                getEst(s0_n6_B24)$me,
                getEst(s0_n12_B3)$me,
                getEst(s0_n12_B6)$me,
                getEst(s0_n12_B12)$me,                
                getEst(s0_n12_B18)$me,
                getEst(s0_n12_B24)$me,
                getEst(s05_n6_B3)$me,
                getEst(s05_n6_B6)$me,
                getEst(s05_n6_B12)$me,                
                getEst(s05_n6_B18)$me,
                getEst(s05_n6_B24)$me,
                getEst(s05_n12_B3)$me,
                getEst(s05_n12_B6)$me,
                getEst(s05_n12_B12)$me,                
                getEst(s05_n12_B18)$me,
                getEst(s05_n12_B24)$me,
                getEst(s10_n6_B3)$me,
                getEst(s10_n6_B6)$me,
                getEst(s10_n6_B12)$me,                
                getEst(s10_n6_B18)$me,
                getEst(s10_n6_B24)$me,
                getEst(s10_n12_B3)$me,
                getEst(s10_n12_B6)$me,
                getEst(s10_n12_B12)$me,                
                getEst(s10_n12_B18)$me,
                getEst(s10_n12_B24)$me)
#word.tab(tabMean, dest="/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/Mean_2025-03-19.docx", help=TRUE)


##Table of max values (ignore CI as doesn't make sense)
tabMax = rbind(getEst(s0_n6_B3)$mx,
               getEst(s0_n6_B6)$mx,
               getEst(s0_n6_B12)$mx,                
               getEst(s0_n6_B18)$mx,
               getEst(s0_n6_B24)$mx,
               getEst(s0_n12_B3)$mx,
               getEst(s0_n12_B6)$mx,
               getEst(s0_n12_B12)$mx,                
               getEst(s0_n12_B18)$mx,
               getEst(s0_n12_B24)$mx,
               getEst(s05_n6_B3)$mx,
               getEst(s05_n6_B6)$mx,
               getEst(s05_n6_B12)$mx,                
               getEst(s05_n6_B18)$mx,
               getEst(s05_n6_B24)$mx,
               getEst(s05_n12_B3)$mx,
               getEst(s05_n12_B6)$mx,
               getEst(s05_n12_B12)$mx,                
               getEst(s05_n12_B18)$mx,
               getEst(s05_n12_B24)$mx,
               getEst(s10_n6_B3)$mx,
               getEst(s10_n6_B6)$mx,
               getEst(s10_n6_B12)$mx,                
               getEst(s10_n6_B18)$mx,
               getEst(s10_n6_B24)$mx,
               getEst(s10_n12_B3)$mx,
               getEst(s10_n12_B6)$mx,
               getEst(s10_n12_B12)$mx,                
               getEst(s10_n12_B18)$mx,
               getEst(s10_n12_B24)$mx)
#word.tab(tabMax, dest="/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/Max_2025-03-19.docx", help=TRUE)


##Full sim results
##Make 30x28 data frame of results
res = rbind(c(0, 6, 3, apply(s0_n6_B3, 2, mean)),
            c(0, 6, 6, apply(s0_n6_B6, 2, mean)),
            c(0, 6, 12, apply(s0_n6_B12, 2, mean)),
            c(0, 6, 18, apply(s0_n6_B18, 2, mean)),
            c(0, 6, 24, apply(s0_n6_B24, 2, mean)),
            c(0, 12, 3, apply(s0_n12_B3, 2, mean)),
            c(0, 12, 6, apply(s0_n12_B6, 2, mean)),
            c(0, 12, 12, apply(s0_n12_B12, 2, mean)),
            c(0, 12, 18, apply(s0_n12_B18, 2, mean)),
            c(0, 12, 24, apply(s0_n12_B24, 2, mean)),
            c(0.05, 6, 3, apply(s05_n6_B3, 2, mean)),
            c(0.05, 6, 6, apply(s05_n6_B6, 2, mean)),
            c(0.05, 6, 12, apply(s05_n6_B12, 2, mean)),
            c(0.05, 6, 18, apply(s05_n6_B18, 2, mean)),
            c(0.05, 6, 24, apply(s05_n6_B24, 2, mean)),
            c(0.05, 12, 3, apply(s05_n12_B3, 2, mean)),
            c(0.05, 12, 6, apply(s05_n12_B6, 2, mean)),
            c(0.05, 12, 12, apply(s05_n12_B12, 2, mean)),
            c(0.05, 12, 18, apply(s05_n12_B18, 2, mean)),
            c(0.05, 12, 24, apply(s05_n12_B24, 2, mean)),
            c(0.1, 6, 3, apply(s10_n6_B3, 2, mean)),
            c(0.1, 6, 6, apply(s10_n6_B6, 2, mean)),
            c(0.1, 6, 12, apply(s10_n6_B12, 2, mean)),
            c(0.1, 6, 18, apply(s10_n6_B18, 2, mean)),
            c(0.1, 6, 24, apply(s10_n6_B24, 2, mean)),
            c(0.1, 12, 3, apply(s10_n12_B3, 2, mean)),
            c(0.1, 12, 6, apply(s10_n12_B6, 2, mean)),
            c(0.1, 12, 12, apply(s10_n12_B12, 2, mean)),
            c(0.1, 12, 18, apply(s10_n12_B18, 2, mean)),
            c(0.1, 12, 24, apply(s10_n12_B24, 2, mean))
            )

res2 = data.frame(res)
names(res2)[1:3] = c("ICC", "nb", "B")

##Append MAD, ASE, ESE, Cov, Left, Right, Power to res2 for all 8 methods
##2, 3, 5, 6, 8, 9, 11, 12
mkDF2 = function(tab) {
    mad = tab[, "MAD"]
    mad2 = mad[!is.na(mad)]
    names(mad2) = c("s_mad", "s_madA", "h_mad", "h_madA", "p_mad", "p_madA", "d_mad", "d_madA")
    ase = tab[, "ASE"]
    ase2 = ase[!is.na(ase)]
    names(ase2) = c("s_ase", "s_aseA", "h_ase", "h_aseA", "p_ase", "p_aseA", "d_ase", "d_aseA")   
    ese = tab[, "ESE"]
    ese2 = ese[!is.na(ese)]
    names(ese2) = c("s_ese", "s_eseA", "h_ese", "h_eseA", "p_ese", "p_eseA", "d_ese", "d_eseA")
    cov = tab[, "Cov"]
    cov2 = cov[!is.na(cov)]
    names(cov2) = c("s_cov", "s_covA", "h_cov", "h_covA", "p_cov", "p_covA", "d_cov", "d_covA")    
    left = tab[, "Left"]
    left2 = left[!is.na(left)]
    names(left2) = c("s_left", "s_leftA", "h_left", "h_leftA", "p_left", "p_leftA", "d_left", "d_leftA")    
    right = tab[, "Right"]
    right2 = right[!is.na(right)]
    names(right2) = c("s_right", "s_rightA", "h_right", "h_rightA", "p_right", "p_rightA", "d_right", "d_rightA") 
    power = tab[, "Power"]
    power2 = power[!is.na(power)]
    names(power2) = c("s_power", "s_powerA", "h_power", "h_powerA", "p_power", "p_powerA", "d_power", "d_powerA")
    c(mad2, ase2, ese2, cov2, left2, right2, power2)
}

temp = rbind(mkDF2(sim1$tab),
             mkDF2(sim1b$tab),
             mkDF2(sim1c$tab),
             mkDF2(sim1d$tab),
             mkDF2(sim1e$tab),
             mkDF2(sim2$tab),
             mkDF2(sim2b$tab),
             mkDF2(sim2c$tab),
             mkDF2(sim2d$tab),
             mkDF2(sim2e$tab),
             mkDF2(sim3$tab),
             mkDF2(sim3b$tab),
             mkDF2(sim3c$tab),
             mkDF2(sim3d$tab),
             mkDF2(sim3e$tab),
             mkDF2(sim4$tab),
             mkDF2(sim4b$tab),
             mkDF2(sim4c$tab),
             mkDF2(sim4d$tab),
             mkDF2(sim4e$tab),            
             mkDF2(sim5$tab),
             mkDF2(sim5b$tab),
             mkDF2(sim5c$tab),
             mkDF2(sim5d$tab),
             mkDF2(sim5e$tab),
             mkDF2(sim6$tab),
             mkDF2(sim6b$tab),
             mkDF2(sim6c$tab),
             mkDF2(sim6d$tab),
             mkDF2(sim6e$tab))             

res3 = cbind(res2, temp)

#write.csv(res3, "/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/Results_2025-03-19.csv", row.names=FALSE)
res3 = read.csv("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Data/Results_2025-03-19.csv", header=TRUE)

##Add balance metrics to res3 (by hand from word doc)
res3$s_sdiff = c(0.3993, 0.2694, 0.1843, 0.1527, 0.1283, 0.2754, 0.1855, 0.1214, 0.1048, 0.0959, 0.3971, 0.2582, 0.1854, 0.1551, 0.1327, 0.2695, 0.2020, 0.1379, 0.1138, 0.0933, 0.3900, 0.2647, 0.1970, 0.1497, 0.1286, 0.2858, 0.1962, 0.1361, 0.1106, 0.0975)
res3$s_dm = c(0.2730, 0.1876, 0.1295, 0.1074, 0.0904, 0.1915, 0.1303, 0.0857, 0.0740, 0.0677, 0.2716, 0.1796, 0.1301, 0.1091, 0.0935, 0.1877, 0.1416, 0.0970, 0.0803, 0.0659, 0.2668, 0.1842, 0.1382, 0.1054, 0.0906, 0.1988, 0.1375, 0.0959, 0.0780, 0.0688)
res3$s_ap = c(0.046, 0.048, 0.042, 0.058, 0.044, 0.046, 0.05, 0.036, 0.048, 0.054, 0.040, 0.05, 0.056, 0.056, 0.052, 0.054, 0.068, 0.056, 0.05, 0.038, 0.048, 0.048, 0.044, 0.042, 0.042, 0.06, 0.064, 0.054, 0.046, 0.058)

res3$h_sdiff = c(0.2024, 0.1389, 0.0965, 0.0784, 0.0679, 0.0825, 0.0616, 0.0438, 0.0353, 0.0323, 0.1988, 0.1419, 0.0995, 0.0815, 0.0689, 0.0912, 0.0639, 0.0450, 0.0354, 0.0330, 0.1918, 0.1342, 0.1017, 0.0806, 0.0709, 0.0905, 0.0649, 0.0437, 0.0384, 0.0336)
res3$h_dm = c(0.1452, 0.0989, 0.0684, 0.0555, 0.0481, 0.0590, 0.0438, 0.0311, 0.0250, 0.0229, 0.1425, 0.1009, 0.0706, 0.0577, 0.0488, 0.0652, 0.0454, 0.0319, 0.0251, 0.0234, 0.1377, 0.0956, 0.0722, 0.0571, 0.0502, 0.0647, 0.0461, 0.0310, 0.0272, 0.0238)
res3$h_ap = c(0, 0, 0, 0, 0.002, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.002, 0, 0, 0, 0, 0)

res3$a_sdiff = c(0.2121, 0.1385, 0.0902, 0.0701, 0.0609, 0.1330, 0.0940, 0.0640, 0.0493, 0.0430, 0.2253, 0.1342, 0.0899, 0.0680, 0.0598, 0.1478, 0.0970, 0.0628, 0.0506, 0.0426, 0.2056, 0.1317, 0.0887, 0.0722, 0.0620, 0.1456, 0.0876, 0.0633, 0.0453, 0.0433)
res3$a_dm = c(0.1512, 0.0985, 0.0639, 0.0497, 0.0431, 0.0945, 0.0667, 0.0453, 0.0349, 0.0305, 0.1599, 0.0955, 0.0638, 0.0482, 0.0424, 0.1049, 0.0688, 0.0445, 0.0358, 0.0301, 0.1467, 0.0936, 0.0629, 0.0511, 0.0439, 0.1035, 0.0622, 0.0448, 0.0321, 0.0307)
res3$a_ap = c(0, 0.002, 0.002, 0, 0, 0.004, 0, 0.002, 0, 0, 0.002, 0, 0, 0, 0, 0.002, 0, 0, 0, 0, 0, 0, 0.002, 0, 0, 0, 0, 0, 0, 0)

res3$d_sdiff = c(0.0559, 0.0276, 0.0139, 0.0087, 0.0075, 0.0354, 0.0164, 0.0075, 0.0057, 0.0041, 0.0580, 0.0274, 0.0138, 0.0094, 0.0070, 0.0363, 0.0165, 0.0081, 0.0056, 0.0042, 0.0607, 0.0281, 0.0136, 0.0097, 0.0067, 0.0363, 0.0163, 0.0083, 0.0057, 0.0042)
res3$d_dm = c(0.0407, 0.0198, 0.0099, 0.0062, 0.0053, 0.0254, 0.0117, 0.0054, 0.0040, 0.0029, 0.0421, 0.0196, 0.0098, 0.0067, 0.0050, 0.0260, 0.0118, 0.0057, 0.0040, 0.0030, 0.0441, 0.0202, 0.0097, 0.0069, 0.0047, 0.0260, 0.0116, 0.0059, 0.0041, 0.0030)
res3$d_ap = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


##Supplementary figures

##Bias
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/Bias_2025-03-19.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$simple[1:5]+0.05, type="n", ylim=c(-0.005, 0.005), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$simple[1:5]+0.05, type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$simpleA[1:5]+0.05, type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$hamlet[1:5]+0.05, type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$hamletA[1:5]+0.05, type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$anova[1:5]+0.05, type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$anovaA[1:5]+0.05, type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d[1:5]+0.05, type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$dA[1:5]+0.05, type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$simple[11:15]+0.05, type="n", ylim=c(-0.005, 0.005), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$simple[11:15]+0.05, type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$simpleA[11:15]+0.05, type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$hamlet[11:15]+0.05, type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$hamletA[11:15]+0.05, type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$anova[11:15]+0.05, type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$anovaA[11:15]+0.05, type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d[11:15]+0.05, type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$dA[11:15]+0.05, type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$simple[21:25]+0.05, type="n", ylim=c(-0.005, 0.005), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$simple[21:25]+0.05, type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$simpleA[21:25]+0.05, type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$hamlet[21:25]+0.05, type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$hamletA[21:25]+0.05, type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$anova[21:25]+0.05, type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$anovaA[21:25]+0.05, type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d[21:25]+0.05, type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$dA[21:25]+0.05, type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)


plot(res3$B[1:5], res3$simple[6:10]+0.05, type="n", ylim=c(-0.005, 0.005), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$simple[6:10]+0.05, type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$simpleA[6:10]+0.05, type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$hamlet[6:10]+0.05, type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$hamletA[6:10]+0.05, type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$anova[6:10]+0.05, type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$anovaA[6:10]+0.05, type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d[6:10]+0.05, type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$dA[6:10]+0.05, type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$simple[16:20]+0.05, type="n", ylim=c(-0.005, 0.005), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$simple[16:20]+0.05, type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$simpleA[16:20]+0.05, type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$hamlet[16:20]+0.05, type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$hamletA[16:20]+0.05, type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$anova[16:20]+0.05, type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$anovaA[16:20]+0.05, type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d[16:20]+0.05, type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$dA[16:20]+0.05, type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$simple[26:30]+0.05, type="n", ylim=c(-0.005, 0.005), xaxt="n", xlab="Number of batches", ylab="Bias", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$simple[26:30]+0.05, type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$simpleA[26:30]+0.05, type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$hamlet[26:30]+0.05, type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$hamletA[26:30]+0.05, type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$anova[26:30]+0.05, type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$anovaA[26:30]+0.05, type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d[26:30]+0.05, type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$dA[26:30]+0.05, type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()


##Mean absolute difference
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/MAD_2025-03-20.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_mad[1:5], type="n", ylim=c(0, 0.03), xaxt="n", xlab="Number of batches", ylab="Mean Abs Dev", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$s_mad[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_madA[1:5], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_mad[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_madA[1:5], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_mad[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_madA[1:5], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_mad[1:5], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_madA[1:5], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_mad[11:15], type="n", ylim=c(0, 0.03), xaxt="n", xlab="Number of batches", ylab="Mean Abs Dev", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$s_mad[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_madA[11:15], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_mad[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_madA[11:15], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_mad[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_madA[11:15], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_mad[11:15], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_madA[11:15], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_mad[21:25], type="n", ylim=c(0, 0.03), xaxt="n", xlab="Number of batches", ylab="Mean Abs Dev", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$s_mad[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_madA[21:25], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_mad[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_madA[21:25], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_mad[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_madA[21:25], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_mad[21:25], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_madA[21:25], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)


plot(res3$B[1:5], res3$s_mad[6:10], type="n", ylim=c(0, 0.03), xaxt="n", xlab="Number of batches", ylab="Mean Abs Dev", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$s_mad[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_madA[6:10], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_mad[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_madA[6:10], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_mad[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_madA[6:10], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_mad[6:10], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_madA[6:10], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_mad[16:20], type="n", ylim=c(0, 0.03), xaxt="n", xlab="Number of batches", ylab="Mean Abs Dev", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$s_mad[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_madA[16:20], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_mad[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_madA[16:20], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_mad[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_madA[16:20], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_mad[16:20], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_madA[16:20], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_mad[26:30], type="n", ylim=c(0, 0.03), xaxt="n", xlab="Number of batches", ylab="Mean Abs Dev", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
abline(h=0, col="lightgrey")
lines(res3$B[1:5], res3$s_mad[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_madA[26:30], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_mad[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_madA[26:30], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_mad[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_madA[26:30], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_mad[26:30], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_madA[26:30], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()


##Average standard error
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/AvgStdErr_2025-03-19.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_ase[1:5], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Average standard error", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ase[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_aseA[1:5], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ase[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_aseA[1:5], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ase[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_aseA[1:5], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ase[1:5], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_aseA[1:5], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ase[11:15], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Average standard error", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ase[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_aseA[11:15], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ase[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_aseA[11:15], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ase[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_aseA[11:15], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ase[11:15], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_aseA[11:15], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ase[21:25], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Average standard error", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ase[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_aseA[21:25], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ase[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_aseA[21:25], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ase[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_aseA[21:25], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ase[21:25], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_aseA[21:25], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)


plot(res3$B[1:5], res3$s_ase[6:10], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Average standard error", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ase[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_aseA[6:10], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ase[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_aseA[6:10], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ase[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_aseA[6:10], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ase[6:10], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_aseA[6:10], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ase[16:20], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Average standard error", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ase[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_aseA[16:20], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ase[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_aseA[16:20], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ase[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_aseA[16:20], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ase[16:20], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_aseA[16:20], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ase[26:30], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Average standard error", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ase[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_aseA[26:30], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ase[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_aseA[26:30], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ase[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_aseA[26:30], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ase[26:30], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_aseA[26:30], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()



##Empirical standard error
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/EmpStdErr_2025-03-20.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_ese[1:5], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Empirical standard error", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ese[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_eseA[1:5], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ese[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_eseA[1:5], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ese[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_eseA[1:5], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ese[1:5], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_eseA[1:5], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ese[11:15], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Empirical standard error", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ese[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_eseA[11:15], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ese[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_eseA[11:15], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ese[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_eseA[11:15], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ese[11:15], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_eseA[11:15], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ese[21:25], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Empirical standard error", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ese[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_eseA[21:25], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ese[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_eseA[21:25], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ese[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_eseA[21:25], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ese[21:25], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_eseA[21:25], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)


plot(res3$B[1:5], res3$s_ese[6:10], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Empirical standard error", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ese[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_eseA[6:10], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ese[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_eseA[6:10], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ese[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_eseA[6:10], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ese[6:10], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_eseA[6:10], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ese[16:20], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Empirical standard error", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ese[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_eseA[16:20], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ese[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_eseA[16:20], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ese[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_eseA[16:20], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ese[16:20], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_eseA[16:20], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_ese[26:30], type="n", ylim=c(0, 0.035), xaxt="n", xlab="Number of batches", ylab="Empirical standard error", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ese[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_eseA[26:30], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_ese[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_eseA[26:30], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_ese[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_eseA[26:30], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_ese[26:30], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_eseA[26:30], type="o", pch=16, lty=3, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()





##Power
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/Power_2025-03-19.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_power[1:5], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_power[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_powerA[1:5], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_power[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_powerA[1:5], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_power[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_powerA[1:5], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_power[1:5], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_powerA[1:5], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_power[11:15], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_power[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_powerA[11:15], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_power[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_powerA[11:15], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_power[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_powerA[11:15], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_power[11:15], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_powerA[11:15], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_power[21:25], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_power[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_powerA[21:25], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_power[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_powerA[21:25], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_power[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_powerA[21:25], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_power[21:25], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_powerA[21:25], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)


plot(res3$B[1:5], res3$s_power[6:10], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_power[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_powerA[6:10], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_power[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_powerA[6:10], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_power[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_powerA[6:10], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_power[6:10], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_powerA[6:10], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_power[16:20], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_power[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_powerA[16:20], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_power[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_powerA[16:20], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_power[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_powerA[16:20], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_power[16:20], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_powerA[16:20], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_power[26:30], type="n", ylim=c(0, 1), xaxt="n", xlab="Number of batches", ylab="Power", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_power[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_powerA[26:30], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_power[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_powerA[26:30], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_power[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_powerA[26:30], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_power[26:30], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_powerA[26:30], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()



##Coverage
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/Coverage_2025-03-19.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_cov[1:5], type="n", ylim=c(0.8, 1), xaxt="n", xlab="Number of batches", ylab="Coverage", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
abline(h=0.95, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_cov[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_covA[1:5], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_cov[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_covA[1:5], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_cov[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_covA[1:5], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_cov[1:5], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_covA[1:5], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_cov[11:15], type="n", ylim=c(0.8, 1), xaxt="n", xlab="Number of batches", ylab="Coverage", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
abline(h=0.95, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_cov[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_covA[11:15], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_cov[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_covA[11:15], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_cov[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_covA[11:15], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_cov[11:15], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_covA[11:15], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_cov[21:25], type="n", ylim=c(0.8, 1), xaxt="n", xlab="Number of batches", ylab="Coverage", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
abline(h=0.95, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_cov[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_covA[21:25], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_cov[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_covA[21:25], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_cov[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_covA[21:25], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_cov[21:25], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_covA[21:25], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)


plot(res3$B[1:5], res3$s_cov[6:10], type="n", ylim=c(0.8, 1), xaxt="n", xlab="Number of batches", ylab="Coverage", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
abline(h=0.95, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_cov[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_covA[6:10], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_cov[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_covA[6:10], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_cov[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_covA[6:10], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_cov[6:10], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_covA[6:10], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_cov[16:20], type="n", ylim=c(0.8, 1), xaxt="n", xlab="Number of batches", ylab="Coverage", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
abline(h=0.95, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_cov[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_covA[16:20], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_cov[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_covA[16:20], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_cov[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_covA[16:20], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_cov[16:20], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_covA[16:20], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)

plot(res3$B[1:5], res3$s_cov[26:30], type="n", ylim=c(0.8, 1), xaxt="n", xlab="Number of batches", ylab="Coverage", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
abline(h=0.95, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_cov[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$s_covA[26:30], type="o", pch=16, lty=3)
lines(res3$B[1:5], res3$h_cov[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$h_covA[26:30], type="o", pch=16, lty=3, col=2)
lines(res3$B[1:5], res3$p_cov[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$p_covA[26:30], type="o", pch=16, lty=3, col=3)
lines(res3$B[1:5], res3$d_cov[26:30], type="o", pch=16, lty=1, col=4)
lines(res3$B[1:5], res3$d_covA[26:30], type="o", pch=16, lty=3, col=4)
legend('bottomright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D", "", "Unadjusted", "Adjusted"), lty=c(1, 1, 1, 1, 0, 1, 3), col=c(1, 2, 3, 4, 1, 1, 1), cex=0.7)
#dev.off()


##Balance metrics
##Standardized difference
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/StdDiff_2025-03-21.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_sdiff[1:5], type="n", ylim=c(0, 0.42), xaxt="n", xlab="Number of batches", ylab="Standardized Difference", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
abline(h=0.2, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_sdiff[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_sdiff[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_sdiff[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_sdiff[1:5], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_sdiff[11:15], type="n", ylim=c(0, 0.42), xaxt="n", xlab="Number of batches", ylab="Standardized Difference", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
abline(h=0.2, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_sdiff[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_sdiff[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_sdiff[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_sdiff[11:15], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_sdiff[21:25], type="n", ylim=c(0, 0.42), xaxt="n", xlab="Number of batches", ylab="Standardized Difference", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
abline(h=0.2, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_sdiff[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_sdiff[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_sdiff[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_sdiff[21:25], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)


plot(res3$B[1:5], res3$s_sdiff[6:10], type="n", ylim=c(0, 0.42), xaxt="n", xlab="Number of batches", ylab="Standardized Difference", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
abline(h=0.2, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_sdiff[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_sdiff[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_sdiff[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_sdiff[6:10], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_sdiff[16:20], type="n", ylim=c(0, 0.42), xaxt="n", xlab="Number of batches", ylab="Standardized Difference", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
abline(h=0.2, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_sdiff[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_sdiff[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_sdiff[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_sdiff[16:20], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_sdiff[26:30], type="n", ylim=c(0, 0.42), xaxt="n", xlab="Number of batches", ylab="Standardized Difference", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
abline(h=0.2, col="lightgrey")
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_sdiff[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_sdiff[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_sdiff[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_sdiff[26:30], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)
#dev.off()


##D-metric
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/Dmetric_2025-03-21.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_dm[1:5], type="n", ylim=c(0, 0.3), xaxt="n", xlab="Number of batches", ylab="D-metric", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_dm[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_dm[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_dm[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_dm[1:5], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_dm[11:15], type="n", ylim=c(0, 0.3), xaxt="n", xlab="Number of batches", ylab="D-metric", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_dm[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_dm[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_dm[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_dm[11:15], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_dm[21:25], type="n", ylim=c(0, 0.3), xaxt="n", xlab="Number of batches", ylab="D-metric", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_dm[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_dm[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_dm[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_dm[21:25], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)


plot(res3$B[1:5], res3$s_dm[6:10], type="n", ylim=c(0, 0.3), xaxt="n", xlab="Number of batches", ylab="D-metric", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_dm[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_dm[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_dm[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_dm[6:10], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_dm[16:20], type="n", ylim=c(0, 0.3), xaxt="n", xlab="Number of batches", ylab="D-metric", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_dm[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_dm[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_dm[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_dm[16:20], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_dm[26:30], type="n", ylim=c(0, 0.3), xaxt="n", xlab="Number of batches", ylab="D-metric", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_dm[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_dm[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_dm[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_dm[26:30], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)
#dev.off()


##ANOVA p-value
#pdf("/Users/jrigdon/OneDrive - Wake Forest Baptist Health/MoTrPAC/Randomization/Output/ANOVAp_2025-03-21.pdf")
par(mfrow=c(2, 3))
plot(res3$B[1:5], res3$s_ap[1:5], type="n", ylim=c(0, 0.08), xaxt="n", xlab="Number of batches", ylab="ANOVA p<0.05", pch=16, main=expression('ICC=0, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ap[1:5], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_ap[1:5], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_ap[1:5], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_ap[1:5], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_ap[11:15], type="n", ylim=c(0, 0.08), xaxt="n", xlab="Number of batches", ylab="ANOVA p<0.05", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ap[11:15], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_ap[11:15], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_ap[11:15], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_ap[11:15], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_ap[21:25], type="n", ylim=c(0, 0.08), xaxt="n", xlab="Number of batches", ylab="ANOVA p<0.05", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=6'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ap[21:25], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_ap[21:25], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_ap[21:25], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_ap[21:25], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)


plot(res3$B[1:5], res3$s_ap[6:10], type="n", ylim=c(0, 0.08), xaxt="n", xlab="Number of batches", ylab="ANOVA p<0.05", pch=16, main=expression('ICC=0, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ap[6:10], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_ap[6:10], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_ap[6:10], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_ap[6:10], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_ap[16:20], type="n", ylim=c(0, 0.08), xaxt="n", xlab="Number of batches", ylab="ANOVA p<0.05", pch=16, main=expression('ICC=0.05, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ap[16:20], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_ap[16:20], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_ap[16:20], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_ap[16:20], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)

plot(res3$B[1:5], res3$s_ap[26:30], type="n", ylim=c(0, 0.08), xaxt="n", xlab="Number of batches", ylab="ANOVA p<0.05", pch=16, main=expression('ICC=0.1, '*'n'[b]*'=12'))
axis(1, at=c(3, 6, 12, 18, 24), c(3, 6, 12, 18, 24))
lines(res3$B[1:5], res3$s_ap[26:30], type="o", pch=16, lty=1)
lines(res3$B[1:5], res3$h_ap[26:30], type="o", pch=16, lty=1, col=2)
lines(res3$B[1:5], res3$a_ap[26:30], type="o", pch=16, lty=1, col=3)
lines(res3$B[1:5], res3$d_ap[26:30], type="o", pch=16, lty=1, col=4)
legend('topright', inset=0.05, bty="n", c("Simple", "Hamlet", "Seq-ANOVA", "Seq-D"), lty=c(1, 1, 1, 1), col=c(1, 2, 3, 4), cex=0.7)
#dev.off()


