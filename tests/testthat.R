# initialize all
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(tictoc)
Rcpp::sourceCpp('src/coxph.cpp')
Rcpp::sourceCpp('src/linear.cpp')
# Rcpp::sourceCpp('src/linearMEXY.cpp')
# Rcpp::sourceCpp('src/lmm.cpp')
Rcpp::sourceCpp('src/logistic.cpp')
# Rcpp::sourceCpp('src/utility.cpp')
source('R/smle.R')
source('R/smle_lmm.R')
source('R/smle_MEXY.R')

# library(TwoPhaseReg)
n = 2000
n2 = 600
true_beta = 0.3
true_gamma = 0.4
true_eta = 0.5
seed = 12345
r = 0.3
N_SIEVE = 8

#### Sieve with histogram bases
set.seed(12345)
U2 = runif(n)
simX = runif(n)
simZ = r*simX+U2
simY = true_beta*simX+true_gamma*simZ+rnorm(n)
order.simY = order(simY)
phase2.id = c(order.simY[1:(n2/2)], order.simY[(n-(n2/2)+1):n])
Bspline_Z = matrix(NA, nrow=n, ncol=N_SIEVE)
cut_z = cut(simZ, breaks=quantile(simZ, probs=seq(0, 1, 1/N_SIEVE)), include.lowest = TRUE)
for (i in 1:N_SIEVE) {
    Bspline_Z[,i] = as.numeric(cut_z == names(table(cut_z))[i])
}
colnames(Bspline_Z) = paste("bs", 1:N_SIEVE, sep="")
dat = data.frame(Y=simY, X=simX, Z=simZ, Bspline_Z)
dat[-phase2.id,"X"] = NA
#
# res = smle(Y="Y", X="X", Z="Z", Bspline_Z=colnames(Bspline_Z), data=dat)
# res

#### Sieve with linear bases
library(splines)
set.seed(12345)
U1 = runif(n)
U2 = runif(n)
simX = runif(n)
simZ_1 = r*simX+U1
simZ_2 = r*simX+U2
simY = true_beta*simX+true_gamma*simZ_1+true_eta*simZ_2+rnorm(n)
order.simY = order(simY)
phase2.id = c(order.simY[1:(n2/2)], order.simY[(n-(n2/2)+1):n])
bs1 = bs(simZ_1, df=N_SIEVE, degree=1, Boundary.knots=range(simZ_1), intercept=TRUE)
bs2 = bs(simZ_2, df=N_SIEVE, degree=1, Boundary.knots=range(simZ_2), intercept=TRUE)
Bspline_Z = matrix(NA, ncol=N_SIEVE^2, nrow=n)
for (i in 1:ncol(bs1)) {
    for (j in 1:ncol(bs2)) {
        idx = i*N_SIEVE+j-N_SIEVE
        Bspline_Z[,idx] = bs1[,i]*bs2[,j]
    }
}
colnames(Bspline_Z) = paste("bs", 1:ncol(Bspline_Z), sep="")
dat = data.frame(Y=simY, X=simX, Z1=simZ_1, Z2=simZ_2, Bspline_Z)
dat[-phase2.id,"X"] = NA
dat["Delta"] = c(rep(1, 546), rep(0, 2000-546))

tic("smle")
# TwoPhase_GeneralSpline
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), data=dat)
# 177.35 sec elapsed
# $coefficients
#            Estimate         SE Statistic      p-value
# Intercept 0.0538190 0.07287324 0.7385291 4.601930e-01
# X         0.3367835 0.09518781 3.5380951 4.030249e-04
# Z1        0.3808042 0.07657005 4.9732790 6.582983e-07
# Z2        0.4135827 0.07696735 5.3734818 7.723062e-08
# Aug 24 : 133.49 sec elapsed
# Sep 13 : 125.75 sec elapsed



# TwoPhase_MLE0
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=NULL, data=dat)
# $coefficients
# Estimate         SE  Statistic      p-value
# Intercept -0.06193395 0.08597975 -0.7203319 4.713207e-01
# X          0.46180716 0.13255814  3.4838084 4.943332e-04
# Z1         0.42190671 0.07844466  5.3783993 7.515098e-08
# Z2         0.44886314 0.07999735  5.6109751 2.011897e-08
#  1724.11 sec elapsed



# TwoPhase_GeneralSpline_logistic
# dat["Y"] = sample(c(0,1), nrow(dat), TRUE)
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), data=dat, model="logistic")
# $coefficients
#              Estimate        SE  Statistic   p-value
# Intercept -0.03133511 0.1527562 -0.2051314 0.8374694
# X          0.35195573 0.3071224  1.1459786 0.2518040
# Z1        -0.09435242 0.1658927 -0.5687557 0.5695219
# Z2        -0.10004598 0.1647833 -0.6071365 0.5437603
# Sep 2: 148.69 sec elapsed
# Sep 8: 116.66 sec elapsed


# TwoPhase_MLE0_logistic
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=NULL, data=dat, model="logistic")
# $coefficients
#              Estimate        SE   Statistic   p-value
# Intercept -0.01363803 0.1830567 -0.07450169 0.9406112
# X          0.13783433 0.2940132  0.46880314 0.6392104
# Z1        -0.02346967 0.1507358 -0.15570075 0.8762689jk
# Z2        -0.03550027 0.1522553 -0.23316271 0.8156351
# 2539.78 sec elapsed

# 100 iterations - abbr
#              Estimate        SE   Statistic   p-value
# Intercept -0.01363803 0.1535751 -0.08880363 0.9292380
# X          0.13783433 0.1623791  0.84884282 0.3959688
# Z1        -0.02346967 0.1495529 -0.15693225 0.8752982
# Z2        -0.03550027 0.1510185 -0.23507233 0.8141526
# 139.8 sec

res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), Delta="Y", data=dat, model="coxph")



res
toc()
# add log likelihood
# add aic calculation

