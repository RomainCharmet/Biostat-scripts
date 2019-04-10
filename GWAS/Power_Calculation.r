set.seed(1)

repetitions = 10000
N = 1370
a = read.table("Score_table.txt",header=TRUE)
SCORE = a$SCORE
AGE = a$AGE
GENDER = a$GENDER
PC1 = a$PC1
PC2 = a$PC2
PC3 = a$PC3
PC4 = a$PC4
BETA = c(0.054925,-0.221521,0.071366,3.304301,5.296428,-10.082558,-15.927069)
INTERCEP = -8.301458
Logit = INTERCEP+BETA[1]*SCORE+BETA[2]*GENDER+BETA[3]*AGE+BETA[4]*PC1+BETA[5]*PC2+BETA[6]*PC3+BETA[7]*PC4
Prob = (1/(1+exp(-Logit)))

significant = matrix(nrow=repetitions, ncol=3)

startT = proc.time()[3]
for(i in 1:repetitions){
  responses          = rbinom(n=N, size=1, prob=Prob)
  model              = glm(responses~SCORE, 
                           family=binomial(link="logit"))
  significant[i,1] = (summary(model)$coefficients[2,4]<.05)
  significant[i,2]   = sum(significant[i,1])
  modelDev           = model$null.deviance-model$deviance
  significant[i,3]   = (1-pchisq(modelDev, 5))<.05
}
endT = proc.time()[3]
endT-startT

sum(significant[,1])/repetitions      # pre-specified effect power for var1
[1] 0.042
sum(significant[,2])/repetitions      # pre-specified effect power for var2
[1] 0.017
sum(significant[,3])/repetitions      # pre-specified effect power for var12
[1] 0.035
sum(significant[,4])/repetitions      # pre-specified effect power for var1X2
[1] 0.019
sum(significant[,5])/repetitions      # pre-specified effect power for var12X2
[1] 0.022
sum(significant[,7])/repetitions      # power for likelihood ratio test of model
[1] 0.168
sum(significant[,6]==5)/repetitions   # all effects power
[1] 0.001
sum(significant[,6]>0)/repetitions    # any effect power
[1] 0.065
sum(significant[,4]&significant[,5])/repetitions
