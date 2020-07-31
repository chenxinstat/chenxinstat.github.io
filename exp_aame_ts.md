Here is an example implementations of Auto Adaptive M-Estimation (AAME) for in R, finding the mean and the autoregressive coefficients of a data, and the innovation density, with a heavy-tailed and autocorrelated error.

The example will have AR(3) error with order p = 3.

1. Run this in R: [AAME_AR functions](functions_osar.r)

2. Load necessary library:
```markdown
library(quadprog)
```

3. Generate a data and fit the model:
```markdown
set.seed(1)
n <- 500
mu <- 100
phi.true <- c(0.8,-0.6,0.4)
p <- length(phi.true)

error <- rep(0,n+p)
innovation <- rep(NA,n)
for (i in 1:n) {
  temp <- 0
  for (j in 1:p) {
    temp <- temp+error[i-j+p]*phi.true[j]
  }
  innovation[i] <- rt(1,df=2)
  error[i+p] <- temp+innovation[i]
}
error <- error[-seq(1:p)]

y <- mu+error
```

4. Estimated mean of the data, and the autoregressive coefficients:
```markdown
fit.osar <- osar(y,p=3)
fit.osar$muhat
fit.osar$phihat
```

5. Plot the estimated innovation density function:
```markdown
par(mfrow=c(1,1),mar=c(3,3,1,1),cex=1,cex.main=1,las=0,mgp=c(1.5,0.5,0))
hist(innovation,xaxt='n',yaxt='n',xlab='',ylab = '',main="",probability = T,nclass = 50)
xp <- seq(-max(fit.osar$knots),max(fit.osar$knots),length=200)
lines(fit.osar$fhat(xp)~xp,type="l",lwd=2)
```
