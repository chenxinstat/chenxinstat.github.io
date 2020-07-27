# Xin Chen
My name is Xin, and I am studying for a PhD：statistics at Colorado State Unversity.

## RESEARCH INTERESTS:
- Nonparametric density estimation
- Robust statistics
- Penalized methods

### Submitted papers

### Code for Auto-Adaptive M-Estimation

```markdown
Here is an example implementation, finding the mean of a heavy-tailed data, and the error density.

library(quadprog)

# generating the data
set.seed(1)
y=rt(500,df=2)

# fit
fit=onesamp(y)

# estimated mean of the data
fit$muhat 

#confidence interval of the mean
fit$confidence.intervals 

# plot the estimated error density function
xp=seq(-max(fit$knots),max(fit$knots),length=100)
hist(y,xaxt='n',yaxt='n',xlab='',ylab = '',main="",probability = T,nclass = 100)
lines(fit$fhat(xp)~xp,lwd=2)
```

### Support or Contact
Email: xin.chen@colostate.edu
