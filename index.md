# Xin Chen
My name is Xin, and I am studying for a PhD：statistics at Colorado State Unversity.

## RESEARCH INTERESTS:
- Nonparametric density estimation
- Robust statistics
- Penalized methods

### Submitted papers

### Code for Auto-Adaptive M-Estimation
Here is an example implementation, finding the mean of a heavy-tailed data, and the error density.

Load necessary library
```markdown
library(quadprog)
```
Generate the data
```markdown
set.seed(1)
y=rt(500,df=2)
```
Fit the model
```markdown
fit=onesamp(y)
```
Estimated mean of the data
```markdown
fit$muhat 
```
Confidence interval of the mean
```markdown
fit$confidence.intervals 
```
Plot the estimated error density function
```markdown
xp=seq(-max(fit$knots),max(fit$knots),length=100)
hist(y,xaxt='n',yaxt='n',xlab='',ylab = '',main="",probability = T,nclass = 100)
lines(fit$fhat(xp)~xp,lwd=2)
```

### Support or Contact
Email: xin.chen@colostate.edu
