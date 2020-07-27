Here is an example implementation of Auto Adaptive M-Estimation, finding the mean of a heavy-tailed data, and the error density.

Load necessary library:
```markdown
library(quadprog)
```
Generate a data and fit the model:
```markdown
set.seed(1)
y <- rt(500,df=2)
fit <- onesamp(y)
```
Estimated mean of the data, and its confidence interval:
```markdown
fit$muhat 
fit$confidence.intervals 
```
Plot the estimated error density function:
```markdown
xp <- seq(-max(fit$knots),max(fit$knots),length=100)
hist(y,xaxt='n',yaxt='n',xlab='',ylab = '',main="",probability = T,nclass = 100)
lines(fit$fhat(xp)~xp,lwd=2)
```
