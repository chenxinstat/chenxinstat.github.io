Here is an example implementation of Auto Adaptive M-Estimation (AAME) for in R, finding the mean of a data and the innovation density, with a heavy-tailed and autocorrelated error.

1. Run this in R: [AAME_AR functions](functions_aame_ar.r)

2. Load necessary library:
```markdown
library(quadprog)
```

3. Generate a data and fit the model:
```markdown
set.seed(1)
y <- rt(500,df=2)
fit <- onesamp(y)
```

4. Estimated mean of the data, and its confidence interval:
```markdown
fit$muhat 
fit$confidence.intervals 
```

5. Plot the estimated error density function:
```markdown
xp <- seq(-max(fit$knots),max(fit$knots),length=100)
hist(y,xaxt='n',yaxt='n',xlab='',ylab = '',main="",probability = T,nclass = 100)
lines(fit$fhat(xp)~xp,lwd=2)
```
