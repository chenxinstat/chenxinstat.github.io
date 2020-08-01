Here is an example implementation of Auto Adaptive M-Estimation (AAME) in R, finding the mean of a data and the error density. The data could be with a heavy-tailed and/or comtaminated error.

1. Run this in R: [AAME functions](r_codes/functions_aame.r)

2. Load necessary library:
```markdown
library(quadprog)
```

3. Generate a data, then fit the model:
```markdown
set.seed(1)
mu <- 100
y <- mu + rt(500,df=2)
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
lines(fit$fhat(xp)~c(xp+fit$muhat),lwd=2)
```
