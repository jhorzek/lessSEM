$$
\begin{array}{l|llll}
\text{penalty} & \text{function} & \text{smooth version} & \text{optimizer} & \text{reference}\\
\hline
\text{ridge} & p( x_j) = \lambda x_j^2 & \lambda x_j^2 & \text{bfgs, glmnet, ista} & \text{(Hoerl \& Kennard, 1970)}\\
\text{lasso} & p( x_j) = \lambda| x_j| & \lambda\sqrt{ x_j^2 + \varepsilon}; \varepsilon > 0 & \text{bfgs, glmnet, ista} & \text{(Tibshirani, 1996)}\\
\text{adaptiveLasso} & p( x_j) = \frac{1}{w_j}\lambda| x_j| & \frac{1}{w_j}\lambda\sqrt{ x_j^2 + \varepsilon}; \varepsilon > 0 & \text{bfgs, glmnet, ista} & \text{(Zou, 2006)}\\
\text{elasticNet} & p( x_j) = \alpha\lambda| x_j| + (1-\alpha)\lambda x_j^2 & \alpha\lambda\sqrt{ x_j^2 + \varepsilon} + (1-\alpha)\lambda x_j^2; \varepsilon > 0 & \text{bfgs, glmnet, ista} & \text{(Zou \& Hastie, 2005)}\\
\text{cappedL1} & p( x_j) = \lambda \min(| x_j|, \theta); \theta > 0 & -- & \text{ista} & \text{(Zhang, 2010)}\\
\text{lsp} & p( x_j) = \lambda \log(1 + |x_j|/\theta); \theta > 0 & -- & \text{ista} & \text{(Candès et al., 2008)} \\
\text{scad} & p( x_j) = \begin{cases}
\lambda |x_j| & \text{if } |x_j| \leq \lambda\\
\frac{-x_j^2 + 2\theta\lambda |x_j| - \lambda^2}{2(\theta -1)} & \text{if } \lambda < |x_j| \leq \lambda\theta \\
(\theta + 1) \lambda^2/2 & \text{if } |x_j| \geq \theta\lambda\\
\end{cases}; \theta > 2 & -- & \text{ista} & \text{(Fan \& Li, 2001)} \\
\text{mcp} & p( x_j) = 
\begin{cases}
\lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
\theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
\end{cases}; \theta > 0 & -- & \text{ista} & \text{(Zhang, 2010)}
\end{array}
$$
