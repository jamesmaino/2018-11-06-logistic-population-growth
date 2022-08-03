---
title: Modelling density dependent (logistic) population growth
excerpt: When population growth can't go on forever
date: '2018-11-06'
isFeatured: true
tags: ["ecology", "population", "exponential", "logistic", "derivation", "R"]
---

import MarkdownWrapper from '../../components/markdown-wrapper' 
export default ({ children }) => <MarkdownWrapper>{children}</MarkdownWrapper>


# Modelling density dependent (logistic) population growth

Let's derive some more population growth functions!

In [a previous post](/posts/20220716-new-post) we derived a function for population growth based on the vital rates of reproduction and mortality. We assumed that the growth rate was constant with respect to the number of individuals in the population or $\frac{{dt}}{{dN}} = rN$. This led to the unrealistic prediction that populations will grow indefinitely. Of course, populations will eventually run into resource problems (e.g. no food) or other density dependent issues (e.g. over-crowding and disease). 

One way to set un upper limit on population growth is to scale the growth rate $rN$ by a function that is 1 when $N$ is small and 0 when $N$ is big. Here is one such function: 

$$\frac{{dN}}{{dt}} = (1 - N/K)rN$$

When $N$ is small the growth rate is approximately $rN$ but as $N \to K$ the growth rate scales to zero. Thus $K$ can be thought of as the carrying capacity of the population. It is important to note that unlike the equation for unconstrained exponential growth, this is a function of convenience rather than having mechanistic interpretation based on how reproduction or mortality depend on population density - it is the simplest function that satisfies our requirements for density dependent growth. 

## Numerical solution
Before we solve this equation with maths (analytically), let's solve it using computers (numerically) for $r=0.3$, and $K=100$ assuming the population started with just 2 individuals. 


```r
library(tidyverse)
library(deSolve)
parameters <- c(r = 0.3, K = 100)
state <- c(N = 2)
logistic <- function(t, state, parameters) {{
  with(as.list(c(state, parameters)), {{
    dN <-  (1 - N / K) * (r * N )
    list(dN)
  }})
}}
times <- seq(0, 50, by = 0.2)
sol = as.data.frame(
  ode(y = state, times = times, func = logistic, parms = parameters)
)
 
```

```r
ggplot(sol) + 
  geom_line(aes(time, N), col = 'red')+
   theme_minimal() +
  theme(
    plot.background = element_rect(fill = rgb(.2,.21,.27)),
    text = element_text(colour = 'grey'), 
    axis.text = element_text(colour = 'grey'), 
    panel.grid = element_line(colour = 'grey')
  )
```
![Logistic equation solved numerically](https://github.com/jamesmaino/2018-11-06-logistic-population-growth/raw/main/plots/ode1.png)

Wow! That was easy. 

## Analytical solution
Now that we've seen the pretty curve, let's use our excitement to try and describe the curve with a single equation without derivatives (solve the differential equation for $N$ as a function of $t$ and get rid of the $\frac{{dN}}{{dT}}$). 

$$\frac{{dN}}{{dt}} = (1 - N/K)rN$$

First, separate the variables. 

$${{rdt}} = \frac{{1}}{{(1 - N/K)N}} {{dN}}$$

Now we want to split up the fraction into something more managable so we assume there exists some $A$ and $B$ such that: 

$$\frac{{1}}{{(1 - N/K)N}} = \frac{{A}}{{(1 - N/K)}} + \frac{{B}}{{N}}$$

$$\frac{{1}}{{(1 - N/K)N}} = \frac{{AN+B(1 - N/K)}}{{(1 - N/K)N}}$$

$$1 = AN + B(1 - N/K)$$

When $N=K$:

$$1 = AK + B(1 - K/K)$$

$$A = 1/K$$

Now substituting $A$ and N = 0:

$$1 = (0)/K + B(1 - (0)/K)$$

$$B = 1, A = 1/K$$

Thus we have:

$${{rdt}} = \frac{{1/K}}{{(1 - N/K)}} + \frac{{1}}{{N}} {{dN}}$$

Now integrate both sides:

$$\int {{rdt}} = \int \frac{{1/K}}{{(1 - N/K)}} + \frac{{1}}{{N}} {{dN}}$$

$$rt = \int \frac{{1/K}}{{(1 - N/K)}}dN + \int \frac{{1}}{{N}} {{dN}}$$

$$rt = \int \frac{{1/K}}{{(1 - N/K)}}dN + \ln|N|$$

The part $\int \frac{{1}}{{(K - N)}}dN$ is a little tricky, but substituting $u = 1 - N/K$ and $du = -1/K dN$ will help. 


From here we can finish the original integration.  

$$rt = \int \frac{{1/K}}{{(1 - N/K)}}dN + \ln|N| + c$$

$$rt =  -\ln|1 - N/K| + \ln|N| + c$$

As N is strictly positive we can remove the absolute operators. Remembering your logarithm rules $\ln(ab) = \ln(a) + \ln(b)$ we get:

$$rt = \ln\frac{{cN}}{{1 - N/K}}$$

$$ce^{{rt}} = \frac{{N}}{{1 - N/K}}$$

We really just want one $N$ so that we can convieniently write the function but for now let's solve the constant $c$ at $t=0$ and $N=N_0$

$$ce^{{r(0)}} = c = \frac{{N_0}}{{1 - N_0/K}} = \frac{{KN_0}}{{K - N_0}}$$

Back to solving for $N$. 

$$ce^{{rt}} = \frac{{N}}{{1 - N/K}}$$

$$(1 - N/K)ce^{{rt}} ={{N}}$$

$$ce^{{rt}} - ce^{{rt}}N/K ={{N}}$$

$$ce^{{rt}} = {{N}} + ce^{{rt}}N/K$$

$$ce^{{rt}} = N({{1}} + ce^{{rt}}/K)$$

$$N = \frac{{ce^{{rt}}}} {{({{1}} + ce^{{rt}}/K)}}$$

Substitute in the constant we solved for earlier.

$$N = \frac{{\frac{{KN_0}}{{K - N_0}}e^{{rt}} }} {{ ( {{1}} + \frac{{KN_0}}{{K - N_0}}e^{{rt}}/K)}}$$

And simplify by multiplying through by $(K - N_0)e^{{rt}}$

$$N = \frac{{KN_0}} {{ {{N_0}} + (K - N_0)e^{{-rt}} }}$$

Let plot it!

```r
r = 0.3
K = 100
N0 = 2
time <- seq(0, 50, by = 0.2)
sol = data.frame(
  time,
  N = K*N0/(N0 + (K - N0)*exp(-r*time))
)
 
```

```r
ggplot(sol) + 
  geom_line(aes(time, N), col = 'red')+
   theme_minimal() +
  theme(
    plot.background = element_rect(fill = rgb(.2,.21,.27)),
    text = element_text(colour = 'grey'), 
    axis.text = element_text(colour = 'grey'), 
    panel.grid = element_line(colour = 'grey')
  )
```

![Logistic equation solved analytically](https://github.com/jamesmaino/2018-11-06-logistic-population-growth/raw/main/plots/ode2.png)

All that calculus for the same line `deSolve` gave us...

[Source code.](https://github.com/jamesmaino/2018-11-06-logistic-population-growth/)
