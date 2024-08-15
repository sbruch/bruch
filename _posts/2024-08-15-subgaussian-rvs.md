---
title: "Concentration of Sub-Gaussian Random Variables"
subtitle: "This post reviews basic tail bounds (Markov, Chebyshev, Chernoff, Hoeffding)
and explains the concept of sub-Gaussian random variables."
layout: default
date: 2024-08-15
keywords: concentration
published: true
category: "Concentration"
---

{% katexmm %}

### Basic Tail Bounds

One of the most basic inequalities governing the tail bounds of a distribution is **Markov's inequality**. It says that, if $X$ is a *non-negative* random variable, then its tail probability has the following upper-bound:

$$
\mathbb{P}[X \geq t] \leq \frac{\mathbb{E}[X]}{t}. \tag{1}
$$

Clearly, this only makes sense if $t > 0$ and $\mathbb{E}[X]$ is finite.


The Markov inequality is tight, in the sense that in cannot be improved in general. To see this, consider the random variable $X$ taking values in $\{ 0, \alpha \}$ with probability $p$ and $1-p$. A straightforward calculation shows that the inequality in (1) becomes equality.

More generally, notice that:

$$
\begin{aligned}
\mathbb{P}[X \geq t] = \frac{\mathbb{E}[X]}{t} &\Rightarrow \int_{t}^\infty f(x)\;dx = \int_{0}^\infty \frac{1}{t} x f(x)\; dx \\
 &\Rightarrow \int_{0}^\infty 1_{x \geq t} f(x)\;dx = \int_{0}^\infty \frac{1}{t} x f(x)\; dx\\
 &\Rightarrow \int_0^\infty (t 1_{x \geq t} - x) f(x) \; dx = 0 \\
 &\Rightarrow x=t1_{x\geq t}.
\end{aligned}
$$

As such, any random variable with the property that $X=t 1_{X \geq t}$ almost surely leads to equality in Markov's bound.



#### From Markov to Chebyshev to Chernoff

A fairly trivial extension of Markov's inequality is the **Chebyshev's inequality.** It is arrived at by a simple application of Equation (1) to the non-negative random variable $X := (Y - \mathbb{E}[Y])^2$ for an arbitrary random variable $Y$:

$$
\mathbb{P}\big[ \lvert Y - \mathbb{E}[Y] \rvert \geq t \big] \leq \frac{\mathop{Var}[Y]}{t^2}. \tag{2}
$$



Similar to the Markov inequality, we can show that Chebyshev's inequality is tight in general. It is easy to construct a random variable that results in equality if we notice that we must have $X = (Y - \mathbb{E}[Y])^2 \in \{0, \alpha\}$. For example:

$$
Y = \begin{cases}
-\sqrt{\alpha} & \text{with probability } p/2 \\
0 & 1-p \\
\sqrt{\alpha} & p/2
\end{cases},
$$



One can just as easily extend Chebyshev's inequality to higher moments, so that:

$$
\mathbb{P}\big[ \lvert Y - \mathbb{E}[Y] \rvert \geq t \big] \leq \frac{\mathbb{E}\big[ \big\lvert Y - \mathbb{E}[Y] \big\rvert^k \big]}{t^k}. \tag{3}
$$

Even more interesting is the fact that we can extend this formulation to functions other than polynomials! For example, consider the exponential function:

$$
\mathbb{P}[ X - \mathbb{E}[X] \geq t ] = \mathbb{P}[ e^{\lambda (X - \mathbb{E}[X])} \geq e^{\lambda t} ] \leq \frac{\mathbb{E}\big[ e^{\lambda (X - \mathbb{E}[X])} \big]}{e^{\lambda t}}. \tag{4}
$$

Optimizing our choice of $\lambda$ so as to obtain the tightest result gives us the **Chernoff bound**!


<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
We can show that, Equation (3) with the optimal choice of $k$ is never worse than Equation (4) with the optimal choice of $\lambda$. Concretely, suppose $X \geq 0$ and that the moment generating function of $X$ exists in a neighborhood around $0$. Then for any given $\delta > 0$, we can see that:

$$
\begin{aligned}
\mathbb{E}[e^{\lambda X}] &= \sum_i \frac{\lambda^i \mathbb{E}[X^i]}{i!} = \sum_i \frac{(\lambda \delta)^i}{i!} \frac{\mathbb{E}[X^i]}{\delta^i} \\
&\geq \sum_i \frac{(\lambda \delta)^i}{i!} \inf_{k=1,2,3,\ldots} \frac{\mathbb{E}[X^k]}{\delta^k} = e^{\lambda \delta} \inf_{k=1,2,3,\ldots} \frac{\mathbb{E}[X^k]}{\delta^k},
\end{aligned}
$$

so that:

$$
\inf_{k=1,2,3,\ldots} \frac{\mathbb{E}[X^k]}{\delta^k} \leq \inf_{\lambda > 0} \frac{\mathbb{E}[e^{\lambda X}]}{e^{\lambda \delta}}.
$$
</div>
</div>

#### Application to Gaussian variables

Consider $X \sim \mathcal{N}(\mu, \sigma^2).$ Let's derive the Chernoff bound for $X$ by optimizing our choice of $\lambda$. First, let's derive the right-hand-side in Equation (4):

$$
\begin{aligned}
e^{-\lambda t} \; \mathbb{E}\big[ e^{\lambda (X - \mu)} \big] &= e^{-\lambda t} \; \int \frac{1}{\sqrt{2\pi \sigma^2}} \exp(-\frac{(x - \mu)^2}{2\sigma^2}) \exp(\lambda (x - \mu))\; dx \\
&\stackrel{y:=x-\mu}{=} e^{-\lambda t} \; \int \frac{1}{\sqrt{2\pi\sigma^2}} \exp(-\frac{y^2 - 2\sigma^2\lambda y + \sigma^4\lambda^2}{2\sigma^2}) \exp(\frac{\sigma^2\lambda^2}{2}) \; dx \\
&= e^{-\lambda t} \; \exp(\frac{\sigma^2\lambda^2}{2}) \int \frac{1}{\sqrt{2\pi\sigma^2}} \exp(-\frac{(y - \lambda)^2}{2}) \; dx \\
&= e^{-\lambda t} \; e^{\sigma^2\lambda^2/2}.
\end{aligned}
$$

Differentiating with respect to $\lambda$ and setting it to $0$:

$$
\begin{aligned}
\frac{\partial}{\partial \lambda} e^{-\lambda t} e^{\sigma^2\lambda^2/2} =(\sigma^2 \lambda - t) e^{-\lambda t} e^{\sigma^2\lambda^2/2} = 0,
\end{aligned}
$$



gives us the optimal choice $\lambda = t/\sigma^2$. Plugging this value back in Equation (4) gives the Chernoff bound for our Gaussian random variable:

$$
\mathbb{P}[ X - \mu \geq t ] \leq \frac{e^{t^2/2\sigma^2}}{e^{t^2/\sigma^2}} = e^{-t^2/2\sigma^2}. \tag{5}
$$



### Sub-Gaussian random variables

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
It is easy to see that, if the following holds for a random variable $X$:

$$
\mathbb{E}\big[ e^{\lambda (X - \mu)}  \big] \leq e^{\sigma^2 \lambda^2/2} \quad \forall \; \lambda \in \mathbb{R}, \tag{6}
$$

for some constant $\sigma$, then inequality (5) holds for $X$. We call such random variables *Sub-Gaussian with parameter* $\sigma$*.*
</div>
</div>


**Example (1)**: Rademacher random variables are sub-Gaussian with parameter $\sigma=1$. That is because:


$$
\begin{aligned}
\mathbb{E}[e^{\lambda X}] &= \frac{1}{2}\big( e^{\lambda} + e^{-\lambda} \big) \\
&= \frac{1}{2} \big( \sum_{k=0} \frac{\lambda^k}{k!} + \frac{(-\lambda)^k}{k!} \big) \\
&= \sum_{k=0} \frac{\lambda^{2k}}{(2k)!} \\
&\leq \sum_{k=0} \frac{\lambda^{2k}}{2^k k!} = e^{\lambda^2/2}.
\end{aligned}
$$



**Example (2)**: Any zero-mean bounded variable $X$ supported on the interval $[a, b]$ is sub-Gaussian with parameter at most $\sigma=b-a$.



Proving this requires the *symmetrization* *argument.* We first introduce an independent copy of $X$ and call it $X^\prime$, and let $\gamma$ be a Rademacher random variable. By the convexity of the exponential and Jensen's inequality, we can write:

$$
\mathbb{E}[e^{\lambda X}] = \mathbb{E}\big[ e^{\lambda X - \mathbb{E}_{X^\prime}[X^\prime]} \big] \leq \mathbb{E}_{X, X^\prime} [e^{\lambda(X - X^\prime)}].
$$

Using the fact that the distribution of $X - X^\prime$ is the same as the distribution of $\gamma (X - X^\prime)$, and applying the result from Example (1), we can write:

$$
\mathbb{E}[e^{\lambda X}] \leq \mathbb{E}_{X, X^\prime, \gamma} \big[ 
e^{\lambda \gamma (X - X^\prime)} \big] \leq \mathbb{E}_{X, X^\prime} \big[ e^{\lambda^2 (X - X^\prime)^2 /2} \big] \leq e^{\lambda^2 (b - a)^2/2}.
$$



Example (2) essentially gives us **Hoeffding's inequality**. Concretely, because the sum of two sub-Gaussian random variables with constants $\sigma_1^2$ and $\sigma_2^2$ is itself a sub-Gaussian random variable with parameter $\sigma_1^2 + \sigma_2^2$, we can show that if $X_i$'s are independent sub-Gaussian random variables with parameter $\sigma_i$, then we have that:

$$
\mathbb{E}[\sum (X_i - \mu_i) \geq t] \leq \exp\Big({-\frac{t^2}{2\sum \sigma_i^2}}\Big).
$$

For bounded random variables, we can derive the familiar form:

$$
\mathbb{E}[\sum (X_i - \mu_i) \geq t] \leq \exp{-\frac{t^2}{2 n (b-a)^2}}.
$$

The result above can be further tightened by showing that, bounded random variable $X \in [a, b]$ is indeed sub-Gaussian with parameter $(b - a)/2$.


<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
Proving that last statement is a bit involved. We'd have to define the following function:

$$
\psi(\lambda) = \log \mathbb{E}[e^{\lambda X}].
$$

The idea here is expand $\psi$ around $0$:

$$
\psi(\lambda) = \lambda \psi^\prime(0) + \frac{\lambda^2}{2} \psi^{\prime\prime}(0) + \mathcal{O}(\lambda^3).
$$

Notice how the right-hand-side has the form we hope to get to in Equation (6): $\exp(\psi) \leq \exp(\lambda \mu + \lambda^2 \sigma^2/2)$? So if we showed that $\psi^\prime(0) \leq \mu$ and $\psi^{\prime\prime}(0) \leq (b-a)^2/4$, then we have proven our claim.



To that end, consider the derivative of $\psi$ with respect to $\lambda$ which is simply:

$$
\psi^\prime(\lambda) = \mathbb{E}[Xe^{\lambda X}] / \mathbb{E}[e^{\lambda X}].
$$

It's clear that $\psi^\prime(0)=\mu$, so that takes care of the first part. Its second derivative then is:

$$
\psi^{\prime\prime}(\lambda) = \frac{\mathbb{E}[X^2 e^{\lambda X}]}{\mathbb{E}[e^{\lambda X}]} - \Big( \frac{\mathbb{E}[X e^{\lambda X}]}{\mathbb{E}[e^{\lambda X}]} \Big)^2.
$$

Now, we have to bound this to show that, when evaluated at $0$, it is no greater than $(b-a)^2/4$. Showing this is a bit tricky and involves a bit of measure theory. Let's rearrange the terms in the last expression so it's visually clearer what's going on:

$$
\psi^{\prime\prime}(\lambda) = \mathbb{E}_X\Bigg[ X^2 \frac{e^{\lambda X}}{\mathbb{E}[e^{\lambda X}]} \Bigg]  - \mathbb{E}_X\Bigg[ X \frac{e^{\lambda X}}{\mathbb{E}[e^{\lambda X}]} \Bigg]^2.
$$

This looks almost like the expression for the variance of a random variable, but with an extra term $e^{\lambda X} / \mathbb{E}[e^{\lambda X}]$. Here, we can appeal to the [Radon-Nikodym theorem](https://en.wikipedia.org/wiki/Radon-Nikodym_theorem). Because the conditions of the theorem hold (i.e., measures are absolutely continuous), we can define a new distribution $Q_\lambda$ so that $e^{\lambda X} / \mathbb{E}[e^{\lambda X}]$ is its Radon-Nikodym derivative with respect to the distribution over $X$. What that helps us accomplish is that, $\psi^{\prime\prime}$ can now be expressed as the variance of a random variable drawn from $Q_\lambda$:

$$
\psi^{\prime\prime}(\lambda) = \mathbb{E}_{X_\lambda \sim Q_\lambda}[X_\lambda^2] - \mathbb{E}_{X_\lambda \sim Q_\lambda}[X_\lambda]^2.
$$

 So, using the fact that the mean is the minimizer of the mean squared error, and knowing that $X_\lambda \in [a, b]$, we can write:

$$
\psi^{\prime\prime}(\lambda) = \mathbb{E}[(X_\lambda - \mathbb{E}[X_\lambda])^2] \leq \mathbb{E}[(X_\lambda - \frac{a+b}{2})^2] \leq (b - \frac{a+b}{2})^2 = (b - a)^2/4.
$$

Putting all the pieces together, we have shown that $\psi(\lambda) \leq \lambda \mu + \lambda^2(b-a)^2/8$, which shows that a random variable bounded to the interval $[a, b]$ is sub-Gaussian with constant $(b - a)/2$.
</div>
</div>

{% endkatexmm %}