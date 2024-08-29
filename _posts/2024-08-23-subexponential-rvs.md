---
title: "Concentration of Sub-Exponential Random Variables"
subtitle: "This post reviews the definition of sub-exponential random variables
 and discusses useful tail bounds (Bernstein) for this family."
layout: default
date: 2024-08-23
keywords: concentration
published: true
category: "Concentration"
---

{% katexmm %}

### From Sub-Gaussian to Sub-Exponential

Recall that a random variable $X$ with mean $\mu$ is [sub-Gaussian](/blog/concentration/2024/08/15/subgaussian-rvs) with constant $\sigma$ if its moment generating function is bounded above as follows:

$$
\mathbb{E}[e^{\lambda (X - \mu)}] \leq e^{\lambda^2 \sigma^2/2} \quad \forall \; \lambda \in \mathbb{R}. \tag{1}
$$

That proved useful because, if $X$ is sub-Gaussian, then we showed that the tail of its distribution is dominated by a Gaussian random variable's with parameters $\mu$ and $\sigma^2$. In particular, we are able to apply the Chernoff bound to sub-Gaussian random variables and derive Hoeffding inequality-type bounds on the tail distribution.



Sub-Gaussianness, however, is somewhat restrictive. It requires that the bound on the moment generating function hold for *all* values of $\lambda$! We can relax that requirement just a little bit and demand that the inequality hold for small values of $\lambda$ around $0$ instead:

$$
\mathbb{E}[e^{\lambda (X - \mu)}] \leq e^{\lambda^2 \nu^2/2} \quad \forall \; \lvert \lambda \rvert \leq \frac{1}{b}. \tag{2}
$$

A random variable $X$ with mean $\mu$ for which inequality (2) holds is said to be ***sub-exponential*** with parameters $(\nu, b)$ for some non-negative $b$ (and defining $1/0=\infty$).



Before we continue, it's worth mentioning that sub-exponentialness is preserved under summation of independent random variables, just as sub-Gaussiannes is. In other words, if we have $n$ sub-exponential random variables $X_i$'s with parameters $(\nu_i, b_i)$, then their sum $X=\sum_i X_i$ is sub-exponential:

$$
\mathbb{E}[e^{\lambda (X - \mathbb{E}[X])}] = \mathbb{E}[\prod_i e^{\lambda (X_i - \mu_i)}] = \prod_i \mathbb{E}[e^{\lambda (X_i - \mu_i)}] \leq \prod_i e^{\lambda^2\nu_i^2/2} = e^{\lambda^2 \sum_i \nu_i^2/2},
$$

which is valid when $\lvert \lambda \rvert \leq 1/\max_i b_i$. So its parameters are $((\sum_i \nu_i^2)^{1/2}, \max_i b_i)$.

#### Equivalent Definitions

There are other ways of characterizing a sub-exponential random variable that are equivalent to the definition in (2). For example, if the moment generating function of $X$ is finite for small values of $\lambda$, then $X$ is indeed a sub-exponential variable! This, and other definitions, are presented in the result below.


<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
The following definitions are equivalent for a centered random variable $X$:

i. There exist constants $\nu$ and $b$ such that inequality (2) holds.

ii. There exists a constant $c_0$ such that $\mathbb{E}[e^{\lambda X}]$ is finite for $\lvert \lambda \rvert \leq c_0$.

iii. There exist constants $c_1$ and $c_2$ such that $\mathbb{P}[\lvert X \rvert \geq t] \leq c_1 e^{-c_2t}$ for all $t > 0$.

iv. The quantity $\sup_{\lambda \geq 2}\Big( \frac{\mathbb{E}[X^k]}{k!} \Big)^{1/k}$ is finite.
</div>
</div>

<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**. (i) $\Leftrightarrow$ (ii) is straightforward.



(ii) $\Rightarrow$ (iii): Applying the Chernoff bound to $X$ gives:

$$
\mathbb{P}[X \geq t] \leq \frac{\mathbb{E}[e^{\lambda X}]}{e^{\lambda t}} \Bigg\rvert_{\lambda = c_0/2} = \mathbb{E}[e^{c_0 X/2}] e^{-c_0 t/2}.
$$

So setting $c_1 = \mathbb{E}[e^{c_0 X/2}] + \mathbb{E}[e^{-c_0 X/2}]$ and $c_2=c_0/2$ gives the desired result.



(iii) $\Rightarrow$ (ii): We start from the definition of expectation of a positive random variable, for a positive $\lambda$:

$$
\begin{aligned}
\mathbb{E}[e^{\lambda \lvert X \rvert}] &= \int_0^\infty \mathbb{P[e^{ \lambda \lvert X \rvert} \geq t]} \; dt \leq 1 + \int_1^{\infty} \mathbb{P}[\lvert X \rvert \geq \frac{\log t}{\lambda}] \; dt \\
&\leq 1 + c_1 \int_1^\infty e^{-c_2 \log t / \lambda} \; dt = 1 + c_1 \int_1^\infty t^{-c_2 / \lambda} \; dt \leq 1 + c_1,
\end{aligned}
$$

when $c_2 / \lambda \geq 2$. Because $\mathbb{E}[e^{-\lambda X}] \leq \mathbb{E}[e^{\lvert \lambda \rvert \lvert X \rvert}]$ and $\mathbb{E}[e^{\lambda X}] \leq \mathbb{E}[e^{\lvert \lambda \rvert \lvert X \rvert}]$, we can conclude that $\mathbb{E}[e^{\lambda X}]$ is finite for $\lvert \lambda \rvert \leq c_2/2$.



(ii) $\Leftrightarrow$ (iv): We can argue using the power series expansion of the moment generating function:

$$
\mathbb{E}[e^{\lambda X}] = 1 + \sum_{k=2}^\infty \frac{\mathbb{E}[X^k]}{k!} \lambda^k.
$$

Notice that the inverse of $\limsup_{\lambda \geq 2} \Big\lvert \frac{\mathbb{E}[X^k]}{k!} \Big\rvert^{1/k}$is the convergence radius of the series. As a result, if (ii) holds, then the series must converge, implying that the quantity in (iv) must be finite. If (iv) holds, then the series must converge, implying that the quantity in (ii) must be finite for within the convergence radius.
</div>
</div>

#### Bernstein's Condition

An alternative way to characterize sub-exponential random variables is by bounding the polynomial moments of $X$, instead of its moment generating function. Concretely, if $X$ is a random variable with mean $\mu$ and variance $\sigma^2$, then we can show that, if the following condition, known as *Bernstein's condition with parameter* $b$*,* holds, then $X$ is a sub-exponential random variable:

$$
\Big\lvert \mathbb{E}\big[ (X - \mu)^k \big] \Big\rvert \leq \frac{1}{2} k! \sigma^2 b^{k-2} \quad \text{for } k = 3, 4, \ldots
$$

<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**. Let's start from the definition of a sub-exponential random variable and expand the left-hand-side of (2):

$$
\mathbb{E}[e^{\lambda (X - \mu)}] = \sum_{k=0}^\infty \frac{\lambda^k \mathbb{E}[(X - \mu)^k]}{k!} = 1 + \frac{\lambda^2 \sigma^2}{2} + \sum_{k=3}^\infty \frac{\lambda^{k} \mathbb{E}[(X - \mu)^{k}]}{k!}.
$$

By Bernstein's condition:

$$
\mathbb{E}[e^{\lambda (X - \mu)}] \leq 1 + \frac{\lambda^2\sigma^2}{2} + \sum_{k=3}^\infty \frac{1}{2} \lvert \lambda \rvert^k \sigma^2 b^{k-2} = 1 + \frac{\lambda^2\sigma^2}{2} + \frac{\lambda^2\sigma^2}{2} \sum_{k=3}^\infty (\lvert \lambda \rvert b)^{k-2},
$$

which, when $\lvert \lambda \rvert b$ is smaller than $1$, is equal to:

$$
1 + \frac{\lambda^2\sigma^2}{2} + \frac{\lambda^2\sigma^2}{2} \Big(\frac{1}{1-\lvert \lambda \rvert b} -1 \Big) \leq \exp\Big( \frac{\lambda^2 \sigma^2}{2(1 - \lvert \lambda \rvert b)} \Big) \leq e^{\lambda^2 (\sqrt{2}\sigma)^2 /2},
$$

when $\lvert \lambda \rvert < 1/2b$. That implies that $X$ is sub-exponential with parameters $(\sqrt{2}\sigma, 2b)$.
</div>
</div>


This is a neat result! For example, it is straightforward to see why Bernstein's condition holds for bounded random variables with the property that $\lvert X - \mu \rvert \leq b$. We proved previously that such random variables are sub-Gaussian. It turns out, they are sub-exponential too!

### Concentration Inequalities

Why is this notion of sub-exponentialness useful anyway? Well, for one thing, the family of sub-exponential random variables includes sub-Gaussian random variables. But more importantly, limiting our focus to a neighborhood around $0$ enables us to derive more expressive bounds on the tail distribution of a sub-exponential $X$. To see why, let's first derive the basic tail bound and then discuss a more elaborate bound that contrasts the concentration of sub-Gaussian and sub-exponential variables.

#### Sub-Exponential Tail Bound

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
For a sub-exponential random variable $X$ with parameters $(\nu, b)$, the following holds:

$$
\mathbb{P}[X \geq \mu + t] \leq \begin{cases}
e^{-t^2/2\nu^2} & \text{if } 0 \leq t \leq \frac{\nu^2}{b},\\
e^{ -t/2b } & \text{if } t > \frac{\nu^2}{b}
\end{cases}. \tag{3}
$$
</div>
</div>

<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**. For valid values of $\lambda$ (i.e., $\lvert \lambda \rvert \leq 1/b$), we have, by definition, that:

$$
\mathbb{P}[X - \mu \geq t] \leq \inf_{\lvert \lambda \rvert \in [0, 1/b]} \frac{\mathbb{E}[e^{\lambda (X - \mu)}]}{e^{\lambda t}} = \inf_{\lvert \lambda \rvert \in [0, 1/b]} \underbrace{e^{\lambda^2 \nu^2/2 - \lambda t}}_{f(\lambda, t)}.
$$

When $t \leq \nu^2/b$, we can solve the unconstrained minimization problem and obtain:

$$
\frac{\partial}{\partial \lambda} f(\lambda, t) = 0 \Rightarrow \nu^2 \lambda - t = 0 \Rightarrow \lambda = \frac{t}{\nu^2}.
$$

Substituting this value of $\lambda$ into $f(\lambda, t)$ gives the desired result, as in the case of a sub-Gaussian random variable.



On the other hand, when $t > \nu^2/b$, then $f(\lambda, t)$ reaches its minimum at the boundary point when $\lambda = 1/b$. That's because $f(\cdot, t)$ is monotonically decreasing in its first argument. Evaluating $f(1/b, t)$ then leads to:

$$
f(\frac{1}{b}, t) = \exp\Big( \frac{\nu^2}{2b^2} - \frac{t}{b} \Big) \leq \exp\Big( -\frac{t}{2b} - \frac{t}{b} \Big) = \exp\Big(-\frac{t}{2b} \Big).
$$
</div>
</div>

#### Bernstein Tail Bound

We can do something even more interesting if we consider Bernstein's condition! If we stare at the proof of the statement above that Bernstein's condition is sufficient to show that a random variable is sub-exponential, we notice an interesting result hiding in there. In particular, we derived, in that proof, the following bound when $\lvert \lambda \rvert b < 1$:

$$
\mathbb{E}[e^{\lambda (X - \mu)}] \leq e^{\lambda^2\sigma^2/2(1 - \lvert \lambda \rvert b)}.
$$

Setting $\lambda = t/(bt + \sigma^2) \in [0, \frac{1}{b})$ results in the following concentration inequality:

$$
\mathbb{P}[X - \mu \geq t] \leq e^{-t^2/2(\sigma^2 + bt)}.\tag{4}
$$

Why is that significant? Consider a bounded random variable $X$ such that $\lvert X - \mu \rvert \leq b$. We argued above that $X$ is sub-exponential as it meets Bernstein's condition. The tail bound in (4) then introduces a dependence on, not just $b$, as the Hoeffding inequality would (as a consequence of $X$ being sub-Gaussian), but also on the variance of $X$!

Because $\sigma^2 \leq b^2$ for a bounded random variable, and for suitably small values of $t$, (4) is never worse than the Hoeffding bound. In particular, when $\sigma^2 \ll b^2$, then (4) offers a much sharper bound on the concentration of $X$.

### Example Sub-Exponential Variables

We have already discussed the case of bounded random variables. Indeed, such variables are sub-Gaussian *and* sub-exponential. Let's now look at a random variable that's sub-exponential but *not* sub-Gaussian: the $\chi^2$ random variable.



Recall that if $Y \sim \mathcal{N}(0, 1)$ is a Gaussian random variable, then $X=Y^2$ is a $\chi^2$ random varaible whose mean is $1$. Let's inspect its moment generating function:

$$
\mathbb{E}[e^{\lambda (X - 1)}] = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} e^{\lambda (y^2 - 1)} e^{-y^2/2} dy = \frac{e^{-\lambda}}{\sqrt{1 - 2\lambda}}.
$$

Obviously, when $\lambda > 1/2$ the moment generating function does not exist! When $\lvert \lambda \rvert \leq 1/4$, then we can show that the moment generating function is bounded above by $e^{2\lambda^2}$. As such, $X$ is sub-exponential with parameters $(2, 4)$.



As another example, consider the sum of $n$ iid $\chi^2$ random variables. Because each variable is itself a sub-exponential with parameters $(2, 4)$, we know from before that their sum, too, is sub-exponential with parameters $(2\sqrt{n}, 4)$. As a result, we can derive the following bound:

$$
\mathbb{P}[\frac{1}{n}\sum X_i - 1 \geq t] \leq e^{-nt^2/8}, \tag{5}
$$

for all $t \in (0, 1)$.

{% endkatexmm %}