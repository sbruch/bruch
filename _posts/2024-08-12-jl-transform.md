---
title: "Johnson-Lindenstrauss Transform"
subtitle: "The JL Lemma offers a simple and elegant machinery to reduce the dimensionality of
a fixed set of high-dimensional points while approximately preserving their pairwise Euclidean
distance. This post explains the concept and describes simple randomized constructions of
such a transformation."
layout: default
date: 2024-08-12
keywords: rnla
published: true
category: "Randomized Numerical Linear Algebra"
---

{% katexmm %}

Say you have a set $\mathcal{V}$ of $m$ points in some high-dimensional space $\mathbb{R}^n$. Let's say $n$ is really large, making any subsequent computation on the $m \times n$ data matrix (let's call it $X$) rather expensive. This is where *dimensionality reduction* can prove helpful.



We know how to apply a number of tricks from basic linear algebra to do just that. For example, we can obtain a low-rank approximation of $X$ by the simple act of projecting it onto the subspace spanned by its top $k$ singular vectors. That gives a matrix $\tilde{X}$ with rank at most $k$ such that $\lVert X - \tilde{X} \rVert_F$ is small.



In some applications, it's not enough to approximately preserve the Frobenius norm of $X$. For example, we may be interested, instead, in approximately preserving the *pairwise Euclidean distance* between the $m$ points---a property that is important in Nearest Neighbor Search. That is, if $\tilde{x}_i \in \mathbb{R}^k$ is the image of $x \in \mathcal{V}$ under our dimensionality reducing transformation (where presumably $k \ll n$), then we want the following to hold for all $1 \leq i \leq j \leq m$: 

$$
(1-\epsilon) \lVert x_i - x_j \rVert_2^2 \leq \lVert \tilde{x}_i - \tilde{x}_j \rVert_2^2 \leq (1+\epsilon) \lVert x_i - x_j \rVert_2^2, \tag{1}
$$

for some $\epsilon \in (0, 1)$. Occasionally I'll write the above inequalities more compactly as:

$$
\lVert \tilde{x}_i - \tilde{x}_j\lVert_2^2 = (1\pm \epsilon) \lVert x_i - x_j \rVert_2^2. \tag{2}
$$



It's fairly easy to see that a bound on the Frobenius norm does not necessarily translate to this $(1 \pm \epsilon)$-approximation of pairwise Euclidean distances! We might in fact be able to lower the error in the Frobenius norm by moving some points farther from each other and others closer to each other. In other words, approximately preserving a *global* notion of distance is not equivalent to approximately preserving our *local* distances required by Equation (2).



That's where the seminal result due to Johnson and Lindenstrauss comes into the picture. It's led to a series of works that offer a simple and elegant, *linear* and *randomized* machinery to reduce dimensionality and yet approximately preserve local distances. These probabilistic algorithms play an important role in randomized numerical linear algebra such as subspace embedding. So it's worth reviewing the basics. Let's dive into it.

### The Johnson-Lindenstrauss Transform

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
**Distributional Johnson-Lindenstrauss (JL) Lemma**. For every positive $n$ and $\epsilon, \delta \in (0, 1)$, there is a distribution $\Pi$ over linear maps in $\mathbb{R}^{k \times n}$ with $k = \Theta(\epsilon^{-2} \log \frac{1}{\delta})$, such that, for any $u \in \mathbb{R}^n$:

$$
\mathbb{P}_{S \sim \Pi}\big[ \lVert Su \rVert_2^2 = (1 \pm \epsilon) \lVert u \rVert_2^2 \big] \geq 1 - \delta. \tag{3}
$$

</div>
</div>

I'll refer to Equation (3) as the *JL Property* and say that a random matrix $S$ *has the JL property* if Equation (2) holds for the distribution over $S$'s.



How does this lemma relate to Equation (2) and our search for a transformation that preserves pairwise Euclidean distances? By setting $u = (x_i - x_j) / \lVert x_i - x_j \rVert_2$ for all $x_i, x_j \in \mathcal{V}$, and choosing $\delta$ to be at most $\delta/{m \choose 2}$, we can apply the union bound and obtain that $\lVert Sx_i - Sx_j \lVert = (1\pm \epsilon) \lVert x_i - x_j \rVert$ with probability at least $1 - \delta$. In fact, that procedure gives the more familiar version of the Johnson-Lindenstrauss lemma, stated below:

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
**Johnson-Lindenstrauss Lemma**. Given a set $\mathcal{V}$ of $m$ points in $\mathbb{R}^n$ and $\epsilon \in (0, 1)$, there exists a mapping $f: \mathbb{R}^n \rightarrow \mathbb{R}^k$ where $k = \Omega(\epsilon^{-2} \log \frac{m}{\delta})$, such that for all pairs $x_i, x_j \in \mathcal{V}$, $\lVert f(x_i) - f(x_j) \lVert_2^2 = (1 \pm \epsilon) \lVert x_i - x_j \lVert_2^2$.
</div>
</div>


Before we continue, I should mention a neat consequence of the distributional JL lemma. In particular, inner products too are approximately preserved! Concretely, if points $u$ and $v$ on the unit sphere $\mathbb{S}^{n-1}$ then:

$$
\mathbb{P} \Big[ \big\lvert \langle Su, Sv \rangle - \langle u, v \rangle \big\rvert \leq \epsilon \big] \geq 1 - \delta. \tag{4}
$$

That's because, we can write:

$$
4\underbrace{\big( \langle Su, Sv \rangle - \langle u, v \rangle \big)}_{\alpha} = \lVert Su + Sv \rVert_2^2 -\lVert Su - Sv \rVert_2^2 - 4\langle u, v \rangle,
$$

so that, with probability at least $1 - \delta$, we have that:

$$
\begin{aligned}
(1 - \epsilon) &\lVert u + v \lVert_2^2 - (1 + \epsilon) \lVert u - v \rVert_2^2 - 4\langle u, v \rangle \\
&\leq 4\alpha \leq (1 + \epsilon) \lVert u+v \lVert_2^2 - (1 - \epsilon) \lVert u - v \rVert_2^2 - 4\langle u, v \rangle \Rightarrow \\
-\epsilon &(2 \lVert u \rVert_2^2 + 2 \lVert v \rVert_2^2) \leq 4\alpha \leq\epsilon (2 \lVert u \rVert_2^2 + 2 \lVert v \rVert_2^2) \Rightarrow
-\epsilon \leq \alpha \leq \epsilon.
\end{aligned}
$$



More generally, for arbitrary $x, y \in \mathbb{R}^n$, the following statement holds:

$$
\mathbb{P} \Big[ \big\lvert \langle Sx, Sy \rangle - \langle x, y \rangle \big\rvert \leq \epsilon \lVert x \rVert_2 \lVert y \rVert_2  \big] \geq 1 - \delta,
$$

which follows by defining $u=x/\lVert x \rVert_2$ and $v=y/\lVert y\rVert_2$ in Equation (4).

### Example Construction

Let's review a few popular constructions of the random matrix $S$.

#### Gaussians

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
Let $S=\frac{1}{\sqrt{k}}R \in \mathbb{R}^{k \times n}$ where the entries of $R$ are iid standard normal random variables. Then $S$ has the JL property.
</div>
</div>

<br>

<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**: First consider the case of $u \in \mathbb{S}^{n-1}$. Notice that $Su$ is a $k$-dimensional Gaussian random vector whose entries are distributed as $\mathcal{N}(0, 1/k)$. That implies that $\lVert Su \rVert_2^2$ is distributed as $X/k$ where $X$ is a $\chi^2$ random variable with $k$ degrees of freedom. We can then apply [known tail bounds for a $\chi^2$ random variable](/blog/concentration/2024/08/23/subexponential-rvs/#example-sub-exponential-variables) with $k$ degrees of freedom, such as:

$$
\mathbb{P} \Big[ \big\lvert \frac{X}{k} - 1 \big\rvert \geq \epsilon \Big] \leq 2\exp(-\frac{k\epsilon^2}{8}), \tag{5}
$$

so that for our specific case, we have that:

$$
\mathbb{P} \Big[ \big\lvert \lVert Su \rVert_2^2 - 1 \big\rvert \geq \epsilon \Big] \leq 2\exp(-\frac{k\epsilon^2}{8}).
$$

Setting the right-hand-side to $\delta$, we get that:

$$
2\exp(-\frac{k\epsilon^2}{8}) = \delta \Rightarrow k=\Omega(\epsilon^{-2}\log\frac{1}{\delta}).
$$

Extending this to $x \in \mathbb{R}^n$ is simply a matter of defining $u = x/\lVert x \rVert_2$, which gives the JL property.
</div>
</div>


That's so elegant, isn't it? All we need to do is to sample a $k \times n$ matrix whose entries are standard normal random variables, then scale it by $1/\sqrt{k}$, and apply the transformation to our points! That's it! We get a whole new set of $k$-dimensional points where the worst distortion of pairwise distances is bounded.

#### Sparse Transformations

It turns out, while the result is simple to prove, it can be a bit expensive to apply in practice. There are other constructions that give much simpler functions $S$ that are faster to use. Sampling a matrix whose entries are Rademacher random variables, for example, is one example:


<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
Let $S=\frac{1}{\sqrt{k}}R \in \mathbb{R}^{k \times n}$ where the entries of $R$ are iid Rademacher random variables. Then $S$ has the JL property.
</div>
</div>

<br>

<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**: As in the previous proof, it is enough to consider points on the unit sphere. As before, we need to argue about the concentration of $\lVert Su \rVert_2^2$ around its mean, which turns out to be $\lVert u \rVert_2^2=1$:

$$
\mathbb{E}\Big[\lVert Su \rVert_2^2\Big] = \mathbb{E}\Big[ \sum_{i=1}^k \big( \sum_{j=1}^n \frac{1}{\sqrt{k}} r_j u_j \big)^2 \Big] = \frac{1}{k} \sum_{i=1}^k \Big( \lVert u \rVert_2^2 + \sum_{j \neq l} \mathbb{E}[r_{ij} r_{il}] u_j u_l \Big) = \lVert u \rVert_2^2=1,
$$

where $r_i$'s are Rademacher random variables. The final step of bounding the tail probability (as in Equation (5)) is a bit involved, so I'll just note that {% cite Achlioptas2001rademacherJL %} (Lemma 5) proves a bound of the following form:

$$
\mathbb{P}\Big[ \big\lvert \lVert Su \rVert_2^2 - 1 \big\rvert \geq \epsilon \Big] \leq \exp(-\frac{k}{2} (\epsilon^2/2 - \epsilon^3/3) ).
$$

Again setting the right-hand-side to $\delta$ gives the desired result: $k=\Omega(\epsilon^{-2} \log \frac{1}{\delta})$.
</div>
</div>


Finally, we can even make $S$ sparse by drawing the entries of $R$ from the following distribution:

$$
\begin{cases}
\sqrt{3} & \text{with probability } \frac{1}{6}\\
0 & \text{with probability } \frac{2}{3} \\
-\sqrt{3} & \text{with probability } \frac{1}{6}\\ 
\end{cases}.
$$

#### Fast JL Transform (FJLT)

There is even a faster-to-apply family of JL transformations that's simply called Fast JL Transform or FJLT for short. Let's see how it works!

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
Let $S = \frac{1}{\sqrt{k}} PHD$ where: $D$ is a $n\times n$ diagonal matrix whose entries are Rademacher random variables; $H$ is an $n \times n$ matrix with the property that $H^\top H = nI$ and $\max \lvert H_{ij} \rvert \leq 1$; and $P$ is a $k \times n$ matrix where each row is $0$ everywhere except in one uniformly random column. Then $S$ has the JL property with $k = \Omega(\epsilon^{-2} n \log (1/\delta) \log (n/\delta))$.
</div>
</div>


In general, a matrix $H$ with that property is called an *unnormalized bounded orthogonal system*. A concrete example is the Hadamard matrix.


<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**: As before, we do this analysis for a vector $u \in \mathbb{R}^n$ with $\lVert u \rVert_2 = 1$.



**Part I**: First, consider $Y=HD$. We can immediately notice a few things about the product $v = Yu$. First, that its $\ell_2$ norm is simply:

$$
\lVert v \rVert_2^2 = \lVert Yu \rVert_2^2 = u^\top D \underbrace{H^\top H}_{nI} D u = n\lVert u \rVert_2^2 = n.
$$

Second, its $i$-th entry is: $v_i = \sum_{j} H_{ij} D_{jj} u_j$. By following the same method to prove that a Rademacher random variable is sub-Gaussian, we can show that each term $H_{ij} D_{jj} u_j$ is a [sub-Gaussian random variable](/blog/concentration/2024/08/15/subgaussian-rvs/) with constant $H_{ij} u_j$. As a result, $v_i$, being the sum of $n$ sub-Gaussian random variables, is itself sub-Gaussian with constant $\sigma^2 = \sum_j H_{ij}^2 u_j^2 \leq \lVert u \rVert_2^2 = 1$.

Considering $v_i$'s are sub-Gaussian with mean $0$, we know the following Chernoff bound holds:

$$
\mathbb{P}[\lvert v_i \rvert \geq t] \leq 2 \exp \Big( -\frac{t^2}{2\sigma^2} \Big) \leq 2 e^{-t^2/2},
$$

so that if $t=\sqrt{2\log (4n / \delta)}$, then $\mathbb{P}[\lvert v_i \rvert \geq t] \leq \delta/2n$. By the application of the union bound, we can deduce the following, which is known as the *flattening lemma*:

$$
\mathbb{P}[\lVert v \rVert_\infty \geq t] \leq \frac{\delta}{2}.
$$

From here on, let's condition on the flattening lemma being true with probability $1 - \delta/2$. This also implies that, with probability at least $1 - \delta/2$, the $\ell_2$ norm of the $i$-th row of $HD$, $Y_i$, is bounded as $\lVert Y_i \rVert_2^2 \leq 2n\log(4n/\delta)$.



**Part II**: Now, let's consider $\frac{1}{\sqrt{k}} Pv$. We have to show that the deviation of the $\ell_2$ norm of this product from the norm of $u$ is bounded. In other words, we want to assert that:

$$
\mathbb{P}\Bigg[\Big\lvert \lVert \frac{1}{\sqrt{k}} Pv\rVert_2^2 - 1 \Big\rvert \geq \epsilon \Bigg] \leq \delta. \tag{P1}
$$

Let $X=\lVert \frac{1}{\sqrt{k}} Pv \rVert_2^2$. This can be expressed as the sum of random variables $X_i$'s where $X_i = v_i^2/k$ for some uniformly random $i \in [n]$. As a result, $\mathbb{E}[X_i] = 1/k$ and therefore $\mathbb{E}[X] = \sum_{i=1}^k \mathbb{E}[X_i] = 1$.



Let us now bound $\mathop{Var}[X_i]$. Notice that:

$$
\begin{aligned}
\mathop{Var}[X_i] &\leq \mathbb{E}[X_i^2] = \frac{1}{k^2} \mathbb{E}[(u^\top (HD)_i^\top \underbrace{(HD)_i}_{Y_i} u)^2] \\
&= \frac{1}{k^2} \mathbb{E}[u^\top Y_i^\top Y_i Y_i^\top Y_i u]
\leq \frac{1}{k^2} \Big\lVert \mathbb{E}[Y_i^\top Y_i Y_i^\top Y_i] \Big\lVert_2,
\end{aligned}
$$

where $Y_i$ is the $i$-th row of $HD$, $\lVert \cdot \rVert_2$ is the operator norm, and the expectation is with respect to the choice of $i$. We can then write:

$$
\begin{aligned}
\mathbb{E}[Y_i^\top Y_i Y_i^\top Y_i] &= \sum_{i=1}^n \frac{1}{n} Y_i^\top (Y_iY_i^\top) Y_i = \sum \frac{1}{n} Y_i^\top Y_i \lVert Y_i \rVert_2^2 \\
&\leq 2 n \log(\frac{4n}{\delta}) \sum \frac{1}{n} Y_i^\top Y_i = 2 n\log(\frac{4n}{\delta}) I.
\end{aligned}
$$

As a result, $\mathop{Var}[X_i] \leq 2n \log(4n/\delta) / k^2$. The variance of $X$ is therefore bounded by $2n\log(4n/\delta)/k$.



Conditioned on the flattening lemma, $X_i \leq 2 \log (4n / \delta)/k = \tau$, making $X_i$ a bounded random variable, hence sub-exponential. Their sum, $X$, is therefore sub-exponential as well. We can now apply Bernstein's inequality as follows:

$$
\begin{aligned}
\mathbb{P}\Big[ \big\lvert X - 1 \big\rvert \geq \epsilon \Big] &\leq 2\exp(-\frac{\epsilon^2}{2(\mathop{Var}[X] + \tau \epsilon)}) \\
&\leq 2 \exp \Big( - \frac{\epsilon^2}{2 \big( 2n\log(4n/\delta)/k + 2\epsilon \log(4n/\delta)/k \big) } \Big) \\
&=2 \exp(-\frac{k\epsilon^2}{4\log(4n/\delta) (n + \epsilon)} ).
\end{aligned}
$$

When $k = \Omega(\epsilon^{-2} n\log(1/\delta) \log(n/\delta))$, the above probability is at most $\delta/2$.



**Part III**: We have thus far conditioned on the event $\mathcal{E}$ where $\lVert v \rVert_\infty \leq \sqrt{2 \log (4n/\delta)}$ with probability $1 - \delta/2$, and showed that inequality (P1) holds. We must now remove our condition:

$$
\mathbb{P}[\mathcal{E}] \cdot \mathbb{P}\Bigg[ \Big\lvert \lVert \frac{1}{\sqrt{k}} Pv \lVert_2^2 - 1 \Big\rvert \leq \epsilon \; \Big\lvert \mathcal{E} \Bigg] \geq (1- \frac{\delta}{2})^2 \geq 1 - \delta.
$$
</div>
</div>

{% endkatexmm %}