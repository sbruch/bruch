---
title: "Schatten Norm Approximation"
subtitle: "Randomized Numerical Linear Algebra is a gorgeous field
that embraces chaos and whose results often seem wonderfully magical.
My very first post on this subject explains a beautiful algorithm to
approximate Schatten norms of a matrix."
layout: default
date: 2024-07-29
keywords: rnla
published: true
category: "Randomized Numerical Linear Algebra"
---

{% katexmm %}

In Numerical Linear Algebra, we are often interested in understanding the effect a linear map
$A \in \mathbb{R}^{m \times n}$ has on a (real) vector space. A map can rotate, shear,
or dilate the space; flip our orientation; expand the space or compress it.
Whatever $A$ does, to say anything intelligent about its behavior, 
we need to be able to *quantify* its effect too!

That's where **matrix norms** come into the picture. You may be familiar with
induced norms such as $\lVert A \rVert_2$ (maximum singular value),
$\lVert A \rVert_1$ (maximum absolute column sum), $\lVert A \rVert_\infty$
(maximum absolute row sum), which capture a lot of information.
$\lVert A \rVert_2$, for example, tells you the largest amount of "stretching"
$A$ is capable of as it transforms your space---to see this, it might help
to go back to the definition of these norms:
$\lVert A \rVert_p = \sup_{\lVert u \rVert_p = 1} \lVert Au \rVert_p$.

In this post, we'll work with a different kind of matrix norm: the Schatten $p$-norms,
defined as follows:
$$\lVert A \rVert_p^p = \sum_{1 \leq i \leq \rho} \sigma^p_i(A),$$
where $\rho$ is the rank of $A$ and $\sigma_i(\cdot)$ is its $i$-th  singular value.
When $p=2$, for example, this is the familiar Frobenius norm; $p=\infty$
gives the spectral norm, which coincides with the induced $2$-norm above.
From here on, when I write $\lVert A \rVert_p$, I'm referring to the
Schatten $p$-norm of $A$, and not the induced norm.

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
One of the interesting facts about Schatten $p$-norm is that they are invariant under
rotation! That is, if I applied an orthogonal transformation
$R \in \mathbb{R}^{m \times m}$ to $A$, $\lVert A \rVert_p$ wouldn't change.
Let's see why:
$$\begin{aligned}
\lVert RA \rVert_p^p &= \lVert R U \Sigma V^\ast \rVert_p^p
    = \lVert \underbrace{Q \Sigma V^\ast}_{B} \rVert_p^p
    = \sum_{1 \leq i \leq \rho} \sigma^p_i = \lVert A \rVert_p^p,
\end{aligned}$$
where $\cdot^\ast$ is the transpose of its argument, and $A=U \Sigma V^\ast$
is the Singular Value Decomposition of $A$. Indeed, because $R$ and $U$ are
orthonormal matrices, $Q=RU$ too will be orthonormal, rendering $Q\Sigma V^\ast$
the SVD of another matrix $B$ whose singular values match those of $A$'s.
You can use a similar argument for transformations of $A$ from the right
(i.e., $AR$ where $R \in \mathbb{R}^{n \times n}$).
</div>
</div>

We see why Schatten $p$-norms are useful. Computing $\lVert A \rVert_p$ is also easy, if we are
allowed to apply SVD to $A$ and extract its singular values. But, let's say,
we don't have access to $A$; it's being kept a secret by someone---let's call
this person the "oracle"---and they won't let us apply SVD to the matrix!
(It doesn't have to be for nefarious reasons, it could
simply be because $A$ is just way too large!)
But we *are* allowed to transform any vector $u$ with $A$ and obtain the output, $Au$.

In what seems so magical, we can, in fact, approximate $\lVert A \rVert_p$ arbitrarily
accurately by repeatedly asking for $Au$ for random vectors $u$! The rest of this
post explains how.

### Approximating $\lVert A \rVert_p$ for symmetric $A$

Suppose for now that $A \in \mathbb{R}^{d \times d}$ and that it is
symmetric---we will come back to the more general case in a bit.
We want to obtain, with probability at least $1 - \delta$,
a $(1+\epsilon)$-approximation of $\lVert A \rVert_p$
for some valid choices of $\delta$ and $\epsilon$.
In other words, we want to find a quantity $X$ such that:
$$\mathbb{P}\Big[ (1 - \epsilon) \lVert A \rVert_p^p \leq X
   \leq (1+ \epsilon) \lVert A \rVert_p^p \Big] \geq 1 - \delta.$$
Also, we'll just focus on $\lVert A \rVert_p^p$, rather than $\lVert A \rVert_p$
because it makes the math less cluttered.

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
Now to the main steps of the algorithm proposed by  {% cite li2014soda %}.
Repeat the following process $T = C / \epsilon^2$ times (for some constant
$C$ to be defined later), indexed by $t$:

* Initialize $u := u^{(0)} \sim \mathcal{N}(\mathbf{0}, I_{d \times d})$,
  drawn from an isotropic Gaussian distribution (whose covariance is
  the identity). This amounts to a vector whose $d$ entries are independent random
  variables sampled from the standard normal distribution.
* Let $u^{(i+1)} \leftarrow Au^{(i)}$ by asking the oracle for the computation,
and repeat that process $\lceil p/2 \rceil$ times.
* Take note of $X^{(t)} \leftarrow \langle
  u^{(\lfloor p/2 \rfloor)}, u^{(\lceil p/2 \rceil)} \rangle$, 
  where $\langle \cdot, \cdot \rangle$ is inner product.

In the end, output the following quantity:
$$
\frac{1}{T} \sum_{1 \leq t \leq T} X^{(t)},
\tag{1}
$$
as the desired approximation of $\lVert A \rVert_p^p$.
</div>
</div>

That was just the (sweet and simple) algorithm. We left out what $C$ should be
and, more importantly, why this procedure even works! So let's get to it!

Before we do though, let's quickly consider the complexity of the algorithm.
If the number of non-zero entries in $A$ is $\mathop{nnz}(A)$, then
each iteration of the algorithm simply runs in $\mathcal{O}(p \cdot \mathop{nnz}(A))$ time.
We repeat the algorithm $T$ times, so that its complexity adds up to
$\mathcal{O}(p \cdot \mathop{nnz}(A) \cdot \epsilon^{-2})$.

### Why does that even work?

We assumed that $A$ was symmetric. That implies that it can be decomposed
into the following form: $U \Sigma U^\ast$ for some orthogonal matrix $U$
(so that $UU^\ast=I$)
and diagonal matrix $\Sigma$. What happens when we take the inner product
inside the sum in Equation (1)? It reduces to:
$$\begin{aligned}
X := \langle u^{(\lfloor p/2 \rfloor)}, u^{(\lceil p/2 \rceil)} \rangle &=
\langle A^{\lfloor p /2 \rfloor} u, A^{\lceil p /2 \rceil} u \rangle \\
&= u^\ast A^p u \\ 
&= u^\ast \Sigma^p u \\
&= \sum_{1 \leq k \leq \rho} \sigma_k^p u_k^2. \tag{2}
\end{aligned}$$
where we used the definition of $u^{(i)}$ to get the first equality,
the definition of inner product for the second,
and the decomposition of $A$ for the third.
In the last expression, $\rho$ is the rank of $A$ and subscripts indicate
specific coordinates.

The output of the algorithm is the mean of $X^{(t)}$'s.
Considering each sample $X^{(t)}$ gives Equation (2), and that
$\mathbb{E}[u_k^2] = \mathop{Var}[u_k] = 1$,
the expected value of $X$ is
an unbiased estimate of $\sum \sigma_k^p$, which is exactly $\lVert A \rVert_p^p$.

#### Analysis of error

I promised we'll get to $C$ at some point. Now that we have verified that
the result of the algorithm is, in fact, an unbiased estimate of $\lVert A \rVert_p^p$,
it's time to fill that gap.

Error analysis often begins with a derivation of or a bound on the variance.
In this case, we should be looking at our random variable $X$ from Equation (2).
It shouldn't be too hard considering we're dealing with *iid* Gaussian variables.
So let's do it:
$$\begin{aligned}
\mathop{Var}[X] &\leq \mathbb{E}[X^2] = \mathbb{E}_u[(u^\ast \Sigma^p u)^2] \\
&= \sum_{j, k} \sigma^p_j \sigma^p_k \; \mathbb{E}[u_j^2 u_k^2] \\
&= \sum_{j} \sigma^{2p}_j \; \underbrace{\mathbb{E}[u_j^4]}_3 +
    \sum_{j \neq k} \sigma^p_j \sigma^p_k \; \underbrace{\mathbb{E}[u_j^2] \mathbb{E}[u_k^2]}_1 \\
&= 3\sum_{j} \sigma^{2p}_j + \sum_{j \neq k} \sigma^p_j \sigma^p_k \\
&\leq 4 \lVert A \rVert_p^{2p}.
\end{aligned}$$

So we know the variance of each sample of $X$ is less than that quantity. Great.
We sample $X$ some $T$ times. Because the $T$ samples are independent,
the variance scales down by a factor of $T$.

Now we can finally get to our approximation error.
We need to understand the maximum probability by which $X$ deviates from
$\lVert A \rVert_p^p$ by a factor of at most $(1 + \epsilon)$.
For that, we will appeal to standard results from concentration of measure.
In particular, we simply apply Chebyshev's inequality:
$$
\mathbb{P}\Big[ \big\lvert X - \mathbb{E}[X] \big\rvert \geq \epsilon \mathbb{E}[X] \Big]
\leq \frac{\mathop{Var}[X]}{\epsilon^2 \mathbb{E}^2[X]},
$$
to the random variable $Z = \sum X^{(t)} / T$. That gives us:
$$
\mathbb{P}\Big[ \big\lvert Z - \lVert A \rVert_p^p \big\rvert \geq
\epsilon \lVert A \rVert_p^p \Big]
\leq \frac{4 \lVert A \rVert_p^{2p}}{T \epsilon^2 \lVert A \rVert_p^{2p}} =
\frac{4}{T \epsilon^2} = \frac{4}{C}.
$$

We wanted the probability above to be at most $\delta$. Setting
the right-hand-side to $\delta$, we get that $C = 4/\delta$.
And there we have it!

### The general case

We saw what we need to ask of the oracle to obtain an approximation of
$\lVert A \rVert_p$ when $A$ is *symmetric*: We initialize a Gaussian vector,
and repeatedly ask the oracle to return the matrix-vector product, which
we then feed back to the oracle. Pretty straightforward. But what about the
more general case where $A \in \mathbb{R}^{m \times n}$?

It turns out, if we formed the following block matrix:
$$
A^\prime = \begin{bmatrix}
0 & A^\ast\\
A & 0
\end{bmatrix},
$$
then $\lVert A^\prime \rVert_p = 2^{1/p} \lVert A \rVert_p$. It's fairly
easy to show this if you rewrite $A^\prime$ in terms of
the SVD of $A=U\Sigma V^\ast$:
$$
A^\prime = \begin{bmatrix}
0 & V \Sigma U^\ast \\
U \Sigma V^\ast & 0
\end{bmatrix} =
\begin{bmatrix}
V & 0 \\
0 & U
\end{bmatrix}
\begin{bmatrix}
0 & \Sigma U^\ast \\
\Sigma V^\ast & 0
\end{bmatrix} =
\begin{bmatrix}
V & 0 \\
0 & U
\end{bmatrix}
\begin{bmatrix}
\Sigma & 0 \\
0 & \Sigma
\end{bmatrix}
\begin{bmatrix}
0 & U^\ast \\
V^\ast & 0
\end{bmatrix}.
$$

We would have to update the algorithm with this logic to ensure
it works for a general matrix $A$.

### Closing Remarks

I hope you found the elegance and simplicity of this randomized algorithm as
compelling as I did and still do---ten years since its publication.
There are numerous other algorithms that feel magical for norm approximation and
other operations on vectors and matrices in this rich literature
that I encourage you to study. There are a number of excellent monographs on the subject
such as {% cite Woodruff_2014 %},
{% cite murray2023randomizednumericallinearalgebra %},
{% cite martinsson2021randomizednumericallinearalgebra %}
and references therein. I hope to highlight some of the more representative
algorithms in future posts.

{% endkatexmm %}