---
title: "Oblivious Subspace Embedding"
subtitle: "In an earlier post, we saw how <a href='/blog/randomized numerical linear algebra/2024/08/12/jl-transform'>the Johnson-Lindenstrauss Transform</a> can reduce the dimensionality of a set of m points while preserving their pairwise distances. That's neat, but what if we need to preserve the pairwise distance between all (infinitely many) points in an entire subspace? This post introduces the concept of subspace embedding, which answers that question."
layout: default
date: 2024-08-25
keywords: rnla
published: true
category: "Randomized Numerical Linear Algebra"
---

{% katexmm %}

Say we have a matrix $A \in \mathbb{R}^{n \times d}$ where $n \gg d$.  If a linear map $S \in \mathbb{R}^{k \times n}$, where $k \ll n$, has the property that, for all $x \in \mathbb{R}^d$ and some $\epsilon \in (0, 1)$, the following holds:

$$
(1 - \epsilon) \lVert Ax \rVert_2 \leq \lVert SAx \rVert_2 \leq (1+\epsilon) \lVert Ax \rVert_2,
$$

then we call $S$ an $\ell_2$-subspace embedding of the column space of $A$, or simply an $\ell_2$-subspace embedding of $A$. Essentially, $S$ preserves norms of the column space of $A$ up to $\epsilon$ amount of error. To simplify notation, I'll write the above expression as: $\lVert SAx \rVert_2 = (1 \pm \epsilon) \lVert Ax \rVert_2$.


<div class="callout-purple" markdown="1">
<div class="with-margin" markdown="1">
We can, without loss of generality, assume that $A$ has orthonormal columns. That is because, if that's not the case, we can redefine $A:=U$ where $U$ comes from the SVD of $A=U\Sigma V^\ast$. Indeed, the sets $\{ y \;|\; y = Ax,\; x \in \mathbb{R}^d \}$ and $\{ y \;|\; y=Uz, \; z \in \mathbb{R}^\rho \}$ where $\rho=\mathop{rank}(A)$ are the same.

What that means is that, we can simplify the condition to: $\lVert SAx \rVert_2 = (1 \pm \epsilon) \lVert x \rVert_2$.  We can also require that the above hold for unit vectors $x$ due to the linearity of the maps involved. That is in turn equivalent to $\lVert I - A^\ast S^\ast S A \lVert_2 \leq \epsilon$ where $\lVert \cdot \rVert_2$ is the *operator norm* of its argument. To see this:


$$
\begin{aligned}
\lVert S&A x \rVert_2 = (1 \pm \epsilon) \lVert x \rVert_2 \Rightarrow \\
& \big\lvert x^\ast \big( A^\ast S^\ast S A - I \big) x \big\rvert \leq \epsilon \lVert x \rVert_2 \Rightarrow \\
& \sup_{\lVert x \rVert_2 = 1} \big\lvert x^\ast \big( A^\ast S^\ast S A - I \big) x \big\rvert \leq \epsilon \Rightarrow \\
& \sup_{\lVert x \rVert_2 = 1} \big\lvert x^\ast \Sigma x \big\rvert \leq \epsilon \Rightarrow \\
& \sigma_1(A^\ast S^\ast S A - I) \leq \epsilon
\end{aligned}
$$

where we used the SVD $(A^\ast S^\ast SA - I)=U \Sigma V^\ast$ and the fact that $U$ and $V$ have orthonormal columns.
</div>
</div>


What properties should $S$ have? Well, first and foremost, we want $k \ll n$ so that $SA$ has as few rows as possible. That is, we want the smallest subspace of the column space of $A$ that offers a $(1 \pm \epsilon)$-approximation of the Euclidean norm. Additionally, we want to be able to compute $SA$ rather quickly. That often means we prefer a sparser $S$.

### Definition

Suppose that we have a distribution $\Pi$ over the space of all $k \times n$ linear maps, where $k$ is a function of $n, d, \epsilon, \delta$, with the following property: For every $S \sim \Pi$ and any fixed matrix $A \in \mathbb{R}^{n \times d}$, with probability at least $1 - \delta$, $S$ is a $(1 \pm \epsilon)$ $\ell_2$-subspace embedding of $A$. Then we call $\Pi$ an $(\epsilon, \delta)$ oblivious $\ell_2$-subspace embedding.



As a side-note, in addition to oblivious embeddings, we can design sampling-based sketching algorithms that are called **Leverage Score Sampling**. I'll discuss these in a later post.

### Construction

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
**Theorem**: Let $S=\frac{1}{\sqrt{k}}R \in \mathbb{R}^{k \times n}$ where the entries of $R$ are iid standard normal random variables. If $k = \Theta\big(\epsilon^{-2} (d + \log(1/\delta)) \big)$, then $S$ is a $(1 \pm \epsilon)$ $\ell_2$-subspace embedding for any matrix $A \in \mathbb{R}^{n \times d}$ with probability at least $1 - \delta$.
</div>
</div>

<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**: Notice that, we are trying to find an $\ell_2$-subspace embedding for the column space of $A$. We can also focus only on unit vectors in the column space, due to linearity. The idea behind the proof is to first lay an [$\epsilon$-net](/blog/randomized numerical linear algebra/2024/08/25/epsilon-nets) over the column space of $A$ with the property that, a [Johnson-Lindenstrauss Transform](/blog/randomized numerical linear algebra/2024/08/12/jl-transform) for the $\epsilon$-net is a $\ell_2$-subspace embedding for the entire space. Let's get to it.



Consider the set $\mathcal{S} = \{ y \;|\; y = Ax,\; \lVert y \rVert_2 = 1 \}$ and its $1/2$-net $\mathcal{N}$. We can write every $y \in S$ as the series: $y = \sum_{i=0}^{\infty} w_i$ where $w_i$ is a scaled member of $\mathcal{N}$ with the property that $\lVert w_i \rVert \leq 2^{-i}$. To see that that's always possible, take $w_0 := v_0 \in \mathcal{N}$ and write $y = w_0 + (y - w_0)$. Clearly $\lVert y - w_0 \rVert_2 \leq 1/2$ because $w_0 \in \mathcal{N}$. Now write $y = w_0 + w_1 + (y - w_0 - w_1)$ where $w_1 = \lVert y - w_0 \rVert_2 v_1$ for some $v_1 \in \mathcal{N}$. Then:

$$
\lVert y - w_0 - w_1 \rVert_2 = \lVert y - w_0 \rVert_2 \cdot \lVert \frac{y - w_0}{\lVert y - w_0 \rVert_2} - v_1 \rVert_2 \leq \frac{\lVert y - w_0 \rVert_2}{2} \leq \frac{1}{4},
$$

as desired. We can continue to expand the residual as above and reason by induction.



Now consider a JLT for $\mathcal{N}$, $S$, so that with probability at least $1 - \delta$, we have that $\langle Sv, Sv^\prime \rangle = \langle v, v^\prime \rangle \pm \epsilon$ for all pairs of $v, v^\prime \in \mathcal{N}$. If we applied the same map to any $y \in \mathcal{S}$, then we'd get, with probability at least $1 - \delta$:

$$
\begin{aligned}
\langle Sy, Sy \rangle &= \langle \sum S w_i, \sum Sw_i \rangle \\
&= \sum_i \langle Sw_i, Sw_i \rangle + \sum_{i \neq j} \langle Sw_i, Sw_j \rangle \\
&= \sum_i \lVert w_i \rVert_2^2 + \sum_{i \neq j} \langle w_i, w_j \rangle \pm \epsilon \sum_{i, j} \lVert w_i \rVert_2 \lVert w_j \rVert_2 \\
&= \langle \sum w_i, \sum w_i \rangle + \mathcal{O}(\epsilon) \\
&= \langle y, y \rangle \pm \mathcal{O}(\epsilon).
\end{aligned}
$$

So, we have established that a JLT for the $1/2$-net over the unit sphere in the column space of $A$ is an $\ell_2$-subspace embedding for the column space of $A$. We know that $\lvert \mathcal{N} \rvert \leq 5^d$, so that $k = \Omega(\epsilon^{-2} \log (5^d/\delta)) = \Omega(\epsilon^{-2} (d + \log (1/\delta)))$ as desired.
</div>
</div>

{% endkatexmm %}