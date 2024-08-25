---
title: "Epsilon Nets"
subtitle: "A common technique in arguments in sketching is to appeal to the concept of an Epsilon Net. This article explains what that is and why it always exists."
layout: default
date: 2024-08-25
keywords: rnla
published: true
category: "Randomized Numerical Linear Algebra"
---

{% katexmm %}
We will be working in the column space of a matrix $A \in \mathbb{R}^{n \times d}$. Without loss of generality, we will assume that $A$ has orthonormal columns; if it doesn't, we can always find an orthonormal basis for its column space by considering its SVD, $A=U\Sigma V^\ast$, and redefining $A:=U$. That is enough because the set $\{ Ax \;|\; x \in \mathbb{R}^d\}$ is the same as the set $\{ Uz \;|\; z \in \mathbb{R}^\rho \}$ where $\rho$ is the rank of $A$.



**Definition**: Define the set $\mathcal{S} = \{ y \;|\; y=Ax, \; \lVert y \rVert_2 = 1 \}$. An $\epsilon$-net is a subset $\mathcal{N} \subset \mathcal{S}$ such that, for every $y \in \mathcal{S}$, there exists a $w \in \mathcal{N}$such that $\lVert y - w \rVert_2 \leq \epsilon$.



### Existence

<div class="callout-yellow" markdown="1">
<div class="with-margin" markdown="1">
**Lemma**: For $\epsilon \in (0, 1)$, there exists an $\epsilon$-net $\mathcal{N}$ for the set $\mathcal{S}$ (defined above) whose size is at most $(1 + 2/\epsilon)^d$.
</div>
</div>


<div class="callout-gray" markdown="1">
<div class="with-margin" markdown="1">
**Proof**: Let's consider the unit sphere in $\mathbb{R}^d: \mathbb{S}^{d-1}$. Instead of finding a net for $\mathcal{S}$, we will argue by finding a net for $\mathbb{S}^{d - 1}$. That's because $A$ provides an isometry when operating on the unit sphere, so that a net for $\mathbb{S}^{d - 1}$ gives us a net for the image of $\mathbb{S}^{d - 1}$ under $A$.



Finding a net over the unit sphere is trivial. We can choose a maximal set $\mathcal{N}^\prime$ such that no two points in the set are within distance $\epsilon$ from each other. This means that balls of radius $\epsilon/2$ centered at points in $\mathcal{N}^\prime$ are disjoint. These balls all fit in the ball with radius $(1 + \epsilon/2)$ centered at the origin. The volume of this enlarged ball is $(1 + \epsilon/2)^d / (\epsilon/2)^d$ times larger than the small balls. As such, there can be at most $(1 + 2/\epsilon)^d$ small balls, implying that $\lvert \mathcal{N}^\prime \rvert \leq (1 + 2/\epsilon)^d$. [This is the covering number of the unit sphere, which is a compact set.]



Now let $\mathcal{N} = \{ y \;|\; y = Ax, \; x \in \mathcal{N}^\prime \}$. We can show by contradiction that, for any $y \in \mathcal{S}$ there must exist a point $w$ in $\mathcal{N}$ such that $\lVert y - w \rVert_2 \leq \epsilon$.$\mathbb{}$ Suppose there exists no such $w$ for some $y=Ax$ for some $x \in \mathbb{S}^{d - 1}$. That implies that, there exists no point $v$ in $\mathcal{N}^\prime$ for which $\lVert x - v \rVert_2 \leq \epsilon$. But that contradicts our construction of $\mathcal{N}^\prime$.
</div>
</div>

{% endkatexmm %}