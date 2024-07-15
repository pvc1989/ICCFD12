<!-- 
layout: home
title: ICCFD12
 -->

# 1. Introduction

## 1.1. High-order schemes for unstructured grids

## 1.2. Shock capturing techniques for DG and FR

## 1.3. Outline of the proposed artificial viscosity

# 2. Methods

## 2.1. DG and FR — a quick review and some remarks

### 2.1.1. Element-wise polynomial approximations

To solve a 2d conservation law (system)
$$
\partial_{t}\,u+\partial_{\vec{r}}\vdot\vec{f}=0,\quad\partial_{\vec{r}}\vdot\vec{f}=\partial_{x}\,f^{x}+\partial_{y}\,f^{y},
$$

one first introduce a coordiante map

$$
\underbrace{(x,y)}_{\vec{r}}\mapsto\underbrace{(\xi,\eta)}_{\vec{\rho}}
\implies
\begin{bmatrix}\partial_{\xi}\,\phi\\
\partial_{\eta}\,\phi
\end{bmatrix}=\underbrace{\begin{bmatrix}\partial_{\xi}\,x & \partial_{\xi}\,y\\
\partial_{\eta}\,x & \partial_{\eta}\,y
\end{bmatrix}}_{\underline{J}}\begin{bmatrix}\partial_{x}\,\phi\\
\partial_{y}\,\phi
\end{bmatrix}=\begin{bmatrix}\partial_{\xi}\,\vec{r}\\
\partial_{\eta}\,\vec{r}
\end{bmatrix}\vdot\partial_{\vec{r}}\,\phi,
$$

in which

$$
U=u\,\underbrace{\det(\underline{J})}_{J},\quad\begin{bmatrix}F^{\xi}\\
F^{\eta}
\end{bmatrix}=J\,\underbrace{\begin{bmatrix}\partial_{x}\,\xi & \partial_{y}\,\xi\\
\partial_{x}\,\eta & \partial_{y}\,\eta
\end{bmatrix}}_{\underline{J}^{-1}}\begin{bmatrix}f^{x}\\
f^{y}
\end{bmatrix}=J\begin{bmatrix}\partial_{\vec{r}}\,\xi\\
\partial_{\vec{r}}\,\eta
\end{bmatrix}\vdot\vec{f}.
$$

Then, one might use either an orthonormal (modal) expansion or a Lagrange (nodal) interpolation for either $$u$$ or $$U\equiv Ju$$ on the $$j$$th element:

$$
u(\vec{r},t)\approx u_{j}^{h}(\vec{r},t)=\begin{cases}
\sum_{n=1}^{N}\hat{u}_{j,n}(t)\,\phi_{j,n}(\vec{\rho}), & \text{Lagrange interpolation},\\
\sum_{m=1}^{M}\tilde{u}_{j,m}(t)\,\psi_{j,m}(\vec{r}), & \text{orthonormal expansion},
\end{cases}
$$

### 2.1.2. Semi-discretized systems from DG and FR

A DG scheme might be formulated either in physical coordinates

$$
\int_{E_j}\phi_n\pdv{u_j^h}{t}
=\int_{E_j}\vec{f_j}^\mathrm{D}\vdot\grad\phi_n
-\oint_{\partial E_j}f^{I}\,\phi_n,\quad \forall n\in\{1,\dots,N\},
$$

or in parametric coordinates

$$
\int_{\mathcal{E}_j}\phi_{n}\pdv{U_j^{h}}{t}=\int_{\mathcal{E}_j}\vec{F}_j^{D}\vdot\grad\phi_{n}-\oint_{\partial\mathcal{E}_j}F^{I}\,\phi_{n},\quad\forall n\in\{1,\dots,N\},
$$

A FR scheme is usually given in parametric coordinates

$$
\begin{aligned}\dv{\hat{U}_{j,m,n}}{t} & =\frac{\partial F_{j}^{\xi,D}(\xi_{m},\eta_{n})}{\partial\xi}+\sum_{a=\pm1}\left[F^{I}-F_{j}^{\xi,D}\right]_{\xi=a,\eta_{n}}\dv{g_{a}(\xi_{m})}{\xi}\\
 & +\frac{\partial F_{j}^{\eta,D}(\xi_{m},\eta_{n})}{\partial\eta}+\sum_{b=\pm1}\left[F^{I}-F_{j}^{\eta,D}\right]_{\xi_{m},\eta=b}\dv{g_{b}(\eta_{n})}{\eta},
\end{aligned}
$$

in which $$g$$'s are correction functions satisfying

$$
g_{+1}(+1)=1,\quad
g_{+1}(-1)=0,\quad
g_{-1}(-\xi)=g_{+1}(+\xi),\quad \forall\xi\in[-1,1],
$$

and look like

![](./fr/HuynhLumping.svg)

The original flux (disconstinuous on element interfaces) is reconstructed to be $$ C_0 $$ constinuous in global:

![](./fr/FRonLegendreRoots.svg)

### 2.1.3. Direct viscous fluxes on element interfaces

Common fluxes on element interfaces:

- Convection: exact or approximate Riemann solvers.
- Diffusion: BR1/BR2, local DG, ***direction DG***, ...

In this work, we use the last one (DDG):
$$
\begin{bmatrix}\partial_{x}\,u\\
\partial_{y}\,u
\end{bmatrix}_{\partial E}=\beta_{0}\,\Delta^{-1}\begin{bmatrix}n_{x}\,u\\
n_{y}\,u
\end{bmatrix}_{\mathrm{R}-\mathrm{L}}+\begin{bmatrix}\partial_{x}\,u\\
\partial_{y}\,u
\end{bmatrix}_{(\mathrm{R}+\mathrm{L})/2}+\beta_{1}\,\Delta
\begin{bmatrix}
n_x\,\partial^2_{xx}\,u + n_y\,\partial^2_{xy}\,u\\
n_x\,\partial^2_{yx}\,u + n_y\,\partial^2_{yy}\,u
\end{bmatrix}_{\mathrm{R}-\mathrm{L}},
$$

If, the interpolation is applied to $u$, the second-order derivatives are

$$
\begin{aligned}\begin{bmatrix}\partial_{xx}^{2}\,u & \partial_{xy}^{2}\,u\\
\partial_{yx}^{2}\,u & \partial_{yy}^{2}\,u
\end{bmatrix} & =\underbrace{\underline{J}^{-1}\begin{bmatrix}\partial_{\xi}\\
\partial_{\eta}
\end{bmatrix}}_{\begin{bmatrix}\partial_{x}\\
\partial_{y}
\end{bmatrix}}\underbrace{\left(\mathinner{\begin{bmatrix}\partial_{\xi}\,u & \partial_{\eta}\,u\end{bmatrix}}\underline{J}^{-T}\right)}_{\begin{bmatrix}\partial_{x}\,u & \partial_{y}\,u\end{bmatrix}}\\
 & =\underline{J}^{-1}\mathinner{\begin{bmatrix}\partial_{\xi}\mathinner{\begin{bmatrix}\partial_{\xi}\,u & \partial_{\eta}\,u\end{bmatrix}}\\
\partial_{\eta}\mathinner{\begin{bmatrix}\partial_{\xi}\,u & \partial_{\eta}\,u\end{bmatrix}}
\end{bmatrix}}\underline{J}^{-T}+\underline{J}^{-1}\begin{bmatrix}\mathinner{\begin{bmatrix}\partial_{\xi}\,u & \partial_{\eta}\,u\end{bmatrix}}\partial_{\xi}\,\underline{J}^{-T}\\
\mathinner{\begin{bmatrix}\partial_{\xi}\,u & \partial_{\eta}\,u\end{bmatrix}}\partial_{\eta}\,\underline{J}^{-T}
\end{bmatrix},
\end{aligned}
$$

in which the derivatives of the Jacobian matrix

$$
\pdv{\underline{J}^{-T}}{\xi}=\left(-\underline{J}^{-1}\,\pdv{\underline{J}}{\xi}\,\underline{J}^{-1}\right)^{T},\quad
\pdv{\underline{J}^{-T}}{\eta}=\left(-\underline{J}^{-1}\,\pdv{\underline{J}}{\eta}\,\underline{J}^{-1}\right)^{T},
$$

should be precomputed and cached on each flux point.

If the interpolation is applied to $U\equiv Ju$, the second-order derivatives are more complex:

$$
\begin{aligned}\begin{bmatrix}\partial^2_{xx}\,u & \partial^2_{xy}\,u\\
\partial^2_{yx}\,u & \partial^2_{yy}\,u
\end{bmatrix} & =\underbrace{\underline{J}^{-1}\begin{bmatrix}\partial_{\xi}\\
\partial_{\eta}
\end{bmatrix}}_{\begin{bmatrix}\partial_{x}\\
\partial_{y}
\end{bmatrix}}\underbrace{\left(\mathinner{\begin{bmatrix}\partial_{\xi}\,U & \partial_{\eta}\,U\end{bmatrix}}\underline{J}^{-T}\,J^{-1}-\mathinner{\begin{bmatrix}\partial_{\xi}\,J & \partial_{\eta}\,J\end{bmatrix}}\underline{J}^{-T}\,J^{-2}\,U\right)}_{\begin{bmatrix}\partial_{x}\,u & \partial_{y}\,u\end{bmatrix}}\\
 & =\underline{J}^{-1}\mathinner{\begin{bmatrix}\partial_{\xi}\,\partial_{\xi}\,U & \partial_{\xi}\,\partial_{\eta}\,U\\
\partial_{\eta}\,\partial_{\xi}\,U & \partial_{\eta}\,\partial_{\eta}\,U
\end{bmatrix}}\underline{J}^{-T}\,J^{-1}+\cdots,
\end{aligned}
$$

$$
\begin{aligned}\partial_{\xi}\left(\mathinner{\begin{bmatrix}\partial_{\xi}\,U & \partial_{\eta}\,U\end{bmatrix}}\underline{J}^{-T}\,J^{-1}\right) & =\mathinner{\begin{bmatrix}\partial_{\xi}\,\partial_{\xi}\,U & \partial_{\xi}\,\partial_{\eta}\,U\end{bmatrix}}\underline{J}^{-T}\,J^{-1}\\
 & +\mathinner{\begin{bmatrix}\partial_{\xi}\,U & \partial_{\eta}\,U\end{bmatrix}}\left(\pdv{\underline{J}^{-T}}{\xi}\,J^{-1}+\underline{J}^{-T}\left(\partial_{\xi}\,J^{-1}\right)\right),
\end{aligned}
$$

$$
\begin{aligned}\partial_{\xi}\left(\mathinner{\begin{bmatrix}\partial_{\xi}\,J & \partial_{\eta}\,J\end{bmatrix}}\underline{J}^{-T}\,J^{-2}\,U\right) & =\mathinner{\begin{bmatrix}\partial_{\xi}\,\partial_{\xi}\,J & \partial_{\xi}\,\partial_{\eta}\,J\end{bmatrix}}\underline{J}^{-T}\,J^{-2}\,U\\
 & +\mathinner{\begin{bmatrix}\partial_{\xi}\,J & \partial_{\eta}\,J\end{bmatrix}}\left(\pdv{\underline{J}^{-T}}{\xi}\,J^{-2}\,U+\underline{J}^{-T}\,\pdv{J^{-2}}{\xi}\,U+\underline{J}^{-T}\,J^{-2}\,\pdv{U}{\xi}\right),
\end{aligned}
$$

in which, the derivatives of Jacobian determinant

$$
\pdv{J^{-1}}{\xi}=-J^{-2}\,\pdv{J}{\xi},\quad\pdv{J^{-2}}{\xi}=-2J^{-3}\,\pdv{J}{\xi},
$$

$$
\pdv{J}{\xi}=\pdv{\det(\underline{J})}{\xi}=\det(\underline{J})\tr(\underline{J}^{-1}\,\pdv{\underline{J}}{\xi})=J\tr(\underline{J}^{-1}\,\pdv{\underline{J}}{\xi}),
$$

should also be cached.

## 2.2. Artificial viscosity based on energy dissipation

### 2.2.1. Quasi-linear artificial viscosity

Take the scalar case in 1d space as an example:

$$
\pdv{u}{t}+\pdv{f(u)}{x}=\pdv{x}\qty({\color{red}\nu(u)}\pdv{u}{x}),\quad \nu(u) \ge 0.
$$

Generalizing to systems:

- For a one-dimensional system, the procedure could be applied to either each conservative variable or each characteristic variable.

- For a multi-dimensional problem, we just apply the proposed procedure to each conservative variable, i.e.

$$
\pdv{t}\mqty[u_1\\ \vdots \\ u_K]+\grad\vdot\mqty[\vb*{f}_1\\ \vdots \\ \vb*{f}_K]
=\grad\vdot\mqty[{\color{red}\nu_1}\grad{u}_1\\ \vdots \\ {\color{red}\nu_K}\grad{u}_K],
$$

Assume the viscosity value is a constant across the entire domain.

The semi-discretized ODE system can be written as

$$
\dv{t}\ket{\hat{u}_{j}}=(\text{inviscid terms})+\nu\left(\underline{D}_{j}\ket{\hat{u}_{j}}+\underline{E}_{j}\ket{\hat{u}_{j-1}}+\underline{F}_{j}\ket{\hat{u}_{j+1}}\right),
$$

Without loss of generality, we hereforth only consider the nodal interpolation in which solution points are also Gaussian quadrature points.

### 2.2.2. Dissipation rate of kinetic energy

One may now define the ***kinetic energy*** on $E_j$ to be

$$
K_{j}\equiv\int_{E_{j}}\frac{1}{2}\qty(u_{j}^{h})^{2}\approx\frac{1}{2}\sum_{n=1}^{N}w_{j,n}\,\hat{u}_{j,n}\,\hat{u}_{j,n}=\frac{1}{2}\underbrace{\begin{bmatrix}\hat{u}_{j,1} & \cdots & \hat{u}_{j,N}\end{bmatrix}}_{\bra{\hat{u}_{j}}}\underbrace{\begin{bmatrix}w_{j,1}\\
 & \ddots\\
 &  & w_{j,N}
\end{bmatrix}}_{\underline{W}_{j}}\underbrace{\begin{bmatrix}\hat{u}_{j,1}\\
\vdots\\
\hat{u}_{j,N}
\end{bmatrix}}_{\ket{\hat{u}_{j}}},
$$

The ***dissipation rate*** of $K_j$ can easily be derived, which is

$$
\dv{K_j}{t}=\bra{\hat{u}_{j}}\underline{W}_{j}\dv{t}\ket{\hat{u}_{j}}=(\text{inviscid terms})+\nu\bra{\hat{u}_{j}}\underline{W}_{j}\,\underline{D}_{j}\ket{\hat{u}_{j}}+(\text{inter-cell viscous terms}),
$$

By ignoring the inviscid terms and inter-cell viscous terms, one gets

$$
\nu_{j}=\frac{\Delta K_{j}}{-G_{j}\,\tau}
\impliedby\dv{t}K_{j}=(\text{inviscid terms})+\nu\underbrace{\bra{\hat{u}_{j}}\underline{W}_{j}\,\underline{D}_{j}\ket{\hat{u}_{j}}}_{G_{j}}+(\text{inter-cell viscous terms}),
$$

### 2.2.3. Oscillation energy estimations

The last question is how to evaluate the ***oscillation energy*** $\Delta K_j$, which quantitatively measures the numerical oscillation on $E_j$.

In one-dimensional cases, we define it to be

$$
\Delta K_j = \int_{x_{j-1/2}}^{x_{j}}\left(u_{j}^{h}-u_{j-1}^{h}\right)^{2}\dd{x}+\int_{x_{j}}^{x_{j+1/2}}\left(u_{j}^{h}-u_{j+1}^{h}\right)^{2}\dd{x},
$$

which is the square of the $L_2$-norm of difference between the solution on $E_j$ and those on its immediate neighbors.

The integrals are as cheap as weighted summations of nodal values if solution points are also Gaussian quadrature points.

The extrapolations are also cheap, since the coordinate map in this case is usually linear and therefore trivial.

However, when it is generalized to *multi-dimensional* cases, extrapolations using nodal expansions in parametric coordinates are generally expensive, since it incurs solving a non-linear algebraic equation

$$
\mqty[x \\ y]_\text{query} = \mqty[x^h(\xi, \eta) \\ y^h(\xi, \eta)],
$$
for the value of $(\xi,\eta)$ from the given $(x,y)$ on a neighboring cell.

Even worse, the solution of this equation might not exist, even for a 4-node quadrangular element:

![](./coordmap.svg)

It is an open problem to design an oscillation measure requiring no extrapolations.

In the current work, we tried the following surface integral of value jumps

$$
\oint_{\partial E_j} \qty(u^h_{j}-u^h_{j'})^2,\quad
j'\in\text{ neighbors of } E_j,
$$

Further tests on more challenging problems are still in progress.

# Results

1. aaaa
2. bbbb
3. ...



## Test Table

| sod                  | Lax                  |
| -- | -- |
| hhhh | mmmm |



| sod                  | Lax                  |
| -- | -- |
| ![](./animation.gif) | ![](./animation.gif) |



## Test Figure

<img src="./animation.gif" style="zoom:50%;" />



### Riemann

<img src="./assets/riemann.png" style="max-width:1080; max-height:200;"/>

# Conclusions

The proposed artificial viscosity

- passes standard cases for testing shock capturing methods.
- is sufficiently large near physical discontinuities, which succeeds in suppressing numerical oscillations.
- is negligible in smooth regions, which maintains the high-order accuracy of the DG or FR solution.

Further improvements might include

- more accurate estimation of oscillation energy.
- more smooth distribution of the artificial viscosity.
- incorporating positivity-preserving mechanisms.
