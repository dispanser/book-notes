---
title:  Lecture Notes "Matrix Methods"
author: Thomas Peiselt
date: 01.03.2020
geometry: margin=2cm
abstract: |
  Matrix Methods for Data Analysis, Signal Processing, and Machine Learning (Gilbert Strang)

  Taken from course [18.065 MIT Open Course Ware](https://ocw.mit.edu/courses/mathematics/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/index.htm)
---

# Introduction

These are my lecture notes following Gilbert Strangs _Matrix Methods for Data
Analysis, Signal Processing, and Machine Learning_ [18.065 MIT Open Course Ware](
https://ocw.mit.edu/courses/mathematics/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/index.htm)

# Lecture 1

$$\mathbf{A} = \begin{bmatrix} a_{11} & a_{12} & \dots  & a_{1n} \\
		                 a_{21} & a_{22} & \dots  & a_{2n} \\
						 \vdots & \vdots & \ddots  & \vdots \\
		                 a_{m1} & a_{m2} & \dots  & a_{mn}
\end{bmatrix}, \quad a_{ij} \in \mathbb{R}$$

## Interpretations of Matrix Multiplication

Consider the matrix-vector product $Ax, x = [x_1, x_2, \dots x_n]^T$. There are multiple interpretations:

1. A _column of dot products_ of rows of $A$ with $x$:
$$ \mathbf{Ax} = \begin{bmatrix}
  row_1(A)^T x \\
  row_2(A)^T x \\
  \vdots \\
  row_n(A)^T x
\end{bmatrix}$$
2. A _linear combinations of columns_ of $A$ with $x_i$:
$$\mathbf{Ax} = \sum_{i=1}^n x_i col_i(A)$$

Consider the matrix multiplication $AB, A \in \mathbb{R}^{m \times n}, B \in \mathbb{R}^{n \times p}$.
Next to the basic interpretation of this
representing a matrix of dot products from rows of $A$ and columns of $B$, we
can also consider this a combination of column $k$ of $A$ by row $k$ of $B$.
Each of these outer products produces a rank 1 matrix of the proper shape, the
final multiplication result is the sum of all these $n$ outer products.


## Column Space, Rank of a Matrix

$C(A)$ is the _column space_ (also: range or image) of $A$: $C(A) = Ax, \quad x \in \mathbb{R}^n$
   
In words: the column space is the image of applying the linear transformation
defined by $A$ to all possible vectors $x$.

$rank(A) = dim(C(A))$

The rank of a matrix is the number of independent columns and rows of $A$.

$$\mathbf{A}
= \begin{bmatrix}
	1 & 3 & 8 \\
	1 & 3 & 8 \\
	1 & 3 & 8 
\end{bmatrix}
=  \begin{bmatrix} 1 \\ 1 \\ 1\end{bmatrix} \begin{bmatrix}1 & 3 & 8\end{bmatrix}$$ 

Any matrix $\mathbf{A}$ can factored it into its column space $\mathbf{C}
\in \mathbb{R}^{m \times r}$ and some matrix $\mathbf{R} \in \mathbb{R}^{r
\times n}$ 

$$\mathbf{A}
= \begin{bmatrix}
	2 & 1 & 3 \\
	3 & 1 & 4 \\
	5 & 7 & 12 
\end{bmatrix}
=  \begin{bmatrix} 2 & 1 \\ 3 & 1 \\ 5 & 7\end{bmatrix}
	\begin{bmatrix}1 & 0 & 1 \\ 0 & 1 & 1\end{bmatrix}$$ 

$C$ contains columns from $A$, and $R$ is (by construction) in reduced
row-echolon form, containing the identity in the left to produce the columns of
$C$, followed by other columns producing the remaining columns of $A$ as linear
combinations of columns from $C$.

## Fundamental Subspaces of $A$

1. column space $C(A)$, $dim = r$
2. row space $C(A^T)$, $dim = r$
3. null space $N(A) = x: Ax = 0$, $dim = n - r$
4. null space $N(A^T)$, $dim = m - r$

The null space is orthogonal to the associated row or column space.

# Lecture 2

## Matrix Factorizations

1. $A = LU$: gauss elimination
   - $L$: lower triangular
   - $U$: upper triangular
   - $A = col_1 row_1 + \begin{bmatrix} 0 & 0 & \dots & 0 
   \\ \hdots & a_{22}' & \dots & a_{2n}' 
   \\ 0  & a_{n2}' & \dots & a_{nn}' \end{bmatrix}$
   - recursively take care of the first (row $\times$ column) with a rank 1 matrix
2. $A = QR$ Gram-Schmidt process
3. $S = Q\Lambda Q^T$: spectral decomposition, eigendecomposition

   A sum of rank $1$ matrices: $S = (Q\Lambda)Q^T = \sum_{i=1}^n \lambda_i q_i q_i^T$
4. $A = X\Lambda X^{-1}$ 
5. $A = U\Sigma V^T$: singular value decomposition

# Lecture 3

Properties of orthogonal matrices, usually named $Q$:

- it's columns are orthogonal
- columns are unit vectors
- $Q^TQ = I$
- $QQ^T = I$ if $Q$ is square.
- $||Qx|| = ||x||$ (because $||Qx||^2 = (Qx)^T(Qx) = x^TQ^TQ^x = x^Tx = ||x||^2$)

## Interesting orthogonal matrices

1. Rotation matrix: $Q = \begin{bmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{bmatrix}$

1. Reflection at $\theta/2$ matrix: $Q = \begin{bmatrix}\cos\theta & \sin\theta \\ \sin\theta & -\cos\theta \end{bmatrix}$ 

1. Householder reflections: $H = I - 2uu^T$ for a _unit vector_ $u$.
   - $H$ is symmetric, so $H^TH = HH$
   - to proof orthonormality, show $H^TH = I$:

     $H^TH = HH = I - 4uu^T + 4uu^Tuu^T = I$ (because $u^Tu = 1$, unit vector)
1. Hadamard matrices: $H_2 = \begin{bmatrix}1 & 1 \\ 1 & -1 \end{bmatrix}$
   - this extends to $H_4, H_8, H_{12}$ via some kind of recursion scheme
   - $H_4 = \begin{bmatrix}H_2 & H_2 \\ H_2 & -H_2\end{bmatrix}$
   - $H_8 = \begin{bmatrix}H_4 & H_4 \\ H_4 & -H_4\end{bmatrix}$
1. Wavelets
1. Eigenvectors of symmetric $S$ or orthogonal $Q$ are orthogonal

## Lecture 4: Eigenvectors and Eigenvalues

### Eigenvectors: Introduction

$Ax = \lambda x$. $x$ is an eigenvector of $A$, $\lambda$ is an eigenvalue of $A$.
Eigenvectors keep their direction unchanged when linear
transformation $A$ is applied to them.

- there are up to $n$ independent eigenvectors for $A \in \mathbb{R}^{n \times n}$
- $A^2x = A(Ax) = A\lambda x = \lambda Ax = \lambda^2 x$
- $A^{-1} x = \frac{1}{\lambda}x$. This does not work when $A$ is not invertible.
- $e^{At}x = e^{\lambda t}x$
- Let $v = \sum c_i x_i$. Then $A^k v = \sum c_i \lambda_i^k x_i$.
  Rotation in the basis of eigenvectors is extremely convenient and cheap

### Matrix Similarity

$B$ is similar to $A$ if $B = M^{-1}AM$.

-  Similar matrices have the same eigenvalues: Let $x$ be an eigenvector of
   $B$. Then $Bx = M^{-1}AMx = \lambda x \Rightarrow MBx = MM^{-1}AMx = AMx = \lambda Mx$.
   $Mx$ is an eigenvector of $A$.
- $AB$ has same eigenvalues as $BA$: $AB = A^{-1}ABA = BA$, so they are similar
  and have the same eigenvalues.
- $A = X\Lambda X^{-1}$ where $X$ is the matrix of eigenvectors of $A$
  - $A$ is similar to $\Lambda$

### Finding Eigenvalues

Let $A = \begin{bmatrix} 0 & 1 \\ -1 & 0\end{bmatrix}$

Looking for: $Ax = \lambda x$... $(A - \lambda I)x = 0 \Rightarrow det(A-\lambda I) = 0)$

$det(A - \lambda I) = det\left(\begin{bmatrix}-\lambda & 1 \\ -1 &-\lambda \end{bmatrix}\right) = \lambda^2 + 1 = 0
\Rightarrow \lambda_1 = i, \lambda_2 = -i$

You can verify with the following rules:

- $trace(A) = \sum \lambda_i$
- $det(A) = \prod \lambda_i$

### Symmetric Matrices

- eigenvalues of $S$ are real
- eigenvectors of $S$ are orthogonal
- all $n$ eigenvectors exist (some may be repeated)
- $S\begin{bmatrix}x_1 & x_2\end{bmatrix} = \begin{bmatrix}x_1 & x_2\end{bmatrix} \Lambda$ 
  is a diagonalization of $S$ 
  - based on the eigenvectors $x_1, x_2$
  - $S = Q\Lambda Q^T$ is the symmetric version of $A = X\Lambda X^{-1}$
  - spectral theorem

# Lecture 5 - Positive Definite Matrices

A matrix $S$ is positive-definite iff and only if:

1. All $\lambda_i > 0$
2. $x^TSx > 0$ when $x > 0$: the _energy_ of $S$
3. $S = A^TA$
4. All leading determinants $> 0$
   - _leading determinants_ are the determinants of all upper-left square matrices
5. All pivots in elimination $> 0$

Surprisingly enough, these five things are equivalent.

# Lecture 6 - Singular Value Decomposition

Moving from square, symmetric matrices and _spectral decomposition_ to
rectangular matrices, where the concept of eigenvalues can not exist due to
dimensionality mismatch.

Without further ado, singular value decomposition for $A \in \mathbb{R}^{m \times n}$: $A = U\Sigma V^T$

- $U \in \mathbb{R}^{m \times m}$ is a _unitary matrix_ containing the left-singular vectors
- $V \in \mathbb{R}^{n \times n}$ is a _unitary matrix_ containing the right-singular vectors
- $\Sigma \in \mathbb{R}^{m \times n}$ is a _diagonal matrix_ containing the singular values

Basic idea: take $A^TA = V\Lambda V^T \in \mathbb{R}^{m \times m}$. This is a
symmetric, positive-definite matrix. Accordingly, $V$ is orthogonal and the
eigenvalues on the diagonal of $\Lambda$ are all positive. Similarly, consider
$AA^T = U\Lambda U^T \in \mathbb{R}^{n \times n}$ with the _same eigenvalues_
for all non-zero elements (the larger of the two matrices will have zeros
on the additional diagonal elements).

Another view on $A = U\Sigma V^T$ is that we're looking for $Av_i = \sigma_i
u_i$. In other words, we are trying to find _orthogonal_ vectors $v_i$
that are transformed into another set of orthogonal vectors $u_i$: $AV = U\Sigma$.

To find these vectors:

- we know that the _eigendecomposition_ of $A^TA = V\Lambda V^T$ produces
  orthonormal $V$ and a positive diagonal $\Lambda$
- $A^TA = (U\Sigma V^T)^TU\Sigma V^T = V\Sigma^TU^TU\Sigma V^T = V \Sigma^T \Sigma V^T 
  \Rightarrow V$ is orthogonal
- let $u_i = \frac{Av_i \sigma_i}{\sigma_i}$:

  We show that $u_j^Tu_i$ is orthonormal:

  $u_j^Tu_i = \frac{(Av_j)}{\sigma_j}^T\frac{(Av_i)}{\sigma_i} = \frac{v_jA^TAv_i}{\sigma_i\sigma_j}
    = \frac{v_j\sigma_i^2vi_i}{\sigma_i\sigma_j} = \frac{\sigma_iv_jvi_i}{\sigma_j}
    = v_j^Tv_j$. $v_j^Tv_i = 0$ for $j \neq i$. $v_j^Tv_i = 1$ for $j = i$, so $U$ is orthogonal.

  Left to show that $u_i$ are eigenvectors of $AA^T$:

  $AA^T u_i
  = AA^T\frac{Av_i}{\sigma_i}
  = A\frac{A^TAv_i}{\sigma_i}
  = \frac{\sigma_i^2Av_i}{\sigma_i}
  = \sigma_i^2\frac{Av_i}{\sigma_i}
  = \sigma_i^2u_i$

- $\Sigma = \Lambda^T \Lambda$ (less the shape differences)

## Geometric interpretation

The three parts of a SVD can be interpreted as a three consecutive linar transformations:

- $V^Tx$: rotation (isometric)
- $\Sigma (V^Tx)$: stretching $V^Tx$ to an ellipse
- $U (\Sigma V^Tx)$: another rotation


# Lecture 7 - Eckart-Young Theorem / PCA

## Eckart-Young: Approximating a Matrix with rank $k$

Let $A = U\Sigma V^T= \Sigma_{i=1}^r\sigma_i u_i v_i^T$, $r = rk(A)$.

Then $A_K = U\Sigma V^T= \Sigma_{i=1}^K\sigma_i u_i v_i^T$ is the rank-$k$
matrix that minimizes $||A - B||$ for any rank-$k$ matrix $B$.

In other words, the first $k$ singular values and vectors provide the best possible
approximation for $A$.

## Norms

Properties of norms:

- $||cv|| = |c|||v||$
- $||v + w|| \leq ||v|| + ||w||$
- $||v|| \geq 0$
- $||Qv|| = ||v||$ (because a rotation is isometric in euclidian space)

### Vector Norms

\begin{align}
	||v||_1 &= |v_1| + \dots + |v_n| \\
	||v||_2 &= \sqrt{|v_1|^2 + \dots + |v_n|^2} \\
	||v||_p &= \sqrt[^p]{|v_1|^p + \dots + |v_n|^p} \\
	||v||_{\infty} &= max|v_i|
\end{align}

Note that $\lim_{p \to \infty} ||v||_p = ||v||_{\infty}$ because the largest
vector component dominates the sum.

### Matrix Norms

\begin{align*}
	||A||_2 &= \sigma_1 \\
	||A||_F &= \sqrt{|a_{11}|^2 + \dots + |a_{mn}|^2} = \sqrt{|\sigma_1|^2 + \dots + |\sigma_r|^2}\\
	||A||_{\text{nuclear}} &= \Sigma_i\sigma_i
\end{align*}

While not really mentioned explicitly in the lecture, the singular value
decomposition leads to a series of $A_k$, where $A_k$ contains the the first
$k$ principal components $\sigma_i u_i v_i^T$.

# Lecture 8 - Vector and Matrix Norms

## Vector Norms

$||v||_0 = \text{number of non-zero components}$. But that's not a norm,
because it violates some of the norm properties above.

Visualizing norms: plot all $||v|| = 1, v \in \mathbb{R}^2$.

- $l_1$: diamond shape centered at the origin
- $l_2$: circle centered at the origin
- $l_\infty$: square centered at the origin
- $l_0$: the entire x and y axis themselves, without the origin.

Norms should be convex, and any norm $l_p$ with $p < 1$ need not apply.

$||v||_S = \sqrt{x^T S x}$ for a symmetric, positive definite matrix $S$.

We discuss how the $l_1$ norm produces a sparse result, whereas $l_2$
doesn't by solving a constraint system $min||x||_p$ where $Ax = b$.

## Matrix Norms

Deriving $||A||_2$ in terms of the vector norm for $x$: 
$$||A||_2 = \max_{\forall x}\frac{||Ax||_2}{||x||_2}$$

The vector that maximizes this is $v_1$, the first singular vector (why?).
$\frac{||Av_1||_2}{||v||_2} = ||Av_1||_2 = ||\sigma_1 u_1||_2 = \sigma_1$.

Multiplying with an orthogonal matrix does not change the norm (for the norms
discussed here, at least).

# Lecture 9 - 11: Different ways to solving Least Squares

## Pseudoinverse $A^+$

- $A^+ = A^{-1}$ when A is invertible
- otherwise, it's the best possible approximation to what an inverse would do
- $A^+Ax = x$ when $Ax \neq 0$: normal inverse behavior for $C(A)$
- $A^+Ax = 0$ when $Ax = 0$: can't recover the information for $N(A)$
- $A^+v_{r+1} = 0$ when $v_{r+1} \in N(A^T)$
  - $r$ is the rank, so $v_{r+1}$ is the first singular vector that's not reachable via $A$.

\begin{align}
	A^{-1} &=  V\Sigma^{-1}U^T & \text{if $A$ is invertible}\\
	A^+ &=  V\Sigma^+U^T & \text{for all $A$}
\end{align}

- $\Sigma^{-1}$ is a square diagonal matrix, without zero elements on the diagonal
  - the diagonal elements are $\frac{1}{\sigma_{ii}}$
- $\Sigma^+$ has rank $r$, and $\frac{1}{\sigma_{ii}}$ for the _non-zero_ diagonal elements

## Solve $A^TA\hat{x} = A^Tb$

There's a solution to $Ax = b$ if and only if all $b$s lie on the same n-dimensional
hyperplane ($A \in \mathbb{R}^{m \times n}$. If that's not the case, we can only
approximately solve it, using least squares (linear regression): 

$$minimize||Ax - b||^2_2 = (Ax-b)^T(Ax-b) = x^TA^TAx - 2b^TAx + b^Tb$$.

The derivative of that leads to: $$A^TA\hat{x} = A^Tb$$

Geometrically, $A\hat{x}$ will be the projection of $b$ onto the column space
$C(A)$, the nearest point to $b$ that can be represented using some $Ax$.

This works well assuming $A$ has independent columns, as only then $A^TA$
is invertible.

Claim: When $N(A) = 0, rk(A) = n$, then $A^+b = (A^TA)^{-1}A^Tb$.

Using SVD, we get: \begin{align}(A^TA)^{-1}A^T &= ((V\Sigma^T U^T)(U\Sigma V^T))^{-1}(V\Sigma^TU^T) \\
	&= ((V\Sigma^T\Sigma V^T))^{-1}(V\Sigma^TU^T) \\
	&= (V(\Sigma^T\Sigma)^{-1} V^T)(V\Sigma^TU^T) \\
    &= (V(\Sigma^T\Sigma)^{-1}\Sigma^TU^T)
\end{align}

## Orthogonalize first

This is just a minor modification of the previous approach: when the input
columns are nearly singular (almost the same direction), pre-process using
Gram-Schmidt to have a numerically saner process.

One improvement to the regular Gram-Schmidt is to reorder the columns in a way
that picks the numerically largest orthogonal vector first.


## Regularization: $(A^TA + \delta^2I)x = A^Tb$

Consider a near-singular $A = U \Sigma V^T$. Then $A^{-1} = V \Sigma^{-1} U^T$
becomes very large because $\Sigma^{-1}$ contains the reciprocal values of the
very small singular values in $A$.

This can lead to numerically unstable computations, where a small deviation in
$A$ changes $x$ considerably.

In such a case, adding a penalty term, e.g. $\delta^2 ||x||^2$ reduces the
magnitude of $x$, leading to a well-conditioned problem. The problem that is
now being solved is

$$
\begin{bmatrix}
	A \\ \delta I
\end{bmatrix}x =
\begin{bmatrix}
	b \\ 0
\end{bmatrix}
$$

This augmented system $A^*x = b^*$ will solve $min(||Ax - b||^2 + \delta^2||x^2_2)$.
The normal equation leads to $A^TA + \delta^2I)x_\delta = A^Tb$.

What happens for $\delta \to 0$? For any $A$, $A^TA + \delta^2I)^{-1}A^T \to A^+$
Checked in the lecture for $1 \times 1$.

The regularization above is known as _ridge regression_. Using the $l_1$-Norm,
the solution gives sparse solution (variable selection) and is called
_the lasso_.



## Gram-Schmidt

Let $A = [a_1 a_2 \dots a_n]$ where the columns $a_i$ are _barely_ independent,
i.e. they are almost dependent. Geometrically, the lines drawn from the columns
would be almost parallel. The goal is to represent $A = QR$, i.e. by
a _orthogonal_ $Q$ and a linear combinations of the columns of $Q$, $R = Q^TA$.

Basic idea: iteratively work through the columns, for each column removing the
parts of the vector that go into the direction of the previous columns (by
computing projections). The result is a vector that is perpendicular to all
previously computed columns. Finally, we make it a unit vector by dividing by
its norm.

\begin{align}
	q_1 &= \frac{a_1}{||a_1||} \\
	A_2 &= a_2 - (q_1^Ta_2)q_1 \\
	q_2 &= \frac{A_2}{||A_2||} \\
	\vdots \\
	A_i &= a_i - \sum_{j=1}^{i-1}(q_j^Ta_i)q_j \\
	q_i &= \frac{A_i}{||A_i||}
\end{align}

### Improving Numerical Robustness

When $a_1$ and $a_2$ are almost in the same direction, the above calculation
produces a very small residual vector $A_2$, which is then scaled up to unit
length. 

Instead of using $A_2$, we could find the largest possible $A_i$ from all
remaining columns, swap columns and proceed. To not take a performance hit by
having to recompute all the intermediate values, we can proceed in a different
order:

- compute $q_i$, starting at $i = 1$.
- next compute $A_j' = A_j - (q_i^TA_j)q_i$ for $j = i+1 \dots n$
  - subtract $q_i$ from all remaining columns
- find the largest $A_j'$ and produce $q_{i+}$ via column exchange

This changes the order from computing $A_i$ in one go into incrementally
adapting $A_j$ with each new $q_i$ that is finalized.

## Krylov Spaces

Solving $Ax = b$, for a large, sparse $A$. Idea: take a vector $b$ and compute
$Ab, A(Ab), A(A(Ab)), \dots, A^{j-1}b$. These vectors span a space called a
_Krylov Space_ $K_j$. We compute $x_j$ as the best vector in $K_j$ that
is hopefully close to the real solution $x$.

The base of $K_j$ may not be very good, so it seems to be a good a idea to
apply Gram-Schmidt to orthonormalize the basis.

# Lecture 12: computing eigenvalues and singular values

## Computing eigenvalues using Gram-Schmidt

\begin{align*}
	A_0 &= Q_0R_0 &\to R_0Q_0 &= A_1 \\
	A_1 &= Q_1R_1 &\to R_1Q_1 &= A_2
\end{align*}

$A_1$ is similar to $A_0$: $UV \sim U^{-1}UVU = VU$.

When proceeding to compute $A_i$, the matrices start having small values,
off-diagonal, especially below the diagonal. At the same time, the values on
the diagonal converge to the eigenvalues.

## Improvement: introduce a shift

$A_0 - sI = Q_0R_0 = R_0Q_0 + sI = A_1$. Eigenvalues shift by $s$.

$A_1 = R_0 Q_0 + sI = R_0 (A_0 - sI)R_0^{-1} + sI = R_0A_0R_0^{-1}$ so $A_1$
is still similar to $A_0$.

To reduce comutations further, it can be shown that

- a Hessenberg-matrix can be computed cheaply
  - Hessenberg-matrix: one diagonal on top of upper triangle form
- the zeros stay there during QR-factorization
- QR-factorizaton can be done cheaper with that form

## Computing singular values

When $A = U\Sigma V^T$, $\Sigma$ contains the singular values on the diagonal.
With $QA = QU\Sigma V^T$ the singular values remain the same, because a
rotation does not change the singular values. As a consequence, it's safe to
multiply a matrix $A$ with a $Q$ (from left or right) without changing the
singular values.

With this in place we can find $Q$s that transform $A$ into a bidiagonal matrix,
which leads to $A^TA$ being tridiagonal (and symmetric), for which we use the
above method to compute the eigenvalues, which are the squared singular values.

# Lecture 13 - Randomized Linear Algebra

When $A$ is really large, we can't work with $A$ directly. Instead, e.g. when
solving $Ax = b$, we sample columns from $A$ and elements of $b$.

Sampling columns from $A$ should be weighted by the squared size of the column.
Consider a matrix multiplication $AB$, which will be sampled:

$$p_j = \frac{||a_j||||b_j^T||}{C}, \text{ where } C = \sum||a_j||||b_j^T|| $$

$$ AB_\text{approx} = pj\frac{a_jb_j^T}{sp_j}$$

I don't see how the $p_j$ in the denominator comes up. Rest of the lecture is
showing that mean comes out correct (mean is $AB$), and computing the variance.
We also show how the choice of $C$ above actually minimizes the variance using
Lagrange multipliers.

# Lecture 14 - Low-Rank Matrices

Note: this will lead up to the _Matrix Inversion Formula_ in signal processing.

General applications involve _recursive_ least squares:

1. solving $(I - UV^T)x = b$, i.e. a small perturbation to $I$
1. add a new measurement (row) to $Ax = b \to A^TA\hat{x} = A^Tb$ 
   - $[A^Tv][\begin{bmatrix}A\\v^T\end{bmatrix}\hat{x}_\text{new} = [A^T v] \begin{bmatrix}b\\b_\text{new}\end{bmatrix}$
   - don't want to compute from scratch, but integrate existing solution and adapt for new measurement
   - $[A^Tv][\begin{bmatrix}A\\v^T\end{bmatrix} = A^T + vv^T$, a rank-1 perturbation of $A^TA$
   - $k$ new measurements fit accordingly
   - since we need $(A^TA)^{-1}$ to compute $\hat{x}$, incremental computations
     help a lot
1. compes up in Kalman filters

## Rank-1 Perturbation of $I$

Question: perturbating a matrix by a rank one: $(I - uv^T)$, what's the inverse
$(I - uv^T)^{-1}$ doing? 

Answer:$(I - uv^T)^{-1} = I + \frac{uv^T}{1-v^Tu}$. Check by multiplying back to get $I$.

Note that the inverse is just a rank-1 change to $I$, because $uv^T$ has rank
$1$ and $1 - v^Tu$ is scalar.

## Rank-k Perturbation of $I$

Let $U, V \in \mathbb{R}^{n \times k}$. Then $I - UV^T$ is a rank-$k$
perturbation of $I$. Then $I + U(I_k - V^TU)^{-1}V^T$ is the rank-$k$ inverse.

## Rank-k Perturbation of $A$

$(A - UV^T)^{-1} = A^{-1} + A^{-1}U(I_k - V^TA^{-1}U)^{-1}V^TA^{-1}$

# Lecture 15

## Derivative: $\frac{dA^{-1}}{dt}$

Assume $A(t)$, i.e. $A$ depends on $t$. We aim to find $A^{-1}(t)$.

\begin{align*}
	B^{-1} - A^{-1} &= B^{-1}(A - B)A^{-1} \\
	\Delta A^{-1} &= (A + \Delta A)^{-1}(-\Delta A)A^{-1} \\
	\frac{\Delta A^{-1}}{\Delta t} &= (A + \Delta A)^{-1}\frac{-\Delta A}{\Delta t}A^{-1}
\end{align*}

TODO: figure out how we got from the first equation to the second!

Let $\Delta t \to 0$: $\frac{dA^{-1}}{dt} = -A^{-1}\frac{dA}{dt}A^{-1}$
Substituting $1/t$ instead of $A$ for the $1$ by $1$ case yields the formula
$\frac{d 1/t}{dt} = -\frac{1}{t^2}$

## Derivative: $\frac{d\lambda}{dt}$

For eigenvector $x$ of $A$: $A(t)x(t) = \lambda(t)x(t)$. Also, by symmetry and
because the eigenvalues of $A$ are equal to those of $A^T$: $y^T(t)A(t) = \lambda(t)y^T(t)$
($x(t)$ is an eigenvector of $A$, and $y(t)$ one of $A^T$). Also, let $y^Tx = 1$,
as a natural way to normalize the lengths of the eigenvectors.

In matrix notation: $AX = X\Lambda, Y^TA = \Lambda Y^T, Y^TX = I$
(for all eigenvectors at once).

$y(t)A(t)x(t) = \lambda (t) y^T(t) x(t) = \lambda (t)$

\begin{align*}
	\frac{d\lambda}{dt} &= \frac{dy^T}{dt}A(t)x(t) &+ y^T(t)\frac{dA}{dt}x(t) &+ y^T(t)A(t)\frac{dx}{dt} \\
						&= \lambda(t)\frac{dy^T}{dt}x(t)   &+ y^T(t)\frac{dA}{dt}x(t) &+ \lambda(t)y^T(t)\frac{dx}{dt} \\
						&= y^T(t)\frac{dA}{dt}x(t) &+ \lambda(t)(\frac{dy^T}{dt}x(t)   &+ y^T(t)\frac{dx}{dt}) \\
						&= y^T(t)\frac{dA}{dt}x(t) &+ \lambda(t)\frac{d1}{dt}  & \\
						&= y^T(t)\frac{dA}{dt}x(t) 
\end{align*}

because $\frac{dy^T}{dt}x(t) + y^T(t)\frac{dx}{dt}$ is the derivative of $y^Tx$, which is $1$.

Rest of lecture: how do eigenvalues of $S$ change when adding a rank-1 $uu^T$? 
They get larger: $\lambda_1 \geq \gamma_1 \geq \lambda_2 \geq \gamma_2 \dots$,
where $\gamma_i$ are eigenvalues of $S$ and $\lambda_i$ are eigenvalues of $S + uu^T$.
This is called _interlacing_. Not going into details here.

# Lecture 16 - 17

## Derivative $\frac{d}{dt}A^2$

$\lim_{t \to 0}\frac{(A + \Delta A)^2 - A}{\Delta t}
= \lim_{t \to 0}\frac{\Delta A A + A \Delta A + (\Delta A)^2}{\Delta t} = A \frac{dA}{dt} + \frac{dA}{dt}A$.

Because matrix multiplication doesn't commute, $\frac{d}{dt}A^2 \neq 2A\frac{dA}{dt}$,
but the two orders of the multiplicatoin appear separately.

## Derivative $\frac{d\sigma}{dt}$

$\sigma(t) = u^TAv$

$\frac{\sigma}{dt} 
  = \frac{du^T}{dt}Av + u^T\frac{dA}{dt}v + u^TA\frac{dv}{dt}
  = u^T\frac{dA}{dt}v + \sigma(\frac{du^T}{dt}v + v^T\frac{dv}{dt})$

$\frac{d\sigma}{dt} = u^T(t)\frac{dA}{dt}v(t)$. 

# Weyl's Inequality

Goal: estimate $\lambda(S+T)$ in terms of $\lambda(S)$ and $\lambda(T)$ for symmetric $S$, $T$.

Weyl: $\lambda_{i+j-1}(S) \leq \lambda_i(S) + \lambda_j(T)$

No proof given (it's in the class notes we don't have access to).

# Lectuer 17 - Townsend

Given $X = \mathbb{R}^{n \times n}$, singular values 
$\sigma_1 \geq \sigma_2 \geq  \sigma_3 ... \geq \sigma_n$. These singular
values give us a lot of information about the matrix:

- number of non-zeros: $rk(X)$
- $dim(C(A)) = dim(C(A^T)) = k$
- $X = u_1v_1^T + ... + u_kv_k^T$
- how closely can we approximate $X$ with a rank-$k$ matrix (ratio between
  $\sigma_1, ...$ and $\sigma_{k+1}$


Q: what properties must $X$ so it has low rank? Think about sending an $n \times n$
picture: if sending the row / col vector components $(u_i, v_i)$ is less
data than sending the actual $X$, we say $X$ is _low rank_.

D: $X$ is said to be _low rank_ if $2kn < n^2 (k < n/2)$. Often, we require $k << n/2$.

Q: What do low rank matrices look like?

- rank-$1$: matrix highly aligned
- example: flags: austrian flag, ... (stripes) are rank 1
- example: triangular flag: (diagonal matrix $X$ w/ $0$ above the diagonal, rest $1$)
  - the inverse $X^{-1}$ is $1$ one the diag, and $-1$ on the diag below (rest all zeros,
    both below the diag and above the diag)
- $\sigma_k(X) = \left[2 sin \frac{\pi (2k-1)}{2(2n+1)}\right]^{-1}$
- intuitively, we need structure aligned with rows and cols to get low-rank

## the japan flag (circle) is surprisingly low cank, too

The argument goes roughly as follows:

- we split the flag into a sum of two pictures: 
  - one contains the parts of the circle
    that are the "around" the maximum inner-circle square (side length: $2\sqrt{r}$)
  - the other one contains the rest, i.e. the inner-circle square
- the latter has rank $1$
- the former can be further split up (due to symmetry) by cutting in half in
  both dimenions
- then, the remaining quarter-circle-thing-without-inner-square can be split up
  into a bit that's more wide than tall (small column space) and the other one
  which is more tall than wide (small col space).

## Numerical rank of a matrix $A$

For $0 < \epsilon < 1$: $rank_\epsilon(X) = k$ when $\sigma_{k+1}(X) \leq
\epsilon \sigma_1(X)$. In other words, we consider small singular values below _tolerance_
$\epsilon$ as $0$, relative to the first singular value $\sigma_1(X)$.

$rank_0(X) = rank(X)$ is the same as the original definition of rank.

Idea: when compressing an image, we ignore extremely small singular values,
because the received picture on the receiver side is not visually different.

Eckart-Young: $\sigma_{k+1}(X) = ||X - X_k||_2$ tells us how well we can
approximate $X$ with a rank-$k$ matrix.

## Matrices of numerically low rank

- every low rank matrix is naturally of numerically low rank
- $H_{jk} = \frac{1}{j+k-1}$: Hilbert matrix, full rank but low numerical rank
- [Vandermonde matrices](https://en.wikipedia.org/wiki/Vandermonde_matrix)

Numerically low rank is bad for computing inverses (almost singular).

Rest of lecture skipped, because it just shows some reasons on why many matrices
are low (numerical) rank, way above my had and applicability.

# Summary: Matrix Decomposition

Let $A \in \mathbb{R}^{n \times n}$. Then $A$ has $n^2$ free variables.
In contrast, $S$ has $\frac{1}{2}n(n+1)$ parameters, because below the
main diagonal there's no freedom to chose values.

- $A = LU$
  - $L$: triangular, $\frac{1}{2}n(n-1)$ free variables
  - $U$: triangular, $\frac{1}{2}n(n+1)$ free variables
  - overall, $n^2$ free variables
- $A = X\Lambda X^{-1}$
  - every eigenvector in $X$ has $n - 1$ free variables (can be scaled,
  so we just fix the first component to $1$ and have $n-1$ left)
  - so $n^2 - n$ in $X$
  - $n$ in the eigenvalues $\Lambda$
- $A = QR$
  - $n-1$ for the first column (normalized, so one fixed param)
  - $n-2$ for the second: normalized and orthogonal to the one, two conditions
  - ... up to $0$ free variables for the last column
  - leading to $\frac{1}{2}n(n-1)$ for $Q$, and $R$ is triangular again (see $LU$)
- $S = Q\Lambda Q^T$
  - $Q$, just as $X$, has $\frac{1}{2}n(n-1)$ free parameters
  - $\Lambda$ has another $n$, so we're good with $\frac{1}{2}n(n+1)$
- $A = QS$ (polar decomposition, orthogonal times symmetric)
  - $Q$, as previously, has $\frac{1}{2}n(n-1)$ free parameters
  - $S$, as in the beginning, has $\frac{1}{2}n(n+1)$ free parameters
  - this again adds up to $n$, phew
- $A = U\Sigma V^T$, $A \in \mathbb{R}^{m \times n}$, so $A$ has $mn$ params
  - let $m \leq n$, and let $A$ be full rank $m$.
  - $U \in \mathbb{R}^{m \times m}$ -- $\frac{1}{2}m(m-1)$
  - $\Sigma \in \mathbb{R}^{m \times n}$ -- $m$ params
  - $V \in \mathbb{R}^{n \times n}$ -- $(n-1) + (n-2) + ... + (n-m) = mn - \frac{1}{2}m(m+1)$$
- now, let the rank of $A$ be $r$
  - $U \in \mathbb{R}^{m \times r}$ -- $(m-1) + (m-2) + ... + (m-r) = mr - \frac{1}{2}r(r+1)$
  - $\Sigma \in \mathbb{R}^{r \times r}$ -- $r$ params
  - $V \in \mathbb{R}^{r \times n}$ -- $(n-1) + (n-2) + ... + (n-r) = rn - \frac{1}{2}r(r+1)$$
  - so that's $mr + nr - r^2$, which means that's the number of parameters for
    a rank-$r$ matrix  $A$ 

# Saddle points

## Saddle points from constraints

Goal: $min\frac{1}{2}x^TSx$ with constraints $Ax = b$.

Lagrange multiplier: $L(x, \lambda) = \frac{1}{2}x^TSx + \lambda^T(Ax - b)$. 
Taking the derivative with respect to $x$ and $\lambda$:

\begin{align*}
	\frac{\partial}{\partial x}L(x, \lambda)       &= Sx + A^T\lambda \\
	\frac{\partial}{\partial \lambda}L(x, \lambda) &= Ax - b 
\end{align*}

Which leads to the matrix representation:

$\begin{bmatrix} S & A^T \\ A & 0 \end{bmatrix}
\begin{bmatrix} x \\ lambda \end{bmatrix}
= \begin{bmatrix} 0 \\ b\end{bmatrix}$

The solution to this system of equation is a _saddle point_ of the lagrange
function, and not a minimum, because the matrix is indefinite. One can see
that the first $n$ pivots will be positive (courtesy of $S$), but the following
will all be negative:

$(AS^{-1})\begin{bmatrix} S & A^T \\ A & 0 \end{bmatrix} 
  \to \begin{bmatrix} S & A^T \\ 0 & -AS^{-1}A^T \end{bmatrix}$

This result comes from doing _block elimination_ on the parts of the matrix, 
i.e. blockwise elimination of the lower rows (coming from $A$).

## Saddles from Rayleigh-Quotient $R(x) = \frac{x^TSx}{x^Tx}$

Similarly to the derivation of $||A||_2 = \frac{||Av_1||_2}{||v_1||_2}$, the
maximum value of the quotient is $\lambda_1$, with $x = q_1$ (first eigenvector).
The minimum value is $\lambda_n$ (with $x = q_n$, the last eigenvector).

For the _Rayleigh-Quotient_, the derivative $\frac{d}{dx}R(x)$ is $0$ at the
eigenvectors, and those are saddle points for $\lambda_2 ... \lambda_{n-1}$.
$\lambda_2 = max_V min_V \frac{x^TSx}{x^Tx}$; saddle points can be expressed as
a max-min problem (no details here, a bit too much) over $(n-k)$-dimensional
spaces $V$ (huh). The min-max is kind of what the saddle is all about.

## Aside: Vandermonde matrix

$$V = \begin{bmatrix}
1 & \alpha_1  & \alpha_1^2 & ...    & \alpha_1^{n-1} \\
1 & \alpha_2  & \alpha_2^2 & ...    & \alpha_2^{n-1} \\
\vdots        & \vdots     & \ddots & \vdots \\
1 & \alpha_m  & \alpha_m^2 & ...    & \alpha_m^{n-1} \\
\end{bmatrix}$$

Expresses some data points $x_1, x_2, ... x_m$ that is fit by a degree-$n$
polynomial when solving $Vc = b$: the first column represents the degree-0, 
followed by successively increasing polynomials, built from $x^k$ for $k = 1 .. n-1$.

# Statistics

- discrete probabilities: $p_1 + ... + p_n = 1$. 
- continuous probabilities: $\int_{-\infty}^{\infty}p(x)dx = 1$

## Mean

- sample mean: $\mu = \frac{1}{N}\sum_{i=1}^Nx_i$ based on $N$ _observations_
- expected mean: $E[x] = m = \sum_{i=1}^np_ix_i$ for an entire _population_ with $n$
  possible different outputs

## Variance

_Variance_ describes how far the data is away from the mean on average. Larger
deviations are weighted higher.

Sample variance, again, is based on $N$ observations: $\frac{1}{N-1}\sum_{i=1}^N(x_i - \mu)^2$.

Expected variance: $\sigma^2 = E[(x-m)^2] = \sum_{i=1}^np_i(x_i-m)^2$:

$E[(x - m)^2] = \sum_{i=1}^np_i(x_i-m)^2 
  = \sum_{i=1}^np_ix_i^2 - 2\sum_{i=1}^n(p_ix_im) + \sum_{i=1}^np_im 
  = \sum_{i=1}^np_ix_i^2 - 2m^2 + m^2 = E[x^2] - m^2$

## Markov's Inequality

Applies when $\forall x_i \geq 0$: Prob$[x \geq a] \leq \frac{E[x]}{a}$

- walks through some example, but omits the proof anyway

## Chebyshev inequality

Prob$\left[|x - m|\right] \geq a \leq \frac{\sigma^2}{a^2}$

Proof omitted: based on Markov, with $a = |x - m|^2$, and everything falls into its place.

## Covariance

- consider $m$ experiments done at once
- out falls $\Sigma \in \mathbb{R}^{m \times m}$
- consider $m = 2$ with two random variables $x, y$: 
  - $V = \sum_{x_i, y_j}p_{ij}\begin{bmatrix}x_i - m_x\\ y_j - m_y\end{bmatrix}[x_i - m_x y_j - m_y]$
  - $V = \begin{bmatrix}\sigma_x^2 & \sigma_{xy} \\ \sigma_{xy} & \sigma_y^2\end{bmatrix}$
  - $V$ is a linear combination of rank-$1$ positive semi-definite matrices and as such is itself psd
  - $p_{ij}$ is the probability that $X = x_i$ and $Y = y_j$
  - $\sum_i p_{ij} = p_j$
  - $\sum_j p_{ij} = p_i$
  - those are called _marginals_ of the joint probability of $x$ and $y$

# Optimization

## Taylor Series

\begin{align*}
	F(x + \Delta x) &\approx F(x) + \Delta x\frac{dF}{dx} + \frac{1}{2}(\Delta x)^2\frac{d^2F}{dx^2} + ... \\
					&\approx \sum_{n=0}^\infty \frac{F^{(n)}(\Delta x)}{n!}(x + \Delta x)^n
\end{align*}

(second line is still wrong, but I can't figure out how)

In case of $x \in \mathbb{R}^n$, this becomes:
$F(x + \Delta x) \approx F(x) + \Delta x^T\nabla F(x) + \Delta x^T H \Delta x$ where
$\nabla F = \begin{bmatrix}\partial F/\partial x_1 \\ \vdots \\ \partial F \partial x_n \end{bmatrix}$ 
is the gradients of $F$ with resepct to the $x_i$ and
$H = \begin{bmatrix}
	\partial F/\partial x_1x_1 & ... & \partial F/\partial x_1x_n  \\ 
	\vdots & & \vdots \\
	\partial F / \partial x_nx_1  & ... & \partial F/\partial x_nx_n  
\end{bmatrix}$ is the symmetric _Hessian_ matrix containing the second derivatives.

There's some symmetry to the _Jacobian_ $J$ with $J_{jk} = \frac{\partial f_j}{\partial x_k}$.
The Jacobian of the gradient is the Hessian, which makes sense because the Jacobian
is the first derivative of the $n$ functions that are the gradients, so it leads to
second derivatives just as the Hessian.

## Minimize $F(x)$ for $x \in \mathbb{R}^n$

We need to solve $f = 0$ (I assume $\nabla F(x) = f)$, which means 
$(f_1 = 0, ... , f_n = 0)$.

### Newton's method for $f = 0$

At the $k$th iteration: $0 = f(x_k) + J(x_k)(x_{k+1}-x_k)$. This is solved by
$x_{k+1} = x_k - J(x_k)^{-1}f(x_k)$ (by multiplying with $J(x_k)^{-1}$
and bringing $x_{k+1}$ to the left).

Let $f(x) = x^2 - 9$. For $x_k = 3$, $x_{k+1} = 3 - \frac{1}{2x_k}f(x_k^2 - 9) = 3$,
so the method is consistent with a solution. However, starting at $x_0 = 0$,
we see that $x_1 = 0 - \frac{1}{0}(-9)$. The reason is that the derivative at
$0$ is $0$, so we are not able to intersect the x-axis.

### Minimizing $F(x)$

1. Steepest descent: $x_{k+1} = x_k - s_k \nabla F$ where $s_k$ is called _learning rate_.
   - $-s_k \nabla F$ is a _direction_ in $\mathbb{R}^n$
   - exact line search: stop when $x_k -s_k \nabla F$ is minimal (bottom of bowl)
   - linear rate of convergence
2. Newton's method: $x_{k+1} = x_k - H^{-1}\nabla F$
   - quadratic rate of convergence (errors are squared in every step)
   - this is great if you start close enough
3. Levenberg-Marquardt
   - attempt to improve on steepest gradient method by using some information
     from the Hessian
   - more robust than _Gauss-Newton_, finding optima in cases where GNA diverges
   - "cheap man's Newton" (Strang)

## Convexity

### Convex set $K$

Given an affine subspace $S$, a _subset_ $K \subseteq S$ is convex if for all
$x, y \in K$ the line connecting $x$ and $y$ is also in $K$.

- the union of two convex sets is usually _not_ convex
- the intersection of two convex sets is always convex

### Convex function $F(x)$

A function $F(x)$ is convex if:

- the _epigraph_ of the function is a convex set.
- $d^2F/d^2x \geq 0$ ($x \in \mathbb{R}$)
- $F(x) = x^TSx$ if $S$ is positive semi-definite

Given two convex functions $F_1(x), F_2(x)$:

- $m(x) = min(F_1, F_2)$ is not convex
- $M(x) = max(F_1, F_2)$ is convex

Connection between the Hessian of a function $f$ and convexity: if the Hessian
is positive semi-definite, $f$ is convex. For strict convexity, positive-definite
is required. 

Example: Given $f(x) = \frac{1}{2}x^TSx - a^Tx - b$, $\nabla f = Sx - a$, and 
the Hessian $H = S$. The minimum of $f(x)$ is where the gradient is zero, at 
$x^* = S^{-1}a = \text{argmin}(f)$.


# Gradient Descent

When the number of variables is too large or the function is too complicated to
take second derivatives, Gauss-Newton cannot be applied, as it depends on the
Hessian. In this case, the natural method relying on only the _gradient_ of
a function $f$ is gradient descent:

$$ x_{k+1} = x_k - S_k\nabla f(x_k)$$

The only thing left to do is to determine the step size $S_k$. $\nabla f(x_k)$
is called the _search direction_.

- _exact line search_: choose $S_k$ to minimize $f(x_{k+1})$ along the search line
- _backtracking_: start with a fixed $S_0$, and iterate: 
  $\frac{1}{2}S_0, \frac{1}{4}S_0, ...$

Example: exact line search on $f(x) = \frac{1}{2}x^TSx$ for 
$S = \begin{bmatrix}1 & 0 \\ 0 & b\end{bmatrix}$, equals 
$f(x) = \frac{1}{2}(x^2 + by^2)$.

- the level lines are ellipses, whose width depends on $b$.
  - level lines: all points $x$ for $f(x) = c)$
- for very small values of $b$, the ellipse is very narrow
- the search directions ($\nabla f$) is perpendicular to level line
- almost parallel to the x-axis, so convergence is very slow (zig-zagging slowly)
- $x_k = b\left(\frac{b-1}{b+1}\right)^k$: $x_k$ flips between positive and negative for adjacent $k$
- $y_k = b\left(\frac{1-b}{1+b}\right)^k$: $y_k$ converges slowly for small $b$ 
  because fraction is close to one
- if $b$ is very large, the roles are flipped: $x_k$ depends on the fraction close to one,
  $y_k$ flips sign.

## Adding momentum

$x_{k+1} = x_k - sz_k$, with $z_k = \nabla f_k + \beta z_{k-1}$. The last term
is called a _memory_ term, providing the momentum by still going into the
direction of the previous step, and transitively also all previous steps.

Skipping: 

- derivation of the results for the introductory example
- showing $x_k = b\left(\frac{\sqrt{b}-1}{\sqrt{b}+1}\right)^k$ 
- $y_k = b\left(\frac{1-\sqrt{b}}{1+\sqrt{b}}\right)^k$
- which is apparently a significant improvement steepest gradient descent 
  with exact line search

Nesterov had a similar idea: $x_{k+1} = x_k + \beta(x_k - x_{k-1}) - s\nabla f(x_k + \gamma(x_k - x_{k-1}))$

$x_{k+1} = y_k - s\nabla f(y_k)$, which brings the same results as the momentum
method, with some different constant.


## SGD: stochastic gradient descent

General problem description: $min_x\frac{1}{n}\sum_{i=1}^nf_i(x)$

- typical form of a ML optimization problem: _finite sum_
- training data: ${(x_1, y_1), ... , (x_n, y_n)} \in mathbb{R}^d \times y$
- both $n$ and $b$ can be huge
  - $n$: number of training data points
  - $d$: number of dimensions (free variables)

Common examples: least Squares: $$\frac{1}{n}||Ax-b||_2^2 
  = \frac{1}{n}\sum_{i=1}^n(a_i^Tx - b_i)^2 
  = \frac{1}{n}\sum_{i=1}^nf_i(x)$$

Deep Neural Netowrks: $$\frac{1}{n}\sum_{i=1}^nloss(y_i, DNN(x; a_i)) = \frac{1}{n}\sum_{i=1}f_i(x)$$

Basic gradient descent iteration: $$x_{k+1} = x_k - \alpha_k \nabla f(x) = x_k - \alpha_k\frac{1}{n}\sum_{i=1}^n\nabla f_i(x_k)$$

Of course, this is the main concern: given the number of samples and the number
of dimensions, computing the gradient wrt all input variables $a_i$ and weights
$x_j$ can take hours just for a single step.

Main idea: use a random sample of inputs to reduce computational effort. Or,
bringing it to the extreme: what about taking a single observation instead?
Can this work? 

There's a nice [simulation](http://fa.bianp.net/teaching/2018/eecs227at/stochastic_gradient.html).
Observations:

- training loss decreases quickly initially, in the _far out zone_
- near the _zone of confusion_, this becomes very erratic
- not finding the optimum exactly is actually not a problem, as we want our
  system to generalize to unseen data

Core principle: replace an expensive subroutine (to compute $\nabla f$ with a
_noisy estimate_ $\nabla f_i$. SGD uses _stochastic gradients_ $g(x)$ such that
$\mathbb{E}[g(x)] = \nabla f(x)$. This is called an _unbiased estimator_.

While an unbiased estimator guarantees that your are stochastically doing the
right thing, the other important aspect of a stochastic method is the _variance_.

Two variants for implementing stochastic process:

1. randomly pick $i$ _with_ replacement
2. randomly pick $i$ _without_ replacement

Even though the first option seems to make more sense intuitively, the second
version is used in every toolkit. It allows to shove all data into the GPU 
after a pre-shuffle phase, so we don't have to jump back and forth between GPU
and CPU to pick the next observation randomly.

Another variant: mini-batch. Instead of picking one individual sample, pick $k$
of them. This reduces variance of the estimate for $\nabla f$, at the expense
of increased computational effort. GPU concurrency can recover the additional
computational effort for the minibatch, as the entire batch can be computed in
a SIMD-style fashion in one pass.

However, too large of a mini batch size reduces the variance too much, leading
to a very small region of confusion, overfitting to the training data and a
lack of generalization.

Selecting the "best" minibatch size is an unsolved problem, just as selecting
the best step size.

# Duality (?)

## Linear Programming

$minxc^Yx = c_1x_1 + c_2x_2 + ... + c_nx_n$ with constraints on x: $Ax = b$ 
and $x \geq 0$. In LP, the constraints are called the _feasible set_ of $x$'s.

- since the optimization target is a linear function, the feasible set is
  always a high-dimensional plane 
- this plane is cut off at the axis's, due to $x \geq 0$
- however, the number of corners grows quickly (exponentially) with $m, n$.

Simplex algorithm (Dantzig):

- start at one corner of the feasible set
- check all the neighbors, go to one that's improving my solution
- a bit like steepest descent along the edges

Karmarkar algorithm:

- travel interior points
- gradient descent but not just along the edges, but along interior points of
  the feasible set
- compute search direction
- stop when you run into the boundaries of the feasible set
- fun fact: LP is in $P$, not $LP$ (not with simplex, though)

### Dual LP

There's an associated twin problem: $max b^Ty$ for $A^Ty \leq c$. Solving this
or its dual solves the other problem as well.

Weak duality: $b^Ty \leq x^Tc$: $ with $Ax = b$, it follows that 
$x^TA^Ty \leq x^Tc$ for $x \geq 0$ ($x \leq 0$ would flip the $\leq$).

"Full" duality: $b^Ty = x^Tc$.

## Max flow = min cut

Given a graph $G = (V, E)$ with a source edge $s \in V$ and a sink edge $t \in V$.
Constraints are given by a maximum capacity for each edge $e$.

Duality: maximizing the cut meets minimizing the flow.

## Two-person games

Payoff matrix: $x \times y$: the amount $x$ is going to pay to $y$; player $x$
is minimizing, $y$ is maximizing return, e.g. 
$C = \begin{bmatrix} 1 & 2 \\ 4 & 8\end{bmatrix}$: $x$ selects the first row,
because those are the smaller values, and $y$ chooses second column, because
that's where the better values are for its goal. However, this _saddle point_
doesn't have to exist; take $C' = \begin{bmatrix} 1 & 8 \\ 4 & 2\end{bmatrix}$:
There's no canonical best choice, so the players will come up with a _mixed
strategy_, where they come up with some probabilities for choosing one of the
rows / columns, to maximize the outcome on average. 
This will lead to a _Nash equilibrium_, when increasing the number of players,
which is a lot more complicated than the two-player version.



# Structure of Neural Networks
o

