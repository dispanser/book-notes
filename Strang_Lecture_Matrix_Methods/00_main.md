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

