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
