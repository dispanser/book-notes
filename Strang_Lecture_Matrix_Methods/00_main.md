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
