<h1 align="center">The problem: solving Ax = b</h1> 

The main goal of this chapter is to understand how to solve the linear system $A\mathbf{x=b}$. 

We first revisit the traditional elimination method, then reformulate it as **Gaussian elimination** and **Gauss-Jordan elimination**, which provide a systematic matrix-based procedure for solving linear systems. After that, we focus on the special case of unique solvability and introduce the **inverse matrix**. Finally, we describe elimination algebraically through **elementary matrices** and **permutation matrices**, which leads naturally to the **LU** and **LDU** factorizations.

<h3 align="center">Gaussian-Jordan Elimination</h3> 

The traditional approach to solve a linear system is **elimination**: 

- We progressively reduce the number of unknowns until one variable can be solved directly, and then use substitution to recover the others.

In this process, we are repeatedly simplifying the system by eliminating variables and transforming it into an equivalent but much easier form. 

For a small system, this can be done equation by equation. As the system grows, such a procedure quickly becomes cumbersome. We want to find a more systematic approach to solve those problems.

```math
\left\{
\begin{aligned}
2x + 4y - 2z &= 2 & (1)\\
4x + 9y - 3z &= 8 & (2)\\
-2x - 3y + 7z &= 10 & (3)
\end{aligned}
\right.
```

This leads naturally to **Gaussian Elimination**. 

```math
\left[
\begin{matrix}
  2 & 4 & -2 \\
  4 & 9 & -3 \\ -2 & -3 & 7
\end{matrix}
\right]\left[
\begin{matrix}
  x \\ y \\ z
\end{matrix}
\right] = \left[
\begin{matrix}
  2 \\ 8 \\ 10
\end{matrix}
\right] \Leftrightarrow A\mathbf{x} = \mathbf{b} \Leftrightarrow
\left[\begin{array}{ccc|c}
  2 & 4 & -2 & 2 \\
  4 & 9 & -3 & 8 \\
  -2 & -3 & 7 & 10
\end{array}\right]
```



Instead of substituting equations one by one, we use row operations on the whole **augmented matrix** to create zeros below the diagonal entries, transforming the system into an upper triangular form.

In the $3\times 3$ case, this means the last row is simplified to involve only $z$, the second row involves $y$ and $z$, and the first row keeps $x,y,z$. Once reaching **Row Echelon Form**, we can solve it by **back substitution**, just as in the traditional method. 
$$
\left[\begin{array}{ccc|c}
  0 & 2 & 3 & 8 \\
  2 & 4 & 1 & 7 \\
  1 & 1 & 2 & 5
\end{array}\right]
\xrightarrow{r_1 \leftrightarrow r_2}
\left[\begin{array}{ccc|c}
  2 & 4 & 1 & 7 \\
  0 & 2 & 3 & 8 \\
  1 & 1 & 2 & 5
\end{array}\right]
\xrightarrow{r_3 = 2r_3 - r_1}
\left[\begin{array}{ccc|c}
  2 & 4 & 1 & 7 \\
  0 & 2 & 3 & 8 \\
  0 & -2 & 3 & 3
\end{array}\right]
\xrightarrow{r_3 = r_3 + r_2}
\left[\begin{array}{ccc|c}
  2 & 4 & 1 & 7 \\
  0 & 2 & 3 & 8 \\
  0 & 0 & 6 & 11
\end{array}\right]
$$

In Row Echelon Form, the system becomes much easier to read. <u>We define that the first non-zero entry in each non-zero row called a **pivot**</u>.

 >But can we simplify the matrix further ?

Instead of stopping after eliminating the entries below the diagonal, **Gauss-Jordan Elimination** continues the reduction by also <u>eliminating the entries above the pivots and normalizing them to $1$.</u> In this way, the matrix is transformed into **Reduced Row Echelon Form (RREF)**, a much cleaner form from which the solution can be read off more directly.

```math
\left[\begin{array}{ccc|c}
  2 & 4 & 1 & 7 \\
  0 & 2 & 3 & 8 \\
  0 & 0 & 6 & 11
\end{array}\right]
\xrightarrow{\frac{1}{6}r_3}
\left[\begin{array}{ccc|c}
  2 & 4 & 1 & 7 \\
  0 & 2 & 3 & 8 \\
  0 & 0 & 1 & \frac{11}{6}
\end{array}\right]
\xrightarrow[r_2 \leftarrow r_2 - 3r_3]{r_1 \leftarrow r_1 - r_3}
\left[\begin{array}{ccc|c}
  2 & 4 & 0 & 7 - \frac{11}{6} \\
  0 & 2 & 0 & 8 - \frac{11}{2} \\
  0 & 0 & 1 & \frac{11}{6}
\end{array}\right]
```

```math
= \left[\begin{array}{ccc|c}
  2 & 4 & 0 & \frac{31}{6} \\
  0 & 2 & 0 & \frac{5}{2} \\
  0 & 0 & 1 & \frac{11}{6}
\end{array}\right]
\xrightarrow{r_1 \leftarrow r_1 - 2r_2}
\left[\begin{array}{ccc|c}
  2 & 0 & 0 & \frac{31}{6} - 5 \\
  0 & 2 & 0 & \frac{5}{2} \\
  0 & 0 & 1 & \frac{11}{6}
\end{array}\right]
= \left[\begin{array}{ccc|c}
  2 & 0 & 0 & \frac{1}{6} \\
  0 & 2 & 0 & \frac{5}{2} \\
  0 & 0 & 1 & \frac{11}{6}
\end{array}\right]\xrightarrow{\frac{1}{2}r_1,\; \frac{1}{2}r_2}
\left[\begin{array}{ccc|c}
  1 & 0 & 0 & \frac{1}{12} \\
  0 & 1 & 0 & \frac{5}{4} \\
  0 & 0 & 1 & \frac{11}{6}
\end{array}\right]
```

---

<h3 align="center">Inverse Matrix</h3> 

After Gaussian-Jordan Elimination, we may ask a more fundamental question:

>When can the system $A\mathbf{x}=\mathbf{b}$ be solved uniquely and directly for every $\mathbf{b}$ ?

This leads to the notion of the **inverse matrix.**  <u>A square matrix $A_{m\times n}$  is called invertible if there exists a matrix</u> $A^{-1}$ <u>such that</u>

```math
AA^{-1} = I\quad \text{and}\quad A^{-1}A = I
```

In that case, multiplying both sides of $A\mathbf{x = b}$ by $A^{-1}$ gives

```math
\mathbf{x} = A^{-1}\mathbf{b}
```

So the inverse matrix provides a direct algebraic way to solve the system.

> [!TIP]
>
> **Def:** For the matrix $A_{m\times n}$ is invertible if there exists a matrix $B$ that satisfies
>
> ```math
> BA = I \quad \text{and} \quad AB = I
> ```
>
> **Claim:** suppose $A$ is invertible, then its inverse is unique.
>
> - Proof:
>
>     - Suppose $A$ has two invertible matrix, then $BA=I$, and $AC=I$. We have
>         
>         ```math
>         B = BI = B(AC) = (BA)C = C
>         ```
>
> **Claim:** suppose there is a non-zero solution $\mathbf{x}$ to $A\mathbf{x}=\mathbf{0}$, then $A$ is not invertible.
>
> - Proof:
>
>     - If $A$ is invertible, then we have the following that contradicts the assumption of non-zero solution.
>         
>         ```math
>         A\mathbf{x}=\mathbf{0}\Rightarrow\mathbf{x}=A^{-1}\mathbf{0} = \mathbf{0}
>         ```
>
> **Claim:** A diagonal matrix has an inverse provided no entries are zero.
>
> - ```math
>     A = 
>     \begin{pmatrix}
>       d_{1} & 0 & \cdots & 0 \\
>       0 & d_{2} & \cdots & 0 \\
>       \vdots & \vdots & \ddots & \vdots \\
>       0 & 0 & \cdots & d_{n}
>     \end{pmatrix} \quad A^{-1} = 
>     \begin{pmatrix}
>       \frac{1}{d_{1}} & 0 & \cdots & 0 \\
>       0 & \frac{1}{d_{2}} & \cdots & 0 \\
>       \vdots & \vdots & \ddots & \vdots \\
>       0 & 0 & \cdots & \frac{1}{d_{n}}
>     \end{pmatrix}
>     ```
>
> **Claim:** If $A$ and $B$ are invertible, then so is $AB$, and $(AB)^{-1}=B^{-1}A^{-1}$.
>
> - Proof: 
>
>     - The idea is to show that $AB$ is invertible by constructing its inverse. It is enough to find a matrix that serves as both a left inverse and a right inverse. Since $A$ and $B$ are invertible, a natrual candidate is 
>         
>         ```math
>         X=B^{-1}A^{-1}
>         ```
>         
>     - Then we can verify that
>         
>         ```math
>         (AB)(B^{-1}A^{-1}) = I \quad (B^{-1}A^{-1})(AB)=I
>         ```
>         
>     - Therefore, $B^{-1}A^{-1}$ is the inverse of $AB$.



The definition of the inverse matrix becomes much more concrete from the viewpoint of Gauss-Jordan Elimination. We can view it as the matrix produced when $A$ is reduced to the identity. To do this, we form the augmented matrix:

```math
\left[ A \text{ }|\text{ }I\right]
```

we perform row operations until the left side becomes the identity matrix. If this is possible, then the right side is exactly $A^{-1}$:

```math
\left[ A \text{ }|\text{ }I\right] \rightarrow \left[ I \text{ }|\text{ }A^{-1}\right]
```

If the left side cannot be reduced to $I$, then $A$ is not invertible.

---

<h3 align="center">Row Operations in Matrix Form</h3> 

>How can we describe the Gauss-Jordan elimination process algebraically ?

Since the whole procedure consists of repeated row operations, we would like to express this as matrix multiplication.

- <u>If we perform a single elementary row operation on the identity matrix $I$, the resulting matrix is called an **elementary matrix**</u> $E$;
- <u>If we exchange two rows of $I$, the resulting matrix is called a **permutation matrix**</u> $P$.

<u>Besides, $E$ and $P$ are invertible matrix, since each such operation can be undone by another row operation of the same type.</u>

> [!NOTE]
>
> ```math
> A=
> \begin{bmatrix}
> 1 & 2 & 3\\
> 4 & 5 & 6\\
> 7 & 8 & 9
> \end{bmatrix}
> ```
>
> Suppose we plan to exchange $R_1$ and $R_3$, 
>
> ```math
> A: R_1\leftrightarrow R_3
> \Leftrightarrow
> PA=
> \begin{bmatrix}
> 0 & 0 & 1\\
> 0 & 1 & 0\\
> 1 & 0 & 0
> \end{bmatrix}
> \begin{bmatrix}
> 7 & 8 & 9\\
> 4 & 5 & 6\\
> 1 & 2 & 3
> \end{bmatrix}
> =
> \begin{bmatrix}
> 7 & 8 & 9\\
> 4 & 5 & 6\\
> 1 & 2 & 3
> \end{bmatrix}
> ```
>
> Suppose we plan to reduce 4 multiples of $R_1$ from $R_2$,
>
> ```math
> A: R_2\rightarrow R_2-4R_1
> \Leftrightarrow
> E_{21}(-4)\cdot A=
> \begin{bmatrix}
> 1 & 0 & 0\\
> -4 & 1 & 0\\
> 0 & 0 & 1
> \end{bmatrix}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 4 & 5 & 6\\
> 7 & 8 & 9
> \end{bmatrix}
> =
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & -3 & -6\\
> 7 & 8 & 9
> \end{bmatrix}
> ```

From these examples, we can see that Gaussian Elimination can be viewed as a sequence of left multiplications by elementary matrices $E$ and permutation matrices $P$, i.e.,

```math
MA = U
```

where $M$ is a product of elementary matrices and permutation matrices.

> [!TIP]
>
> **Def:** $A_{n\times n}$ matrix is non-singular if it has full set of $n$ (non-zero) pivots.
>
> Claim: $A$ is invertible if and only if it is non-singular.
>
> - Proof
>
>     - $\Rightarrow$ Suppose $A$ is singular. Then after Gauss Elimination, $A$ can be reduced to a matrix with at least one zero row. Hence, there exists an invertible matrix $M$ such that $MA$ has a zero row.
>
>         - Now suppose $A$ is invertible, then there would exist $A^{-1}$ such that
>             
>             ```math
>             AA^{-1} = I
>             ```
>             
>         - Multiplying on the left by $M$, we get
>             $$
>             (MA)A^{-1} = M
>             $$
>     
>         - Since $MA$ has a zero row, the product must also have a zero row. Therefore, $M$ has a zero row, which is impossible because $M$ is invertible. This contradiction shows that $A$ is not invertible.
>     
>     - $\Leftarrow$ Suppose $A$ is non-singular
>     
>         - By Gauss-Jordan Elimination, we can find a matrix $B$ such that $AB = I$, since $Ax_i = \mathbf{e_i}$ is solvable for all $i=1,...,n$;
>         - Gauss-Jordan Elimination is really a sequence of multiplication elementary matrices on the left, which can be denoted as $D^{-1}(E...P..E)A=I$ where (1) $E_{ij}$ to subtract a mulitple of $l_{ij}$ of row $j$ from row $i$, and for (2) $P_{ij}$ denotes to exchange row $j$ and row $i$, and for (3) $D^{-1}$ denotes to divide all rows by their pivots.
>         - That is, there is a matrix $C$ that satisfies $CA=I$. Therefore, $B=C=A^{-1}$ and $A$ is invertible.
>
> 

---

<h3 align="center">LU factorization</h3> 

The Gauss elimination process can be written as a sequence of left multiplication by elementary matrices:

```math
MA=U\quad \text{suppose } M = E_kE_{k-1}\cdots E_2E_1A
```

Since each elementary matrix $E$ is invertible, their product $M$ is invertible as well (proved before). Therefore,

```math
A = M^{-1}U
```

Let explore what the $M^{-1}$ looks like. We first look at a single elimination step

```math
R_i\leftarrow R_i - l_{ij}R_j
```

then the elementary matrix is

```math
E_{ij}(-l_{ij}) = I - l_{ij}e_{ij}
```

and its inverse is

```math
E_{ij}(-l_{ij})^{-1}=E_{ij}(l_{ij}) = I + l_{ij}e_{ij}
```

So the inverse keeps the same position, but reverses the sign of the elimination multiplier. Since each $E^{-1}$ is lower triangular, and the product of lower triangular matrices is still lower triangular, $M^{-1}$ is also lower-triangular.

```math
L = M^{-1}
```

and obtain the factorization
```math
A = LU
```

If all diagonal entries of $U$ are non-zero, then we can further refine $LU$ factorization into $LDU$ factorization.

We observe that $L$ already has $1$'s on the diagonal. To make the factorization more symmetric, we can write $U=DU'$ as long as each row of $U$ can be obtained by scaling a row whose diagonal entry is $1$. **[If a diagonal entry is zero, the corresponding row must be a zero row]**

```math
U = DU^{’}
```

where $D$ is diagonal and $U^{’}$ is an upper triangular matrix with $1$'s on its diagonal. We can obtain

```math
A = LDU^{’}
```

> [!NOTE]
>
> Example: 
>
> We now illustrate how Gaussian Elimination leads to LU factorization. Let us starting from
>
> ```math
> A=
> \begin{bmatrix}
> 1 & 2 & 3\\
> 4 & 5 & 6\\
> 7 & 8 & 9
> \end{bmatrix}
> ```
>
> we apply elimination to transform $A$ into an upper triangular matrix $U$.
>
> ```math
> \begin{bmatrix}
> 1 & 2 & 3\\
> 4 & 5 & 6\\
> 7 & 8 & 9
> \end{bmatrix}
> \xrightarrow{\,R_2 \leftarrow R_2 - 4R_1\,}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & -3 & -6\\
> 7 & 8 & 9
> \end{bmatrix}
> \xrightarrow{\,R_3 \leftarrow R_3 - 7R_1\,}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & -3 & -6\\
> 0 & -6 & -12
> \end{bmatrix}
> \xrightarrow{\,R_3 \leftarrow R_3 - 2R_2\,}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & -3 & -6\\
> 0 & 0 & 0
> \end{bmatrix}
> =U
> ```
>
> Each elimination step can be represented by left multiplication by an **elementary matrix**, 
>
> ```math
> E_{21}=
> \begin{bmatrix}
> 1 & 0 & 0\\
> -4 & 1 & 0\\
> 0 & 0 & 1
> \end{bmatrix},
> \qquad
> E_{31}=
> \begin{bmatrix}
> 1 & 0 & 0\\
> 0 & 1 & 0\\
> -7 & 0 & 1
> \end{bmatrix},
> \qquad
> E_{32}=
> \begin{bmatrix}
> 1 & 0 & 0\\
> 0 & 1 & 0\\
> 0 & -2 & 1
> \end{bmatrix}
> ```
>
> So the whole elimination process can be written as
>
> ```math
> E_{32}E_{31}E_{21}A=U
> ```
>
> Let
>
> ```math
> M = E_{32}E_{31}E_{21}
> ```
>
> Then $MA = U$, and since $M$ is invertible, we obtain
>
> ```math
> A = M^{-1}U
> ```
>
> Now each inverse elementary matrix keeps the same lower-triangular position, but reverses the sign of the elimination multiplier. Hence,
>
> ```math
> E_{21}^{-1}=
> \begin{bmatrix}
> 1 & 0 & 0\\
> 4 & 1 & 0\\
> 0 & 0 & 1
> \end{bmatrix},
> \qquad
> E_{31}^{-1}=
> \begin{bmatrix}
> 1 & 0 & 0\\
> 0 & 1 & 0\\
> 7 & 0 & 1
> \end{bmatrix},
> \qquad
> E_{32}^{-1}=
> \begin{bmatrix}
> 1 & 0 & 0\\
> 0 & 1 & 0\\
> 0 & 2 & 1
> \end{bmatrix}
> ```
>
> ```math
> M^{-1}=E_{21}^{-1}E_{31}^{-1}E_{32}^{-1}
> =
> \begin{bmatrix}
> 1 & 0 & 0\\
> 4 & 1 & 0\\
> 7 & 2 & 1
> \end{bmatrix}
> =L
> ```
>
> Therefore,
>
> ```math
> A = LU
> \Leftrightarrow
> \begin{bmatrix}
> 1 & 2 & 3\\
> 4 & 5 & 6\\
> 7 & 8 & 9
> \end{bmatrix}
> =
> \begin{bmatrix}
> 1 & 0 & 0\\
> 4 & 1 & 0\\
> 7 & 2 & 1
> \end{bmatrix}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & -3 & -6\\
> 0 & 0 & 0
> \end{bmatrix}
> ```
>
> ```math
> U=
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & -3 & -6\\
> 0 & 0 & 0
> \end{bmatrix}
> =
> \begin{bmatrix}
> 1 & 0 & 0\\
> 0 & -3 & 0\\
> 0 & 0 & 0
> \end{bmatrix}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & 1 & 2\\
> 0 & 0 & 1
> \end{bmatrix}
> =DU'
> ```
>
> ```math
> A=
> \begin{bmatrix}
> 1 & 0 & 0\\
> 4 & 1 & 0\\
> 7 & 2 & 1
> \end{bmatrix}
> \begin{bmatrix}
> 1 & 0 & 0\\
> 0 & -3 & 0\\
> 0 & 0 & 0
> \end{bmatrix}
> \begin{bmatrix}
> 1 & 2 & 3\\
> 0 & 1 & 2\\
> 0 & 0 & 1
> \end{bmatrix}
> =LDU'
> ```
>
> 
>
> In this example, the last diagonal entry of $U$ is zero. This is still acceptable, since the last row of $U$ is already a zero row.



> [!TIP]
>
> **Claim:**  If $A = L_1D_1U_1$ and $A = L_2D_2U_2$ where the $L_s$ are lower triangular matrix with unit diagonal. The $U_s$ are upper-triangle matrix with unit diagonal and $D_s$ are the diagonal matrix with non-zero on the diagonal, then $L_1 = L_2, D_1=D_2,U_1=U_2$. **[LU factorization uniqueness]** 
>
> - Proof [Homework, I will update later.]
>
> **Claim:** If a symmetric matrix is factored into $LDU$ with no row exchanges, then $U=L^T$
>
> - Proof
>     - $A = LDU = A^T = (LDU)^T = U^TDL^T$
>     - Since the factorization is unique, we have $U = L^T$ and $L = U^T$.
>     - Here the transpose of a matrix means row columns exchange i.e, $(A^T)_{ij} = A_{ji}$. I will explain this notion in the next chapter vector space and subspaces 😃.
>
> **Claim:** If row exchanges are needed, how can we execute the factorization?
>
> - If row exchanges are needed during elimination, we can do them in advance. The product $PA$ will be in right order. So that no exchange needed for $PA$, which means $PA = LU$;
>
> - If we hold row exchanges until after elimination, we then have
>     
>     ```math
>     L^{-1}A=PU\Rightarrow A = LPU
>     ```
>
> ---
>
> <h6 align="center">One Square System = Two Triangular Systems</h6> 
>
> Suppose Gaussian Elimination requires no row exchanges, which means 
>
> ```math
> A = LU \Rightarrow A\mathbf{x} = \mathbf{b} \Rightarrow LU\mathbf{x} = \mathbf{b} \Rightarrow L\mathbf{c}=\mathbf{b}\Rightarrow U\mathbf{x} = \mathbf{c}
> ```
>
> We have $U\mathbf{x} = \mathbf{c}$ where $L\mathbf{c} = \mathbf{b}$. The computation steps is (1) Factor $A=LU$ by Gaussian Elimination; (2) Solve $L\mathbf{c}=\mathbf{b}$ and then $U\mathbf{x} = \mathbf{c}$.
>
> ```math
> A = \left[
> \begin{matrix}
> a_{11} & \cdots & a_{1n}\\ a_{21} & \cdots &  a_{2n} \\ \vdots & \vdots & \vdots \\ a_{n1} & \cdots & a_{nn} \\
> \end{matrix}
> \right] = LU
> ```
>
> To eliminate $a_{21}$, we first compute the multiplier
>
> ```math
> l_{21} = \frac{a_{21}}{a_{11}}
> ```
>
> and then perform the row update
>
> ```math
> R_2\leftarrow R_2 - l_{21}R_1
> ```
>
> Now let us estimate the computation cost for one row update:
>
> - Computing the multiplier $l_{21}$ cost: $1$ division
> - Updating the remaining $(n-1)$ entries costs: $(n-1)$ multiplications, and $(n-1)$ subtractions.
>
> So the total cost for eliminating $a_{21}$ is $1+2(n-1)$. This can be generalized to $a_{i1}$. So to eliminate the entries in the first column, we need to compute $(n-1)(1+2(n-1))$. Similarly, there are $(n-k)$ rows below the pivot $a_{kk}$, so the total cost for step $k$ is
>
> ```math
> (n-k)\bigl(1+2(n-k)\bigr)
> ```
>
> Summing over all pivot columns, we obtain
>
> ```math
> \sum_{k=1}^{n-1}(n-k)\bigl(1+2(n-k)\bigr)\Rightarrow \frac{n(n-1)}{2}
> + \frac{(n-1)n(2n-1)}{3}
> ```
>
> Its leading term is $\frac{2}{3}n^3$. So Gaussian elimination / LU factorization has complexity
>
> ```math
> O(n^3)
> ```
>
> Next, to solve $L\mathbf{c=b}$, we use forward substitution since $L$ is lower triangular matrix.
>
> From the first row, we can immediately get the following, which requires no arithmetic operation:
>
> ```math
> c_1 = b_1
> ```
>
> From the second row, we can compute the following, which requires 1 multiplication and one subtraction:
>
> ```math
> p_{21}c_1 + c_2 = b_2\Rightarrow c_2 = b_2 - p_{21}c_1
> ```
>
> In general, from the $k$-th row, we can obtain the following that requires $k-1$ multiplication and $k-1$ additions.
>
> ```math
> p_{k1}c_1+\cdots+p_{k,k-1}c_{k-1}+c_k=b_k\Rightarrow c_k=b_k-\sum_{j=1}^{k-1}p_{kj}c_j
> ```
>
> Therefore, the total cost of solving $L\mathbf{c=b}$ is
> ```math
> 2\sum_{k=1}^{n-1} k
> =2\cdot \frac{n(n-1)}{2}
> =n(n-1),\quad \text{which is }O(n^2)
> ```
>
> After solving $L\mathbf{c=b}$, we need to solve $U\mathbf{x}=\mathbf{c}$. Since $U$ is an upper triangular matrix, this step is done by back substitution. Its complexity is essentially the same as that of solving $L\mathbf{c=b}$. Indeed, each row requires a number of arithmetic operations proportional to the number of already known entries, so the total cost is again of order
>
> ```math
> O(n^2)
> ```
>
> 
>
> ---
>
> Combining all three steps, we can get
>
> ```math
> \begin{aligned}
> &A = LU      & \text{cost} &\sim O(n^3) \\
> &Lc = b      & \text{cost} &\sim O(n^2) \\
> &Ux = c      & \text{cost} &\sim O(n^2)
> \end{aligned}
> ```
>
> then the total cost is
>
> ```math
> O(n^3) + O(n^2) + O(n^2) = O(n^3)
> ```
>
> Therefore, the overall complexity is dominated by the $LU$ factorization step. 
>
> However, if we need to solve $A\mathbf{x=b}$ for many different right-hand $\mathbf{b}$, then factorization is done only once. After that, each new system only requires solving the latter two steps, which costs only $O(n^2)$.

---

This version mainly aims to clarify the main storyline of the chapter.  More detailed exercises, theorem proofs, and visual summaries will be added in future updates. 
The next chapter will shift the focus from **how to solve \(Ax=b\)** to **how to understand its structure**, leading naturally to **vector spaces** and **subspaces**.



