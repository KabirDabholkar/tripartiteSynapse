There are a total of 54 reactions that involve Ca2+ binding to a state of the sensor and forming a different state.
$$Ca^{2+}+S_1 \rightarrow S_2$$

Let $R_i$ be the event that in period of $\Delta T$, the ith reaction takes place at least once.
$$r_i\Delta T=P(R_i)=k_i[Ca^{2+}][S_i]\Delta T$$

Let $B_N$ be the event that exactly N reactions occur in $\Delta T$.


Probability that no reactions occur in $\Delta T$:
$$P(B_{N=0})=lim_{K\to \infty}\big(1-\sum_i r_i\epsilon\big)^K=e^{-\Delta T R_{tot}}$$


Let $r_{bind}=\sum r_i$, the sum of all the binding reaction rates.

Now,
$$P(R_i|B_{N=0})=\frac{P(R_i \cap B_{N=0})}{P(B_{N=0})}=0$$

Similarly we can calculate $P(B_{N=1})$:
\begin{align*}
	P(B_{N=1})&=\sum_i P(R_i) lim_{K \to \infty} \big[ 1-\sum_{j\neq i}r_j \epsilon  \big]^K \\
		  &=\sum_i r_i \Delta T e^{-\Delta T(r_{bind}-r_i)} \\
	P(R_l|B_{N=1})&=\frac{P(R_l \cap B_{N=1})}{P(B_{N=1})}\\
		  &=\frac{r_l\Delta T e^{-\Delta T(R_{bind}-R_l)}}{\sum_i r_i \Delta T e^{-\Delta T(r_{bind}-r_i)}} \\
		  &=\frac{k_l[S_l]}{\sum_i k_i [S_i] e^{-\Delta T(r_l-r_i)}}
\end{align*}
We can approximate the exponential to one.

This extends to $N=2$ as follows:
\begin{align*}
	P(B_{N=2})&=\sum_{i,j} P(R_i) P(R_j) \lim_{K\to \infty} \big[ 1-\sum_{k\neq i,j}r_k \epsilon  \big]^K \\
		  &=\sum_{i,j} r_i r_j \Delta T e^{-\Delta T(r_{bind}-r_i-r_j)} \\
	P(R_m,R_n|B_{N=2})&=\frac{r_m r_n (\Delta T)^2 e^{-\Delta T(R_{bind}-R_m-R_n)}}{\sum_{i,j} r_i r_j (\Delta T)^2 e^{-\Delta T(R_{bind}-r_i-r_j)}} \\
		  &=\frac{k_m k_n[S_m][S_n]}{\sum_{i,j} k_i k_j[S_i][S_j] e^{-\Delta T((r_m+r_n)-(r_i+r_j))}} \\
		  &\approx \frac{k_m k_n[S_m][S_n]}{\sum_{i,j} k_i k_j[S_i][S_j]}
\end{align*}

We can similarly extend this to higher Ns. In extremely rare circumstances, N=5 can occur.
We are assuming the independence of the $R_i$s.


