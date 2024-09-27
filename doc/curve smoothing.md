# Curve Smoothing

三角网格上的曲线（表示为带拓扑的点集）通常带有锯齿，不平滑。在网格上的想要得到较为平滑的曲线通常有两种做法（我了解的），一种是基于连续曲线差值[^1]，第二种是在曲面上建立平滑的标量场提取等值线[^2]。

## 优化形式

给定网格 $\mathcal M$ 和曲线 $\Omega\subset\mathcal M$，以曲线作为约束，求解网格上的标量场 $\phi$ 是我们的优化目标，即
$$
\begin{equation}
\label{eq:target}
\begin{aligned}
\arg\min_{\phi:\mathcal M\rightarrow\mathbb R} E(\phi)\\
\text{s.t.}\quad \phi(\Omega)=g
\end{aligned}
\end{equation}
$$
在处理曲线优化的问题上，为了保证平滑性，我们逼近原始曲线即可。

故将原始曲线作为软约束带入到二次规划[^2]有
$$
\begin{equation}
\label{eq:soft_qud}
\min_\phi \left\|\Delta\phi-\nabla\cdot\mathbf X\right\|^2+\lambda\left\|\phi(\Omega)-g\right\|^2
\end{equation}
$$
网格上的拉普拉斯矩阵通常是半正定的，我们把上述的问题当作线性问题求解即是求解一个法方程。



## 附录

### 公式 $\ref{eq:soft_qud}$ 推导

这里我们考虑带面积权的最优化问题
$$
\min_\phi \left\|\Delta\phi-\nabla\cdot\mathbf X\right\|^2_M+\lambda\left\|\phi(\Omega)-g\right\|^2_M
$$
写成矩阵形式有
$$
\begin{aligned}
\mathbf{(Lx-b)^\top(Lx-b)+\lambda(Ax-g)^\top(Ax-g)}
\end{aligned}
$$
根据 The Matrix Cookbook[^3]的矩阵求导且令导数为0有
$$
(\mathbf{L^\top L+\lambda A^\top A)x=2L^\top b +2\lambda A^\top g}
$$




[^2]: A Heat Flow Based Relaxation Scheme for *n* Dimensional Discrete Hyper Surfaces
[^3]: The Matrix Cookbook
