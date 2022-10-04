---
title: 数字信号处理Chapter01
date: 2022-10-04 23:01:23
categories:
- [Learn]
tags: 
- [Review]
- [Math]
---
## 第一章 离散时间信号与系统

复习应试
<!-- more -->
### 1.1 符号表示及基础

离散时间信号通常用序列：
$ \{x(n)\}$ ，$ n $ 为 $ 0,1,2 ...$ , $ x(n) $ 表示为序列中第 $ n $ 个样本值。

$ \{·\} $ 表示全部样本值的集合

$ \{x*(n)\}$ 表示复序列的共轭

连续时间序列 $ \{t\} $ 与离散时间序列 $ \{x(n)\}$ 的关系：

$$
x(n) = x_a(t) |_{t = nT} =x_a(nT) \tag {1.1}
$$
其中采样频率$ f_s = \frac{1}{T} $（T为采样周期，即两个样本间的时间间隔）

周期序列表示为 $ \widetilde{x}(n) $ 
其中
$$
\widetilde{x}(n) = x(n+kN) , 0 \leq n \leq N-1 ,k为任意整数 \tag{1.2}
$$

#### 1.1.1 常见典型序列

1. 单位脉冲序列
$$
\delta (n)=\left \{
\begin{aligned}
1, n = 0\\
0, n \neq 0 \\
\end{aligned}
\right. \tag{1.3}
$$

2. 单位阶跃序列
$$
u (n)=\left \{
\begin{aligned}
1, n \geq 0\\
0, n < 0 \\
\end{aligned}
\right. \tag{1.4}
$$

3. 矩形序列
$$
R_N (n)=\left \{
\begin{aligned}
1, n \leq n \leq N-1 \\
0, n < 0,n \geq N \\
\end{aligned}
\right. \tag{1.5}
$$

4. 实指数序列
$$
x(n) = a^n u(n) \tag{1.6}
$$

$ a \neq 0, |a| < 1 $ 时收敛，$|a| \geq 1$ 时发散

5. 正弦序列
$$
x(n) = sin(\omega_0n)
\tag{1.7}
$$
$\omega_0$为数字角频率，单位为弧度 $rad$

6. 复指数序列
$$x(n) = (re^{j\omega_0})^n = r^n[cos(\omega_0n)+jsin(\omega_0n)] 
\tag{1.8}
$$

#### 1.1.2 序列的运算

1. 序列的加法
$$
z(n) = x(n) + y(n) \tag{1.9}
$$

2. 序列的相乘
$$
z(n) = x(n)  y(n) \tag{1.9}
$$

3. 序列的位移
$$
z(n) = x(n-n_0) \tag{1.9}
$$
当 $ n_0 > 0 $ 时 $z(n)$ 是 $ x(n) $ 的延迟；当 $ n_0 < 0 $ 时 $z(n)$ 超前于 $ x(n) $ ；

4. 序列的能量及序列的绝对值
序列的能量定义为序列样本值的平方和
$$
S = \sum^{\infty}_{n = -\infty} |x(n)|^2 
\tag{1.10}
$$
如果序列 $x(n)$ 满足 $S < \infty$ 则为平方可和序列
如果序列满足
$$
\sum^{\infty}_{n = -\infty} |x(n)| < \infty \tag{1.11}
$$
则为绝对可和序列
如果序列的每一个样本值的绝对值均小于某一个有限的正整数 $B_x$ 则 $x(n)$ 为有界序列，即
$$
|x(n)| \leq B_x < \infty 
\tag{1.12}
$$

5. 实序列的偶部和奇部
任何序列均可以分解成偶对成序列和奇对称序列的和的形式，即
$$
x(n) = x_e(n) + x_o(n)
\tag{1.13}
$$
$x_e(n)$ 和 $x_o(n)$ 分别称为 $x(n)$ 的偶部和基部，其分别等于
$$
x_e(n) = \frac{1}{2}[x(n) + x(-n)]
\tag{1.14a}
$$
$$
x_o(n) = \frac{1}{2}[x(n) - x(-n)]
\tag{1.14b}
$$

6. 任意序列的单位脉冲表示
任一序列 $x(n)$ 都可以表示成单位脉冲序列移位的加权和，即
$$
x(n) = \sum^{\infty}_{m = -\infty}x(m)\delta(n-m)
\tag{1.15}
$$

### 1.2 离散时间信号的傅里叶变换与 $ \mathscr{z}$ 变换

#### 1.2.1 离散时间信号的傅里叶变换

离散时间傅里叶变换 $ DTFT $ (discrete-time Fourier tansform) ,序列的 $DTFT$ 定义为：
$$
X(e^{j\omega}) = \sum^{\infty}_{n = -\infty}x(n)e^{-j\omega n},\omega = \frac{2\pi f}{f_S}
\tag{1.16}
$$
式中， $ \omega $ 为数字角频率，它是频率 $f$ 对采样频率 $f_s$ 作归一化后的角频率。
