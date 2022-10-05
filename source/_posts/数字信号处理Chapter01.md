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
<!--more-->
### 1.1 符号表示及基础

离散时间信号通常用序列：
$ \{x(n)\}$ ，$ n $ 为 $ 0,1,2 ...$ , $ x(n) $ 表示为序列中第 $ n $ 个样本值。

$ \{·\} $ 表示全部样本值的集合

$ \{x*(n)\}$ 表示复序列的共轭

连续时间序列 $ x\{t\} $ 与离散时间序列 $ \{x(n)\}$ 的关系：

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
z(n) = x(n)  y(n) \tag{1.10}
$$

3. 序列的位移
$$
z(n) = x(n-n_0) \tag{1.11}
$$
当 $ n_0 > 0 $ 时 $z(n)$ 是 $ x(n) $ 的延迟；当 $ n_0 < 0 $ 时 $z(n)$ 超前于 $ x(n) $ ；

4. 序列的能量及序列的绝对值
序列的能量定义为序列样本值的平方和
$$
S = \sum^{\infty}_{n = -\infty} |x(n)|^2
\tag{1.12}
$$
如果序列 $x(n)$ 满足 $S < \infty$ 则为平方可和序列
如果序列满足
$$
\sum^{\infty}_{n = -\infty} |x(n)| < \infty \tag{1.13}
$$
则为绝对可和序列
如果序列的每一个样本值的绝对值均小于某一个有限的正整数 $B_x$ 则 $x(n)$ 为有界序列，即
$$
|x(n)| \leq B_x < \infty
\tag{1.14}
$$

5. 实序列的偶部和奇部
任何序列均可以分解成偶对成序列和奇对称序列的和的形式，即
$$
x(n) = x_e(n) + x_o(n)
\tag{1.15}
$$
$x_e(n)$ 和 $x_o(n)$ 分别称为 $x(n)$ 的偶部和基部，其分别等于
$$
x_e(n) = \frac{1}{2}[x(n) + x(-n)]
\tag{1.15a}
$$
$$
x_o(n) = \frac{1}{2}[x(n) - x(-n)]
\tag{1.15b}
$$

6. 任意序列的单位脉冲表示
任一序列 $x(n)$ 都可以表示成单位脉冲序列移位的加权和，即
$$
x(n) = \sum^{\infty}_{m = -\infty}x(m)\delta(n-m)
\tag{1.16}
$$

### 1.2 离散时间信号的傅里叶变换与 $ \mathscr{z} $ 变换

#### 1.2.1 离散时间信号的傅里叶变换

离散时间傅里叶变换 $ DTFT $ (discrete-time Fourier tansform) ,序列的 $DTFT$ 定义为：
$$
X(e^{j\omega}) = \sum^{\infty}_{n = -\infty}x(n)e^{-j\omega n},\omega = \frac{2\pi f}{f_S}
\tag{1.17}
$$
式中， $ \omega $ 为数字角频率，它是频率 $f$ 对采样频率 $f_s$ 作归一化后的角频率。
$X(e^{j\omega})$ 时 $\omega$ 的连续函数，且周期为 $2\pi$
式（$1.17$）级数不一定总是收敛的，当 $x(n)$ 绝对可和时，它的 $DTFT$ 一定存在。
离散时间信号的傅里叶逆变换（$IDTFT$）：
$$
x(n) = \frac{1}{2\pi}\int_{-\pi}^{\pi} X(e^{j\omega})e^{j\omega m} d\omega
\tag{1.18}
$$
$x(n)$ 和 $X(e^{j\omega})$ 对应关系可表示为：$X(e^{j\omega}) = DTFT[x(n)]$ ,$x(n)=IDTFT[X(e^{j\omega})]$

$X(e^{j\omega})$ 的几种表示方法：
$$
X(e^{j\omega}) = Re[X(e^{j\omega})]+jIm[X(e^{j\omega})] = |X(e^{j\omega})|e^{j\phi(\omega)}
\tag{1.19}
$$
$Re[·]$ 和 $Im[·]$ 表示取实部和虚部。
$|X(e^{j\omega})|$ 为离散序列 $x(n)$ 的幅度谱，$\phi(\omega)$为离散序列的相位谱。

$DTFT$ 的主要特性
|序列|$DTFT$|
|:---:|:---:|
|$ax(n)+by(n)$|$aX(e^{j\omega})+Y(e^{j\omega})$|
|$x^*(n)$|$X^*(e^{-j\omega})$|
|$x^*(-n)$|$X^*(e^{j\omega})$|
|$x(n-n_0)$|$e^{-jn_0\omega}X(e^{j\omega})$|
|$e^{j\omega_0 n}x(n)$|$X(e^{j(\omega - \omega_0)})$|
|$Re[x(n)]$|$X_e(e^{j\omega})$  [$X(e^{j\omega})$ 的共轭偶对称部分]|
|$jIm[x(n)]$|$X_o(e^{j\omega})$ [$X(e^{j\omega})$ 的共轭奇对称部分]|
|$x(n)$ 为实序列|$X(e^{j\omega}) = X^*(e^{-j\omega})$
||$Re[X(e^{j\omega})] = Re[X(e^{-j\omega})]$|
||$Im[X(e^{j\omega})] = -Im[X(e^{-j\omega})]$|
||$arg[X(e^{j\omega})] = -arg[X(e^{-j\omega})]$|
|$x_e(n)$ [$x(n)$ 的共轭偶对称部分]|$Re[X(e^{j\omega})]$|
|$x_o(n)$  [$x(n)$ 的共轭偶奇称部分]|$jIm[X(e^{j\omega})]$|

#### 1.2.2 $\mathscr{z}$变换

序列 $x(n)$ 的 $\mathscr{z}$ 变换定义为：
$$
X(z) = \sum^{\infty}_{n = -\infty}x(n)z^{-n} ,(n = 0时为单边z变换)
\tag{1.20}
$$
上式中 $z$ 为复变量，也可记为 $\mathscr{Z}[x(n)] = X(z)$
对于所有的序列或所有的 $z$ 值，$z$变换并不总是收敛，使 $z$ 变换收敛的 $z$ 值的集合称作收敛区域，一般为 $z$ 平面上的一个环形区域，该区域为:
$$
R_{x^-} <|z|<R_{x^+}
\tag{1.21}
$$
其中 $R_{x^-}$ 可以小到0，$R_{x^+}$ 可以大到 $\infty$

以下讨论几种序列的收敛域

1. 有限长序列
仅有有限个数的序列值是非零值，从而
