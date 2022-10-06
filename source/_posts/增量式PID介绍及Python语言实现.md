---
title: 增量式PID介绍（Python语言实现）
date: 2022-07-08 15:01:57
categories:
- [Control]
- [Python]
tags: 
- Control
- Python
- PID
# banner_img: /image/pidzengliang/titu.jpg
--- 

## 增量式PID介绍及Python语言实现

### 摘要

本文主要记录在电子设计大赛训练时，做电源时需要对电路中某些输出量进行稳定控制，PID控制算法可以很好的实现该功能，其中增量式更符合对电源的电压电流量的控制。本文主要粗略对其进行介绍再使用Python实现一个电压输出稳定控制
<!-- more -->
### 1.增量式PID算法介绍

以下仅限与对增量式PID控制算法的介绍，对于PID控制算法的概述（General Introduction）见“PID控制器入门概要”

增量式PID控制将当前时刻的控制量和上一时刻的控制量做差，以差值为新的控制量，是一种递推式的算法。

其公式表示为：
$$
\Delta u(k) = K_p (e(k) - e(k-1)) + K_i e(k) + K_d(e(k) - 2e(k-1) + e(k-2)) \tag{1}
$$
其中:

- $K_p$为比例系数
- $K_i$为积分系数
- $K_d$为微分系数

因此若要使用增量式PID算法需要保存$e(k-1),e(k-2)$两个时刻的输出值，在加上输出回馈值$e(k)$的积分
注意：增量式PID的运算结果为$u(k)$同$u(k-1)$之间的差值，因此输出结果需叠加，即为：
$$
u(k) = u(k-1) + \Delta(k)
$$
![增量PID流程框图](/image/pidzengliang/zengliangpid.png)
厘清$u(k),u(k-1),\Delta u,e(k),e(k-1),e(k-2)$

### 2.使用Python语言进行实现

下面使用Python语言进行展示，得到输出曲线
假设如下案例，现有一DCDC升压（BOOST）电路，需要将15V的初始电压提升到30V，控制方式为驱动MOS管周期内导通与关断时间比值即可。但现在只关注电压输入与控制后输出的结果，不考虑PWM占空比，定时器计数等问题。

```Python
# 名称：PID算法（增量式）
# 作者：Shixin
# 更新时间：2022.07.08 21.55
# 版本: 1.0
import numpy
import pandas
import matplotlib.pyplot

def main():
    Setpoint = 30.0
    out = 0.0
    k_p = 0.05
    k_i = 0.2
    k_d = 0.1
    e_k1 = 0.0
    e_k2 = 0.0
    time = range(1, 50, 1)
    output = numpy.arange(1.0, 51.0, 1.0)

    for i in time:
        e_k = Setpoint - out
        outk = k_p * (e_k - e_k1) + k_i * e_k + k_d * (e_k - 2 * e_k1 + e_k2)
        e_k2 = e_k1
        e_k1 = e_k
        out = out + outk
        print(out)
        output[i] = out

    matplotlib.pyplot.plot(range(1, 51, 1), output)
    matplotlib.pyplot.show()

    print(output)

if __name__ == '__main__':
    main()

```

通过不断改变$K_p,K_i,K_d$三个值最终可以得到一个十分合适的输出曲线：
![PID输出曲线](/image/pidzengliang/pidout.png)
观察该曲线可以发现系统可以很达到了稳定状态，且到达稳定点的时间也很快
观察如下数值结果：
![PID输出数值结果](/image/pidzengliang/pidshuzhijieguo.png)
可以看到在数值上达到预期值的效果也十分理想，误差小于0.001

### 3.其他的一些问题

当我在做DC-DC电源转换器直流稳压源时，发现由于系统反应过于迅速以至于转换器来不及反应导致系统无法达到预期目的，因此需要将参数调小并留有一定的系统反应时间。
