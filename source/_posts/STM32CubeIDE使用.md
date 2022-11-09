---
title: STM32CubeIDE使用
date: 2022-10-30 12:58:09
categories:
    - [STM32]
    - [嵌入式]
tags:
    - STM32
    - 嵌入式
    - Eclipse
---

## STM32Cube使用基础（CubeMX代码生成）

STM32Cuben 开发套件有 STM32CubeIDE 和 STM32CubeMX 两个，其前者为基于 Eclipse 开发的新一代 STM32 开发集成环境（不过无代码提示😅）内置了 STM32CubeMX 但只能使用 HAL 库生成 STM32CubeIDE 工具链下的代码；而 STM32CubeMX 则只专注于代码生成，且可以生成 MDK-Keil, STM32CubeIDE, Makefile 等多种工具链下的代码。
