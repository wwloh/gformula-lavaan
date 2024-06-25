# Parametric g-formula using lavaan

This repository contains the R scripts for the applied example used for illustration, and simulation study, in the paper, *Parametric g-Formula for Testing Time-Varying Causal Effects: What It Is, Why It Matters, and How to Implement It in Lavaan*.

The paper has been accepted for publication in *Multivariate Behavioral Research*.

Loh, W. W., Ren, D., & West, S. G. (2024). Parametric g-formula for testing time-varying causal effects: What it is, why it matters, and how to implement it in lavaan. Multivariate Behavioral Research. https://doi.org/10.1080/00273171.2024.2354228

A previous version of this manuscript is available as a preprint on PsyArXiv at https://doi.org/10.31234/osf.io/m37uc.

****

**Abstract**

Psychologists leverage longitudinal designs to examine the causal effects of a focal predictor (i.e., treatment or exposure) over time. But causal inference of naturally observed time-varying treatments is complicated by treatment-dependent confounding in which earlier treatments affect confounders of later treatments. In this tutorial article, we introduce psychologists to an established solution to this problem from the causal inference literature: the parametric g-computation formula. We explain why the g-formula is effective at handling treatment-dependent confounding. We demonstrate that the parametric g-formula is conceptually intuitive, easy to implement, and well-suited for psychological research. We first clarify that the parametric g-formula essentially utilizes a series of statistical models to estimate the joint distribution of all post-treatment variables. These statistical models can be readily specified as standard multiple linear regression functions. We leverage this insight to implement the parametric g-formula using lavaan, a widely adopted R package for structural equation modeling. Moreover, we describe how the parametric g-formula may be used to estimate a marginal structural model whose causal parameters parsimoniously encode time-varying treatment effects. We hope this accessible introduction to the parametric g-formula will equip psychologists with an analytic tool to address their causal inquiries using longitudinal data.

**Keywords**
Causal inference; Longitudinal data; Propensity scores; Post-treatment confounding; Time-varying confounding
