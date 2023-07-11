---
title: 'epiworldR: Fast Agent-Based Epi Models'
authors:
- affiliation: 1
  name: Derek Meyer
  orcid: 0009-0005-1350-6988
- affiliation: 1
  name: George G Vega Yon
  orcid: 0000-0002-3171-0844
date: "03 July 2023"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- ABM
- epiworldR
- parallel computing
affiliations:
- index: 1
  name: Division of Epidemiology, University of Utah
---

# Introduction

Agent-based modeling (ABM) has emerged as a powerful computational approach to studying complex systems across various fields, including social sciences and epidemiology. By simulating the interactions and behaviors of individual entities, known as agents, ABM provides a unique lens through which researchers can analyze and understand the emergent properties and dynamics of these systems. The epiworldR package provides a flexible framework for ABM implementation and methods for prototyping disease outbreaks and transmission models using a C++ backend. It supports multiple epidemiological models, including the Susceptible-Infected-Susceptible (SIS), Susceptible-Infected-Removed (SIR), Susceptible-Exposed-Infected-Removed (SEIR), and others, involving arbitrary mitigation policies and multiple-disease models. Users can specify infectiousness/susceptibility rates as a function of agents' features, providing great complexity for the model dynamics.

# Key Features

Built on the robust foundation of C++, this package combines efficient computation with the flexibility of R, providing a seamless user experience. One of its standout features is its ability to handle multiple viruses, allowing researchers to simulate and analyze the dynamics of various infectious diseases simultaneously. Additionally, the package boasts a multi-tool capability, encompassing a diverse set of epidemiological models, enabling users to choose the most appropriate model for their specific research question or scenario. 'epiworldR' also offers global actions capability, allowing users to implement interventions and policies at a global scale, mimicking real-world scenarios and aiding in the assessment of their impact. Furthermore, the inclusion of popular epidemiological models within the package ensures that researchers have access to well-established and widely recognized frameworks, making it easier to compare and benchmark their results. With its extensive range of features, 'epiworldR' empowers researchers and practitioners in the field of epidemiology to conduct thorough analyses, make informed decisions, and contribute to the advancement of public health.

# Existing Alternatives

Several existing alternatives to 'epiworldR' exist. These alternatives provide researchers with a range of options, each with its own unique features and strengths, enabling the exploration and analysis of infectious disease dynamics through agent-based modeling. The below table is a manually curated list of existing alternatives:

| Package                                                       | Multiple-Viruses | Multiple-Tools | Multiple-Runs | Global Actions | Dependencies                                                                                             | Activity                                                                                                               |
|:--------------------------------------------------------------|:-----------------|:---------------|:--------------|:---------------|:---------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------|
| [**epiworldR**](https://cran.r-project.org/package=epiworldR) | yes              | yes            | yes           | yes            | [![status](https://tinyverse.netlify.com/badge/epiworldR)](https://CRAN.R-project.org/package=epiworldR) | [![Activity](https://img.shields.io/github/last-commit/UofUEpiBio/epiworldR)](https://github.com/UofUEpiBio/epiworldR) |
| [**ABM**](https://cran.r-project.org/package=ABM)             | \-               | \-             | \-            | \-             | [![status](https://tinyverse.netlify.com/badge/ABM)](https://CRAN.R-project.org/package=ABM)             | [![Activity](https://img.shields.io/github/last-commit/junlingm/ABM)](https://github.com/junlingm/ABM)                 |
| [**abmR**](https://cran.r-project.org/package=abmR)           | \-               | \-             | \-            | \-             | [![status](https://tinyverse.netlify.com/badge/abmR)](https://CRAN.R-project.org/package=abmR)           | [![Activity](https://img.shields.io/github/last-commit/bgoch5/abmR)](https://github.com/bgoch5/abmR)                   |
| [**cystiSim**](https://cran.r-project.org/package=cystiSim)   | \-               | \-             | \-            | \-             | [![status](https://tinyverse.netlify.com/badge/cystiSim)](https://CRAN.R-project.org/package=cystiSim)   | [![Activity](https://img.shields.io/github/last-commit/brechtdv/cystiSim)](https://github.com/brechtdv/cystiSim)       |
| [**villager**](https://cran.r-project.org/package=villager)   | \-               | \-             | \-            | \-             | [![status](https://tinyverse.netlify.com/badge/villager)](https://CRAN.R-project.org/package=villager)   | [![Activity](https://img.shields.io/github/last-commit/zizroc/villager)](https://github.com/zizroc/villager)           |

# Conclusion

The development of the 'epiworldR' package has ushered in a new era of agent-based modeling in the field of social science and epidemiology. By harnessing the power of C++ and the flexibility of R, this comprehensive package offers a multitude of features that enhance the modeling and analysis of complex infectious disease dynamics. The package's ability to handle multiple viruses, its diverse set of epidemiological models, and the inclusion of global actions capability empower researchers to explore a wide range of scenarios and make informed decisions regarding public health interventions. The 'epiworldR' package serves as a valuable resource for the social sciences and epidemiological communities, enabling the study of real-world phenomena, prediction of outcomes, and policy analysis. As the field of epidemiology continues to evolve, 'epiworldR' stands at the forefront, providing researchers and practitioners with a powerful tool to navigate the complexities of infectious diseases and contribute to the improvement of global health outcomes.
