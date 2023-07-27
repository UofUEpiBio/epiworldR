---
title: 'epiworldR: Fast Agent-Based Epi Models'
authors:
- affiliation: 1
  name: Derek Meyer
  orcid: 0009-0005-1350-6988
- affiliation: 1
  name: George G Vega Yon
  orcid: 0000-0002-3171-0844
date: "13 July 2023"
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
  name: Division of Epidemiology, University of Utah, Salt Lake City, Utah
---

# Introduction

Agent-based modeling (ABM) has emerged as a powerful computational approach to studying complex systems across various fields, including social sciences and epidemiology. By simulating the interactions and behaviors of individual entities, known as agents, ABM provides a unique lens through which researchers can analyze and understand the emergent properties and dynamics of these systems. The `epiworldR` package provides a flexible framework for ABM implementation and methods for prototyping disease outbreaks and transmission models using a C++ backend. It supports multiple epidemiological models, including the Susceptible-Infected-Susceptible (SIS), Susceptible-Infected-Removed (SIR), Susceptible-Exposed-Infected-Removed (SEIR), and others, involving arbitrary mitigation policies and multiple-disease models. Users can specify transmission/susceptibility rates as a function of agents' features, providing great complexity for the model dynamics.

# Key Features

Built on the robust foundation of C++, this package combines efficient computation with the flexibility of R, providing a seamless user experience with several standout features:

-   Multiple Viruses - This feature allows researchers to simulate and analyze the dynamics of various infectious diseases simultaneously.
-   Multiple Tools - This capability enables users to design model objects that may reflect real-world interventions. Examples include vaccines, mask wearing protocols, social distancing, and much more.
-   Multiple Runs - `epiworldR` has the ability to run simulations multiple times, providing a more complete picture of simulation results and its behaviors.
-   Global Actions - Similar to adding multiple tools, global actions allow users to implement interventions and policies at a global scale at any point during a simulation, mimicking real-world scenarios and aiding in the assessment of their impact.
-   Built-In Epidemiological Models - This feature ensures that researchers have access to well-established and widely recognized frameworks, making it easier to compare and benchmark their results. Included are popular epidemiological models such as: SIR, SIS, SEIR, SIR connected, SIS connected, SEIR connected, and more.

With its extensive range of features, `epiworldR` empowers researchers and practitioners in the field of epidemiology to conduct thorough analyses, make informed decisions, and contribute to the advancement of public health.

# Existing Alternatives

Several alternatives to `epiworldR` exist and provide researchers with a range of options, each with its own unique features and strengths, enabling the exploration and analysis of infectious disease dynamics through agent-based modeling. Below is a manually curated table of existing alternatives including ABM [@ABM], abmR [@abmR], cystiSim [@cystiSim], villager [@villager], and RNetLogo [@RNetLogo].

| Package                                                                     | Multiple Viruses | Multiple Tools | Multiple Runs | Global Actions | Built-In Epi Models | Dependencies                                                                                             | Last Commit                                                                                                               |
|:--------|:--------|:--------|:--------|:--------|---------|:--------|:--------|
| [**epiworldR**](https://cran.r-project.org/package=epiworldR)               | yes              | yes            | yes           | yes            | yes                 | 1 | July 2023 |
| [**ABM**](https://cran.r-project.org/package=ABM)                           | \-               | \-             | \-            | yes            | yes                 | 2             | April 2023                 |
| [**abmR**](https://cran.r-project.org/package=abmR)                         | \-               | \-             | yes           | \-             | \-                  | 156          | June 2023                   |
| [**cystiSim**](https://cran.r-project.org/package=cystiSim)                 | \-               | yes            | yes           | \-             | \-                  | 28   | February 2020      |
| [**villager**](https://cran.r-project.org/package=villager)                 | \-               | \-             | \-            | yes            | \-                  | 26  | December 2022          |
| [**RNetLogo**](https://cran.r-project.org/web/packages/RNetLogo/index.html) | \-               | yes            | yes           | yes            | \-                  | 7   | June 2017               |

# Statement of Need

The `epiworldR` package addresses the need for sophisticated epidemiological research by offering a user-friendly agent-based modeling (ABM) tool for simulating and analyzing infectious disease dynamics. It caters to a diverse audience of researchers, epidemiologists, public health practitioners, social scientists, and policymakers, enabling them to study complex disease spread patterns, evaluate intervention strategies, and prepare for emerging pathogens. Unlike traditional epidemiological models, epiworldR's ABM approach captures individual-level behaviors and interactions, providing a more realistic representation of disease dynamics. This distinction sets it apart from compartmental models and other ABM frameworks, making epiworldR a valuable and accessible tool for advancing infectious disease modeling and enhancing outbreak preparedness.

# Conclusion

The development of the `epiworldR` package has ushered in a new era of agent-based modeling in the field of social science and epidemiology. By harnessing the power of C++ and the flexibility of R, this comprehensive package offers a multitude of features that enhance the modeling and analysis of complex infectious disease dynamics. The package's ability to handle multiple viruses and tools, its diverse set of epidemiological models, its capability to run simulations multiple times, and the inclusion of global actions capability empower researchers to explore a wide range of scenarios and make informed decisions regarding public health interventions. `epiworldR` serves as a valuable resource for the social science and epidemiological communities, enabling the study of real-world phenomena, prediction of outcomes, and policy analysis. As the field of epidemiology continues to evolve, `epiworldR` stands at the forefront, providing researchers and practitioners with a powerful tool to navigate the complexities of infectious diseases and contribute to the improvement of global health outcomes.
