# Project Description

## Background and Motivation

Influenza surveillance today relies mostly on clinical reporting — people visiting hospitals, getting tested, and those results being officially recorded. This works reasonably well in countries like the United States, but even there, many infections go unreported due to mild symptoms, lack of testing, or delays in healthcare access. In countries like India, this gap is even larger because clinical surveillance infrastructure is more limited and inconsistent.

Wastewater surveillance has emerged as a way to measure virus circulation at the community level, independent of individual testing behavior. At the same time, viral genome sequencing provides insight into how conserved or variable different viral genes are during a season. This project explores how these two approaches can be combined to better estimate the *true* disease burden of influenza.

---

## Core Idea of the Project

**Can wastewater surveillance and viral genomics be used to estimate how much influenza is being under-reported clinically?**

The central hypothesis is that wastewater viral concentration reflects the actual number of infected individuals more directly than reported clinical case counts. By learning the relationship between wastewater viral load and clinical illness indicators in a well-instrumented setting (California), we can build a model that estimates true infection burden. This framework can then be extended to countries like India, where wastewater data may eventually become available but clinical reporting remains incomplete.

---

## Role of Genomics

This project focuses on two Influenza A gene segments:

- **Hemagglutinin (HA)**: a surface protein under strong immune pressure and known to be highly variable.
- **Nucleoprotein (NP)**: a conserved internal protein essential for viral replication.

By analyzing within-season genomic diversity of these genes, the project aims to:
- Understand how much viral evolution occurs during a single influenza season.
- Interpret whether observed genetic variability meaningfully affects surveillance signals.
- Provide biological context for why some surveillance targets are more stable than others.

Genomics is not used here to predict individual variants or mutations, but rather to support interpretation of surveillance reliability and biological constraints.

---

## Data Sources and Approach

The project uses publicly available data, including:
- Influenza A wastewater surveillance data from the CDC
- Clinical influenza-like illness (ILI) data from FluView
- Influenza A HA and NP sequences from NCBI, along with collection dates

All analysis is computational and statistical. No laboratory or wet-lab work is involved.

---

## Intended Outcome

By the end of the project, the goal is to:
- Demonstrate how wastewater viral load correlates with clinical influenza indicators.
- Quantify the degree to which clinical surveillance may under-estimate true infection levels.
- Propose a modeling framework that could be applied in settings like India once wastewater surveillance data becomes available.

The project does **not** claim to replace clinical surveillance, predict individual infections, or detect outbreaks before they happen. Instead, it aims to complement existing systems by improving population-level estimates.

---

## Scope and Limitations

This project is intentionally scoped to:
- One influenza season
- Aggregated weekly data
- Population-level trends rather than individual outcomes

The analysis is meant to be exploratory and explanatory, not definitive. All conclusions are framed with biological and data limitations in mind.

---

## Why This Matters

Accurate disease burden estimates are essential for public health planning, resource allocation, and policy decisions. In regions where clinical reporting is incomplete, relying solely on reported case counts can lead to underestimation of disease impact. This project explores how combining wastewater surveillance with viral genomics can help close that gap in a realistic and data-driven way.
