# Influenza Surveillance via Genomics & Wastewater

## Project Status
** – Project initialization & scoping**

This repository documents the early-stage development of a computational research project exploring how **viral genomics** and **wastewater surveillance** could be used to estimate influenza disease burden in settings with limited clinical reporting.

At this stage, the project is focused on **problem formulation, data discovery, and methodological design**. No final conclusions are claimed.

---

## Motivation
Influenza surveillance traditionally relies on **clinical reporting**, which depends on healthcare access, testing availability, and reporting compliance. In many regions (e.g., India), these systems significantly undercount true infections.

Wastewater surveillance offers a population-level signal that does not depend on individual testing behavior. Viral genomics provides additional context about circulating strains and their biological characteristics.

This project asks:

> *If wastewater viral load and viral genomic data were available, could they be used to estimate the true clinical burden of influenza in under-reported settings?*

---

## Initial Research Questions

The guiding questions are:

1. How do influenza viral genomic characteristics (conserved vs variable genes) behave within a single season?
2. How does wastewater influenza viral load relate to reported clinical indicators?
3. Can wastewater + genomics be used to **correct or calibrate** clinical case counts?
4. Can a model trained on U.S. data be conceptually extended to countries with weaker surveillance infrastructure?

---

## Data Identified (Preliminary)

### Viral Genomics
- Source: **NCBI GenBank**
- Genes of interest:
  - **Segment 4: Hemagglutinin (HA)** – antigenically variable
  - **Segment 5: Nucleoprotein (NP)** – conserved structural gene
- Metadata: accession IDs, collection dates

### Wastewater Surveillance
- Source: **CDC National Wastewater Surveillance System (NWSS)**
- Metric: Influenza A viral concentration (PCR-based)
- Aggregation: weekly, state-level (California as pilot)

### Clinical Surveillance
- Source: **CDC FluView / ILINet**
- Metric: weekly influenza-like illness (ILI) percentages
- Aggregation: weekly, state-level
---

## Notes
This README will evolve weekly to reflect the project’s current state. Claims and conclusions will only be added once supported by analysis.
