# Prospectus: Benchmarking DIA Search Softwares for Proteomics Analysis

## Background
Data-Independent Acquisition (DIA) has become a cornerstone of modern proteomics, enabling comprehensive protein quantification. However, the field is evolving rapidly, with vendors shifting towards licensed, pay-to-use models. This creates challenges for accessibility, particularly in academic and resource-limited settings. To address this, we propose a collaborative benchmarking study to evaluate free-to-use and licensed DIA search softwares, identifying alternatives, strengths, weaknesses, and opportunities for future development.

## Objectives
- Evaluate and compare the performance of free-to-use and licensed DIA search softwares.
- Identify weaknesses to inform software development and capability building.
- Foster collaboration among Australian bioinformaticians to tackle DIA analysis challenges.
- Develop expertise in using diverse DIA tools and promote open-source alternatives.

## Software to Benchmark
### Free-to-Use
- DIA-NN (open-source, widely used for DIA data processing).
- OpenSwath (part of OpenMS suite, free and extensible).
- MSFragger-DIA (academic tool for DIA searches).

### Licensed
- Spectronaut (Pulsar engine, vendor-specific).
- Proteome Discoverer (with CHIMERYS for DIA).
- DIA-Umpire (integrated in some commercial pipelines).

(Note: Final list to be refined based on collaborator input and availability.)

## Datasets
- Recent DIA-focused datasets from Australian labs/groups/facilities (e.g., MTI 2-specific datasets, excluding sepsis data).
- Preference for publicly available or shareable data to ensure reproducibility.
- Mix of standards (e.g., controlled spike-ins) and real biological samples for comprehensive evaluation.

## Rationale
- Vendors are increasingly adopting pay-to-use models, limiting access to advanced DIA tools.
- Benchmarking will uncover viable free alternatives and highlight areas for improvement.
- This study will build national capability in DIA analysis and unite bioinformaticians on a shared challenge.
- Relative benchmarking across identical samples will provide actionable insights into software performance.

## Resources Available
- ABLeS Software Accelerator allocation at NCI for open-source software testing.
- On-premises resources for licensed software comparisons (to handle proprietary constraints).
- Collaborative access to computational infrastructure and expertise from participating groups.

## Benchmarking Approach
- **Relative Benchmarking:** Analyze the same samples across all softwares to compare metrics like protein identification rates, quantification accuracy, false discovery rates, and computational efficiency.
- **Standards vs. Real Samples:** Include both controlled standards for precision testing and complex biological samples for real-world applicability.
- **Metrics:** Sensitivity, specificity, runtime, ease of use, and integration with downstream tools.
- **Workflow:** Standardized pipelines for data processing, with version control and reproducible scripts.

## Anticipated Outcomes
- Comprehensive report on software performance, including strengths/weaknesses and recommendations.
- Installation guides for softwares at NCI (where applicable).
- Open-source repository of benchmarking scripts and results.
- Foundation for future developments, such as hybrid tools or improvements to free softwares.
- Strengthened network of Australian bioinformaticians focused on proteomics challenges.

## Call for Collaboration
We invite bioinformaticians, proteomics experts, and facilities across Australia to join this initiative. Contributions could include dataset provision, software expertise, computational resources, or analysis support. 

**Contact:** [Your Name/Email] for expressions of interest or to discuss participation.

*This prospectus is a draft and open to feedback. Estimated timeline: 6-12 months.*