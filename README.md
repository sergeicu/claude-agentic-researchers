# Specialized Research Agents for Claude Code 🔬

AI-powered research assistants for **Computational Biology** and **Clinical Trial** researchers using Claude Code.

[![License: Unlicense](https://img.shields.io/badge/license-Unlicense-blue.svg)](http://unlicense.org/)
[![Claude Code](https://img.shields.io/badge/Claude-Code-blueviolet)](https://docs.claude.com/en/docs/claude-code)

## 🎯 Why Specialized Research Agents?

Research in computational biology and clinical trials requires domain-specific expertise:

- **Computational Biology**: Know genomics databases (TCGA, GEO), bioinformatics journals, analysis tools
- **Clinical Trials**: Understand trial design, regulatory requirements, clinical literature

These agents provide **specialized knowledge** for each research domain, from hypothesis formation through data analysis.

## 🚀 Quick Start

```bash
# Clone this repository to your research project
git clone https://github.com/yourusername/claude-agentic-researchers.git
cd claude-agentic-researchers

# Run setup script to verify agents
chmod +x setup-claude-agents.sh
./setup-claude-agents.sh

# Start Claude Code from this directory
claude
```

That's it! The research agents are now available via `@agent-name`.

## 📦 Available Agents

### 🧬 Computational Biology (4 agents)

For genomics, transcriptomics, proteomics, and systems biology research.

| Agent | Color | Focus |
|-------|-------|-------|
| **compbio-hypothesis-generator** | 🟣 Purple | Generate research questions from omics data, form molecular hypotheses |
| **compbio-dataset-specialist** | 🟢 Green | Expert in TCGA, GTEx, GEO, UniProt, single-cell datasets |
| **compbio-literature-specialist** | 🔵 Blue | Search Bioinformatics, Genome Biology, Nature Methods journals |
| **compbio-infrastructure-specialist** | 🔵 Cyan | Setup SRA Toolkit, NCBI, conda environments, handle FASTQ/BAM files |

### 🏥 Clinical Trial Research (4 agents)

For clinical trial design, patient data analysis, and medical research.

| Agent | Color | Focus |
|-------|-------|-------|
| **clinical-hypothesis-generator** | 🟣 Purple | Translate preclinical → clinical trials, design Phase I/II/III studies |
| **clinical-dataset-specialist** | 🟢 Green | Expert in ClinicalTrials.gov, MIMIC, Medicare claims, patient registries |
| **clinical-literature-specialist** | 🔵 Blue | Search NEJM, Lancet, JAMA, find pivotal trials and guidelines |
| **clinical-infrastructure-specialist** | 🔵 Cyan | Setup ClinicalTrials.gov API, REDCap, HIPAA-compliant systems |

## 💡 Usage Examples

### Computational Biology Workflow

**Starting a Cancer Genomics Project:**

```bash
# 1. Generate research hypotheses
@compbio-hypothesis-generator
I have RNA-seq data from 50 melanoma tumors (25 responders, 25 non-responders to immunotherapy). What interesting questions can I explore?

# 2. Review existing literature
@compbio-literature-specialist
Find recent papers on immunotherapy resistance mechanisms in melanoma

# 3. Find relevant datasets
@compbio-dataset-specialist
Find publicly available melanoma RNA-seq datasets I can use for validation

# 4. Setup analysis environment
@compbio-infrastructure-specialist
Setup conda environment for RNA-seq analysis with STAR, DESeq2, and GSEA tools
```

**Example Output from compbio-hypothesis-generator:**

Based on your data, here are testable hypotheses:

1. **Immune Evasion Hypothesis**: Non-responders have upregulated immune checkpoints beyond PD-L1
   - Analysis: Differential expression + immune deconvolution
   - Validation: IHC for novel checkpoints

2. **Tumor Microenvironment Hypothesis**: Responders have higher CD8+ T cell infiltration
   - Analysis: CIBERSORT deconvolution + correlation with response
   - Validation: Multiplex immunofluorescence

3. **Predictive Signature Hypothesis**: A 20-gene classifier can predict response >85% accuracy
   - Analysis: Feature selection + ML (elastic net, random forest)
   - Validation: Independent cohort

### Clinical Trial Workflow

**Designing a Clinical Trial:**

```bash
# 1. Formulate trial hypothesis
@clinical-hypothesis-generator
We showed 80% tumor shrinkage with our KRAS G12C inhibitor in mouse models. Help me design a Phase I/II trial.

# 2. Review similar trials
@clinical-literature-specialist
Find completed Phase I/II trials for other KRAS G12C inhibitors (sotorasib, adagrasib)

# 3. Get historical control data
@clinical-dataset-specialist
Find historical response rates for standard chemotherapy in KRAS G12C+ NSCLC

# 4. Setup trial infrastructure
@clinical-infrastructure-specialist
Setup REDCap for collecting patient data and configure ClinicalTrials.gov API access
```

**Example Output from clinical-hypothesis-generator:**

**Recommended Trial Design:**

**Phase I** (Dose Escalation):
- Design: 3+3 dose escalation (50mg, 100mg, 200mg, 400mg daily)
- Sample size: 12-24 patients
- Primary endpoint: Safety, MTD
- Duration: 12-18 months

**Phase II** (Efficacy):
- Design: Single-arm expansion at RP2D
- Sample size: 90 patients (Simon two-stage)
- Population: KRAS G12C+ metastatic NSCLC, post-platinum therapy, ECOG 0-1
- Primary endpoint: ORR by RECIST 1.1
- Success criterion: ORR ≥40% (vs historical 15% for docetaxel)
- Secondary: PFS, OS, ctDNA clearance, resistance mechanisms

## 🔬 Agent Capabilities

### compbio-hypothesis-generator

**Specializes in:**
- Generating testable hypotheses from RNA-seq, WGS, proteomics, single-cell data
- Understanding biological mechanisms and pathways
- Designing computational analysis strategies
- Suggesting appropriate validation experiments

**Example questions:**
- "What genes drive drug resistance in my cancer cells?"
- "Design a hypothesis for comparing tumor vs normal tissue"
- "What pathway analysis should I do with my differentially expressed genes?"

### compbio-dataset-specialist

**Knows about:**
- **Genomics**: TCGA, ICGC, 1000 Genomes, gnomAD, UK Biobank
- **Expression**: GEO, GTEx, ArrayExpress, Human Protein Atlas
- **Single-cell**: CELLxGENE, Single Cell Portal, 10x datasets
- **Proteomics**: PRIDE, UniProt, BioGRID, STRING
- **Functional**: DepMap, ENCODE, LINCS

**Example questions:**
- "Find scRNA-seq datasets for immune cells in tumors"
- "Where can I get TCGA breast cancer multi-omics data?"
- "Compare GTEx vs TCGA for normal tissue expression"

### compbio-literature-specialist

**Knows journals:**
- Bioinformatics, Genome Biology, Nature Methods, BMC Bioinformatics
- Nature Genetics, Cell Genomics, JCO Precision Oncology
- PLOS Computational Biology, Nucleic Acids Research

**Searches:**
- PubMed with MeSH terms, bioRxiv for preprints
- Methods papers, benchmarking studies, tool comparisons

**Example questions:**
- "Find best practices for scRNA-seq clustering"
- "Papers about CRISPR screens in cancer"
- "Compare tools for variant calling"

### compbio-infrastructure-specialist

**Sets up:**
- Conda/Bioconda environments
- SRA Toolkit for downloading sequencing data
- NCBI, Ensembl, UniProt API access
- Cloud platforms (AWS, Google Cloud) for genomics
- Workflow managers (Nextflow, Snakemake)

**Handles:**
- FASTQ, BAM, VCF file formats
- Large data downloads and processing
- Reference genome setup
- Docker containers for reproducibility

**Example questions:**
- "Download RNA-seq data from GEO accession GSE12345"
- "Setup STAR for aligning RNA-seq reads"
- "Create conda environment for single-cell analysis"

### clinical-hypothesis-generator

**Specializes in:**
- Translating preclinical findings to clinical trials
- Designing Phase I/II/III trials
- Selecting appropriate patient populations and endpoints
- Understanding regulatory requirements (FDA, EMA)

**Knows:**
- Trial phases, designs (RCT, single-arm, basket, umbrella)
- Endpoints (OS, PFS, ORR, DFS, QoL)
- Sample size calculations
- Biomarker strategies

**Example questions:**
- "Design a Phase II trial for novel immunotherapy combination"
- "What endpoints should I use for adjuvant breast cancer trial?"
- "How do I design a biomarker-driven basket trial?"

### clinical-dataset-specialist

**Knows about:**
- **Trial Data**: ClinicalTrials.gov, AACT database, Vivli, YODA
- **EHR Data**: MIMIC, eICU, UK Biobank, All of Us
- **Claims**: Medicare, Optum, MarketScan
- **Registries**: SEER (cancer), OPTN (transplant), rare disease registries
- **Safety**: FAERS, Sentinel

**Example questions:**
- "Find completed melanoma immunotherapy trials with results"
- "Access MIMIC-IV for ICU outcomes research"
- "Get historical control data for metastatic NSCLC survival"

### clinical-literature-specialist

**Knows journals:**
- NEJM, Lancet, JAMA, BMJ (general medicine)
- JCO, Lancet Oncology, JAMA Oncology (oncology)
- Circulation, JACC (cardiology)
- Cochrane Library (systematic reviews)

**Searches:**
- PubMed Clinical Queries
- ClinicalTrials.gov for trial protocols
- FDA approval documents, EMA EPARs
- NCCN, ASCO, ESMO guidelines

**Critically appraises:**
- Trial design, endpoints, statistical analysis
- Bias, conflicts of interest
- Clinical significance vs statistical significance

**Example questions:**
- "Find pivotal trials that led to pembrolizumab approval"
- "What's the standard first-line treatment for HER2+ breast cancer?"
- "Review systematic reviews of CAR-T therapy in lymphoma"

### clinical-infrastructure-specialist

**Sets up:**
- ClinicalTrials.gov API, AACT database queries
- PubMed E-utilities for literature search
- REDCap for clinical trial data collection
- MIMIC/eICU access via PhysioNet
- UK Biobank, Medicare data applications

**Ensures compliance:**
- HIPAA de-identification
- IRB/ethics requirements
- Data Use Agreements (DUAs)
- 21 CFR Part 11 (electronic records)

**Example questions:**
- "Setup REDCap for Phase II trial with adverse event tracking"
- "Query ClinicalTrials.gov API for all recruiting lung cancer trials"
- "Access MIMIC-IV data on Google BigQuery"

## 📁 Repository Structure

```
claude-agentic-researchers/
├── README.md
├── setup-claude-agents.sh
├── .claude/
│   └── agents/
│       ├── compbio-hypothesis-generator.md
│       ├── compbio-dataset-specialist.md
│       ├── compbio-literature-specialist.md
│       ├── compbio-infrastructure-specialist.md
│       ├── clinical-hypothesis-generator.md
│       ├── clinical-dataset-specialist.md
│       ├── clinical-literature-specialist.md
│       └── clinical-infrastructure-specialist.md
├── research-phase-design.md (design docs)
└── future-analysis-phase.md (future plans)
```

## 🛠️ Requirements

- [Claude Code](https://docs.claude.com/en/docs/claude-code) installed
- Bash (for setup script) - works on macOS, Linux, Windows (Git Bash/WSL)

## 🎓 Complete Research Workflows

### Computational Biology: Cancer Drug Resistance Study

```bash
# 1. Hypothesis Generation
@compbio-hypothesis-generator
I have RNA-seq from 200 cancer cell lines with drug sensitivity data. Generate hypotheses about resistance mechanisms.

# 2. Literature Review
@compbio-literature-specialist
Find papers on EGFR inhibitor resistance mechanisms in lung cancer

# 3. Dataset Acquisition
@compbio-dataset-specialist
Find CCLE or DepMap data with both gene expression and drug response

# 4. Infrastructure Setup
@compbio-infrastructure-specialist
Setup Python environment with scanpy, DESeq2, and GSEA tools

# 5. Analysis (your own work or use future analysis agents)
# - Differential expression analysis
# - Pathway enrichment
# - Machine learning for prediction
# - Validation experiments
```

### Clinical Trial: Immunotherapy Combination Trial

```bash
# 1. Trial Design
@clinical-hypothesis-generator
Design a Phase Ib/II trial combining anti-PD-1 + novel TGF-beta inhibitor in melanoma patients who progressed on anti-PD-1 monotherapy

# 2. Literature & Regulatory Review
@clinical-literature-specialist
Find completed combination immunotherapy trials and their endpoints. What did FDA require for approval?

# 3. Historical Controls
@clinical-dataset-specialist
Find response rates for salvage therapy in anti-PD-1 refractory melanoma

# 4. Trial Setup
@clinical-infrastructure-specialist
Setup REDCap with:
- Eligibility screening forms
- Dose escalation tracking (3+3 design)
- Adverse event reporting (CTCAE grading)
- Response assessment (RECIST 1.1)
- Patient-reported outcomes (PRO-CTCAE)

Register trial on ClinicalTrials.gov

# 5. Conduct trial & publish results
```

## 🆚 Why Domain-Specific Agents?

**Generic Research Agents**:
- Broad knowledge but shallow expertise
- Don't know specialized databases (TCGA, ClinicalTrials.gov)
- Miss domain-specific nuances

**Specialized Agents (this repo)**:
- **CompBio agents**: Know difference between RNA-seq vs scRNA-seq, understand GTEx vs TCGA use cases, speak the language of genomics
- **Clinical agents**: Understand Phase I vs III trial designs, know RECIST criteria, familiar with FDA requirements

**Result**: More accurate recommendations, better research design, time saved.

## 📖 Best Practices

### When to Use Which Agent

**Computational Biology**:
1. **Starting new analysis** → compbio-hypothesis-generator
2. **Finding data** → compbio-dataset-specialist
3. **Understanding methods** → compbio-literature-specialist
4. **Setting up tools** → compbio-infrastructure-specialist

**Clinical Trial Research**:
1. **Designing trial** → clinical-hypothesis-generator
2. **Finding similar trials** → clinical-literature-specialist
3. **Getting trial/patient data** → clinical-dataset-specialist
4. **Setting up infrastructure** → clinical-infrastructure-specialist

### Combine Agents for Complex Tasks

```bash
# Multi-agent collaboration
@compbio-dataset-specialist @compbio-literature-specialist
Find melanoma scRNA-seq datasets AND papers analyzing immune cell populations in these datasets

@clinical-hypothesis-generator @clinical-literature-specialist
Design Phase II trial for drug X AND review similar completed trials to inform the design
```

## 🤝 Contributing

Contributions are welcome! Ideas for new agents:

**More Domain-Specific Agents**:
- **neuroscience-researcher**: Neuroimaging, electrophysiology, behavioral data
- **structural-biologist**: PDB, AlphaFold, cryo-EM databases
- **pharmacologist**: DrugBank, ChEMBL, ADME/Tox prediction
- **epidemiologist**: NHANES, population health data, statistical analysis

**Analysis-Phase Agents** (for data you've already collected):
- **data-cleaner**, **eda-specialist**, **statistician**, **ml-engineer**

See [CONTRIBUTING.md](CONTRIBUTING.md) (if you create one) for guidelines.

## 📄 License

This project is released into the public domain under The Unlicense. You can do whatever you want with it - no attribution required!

See the [LICENSE](LICENSE) file or visit [unlicense.org](https://unlicense.org).

## 🙏 Acknowledgments

- Built for use with [Claude Code](https://www.anthropic.com/claude/code)
- Inspired by [ciign/agentic-engineering](https://github.com/ciign/agentic-engineering)
- Designed for the research community: computational biologists and clinical trial researchers

## 📞 Support & Community

- **Issues**: [GitHub Issues](https://github.com/yourusername/claude-agentic-researchers/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/claude-agentic-researchers/discussions)
- **Claude Code Docs**: [docs.claude.com/claude-code](https://docs.claude.com/en/docs/claude-code)

## 🌟 Use Cases

### Computational Biology
- PhD students analyzing omics data
- Postdocs designing experiments
- Bioinformatics core facilities
- Pharma computational biologists

### Clinical Trial Research
- Clinical trial coordinators
- Medical oncologists designing trials
- Regulatory affairs specialists
- Clinical research organizations (CROs)

---

**Accelerate your research with specialized AI agents! 🚀**

*Made with ❤️ for computational biologists and clinical researchers*
