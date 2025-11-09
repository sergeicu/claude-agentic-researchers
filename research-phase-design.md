# Research Phase Agents - Design Document

## Overview

This document outlines the design for research-focused agents that support the scientific discovery workflow. Unlike traditional data science agents that focus on modeling and analysis, these agents focus on the critical but often-overlooked **research discovery phase**: forming hypotheses, reviewing literature, and setting up research infrastructure.

## Research Workflow Steps

1. **Hypothesis Formation** - What questions should we ask?
2. **Infrastructure Setup** - How do we access research resources (papers, datasets)?
3. **Literature Search** - What's been done before?
4. **Paper Analysis** - What do these papers say?
5. **Synthesis** - What's the state of knowledge? What gaps exist?
6. **Hypothesis Refinement** - What should we test based on prior work?
7. **Data Acquisition** - Download relevant datasets

## Agent Design Options

### Option A: Minimal (3 Agents)

**1. hypothesis-generator**
- Forms initial and refined hypotheses
- Identifies interesting research questions
- Suggests testable hypotheses

**2. research-infrastructure-specialist**
- Sets up MCP tools for paper/data access
- Configures research infrastructure
- Helps download datasets and papers

**3. literature-analyst**
- Searches for relevant papers
- Reads and analyzes papers
- Synthesizes findings across papers
- Identifies research gaps

**Pros**: Simple, clear boundaries
**Cons**: literature-analyst is very broad (search + read + synthesize)

---

### Option B: Expanded (6 Agents)

**1. hypothesis-generator**
- Forms and refines hypotheses
- Identifies research questions

**2. research-infrastructure-specialist**
- Sets up MCP tools
- Configures access to research platforms

**3. dataset-finder**
- Specialized in finding research datasets
- Knows major repositories (Kaggle, UCI, government data, domain-specific)
- Downloads and validates datasets

**4. literature-searcher**
- Searches for relevant papers
- Knows how to query PubMed, arXiv, Google Scholar
- Filters and ranks search results

**5. paper-analyst**
- Deep reading of individual papers
- Extracts methodologies, findings, limitations
- Understands scientific paper structure

**6. research-synthesizer**
- Combines findings across multiple papers
- Identifies consensus and contradictions
- Maps the research landscape
- Identifies gaps

**Pros**: Highly specialized, each agent has focused role
**Cons**: May be too granular, more agents to manage

---

### Option C: Balanced (4-5 Agents) - RECOMMENDED

**1. hypothesis-generator** 🧠
- **Focus**: Research question formation and hypothesis development
- **Responsibilities**:
  - Generate interesting research questions from data
  - Form testable hypotheses
  - Refine hypotheses based on literature review
  - Identify assumptions and predictions
  - Suggest experimental approaches

**2. research-infrastructure-specialist** 🔧
- **Focus**: Setting up research tools and data access
- **Responsibilities**:
  - Recommend which platforms to use (PubMed, arXiv, bioRxiv, etc.)
  - Setup and configure MCP tools for research
  - Help access paywalled content (via institutional access)
  - Configure API keys and authentication
  - Troubleshoot access issues

**3. literature-analyst** 📚
- **Focus**: Finding, reading, and synthesizing research papers
- **Responsibilities**:
  - Search for relevant papers across platforms
  - Read and extract key information from papers
  - Summarize individual papers
  - Compare and contrast methodologies
  - Identify research gaps and opportunities
  - Track citations and influence

**4. dataset-specialist** 💾
- **Focus**: Finding, accessing, and validating research datasets
- **Responsibilities**:
  - Know major dataset repositories
  - Search for relevant datasets
  - Evaluate dataset quality and relevance
  - Download and validate datasets
  - Understand data formats and standards
  - Handle large dataset downloads

**5. [Optional] research-validator** ✅
- **Focus**: Assessing research quality and validity
- **Responsibilities**:
  - Evaluate study design quality
  - Identify potential biases
  - Check statistical validity
  - Assess reproducibility
  - Flag questionable research practices

**Pros**: Good balance of specialization and manageability
**Cons**: Still need to split literature-analyst into search vs analysis?

---

## Recommended Design: Option C (4 Core Agents)

Let's go with **4 core agents** that cover the essential research workflow:

### 1. hypothesis-generator 🧠

**Color**: Purple (creative/strategic thinking)
**Model**: Sonnet (needs reasoning capability)

**Primary Focus**: Generate and refine research questions and hypotheses

**Core Responsibilities**:
- Analyze datasets and suggest interesting research questions
- Form clear, testable hypotheses
- Refine hypotheses based on literature review findings
- Identify key assumptions and predictions
- Suggest appropriate research methodologies
- Frame questions that address research gaps

**Technical Expertise**:
- Scientific method and hypothesis formation
- Understanding of research design
- Domain knowledge across scientific fields
- Statistical hypothesis framing
- Experimental design principles

**When to Use**:
- Starting a new research project
- Have data but don't know what to investigate
- After literature review, need to refine hypothesis
- Want to generate multiple testable hypotheses
- Need to frame research questions clearly

**Example Use Cases**:
- "I have a PubMed dataset - what interesting questions could I explore?"
- "Based on the literature review, help me refine my hypothesis"
- "Generate testable hypotheses about collaboration patterns in COVID research"
- "What assumptions am I making with this hypothesis?"

**Example Workflow**:
```
User: "I have PubMed data on 50K papers from 2019-2024"

@hypothesis-generator analyzes the data structure and suggests:
1. "How did research collaboration patterns change during COVID-19?"
2. "Which countries increased their research output the most?"
3. "Did international co-authorship increase in pandemic-related research?"
4. "Can we predict highly-cited papers from early collaboration metrics?"

User chooses #3, agent refines:
"Hypothesis: International co-authorship rates in COVID-19 research
increased by >50% compared to pre-pandemic baselines, with sustained
elevation through 2024."

Testable predictions:
- Proportion of papers with authors from 2+ countries increased
- This increase is specific to COVID research (control: other diseases)
- Effect persisted beyond initial pandemic period (2020-2021)
```

---

### 2. research-infrastructure-specialist 🔧

**Color**: Cyan (infrastructure/devops analog)
**Model**: Sonnet

**Primary Focus**: Set up tools and infrastructure to access research resources

**Core Responsibilities**:
- Recommend which platforms to search (PubMed, arXiv, bioRxiv, Google Scholar, etc.)
- Setup MCP tools for accessing papers and datasets
- Configure authentication and API keys
- Help navigate institutional access
- Troubleshoot access issues
- Know about research platform APIs and access methods

**Technical Expertise**:
- MCP (Model Context Protocol) server setup
- Research platform APIs (PubMed, arXiv, Semantic Scholar, etc.)
- Authentication methods (API keys, OAuth, institutional access)
- Dataset repositories (Kaggle, UCI ML, Hugging Face datasets, domain-specific repos)
- File formats for research data (CSV, Parquet, HDF5, FITS, etc.)
- Bulk download strategies

**When to Use**:
- Need to access research papers systematically
- Want to setup automated literature monitoring
- Need to find and download large datasets
- Encountering access issues to papers or data
- Setting up reproducible research infrastructure

**Example Use Cases**:
- "Setup tools to search PubMed and download abstracts"
- "I need to access papers from arXiv - what MCP tools are available?"
- "Configure access to bioRxiv for preprint analysis"
- "Find and download the PubMed baseline files"
- "Setup automated alerts for new papers in my research area"

**Example Workflow**:
```
User: "I need to analyze COVID-19 research papers from multiple sources"

@research-infrastructure-specialist recommends:

1. For PubMed:
   - Install 'pubmed-mcp-server' for E-utilities API access
   - Provides NCBI API key setup instructions
   - Explains rate limits and bulk download strategies

2. For arXiv:
   - Install 'arxiv-mcp-server' for arXiv API
   - Shows how to search by category and date

3. For Semantic Scholar:
   - Install 'semantic-scholar-mcp' for citation networks
   - Explains how to get API key

4. Helps configure .claude/config.json to add these MCP servers

5. Provides example queries for each platform
```

**Key Capabilities**:
- **MCP Discovery**: Know which MCP tools exist for research
- **Configuration**: Help edit config files and troubleshoot
- **Best Practices**: Rate limiting, caching, bulk downloads
- **Format Knowledge**: Understand BibTeX, EndNote, RIS, etc.
- **Data Sources**: Know specialized repositories (GenBank, Protein Data Bank, clinical trial registries)

---

### 3. literature-analyst 📚

**Color**: Blue (analytical/knowledge)
**Model**: Sonnet

**Primary Focus**: Search, read, analyze, and synthesize research literature

**Core Responsibilities**:
- Craft effective search queries for research platforms
- Search for relevant papers using MCP tools or web search
- Read and understand research papers (PDFs)
- Extract key information (methods, findings, limitations)
- Summarize individual papers
- Compare and contrast multiple papers
- Identify consensus and contradictions in literature
- Map the research landscape
- Identify research gaps and opportunities
- Track citations and influence

**Technical Expertise**:
- Scientific paper structure (IMRAD format)
- Search query optimization (Boolean operators, MeSH terms, field tags)
- Critical reading and evaluation
- Research synthesis methodologies
- Citation analysis
- Understanding research methodologies across domains
- Identifying study limitations and biases

**When to Use**:
- Need to review what's been done in an area
- Looking for papers on a specific topic
- Want summary of a research paper
- Comparing different approaches in literature
- Identifying what's known vs unknown
- Understanding the state of a research field

**Example Use Cases**:
- "Find papers about collaboration networks in medical research"
- "Read this PDF and summarize the key findings"
- "Compare the methodologies used in these 5 papers"
- "What's the current state of knowledge on this topic?"
- "Synthesize findings from 20 papers on COVID research collaboration"
- "Has anyone already tested this hypothesis?"

**Example Workflow**:
```
User: "Find papers about international collaboration in pandemic research"

@literature-analyst:

1. Crafts search query:
   PubMed: "(COVID-19 OR pandemic) AND (international collaboration OR
   cross-country OR multi-country) AND (research OR publication)"

2. Searches using MCP tools, finds 156 relevant papers

3. Ranks by relevance and citation count, identifies top 20

4. Reads top 5 most relevant papers

5. Summarizes each:
   - Paper 1 (Smith et al., 2021): Found 300% increase in international
     co-authorship on COVID papers. Used bibliometric analysis of 45K papers.
   - Paper 2 (Chen et al., 2022): Network analysis showed US-China
     collaboration decreased but Europe-Asia increased...
   [continues]

6. Synthesizes across papers:
   - Consensus: International collaboration definitely increased
   - Contradiction: Magnitude varies (200%-500% depending on definition)
   - Gap: Most studies end in 2021, lack recent data
   - Methods: Primarily bibliometric, few use network analysis
   - Limitation: Most focus on publications, not preprints

7. Recommendations:
   - Your hypothesis is supported by existing literature
   - Gap exists: extending analysis through 2024
   - Consider adding network analysis perspective
   - Include preprints for more complete picture
```

**Advanced Capabilities**:
- **PDF Reading**: Extract text, figures, tables from PDFs
- **Citation Network**: Understand "cited by" and "references"
- **Meta-Analysis**: Synthesize quantitative findings across studies
- **Quality Assessment**: Identify methodological strengths/weaknesses
- **Trend Identification**: See how research evolved over time

---

### 4. dataset-specialist 💾

**Color**: Green (data/resources)
**Model**: Sonnet

**Primary Focus**: Finding, evaluating, and acquiring research datasets

**Core Responsibilities**:
- Know major dataset repositories and how to search them
- Find datasets relevant to research questions
- Evaluate dataset quality, size, and relevance
- Understand data provenance and licensing
- Download and validate datasets
- Handle large dataset downloads
- Recommend appropriate datasets for research questions

**Technical Expertise**:
- Dataset repositories (Kaggle, UCI ML, Hugging Face, Data.gov, domain-specific)
- Data formats (CSV, Parquet, HDF5, JSON, SQL, domain-specific)
- Data validation and quality checks
- Understanding of dataset documentation
- Data licensing (CC0, CC-BY, research-only, etc.)
- APIs for dataset access
- Large file download strategies

**Dataset Sources Knowledge**:

**General**:
- Kaggle datasets
- UCI Machine Learning Repository
- Hugging Face datasets
- Data.gov
- Google Dataset Search
- OpenML

**Academic/Research**:
- Harvard Dataverse
- Zenodo
- figshare
- Dryad
- OSF (Open Science Framework)

**Domain-Specific**:
- Biomedical: PubMed, ClinicalTrials.gov, TCGA, UK Biobank
- Genomics: GenBank, ENCODE, 1000 Genomes
- Astronomy: SDSS, NASA archives
- Climate: NOAA, NASA Earth Data
- Social Science: ICPSR, Pew Research

**When to Use**:
- Need a dataset for research question
- Looking for benchmark datasets
- Want to validate findings on multiple datasets
- Need specific types of data (medical, financial, etc.)
- Evaluating dataset quality
- Don't know where to find data

**Example Use Cases**:
- "Find datasets about academic publication networks"
- "I need medical imaging data for cancer research"
- "What's the best dataset for testing this hypothesis?"
- "Download the latest PubMed baseline files"
- "Find benchmark datasets for citation prediction"
- "Is this dataset appropriate for my research question?"

**Example Workflow**:
```
User: "I need a dataset of academic papers with citation information"

@dataset-specialist:

1. Identifies relevant datasets:

   Option 1: Microsoft Academic Graph (MAG)
   - Size: 260M+ papers
   - Coverage: Broad, all fields
   - Pros: Comprehensive, free
   - Cons: Discontinued (last update 2021)
   - License: Open Data Commons

   Option 2: OpenAlex
   - Size: 250M+ papers
   - Coverage: All fields, actively maintained
   - Pros: Free, good API, citation network
   - Cons: Newer, less validated
   - License: CC0

   Option 3: Semantic Scholar Datasets
   - Size: 200M+ papers
   - Coverage: CS, biomed, strong
   - Pros: Excellent metadata, embeddings available
   - Cons: Field-biased
   - License: Open Data Commons

   Option 4: PubMed (if biomedical only)
   - Size: 35M+ papers
   - Coverage: Biomedical only
   - Pros: Highest quality metadata, MeSH terms
   - Cons: Limited to biomed
   - License: Public domain

2. Recommendation:
   "For broad academic research: Use OpenAlex (actively maintained)
    For biomedical focus: Use PubMed + PMC
    For CS/AI focus: Semantic Scholar"

3. Provides download instructions:
   - OpenAlex: Snapshot download or API access
   - Shows how to download in chunks
   - Estimates: ~500GB compressed, 2TB uncompressed
   - Suggests filtering to reduce size

4. Validates after download:
   - Checks file integrity
   - Verifies record counts
   - Confirms schema matches documentation
```

**Key Capabilities**:
- **Discovery**: Know where to look for specific types of data
- **Evaluation**: Assess quality, coverage, recency
- **Comparison**: Help choose between similar datasets
- **Access**: Navigate different download methods (bulk, API, streaming)
- **Validation**: Ensure downloaded data is correct and complete

---

## How Research Agents Work Together

**Complete Research Workflow**:

```
Phase 1: Initial Exploration
@hypothesis-generator
"I have PubMed data - what could I investigate?"
→ Suggests 5 interesting research questions

User selects: "International collaboration patterns"

Phase 2: Infrastructure Setup
@research-infrastructure-specialist
"Setup tools to search literature and access datasets"
→ Installs PubMed MCP, Semantic Scholar MCP, configures API keys

Phase 3: Literature Review
@literature-analyst
"Find papers about international research collaboration"
→ Searches, finds 50 relevant papers, reads top 10, synthesizes findings

Phase 4: Hypothesis Refinement
@hypothesis-generator
"Based on the literature, refine my hypothesis"
→ "International collaboration in COVID research increased 3x and
   remained elevated through 2024, driven by urgent need for
   distributed expertise and data sharing"

Phase 5: Data Acquisition
@dataset-specialist
"Find datasets with publication and collaboration data 2019-2024"
→ Recommends OpenAlex, provides download instructions

@research-infrastructure-specialist
"Help download OpenAlex COVID subset"
→ Sets up download script, filters for relevant papers

Phase 6: Ready for Analysis
→ Hand off to analysis-phase agents (data-cleaner, eda-specialist, etc.)
```

---

## Agent Interaction Patterns

### Sequential (Common)
```
@hypothesis-generator → @literature-analyst → @hypothesis-generator → @dataset-specialist
```
Form hypothesis → Review literature → Refine → Get data

### Parallel
```
@literature-analyst + @dataset-specialist (both at same time)
```
Search papers AND find datasets simultaneously

### Iterative
```
@literature-analyst ⇄ @hypothesis-generator
```
Read papers → refine question → search more → refine more

### Infrastructure Support
```
@research-infrastructure-specialist (enables all others)
```
Setup once, then other agents use the infrastructure

---

## Design Principles

1. **Research-First**: Focus on discovery, not just analysis
2. **Infrastructure-Aware**: Recognize need for tool setup
3. **Literature-Centric**: Papers and prior work are first-class citizens
4. **Hypothesis-Driven**: Scientific method at the core
5. **Reproducible**: Setup can be documented and shared
6. **Platform-Agnostic**: Work across different research platforms
7. **Domain-Flexible**: Work for any scientific field

---

## Future Enhancements

**Potential 5th Agent**: research-validator
- Evaluate research quality
- Identify biases and limitations
- Assess reproducibility
- Check statistical validity

**Potential 6th Agent**: citation-analyst
- Analyze citation networks
- Identify influential papers
- Track research impact
- Map research communities

**Integration with Analysis Phase**:
- Research agents output → Clean hypotheses and curated datasets
- Analysis agents input → Ready for EDA, statistics, modeling

---

## Comparison to Web Dev Agents (ciign)

| Web Dev (ciign) | Research (sergeicu) |
|----------------|---------------------|
| full-stack-developer | hypothesis-generator |
| backend-specialist | literature-analyst |
| devops-engineer | research-infrastructure-specialist |
| database-designer | dataset-specialist |

**Key Difference**: Research agents focus on knowledge discovery and literature, not code development

---

## Conclusion

**Recommended: 4 Core Research Agents**

1. **hypothesis-generator** - What to investigate
2. **research-infrastructure-specialist** - How to access resources
3. **literature-analyst** - What's been done
4. **dataset-specialist** - Getting the data

These four agents cover the essential research discovery workflow and provide unique value not found in existing agent collections. They work together to transform a vague research interest into a well-defined, literature-informed hypothesis with curated data ready for analysis.
