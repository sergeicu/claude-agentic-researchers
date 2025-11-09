# Research Discovery Agents 🔬

Specialized AI agents for scientific research discovery and literature analysis with Claude Code. These agents help you form hypotheses, review literature, and set up research infrastructure before diving into data analysis.

**Transform your research workflow with AI-powered research specialists.**

[![License: Unlicense](https://img.shields.io/badge/license-Unlicense-blue.svg)](http://unlicense.org/)
[![Claude Code](https://img.shields.io/badge/Claude-Code-blueviolet)](https://docs.claude.com/en/docs/claude-code)

## 🎯 Why Research Agents?

Most data science tools focus on modeling and analysis, but **research starts earlier**:
- What questions should I ask?
- What's already been discovered?
- Where do I find relevant data and papers?
- How do I access research resources?

These agents fill that gap, covering the **research discovery phase** that comes before traditional data science.

## 🚀 Quick Start

```bash
# Navigate to your project directory
cd /path/to/your/research/project

# Clone or copy this repository
git clone https://github.com/yourusername/research-agents.git
cd research-agents

# Run the setup script
chmod +x setup-claude-agents.sh
./setup-claude-agents.sh

# Start Claude Code from your research directory
cd /path/to/your/research/project
claude
```

That's it! The research agents are now ready to use.

## 📦 What's Included

This repository contains **4 specialized research agents** covering the scientific discovery workflow:

### 🧠 Research Discovery
- **hypothesis-generator** - Form research questions and testable hypotheses
- **research-infrastructure-specialist** - Setup MCP tools for accessing papers and datasets
- **literature-analyst** - Search, read, analyze, and synthesize research literature
- **dataset-specialist** - Find, evaluate, and acquire research datasets

## 🔬 The Research Agents

### hypothesis-generator

**Generate and refine research questions and hypotheses**

This agent helps you:
- Identify interesting research questions from data
- Form clear, testable hypotheses
- Refine hypotheses based on literature review
- Suggest appropriate methodologies
- Frame questions that address research gaps

**Example scenarios**:
- Starting a new research project
- Have data but unsure what to investigate
- Need to refine vague research ideas
- Want to generate multiple hypothesis options

### research-infrastructure-specialist

**Setup tools and infrastructure to access research resources**

This agent helps you:
- Recommend which platforms to use (PubMed, arXiv, Semantic Scholar, etc.)
- Setup MCP (Model Context Protocol) tools for research
- Configure authentication and API keys
- Access datasets from various repositories
- Troubleshoot access and download issues

**Example scenarios**:
- Need to systematically search research papers
- Want to setup automated literature monitoring
- Looking for specific datasets
- Encountering access issues to papers or data

### literature-analyst

**Search, read, analyze, and synthesize research literature**

This agent helps you:
- Craft effective search queries for research platforms
- Find relevant papers using MCP tools
- Read and summarize research papers (including PDFs)
- Compare and contrast multiple papers
- Identify consensus, contradictions, and research gaps
- Map the research landscape

**Example scenarios**:
- Need to understand what's been done in your area
- Want comprehensive literature review
- Comparing different methodologies
- Identifying research opportunities
- Checking if your hypothesis has been tested

### dataset-specialist

**Find, evaluate, and acquire research datasets**

This agent helps you:
- Know major dataset repositories (Kaggle, UCI, Hugging Face, domain-specific)
- Find datasets relevant to research questions
- Evaluate dataset quality and appropriateness
- Navigate download methods (bulk, API, streaming)
- Validate downloaded datasets
- Understand data licensing

**Example scenarios**:
- Need data for your research question
- Looking for benchmark datasets
- Evaluating dataset quality
- Don't know where to find specialized data
- Need to download large research datasets

## 💡 Usage Examples

### Single Agent

```bash
# Generate research hypotheses
@hypothesis-generator I have a PubMed dataset of 50K papers from 2019-2024, what interesting questions could I explore?

# Setup research infrastructure
@research-infrastructure-specialist Setup tools to search PubMed and download paper abstracts

# Review literature
@literature-analyst Find papers about international collaboration in pandemic research

# Find datasets
@dataset-specialist Find datasets about academic publication networks with citation data
```

### Multiple Agents

```bash
# Research infrastructure + literature review
@research-infrastructure-specialist @literature-analyst Setup PubMed access and find papers on COVID research collaboration

# Hypothesis + literature review
@hypothesis-generator @literature-analyst I want to study collaboration patterns - help me form and validate hypotheses

# Infrastructure + dataset acquisition
@research-infrastructure-specialist @dataset-specialist Find and download datasets for studying research impact
```

### Complete Research Workflow

```bash
# 1. Form initial hypothesis
@hypothesis-generator I have a dataset of 50K biomedical papers - what could I investigate about collaboration patterns?

# 2. Setup research infrastructure
@research-infrastructure-specialist Setup tools to access PubMed, arXiv, and Semantic Scholar

# 3. Review existing literature
@literature-analyst What's the current state of research on international scientific collaboration?

# 4. Refine hypothesis based on findings
@hypothesis-generator Based on the literature review, help me refine my hypothesis about COVID-era collaboration

# 5. Find appropriate datasets
@dataset-specialist Find datasets with publication metadata and collaboration data from 2019-2024

# 6. Download and prepare data
@research-infrastructure-specialist Help me download the OpenAlex dataset subset for COVID research

# → Now ready for analysis phase (data cleaning, EDA, modeling, etc.)
```

## 📁 Directory Structure

After running the setup script, you'll have:

```
.claude/
└── agents/
    ├── hypothesis-generator.md
    ├── research-infrastructure-specialist.md
    ├── literature-analyst.md
    └── dataset-specialist.md
```

## 🎨 Agent Visual Guide

Each agent has a unique color in the Claude Code UI:

| Agent | Color | Primary Focus |
|-------|-------|---------------|
| 🟣 hypothesis-generator | Purple | Research question formation |
| 🔵 research-infrastructure-specialist | Cyan | Tool setup & data access |
| 🔵 literature-analyst | Blue | Literature review & synthesis |
| 🟢 dataset-specialist | Green | Dataset discovery & acquisition |

## 🎓 Agent Capabilities

### hypothesis-generator
- Generate interesting research questions from data
- Form testable hypotheses with clear predictions
- Refine hypotheses based on literature findings
- Identify assumptions and research design considerations
- Suggest appropriate methodologies
- Frame questions that address identified gaps

### research-infrastructure-specialist
- Know available MCP tools for research (PubMed, arXiv, Semantic Scholar, etc.)
- Setup and configure MCP servers
- Configure API keys and authentication
- Navigate institutional access
- Understand research platform APIs
- Handle bulk downloads and rate limiting
- Troubleshoot access issues

### literature-analyst
- Craft effective search queries (Boolean operators, MeSH terms, field tags)
- Search across multiple research platforms
- Read and understand scientific papers (including PDFs)
- Extract key information (methods, findings, limitations)
- Summarize individual papers clearly
- Compare and contrast multiple papers
- Synthesize findings across literature
- Identify consensus and contradictions
- Map research landscape and trends
- Identify research gaps and opportunities

### dataset-specialist
- Know major dataset repositories (general and domain-specific)
- Search for relevant datasets
- Evaluate dataset quality, coverage, and recency
- Understand data provenance and licensing
- Compare alternative datasets
- Navigate different download methods
- Validate downloaded data
- Handle large dataset downloads efficiently

## 🔄 Research Workflow

### Starting a New Research Project

**Step 1: Explore & Hypothesize**
```bash
@hypothesis-generator What interesting questions can I explore in this dataset?
```

**Step 2: Setup Infrastructure**
```bash
@research-infrastructure-specialist Setup tools to search PubMed and arXiv
```

**Step 3: Literature Review**
```bash
@literature-analyst Find and synthesize papers on [your topic]
```

**Step 4: Refine Hypothesis**
```bash
@hypothesis-generator Based on the literature, refine my research hypothesis
```

**Step 5: Acquire Data**
```bash
@dataset-specialist Find datasets for testing this hypothesis
@research-infrastructure-specialist Help download the recommended dataset
```

**Step 6: Analysis** (Future: use analysis-phase agents)
- Data cleaning and preprocessing
- Exploratory data analysis
- Statistical testing
- Modeling and visualization

### Investigating a Specific Question

**Example: "How did COVID-19 affect research collaboration?"**

```bash
# Step 1: Refine the question
@hypothesis-generator Help me form a testable hypothesis about COVID's impact on research collaboration

# Step 2: Setup tools
@research-infrastructure-specialist Setup PubMed and Semantic Scholar access for COVID research analysis

# Step 3: Literature review
@literature-analyst Find papers that studied collaboration patterns during COVID-19

# Step 4: Synthesize findings
@literature-analyst Compare methodologies and findings across the top 10 papers

# Step 5: Identify gap
@hypothesis-generator Based on literature gaps, refine my specific research contribution

# Step 6: Get data
@dataset-specialist Find datasets with publication and collaboration data from 2019-2024
```

## 🛠️ Requirements

- [Claude Code](https://docs.claude.com/en/docs/claude-code) installed
- Bash (for the setup script) - works on:
  - macOS
  - Linux
  - Windows (Git Bash, WSL, or Cygwin)

## 📖 Best Practices

### When to Use Each Agent

**Starting Research**:
1. `@hypothesis-generator` - Form initial research questions
2. `@research-infrastructure-specialist` - Setup access to research resources
3. `@literature-analyst` - Understand existing work
4. `@hypothesis-generator` - Refine based on literature
5. `@dataset-specialist` - Acquire relevant data

**Literature Review**:
1. `@research-infrastructure-specialist` - Setup search tools
2. `@literature-analyst` - Search and read papers
3. `@literature-analyst` - Synthesize findings
4. `@literature-analyst` - Identify gaps

**Data Acquisition**:
1. `@dataset-specialist` - Find relevant datasets
2. `@dataset-specialist` - Evaluate quality and appropriateness
3. `@research-infrastructure-specialist` - Help with download and access

### Tips for Best Results

**Be Specific**:
- "Find papers about neural networks" → "Find papers about transformer architectures for biomedical NLP published 2020-2024"
- "I need data" → "I need a dataset of research papers with citation counts and author affiliations"

**Provide Context**:
- Share your research domain
- Mention what you've already tried
- Explain your constraints (time, access, computational resources)

**Iterate**:
- Start broad, then narrow down
- Refine hypotheses as you learn more
- Ask follow-up questions

**Combine Agents**:
- Use multiple agents for comprehensive workflows
- Infrastructure agent enables the others
- Hypothesis agent should bookend your research (start and refine)

## 🔮 Future: Analysis Phase Agents

The current repository focuses on **research discovery**. Future agents will cover the **analysis phase**:

- **data-cleaner** - Data quality and preprocessing
- **eda-specialist** - Exploratory data analysis
- **statistician** - Statistical hypothesis testing
- **feature-engineer** - Feature engineering for modeling
- **ml-engineer** - Machine learning and deep learning
- **nlp-specialist** - Natural language processing (great for papers!)
- **science-communicator** - Visualization and reproducible notebooks

See [future-analysis-phase.md](future-analysis-phase.md) for detailed plans.

## 🆚 How This Differs from Web Dev Agents

**Web Development Agents (e.g., ciign)**:
- Focus: Building software (frontend, backend, databases)
- Tools: React, Node.js, PostgreSQL, Docker
- Workflow: Design → Implement → Test → Deploy

**Research Discovery Agents (this repo)**:
- Focus: Scientific discovery (hypotheses, literature, data)
- Tools: PubMed, arXiv, datasets, research papers
- Workflow: Question → Literature → Hypothesis → Data → Analysis

**Complementary, not competing**: Use web dev agents for building research tools, use research agents for doing research.

## 🤝 Contributing

Contributions are welcome! Here's how you can help:

1. **Fork the repository**
2. **Create a feature branch** (`git checkout -b feature/new-research-agent`)
3. **Make your changes**
4. **Test thoroughly**
5. **Commit your changes** (`git commit -m 'Add new research agent'`)
6. **Push to the branch** (`git push origin feature/new-research-agent`)
7. **Open a Pull Request**

### Ideas for New Research Agents

**Specialized Research**:
- **citation-analyst** - Citation network analysis and impact tracking
- **research-validator** - Evaluate research quality and reproducibility
- **preprint-monitor** - Track and analyze preprints and early research
- **grant-writer** - Help write research proposals and grants

**Domain-Specific**:
- **bioinformatics-specialist** - Genomics, proteomics, clinical data
- **clinical-trial-analyst** - Clinical trial design and analysis
- **meta-analyst** - Systematic reviews and meta-analysis
- **patent-analyst** - Patent search and analysis

**Analysis Phase** (see future-analysis-phase.md):
- **data-cleaner**, **eda-specialist**, **statistician**, **ml-engineer**, etc.

## 📄 License

This project is released into the public domain under The Unlicense. You can do whatever you want with it - no attribution required!

See the [LICENSE](LICENSE) file for details or visit [unlicense.org](https://unlicense.org).

## 🙏 Acknowledgments

- Built for use with [Claude Code](https://www.anthropic.com/claude/code)
- Inspired by [ciign/agentic-engineering](https://github.com/ciign/agentic-engineering) for web development
- Designed to support the scientific research community

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/research-agents/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/research-agents/discussions)
- **Documentation**: [Claude Code Docs](https://docs.claude.com/en/docs/claude-code)

## 🌟 Use Cases

### Academic Research
- PhD students exploring dissertation topics
- Postdocs planning new research directions
- Faculty designing grant proposals
- Research groups coordinating literature reviews

### Industry Research
- R&D teams investigating new technologies
- Data scientists exploring research datasets
- Product teams validating hypotheses with academic literature
- Innovation teams identifying research opportunities

### Meta-Research
- Studying research trends and patterns
- Analyzing collaboration networks
- Understanding citation dynamics
- Mapping research landscapes

### Learning & Education
- Students learning research methodology
- Researchers entering new fields
- Teaching research skills
- Understanding scientific domains

## 💡 Example: Complete Research Project

**Research Question**: "How has international collaboration in AI research evolved?"

```bash
# Phase 1: Hypothesis Formation
@hypothesis-generator
"I'm interested in AI research collaboration patterns. What specific, testable hypotheses could I explore?"

→ Suggests: "International co-authorship in AI research increased 200% from 2015-2024,
with US-China collaboration being the largest bilateral partnership until 2022"

# Phase 2: Infrastructure Setup
@research-infrastructure-specialist
"Setup tools to access arXiv (for AI papers) and Semantic Scholar (for citation data)"

→ Installs and configures MCP servers, provides API keys setup

# Phase 3: Literature Review
@literature-analyst
"Find papers that studied collaboration patterns in AI/ML research"

→ Finds 25 relevant papers, reads top 10, provides synthesis:
- Most studies end in 2020
- Focus on publication counts, not citation impact
- Gap: Recent geopolitical impacts not well studied

# Phase 4: Hypothesis Refinement
@hypothesis-generator
"Based on the literature gaps, refine my hypothesis to make a novel contribution"

→ Refined: "US-China AI collaboration peaked in 2019-2020 and declined 40% by 2023
due to geopolitical tensions, while US-Europe and China-Europe collaborations increased
to compensate, maintaining overall international collaboration rates"

# Phase 5: Data Acquisition
@dataset-specialist
"Find datasets with AI research papers, author affiliations, and citations from 2015-2024"

→ Recommends: Semantic Scholar AI dataset (2M papers) or arXiv CS.AI bulk download

@research-infrastructure-specialist
"Help me download the Semantic Scholar AI dataset"

→ Sets up download script, filters for AI papers, validates data

# Phase 6: Analysis
→ Now use analysis-phase agents or traditional data science workflow
- Clean author affiliations (extract countries)
- EDA: Visualize collaboration trends over time
- Statistical testing: Test significance of decline
- Network analysis: Map collaboration networks
- Visualization: Create publication-quality figures

# Result: Novel research contribution with solid literature foundation
```

## 🎯 Why Use These Agents?

**Save Time**:
- Automate literature search and synthesis
- Quickly setup research infrastructure
- Find datasets faster

**Improve Quality**:
- Form better hypotheses based on literature
- Identify research gaps systematically
- Access more comprehensive data sources

**Learn Research Skills**:
- See how experts approach research
- Learn best practices for literature review
- Understand research methodology

**Stay Organized**:
- Systematic approach to research discovery
- Clear workflow from question to data
- Reproducible research setup

---

**Made with ❤️ for the research community**

*Accelerate your research discovery with Claude Code! 🚀*
# claude-agentic-researchers
