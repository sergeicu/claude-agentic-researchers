# Future: Analysis Phase Agents

## Overview

This document describes the **Analysis Phase** agents that are planned for future implementation. The current repository focuses exclusively on **Research Phase** agents (hypothesis generation, literature review, research infrastructure). The analysis phase agents will complement the research agents by handling the traditional data science workflow once data has been acquired and research context established.

## Why Defer Analysis Phase?

The research phase agents are:
1. **Novel** - No existing agent collections focus on research infrastructure and literature analysis
2. **Enabling** - They set up the foundation for informed analysis
3. **Underserved** - Research discovery tools are rare compared to ML/modeling tools

Analysis phase agents are:
1. **Well-covered** - Many existing tools and workflows already exist
2. **Dependent** - They work best when built on solid research foundation
3. **Can use general-purpose agents** - Claude Code's base capabilities handle many analysis tasks well

## The 7 Analysis Phase Agents

### 1. data-cleaner 🧹

**Primary Focus**: Data quality, preprocessing, and cleaning

**Core Responsibilities**:
- Handle missing values (imputation strategies)
- Detect and handle outliers
- Remove duplicates and inconsistencies
- Standardize data formats
- Validate data quality
- Fix data type issues

**Technical Expertise**:
- pandas, polars for data manipulation
- numpy for numerical operations
- Data validation libraries (great_expectations, pandera)
- Imputation techniques (mean/median, KNN, MICE)
- Outlier detection (IQR, Z-score, isolation forest)

**When to Use**:
- Raw data has quality issues
- Missing values need handling
- Data types are inconsistent
- Outliers need investigation
- Before any analysis or modeling

**Example Use Cases**:
- "This dataset has 30% missing values in key columns, help me clean it"
- "Detect and handle outliers in this medical measurements dataset"
- "Standardize date formats across these merged datasets"

---

### 2. eda-specialist 🔍

**Primary Focus**: Exploratory data analysis and pattern discovery

**Core Responsibilities**:
- Understand data distributions
- Identify correlations and relationships
- Discover patterns and anomalies
- Segment and group data
- Generate initial insights
- Create exploratory visualizations

**Technical Expertise**:
- pandas for data manipulation
- Statistical analysis (descriptive statistics)
- Correlation analysis (Pearson, Spearman, distance correlation)
- Automated EDA tools (pandas-profiling, sweetviz)
- Dimensionality reduction (PCA, t-SNE, UMAP)
- Clustering for exploration (K-means, DBSCAN)

**When to Use**:
- First time exploring a clean dataset
- Looking for patterns and relationships
- Understanding data structure
- Before feature engineering
- Generating analysis ideas

**Example Use Cases**:
- "What are the main patterns and correlations in this clinical trial data?"
- "Explore this dataset and tell me what's interesting"
- "Are there natural groupings in this patient data?"

---

### 3. statistician 📊

**Primary Focus**: Rigorous statistical analysis and hypothesis testing

**Core Responsibilities**:
- Design statistical tests
- Hypothesis testing (t-tests, ANOVA, chi-square)
- A/B test analysis
- Confidence intervals and p-values
- Check statistical assumptions
- Power analysis and sample size calculation

**Technical Expertise**:
- scipy.stats for statistical tests
- statsmodels for advanced statistical modeling
- Hypothesis testing frameworks
- Experimental design principles
- Bayesian statistics (PyMC, Stan)
- Effect size calculation
- Multiple testing correction (Bonferroni, FDR)

**When to Use**:
- Testing specific hypotheses
- Need rigorous statistical inference
- A/B testing results
- Ensuring scientific validity
- Publication-ready analysis

**Example Use Cases**:
- "Test if treatment A is significantly better than treatment B"
- "Is this observed difference statistically significant?"
- "Design a statistical test for this research hypothesis"
- "Calculate required sample size for this experiment"

---

### 4. feature-engineer ⚙️

**Primary Focus**: Feature creation, transformation, and selection

**Core Responsibilities**:
- Create new features from existing data
- Transform features (scaling, encoding, binning)
- Extract features from complex data types (dates, text, nested data)
- Select important features
- Handle categorical variables
- Create interaction terms

**Technical Expertise**:
- scikit-learn preprocessing and feature selection
- Feature engineering techniques (polynomial, interactions)
- Encoding strategies (one-hot, label, target, embeddings)
- Scaling methods (standard, minmax, robust)
- Datetime feature extraction
- Text feature extraction (TF-IDF, n-grams)
- Automated feature engineering (Featuretools)

**When to Use**:
- Preparing data for modeling
- Improving model performance
- Creating domain-specific features
- Handling categorical variables
- Engineering features from dates/text

**Example Use Cases**:
- "Create time-based features from publication dates"
- "Encode these categorical variables appropriately"
- "Generate polynomial and interaction features"
- "Extract features from this JSON column"

---

### 5. ml-engineer 🤖

**Primary Focus**: Machine learning and deep learning model development

**Core Responsibilities**:
- Select appropriate algorithms
- Build and train models (classical ML and deep learning)
- Hyperparameter tuning
- Cross-validation strategies
- Model evaluation and validation
- Prevent overfitting
- Model interpretation

**Technical Expertise**:

**Classical ML**:
- scikit-learn (linear models, trees, ensembles, SVMs)
- xgboost, lightgbm, catboost (gradient boosting)
- Model selection strategies
- Hyperparameter optimization (grid search, random search, Bayesian optimization)

**Deep Learning**:
- PyTorch and/or TensorFlow/Keras
- Neural network architectures (MLPs, CNNs, RNNs, Transformers)
- Transfer learning and fine-tuning
- GPU optimization
- Model checkpointing and versioning

**Model Evaluation**:
- Cross-validation (k-fold, stratified, time-series)
- Metrics selection (accuracy, precision, recall, F1, AUC-ROC, MSE, MAE)
- Learning curves and validation curves
- Model interpretation (SHAP, LIME, feature importance)

**When to Use**:
- Building predictive models
- Classification or regression tasks
- Deep learning for complex patterns
- Model optimization needed
- Comparing multiple algorithms

**Example Use Cases**:
- "Build a gradient boosting model to predict publication impact"
- "Create a neural network for this image classification task"
- "Compare performance of random forest vs XGBoost"
- "Optimize hyperparameters for this model"

---

### 6. nlp-specialist 📝

**Primary Focus**: Natural language processing and text analysis

**Core Responsibilities**:
- Text preprocessing and cleaning
- Named entity recognition (NER)
- Text classification and sentiment analysis
- Topic modeling
- Text embeddings and semantic similarity
- Information extraction
- Language model fine-tuning

**Technical Expertise**:
- transformers (Hugging Face) for modern NLP
- spaCy for efficient NLP pipelines
- NLTK for traditional NLP tasks
- Text preprocessing (tokenization, lemmatization, stopwords)
- Pre-trained models (BERT, GPT, domain-specific models like BioBERT, SciBERT)
- Named entity recognition
- Topic modeling (LDA, BERTopic)
- Text classification and sequence labeling
- Semantic search and embeddings

**When to Use**:
- Working with text data (abstracts, papers, clinical notes)
- Extracting entities from text
- Classifying documents
- Finding similar texts
- Analyzing research papers

**Example Use Cases**:
- "Extract disease and drug names from these PubMed abstracts"
- "Classify research papers into research areas"
- "Find papers similar to this abstract using semantic search"
- "Perform topic modeling on 10,000 paper abstracts"
- "Fine-tune BioBERT for this medical text classification task"

---

### 7. science-communicator 📈

**Primary Focus**: Visualization, documentation, and reproducible research

**Core Responsibilities**:

**Visualization**:
- Create publication-quality figures
- Design effective data visualizations
- Follow visualization best practices
- Interactive dashboards

**Reproducibility**:
- Organize notebooks for clarity
- Document analysis workflow
- Create reproducible environments
- Version control best practices
- Write clear analysis documentation

**Technical Expertise**:

**Visualization**:
- matplotlib for publication figures
- seaborn for statistical visualizations
- plotly for interactive visualizations
- Visualization design principles
- Color theory for accessibility
- Figure composition and layout

**Reproducibility**:
- Jupyter notebook organization
- Documentation best practices
- Environment management (conda, pip, requirements.txt)
- Code organization and modularity
- Git for version control
- Testing analysis code

**When to Use**:
- Creating visualizations for papers/presentations
- Making analysis reproducible
- Organizing messy notebooks
- Documenting findings
- Sharing research results

**Example Use Cases**:
- "Create publication-quality figures for this analysis"
- "Organize this messy notebook into a clean, reproducible workflow"
- "Design an interactive dashboard for exploring this data"
- "Document this analysis so others can reproduce it"

---

## How Analysis Agents Work Together

**Typical Analysis Workflow**:

1. **@data-cleaner** → Clean and preprocess raw data
2. **@eda-specialist** → Explore patterns and generate insights
3. **@statistician** → Test statistical hypotheses rigorously
4. **@feature-engineer** → Create features for modeling
5. **@ml-engineer** OR **@nlp-specialist** → Build predictive models
6. **@ml-engineer** → Evaluate and interpret models
7. **@science-communicator** → Create visualizations and documentation

**Example: Complete Analysis Pipeline**:

```
Research phase complete → Data acquired

@data-cleaner
"Clean this PubMed collaboration dataset"
→ Handles missing author affiliations, standardizes country names

@eda-specialist
"Explore collaboration patterns across countries"
→ Identifies surge in international collaborations during COVID-19

@statistician
"Test if international collaboration rate significantly increased"
→ Runs statistical tests, confirms significance with p<0.001

@feature-engineer
"Create features: collaboration network metrics, temporal features"
→ Generates graph centrality, time-based features

@nlp-specialist
"Extract research topics from abstracts"
→ Uses BERTopic to identify 50 main research topics

@ml-engineer
"Predict which papers will have high collaboration impact"
→ Builds XGBoost model, achieves 0.85 AUC

@science-communicator
"Create publication figures and reproducible notebook"
→ Generates publication-quality visualizations, documents entire workflow
```

---

## Integration with Research Phase

The analysis phase agents are designed to work seamlessly after the research phase:

**Research Phase** (Current repository focus):
1. @hypothesis-generator → "What questions should we ask?"
2. @research-infrastructure-specialist → "Setup tools to access data/papers"
3. @literature-analyst → "What has been done before?"
4. @hypothesis-generator → "Refined hypothesis based on literature"

**Analysis Phase** (This document):
5. @data-cleaner → Clean acquired data
6. @eda-specialist → Explore the data
7. @statistician → Test hypotheses
8. ... (continue with modeling, visualization)

---

## Implementation Priority

When implementing these agents in the future, suggested priority order:

**Tier 1 (Essential)**:
1. **data-cleaner** - Almost always needed first
2. **eda-specialist** - Critical for understanding data
3. **science-communicator** - Needed to share results

**Tier 2 (High Value)**:
4. **statistician** - Important for scientific rigor
5. **nlp-specialist** - Especially valuable for text-heavy research

**Tier 3 (Specialized)**:
6. **feature-engineer** - Needed when building models
7. **ml-engineer** - Needed for predictive modeling tasks

---

## Design Philosophy

These agents follow the same principles as the research phase agents:

1. **Focused Expertise**: Each agent has a clear, specific role
2. **Complementary**: Agents work together without overlapping
3. **Practical**: Solve real problems data scientists face
4. **Tool-Specific Knowledge**: Deep expertise in relevant libraries
5. **Best Practices**: Encode domain knowledge and best practices
6. **Clear Invocation**: Obvious when to use each agent

---

## Future Enhancements

Additional specialized agents could be added:

- **time-series-analyst** - Forecasting, seasonality, temporal patterns
- **experiment-designer** - A/B testing, causal inference, experimental design
- **bioinformatics-specialist** - Genomics, proteomics, clinical trial analysis
- **geospatial-analyst** - Geographic data, mapping, spatial analysis
- **computer-vision-specialist** - Image analysis, medical imaging
- **recommender-specialist** - Recommendation systems, collaborative filtering

---

## Conclusion

The analysis phase agents represent the traditional data science workflow. While deferred for now, they are well-designed and ready for implementation when needed. The current focus on research phase agents provides unique value by addressing the often-overlooked front end of scientific discovery: forming hypotheses, reviewing literature, and setting up research infrastructure.

Together, the research phase and analysis phase agents will provide comprehensive coverage of the entire scientific research workflow, from initial curiosity to publication-ready results.
