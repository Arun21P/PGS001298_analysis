## Requirements
- Python 3.9+
- pandas
- matplotlib

# Polygenic Score (PGS) Analysis – PGS001298

This repository contains a Python pipeline to **preprocess, analyze, and visualize** a Polygenic Score (PGS) scoring file, specifically **PGS001298_hmPOS_GRCh38**.

The workflow includes:

- Unzipping the `.gz` scoring file
- Cleaning and harmonizing the data
- Exploratory data analysis
- Chromosome-specific effect weight visualization

---

## 📁 Repository Structure
```
PGS001298_analysis/
│
├── data/
│   └── PGS001298_hmPOS_GRCh38.txt.gz
│
├── output/
│   ├── PGS001298_hmPOS_GRCh38_cleaned.txt
│   ├── Distribution_of_Effect_Weight.png
│   └── chr21_effect_weight_hist.png
│
├── scripts/
│   ├── A_assignment.py
│   ├── B_assignment.py
│   └── PGS_analysis_google_colab.ipynb
│
├── Coding Assessment - Bioinformatics.docx
├── requirements.txt
├── .gitignore
└── README.md
```

