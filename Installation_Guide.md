# Detailed Installation Guide

## Table of Contents
1. [Quick Start](#quick-start)
2. [System Requirements](#system-requirements)
3. [Step-by-Step Installation](#step-by-step-installation)
4. [Java Configuration](#java-configuration)
5. [Verification](#verification)
6. [Common Issues and Solutions](#common-issues-and-solutions)

---

## Quick Start

For users with working conda and basic familiarity:

```bash
git clone https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline-updated.git
cd Virus-Variant-Calling-Pipeline-updated
conda env create -f environment.yml
conda activate dengue_pipeline
pip install -r requirements.txt
pip install .
java -version  # Verify Java is available
