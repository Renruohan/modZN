# modZN – automated overlap & primer designer for modAssembly
This repository contains the Python implementation of **modZN**, the computational design tool accompanying:

> **“Modular multi-fragment enzymatic assembly based on exonuclease resistant bases”**  
> Rong Chen, Shuang Feng, Ruohan Ren, Zhiguang Wu, Zixuan Bian, Richard Lakerveld, Ninuo Xia, Fang Fang, and Bryan Wei


In the paper, **modAssembly** is a novel DNA assembly method that uses **phosphorothioate-modified bases** to generate precise 3′ sticky ends resistant to exonuclease digestion.  
**modZN** is the computational counterpart — an automated primer and overlap sequence designer ensuring all sticky ends have comparable hybridization free energies (ΔG), thereby maximizing one-pot assembly efficiency.

This repository provides the **GUI-based Python tool** that implements modZN, with versions for both **circular** and **linear** DNA assemblies.

---

## 🧬 Background

Traditional methods like Gibson or Golden Gate assemblies suffer from sequence restrictions or declining efficiency as fragment count increases.  
**modAssembly**, introduced in the cited paper, enables assembly of **up to 20 fragments** with high fidelity and flexibility by combining:

- PCR using **phosphorothioate-modified primers**,
- controlled λ-exonuclease digestion to produce **defined sticky ends**, and
- one-pot **self-assembly** through complementary overhangs.

The **modZN** program automates the design of these overlapping sequences and primers.

---


## ⚙️ What modZN Does

The tool automates the **rational design** of overlaps and primers for modAssembly by:

1. **Scanning** through fragment ends to find candidate overlaps.  
2. **Calculating ΔG values** using the **nearest-neighbor thermodynamic model**.  
3. **Balancing ΔG** across all interfaces (ΔΔG < threshold) to achieve similar melting behavior.  
4. **Generating PCR primers** with phosphorothioate modifications (represented as `*`).  
5. **Calling IDT’s OligoAnalyzer** (via Selenium) to calculate realistic melting temperatures (Tm).  
6. **Displaying results in a GUI** and exporting primer lists as a `.csv` file.

There are two independent GUI tools provided:
- **modzn_circular.py** – for plasmid/circular assemblies.  
- **modzn_linear.py** – for linear assemblies.

---

## 📦 Requirements

- **Python:** 3.8 or higher  
- **Dependencies:**  
  pip install selenium webdriver-manager PySimpleGUI

---

## ▶️ Usage

python modzn_circular.py


