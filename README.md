# modZN – automated overlap & primer designer for modAssembly
This repository contains the Python implementation of **modZN**, the computational design tool accompanying:

> **“Modular multi-fragment enzymatic assembly based on exonuclease resistant bases”**  
> Rong Chen, Shuang Feng, Ruohan Ren, Zhiguang Wu, Zixuan Bian, Richard Lakerveld, Ninuo Xia, Fang Fang, and Bryan Wei


In the paper, **modAssembly** is a novel DNA assembly method that uses **phosphorothioate-modified bases** to generate precise 3′ sticky ends resistant to exonuclease digestion.  
**modZN** is the computational counterpart — an automated primer and overlap sequence designer ensuring all sticky ends have comparable hybridization free energies (ΔG), thereby maximizing one-pot assembly efficiency.

This repository provides the **GUI-based Python tool** that implements modZN, with versions for both **circular** and **linear** DNA assemblies. If you would like to use modZN more conveniently, the web page tool can be accessed at http://8.130.46.70/#/seq.

---

## 🧬 Background

Traditional methods like Gibson or Golden Gate assemblies suffer from sequence restrictions or declining efficiency as fragment count increases.  
**modAssembly**, introduced in the cited paper, enables assembly of **up to 20 fragments** with high fidelity and flexibility by combining:

- PCR using **phosphorothioate-modified primers**,
- controlled λ-exonuclease digestion to produce **defined sticky ends**, and
- one-pot **self-assembly** through complementary overhangs.

The **modZN** program automates the design of these overlapping sequences and primers.


![workflow](https://github.com/Renruohan/modZN/blob/main/workflow.png)
The **Workflow** of modAssembly. (A) The schematic workflow of three-fragment **modAssembly**, including fragment amplification using modified primers (step I), lambda exonuclease digestion (step II) and fragment hybridization (step III). The asterisks denote the phosphorothioate modifications. The green, blue and orange segments denote the overlapping segments of F3/F1, F1/F2 and F2/F3, respectively. (B) The primer design of **modZN**. The process includes the search of overlapping segments with specific ΔG requirement (step 1) and the generation of full-length primer sequences with phosphorothioate modifications (step 2).

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

## 💻 Usage
Since we use IDT for primer design calculations, please first register an IDT account and enter it in **./modZN/IDT_account.txt**.
Then, simply run the following command.

- Python GUI_for_modZN.py

[You can adjust parameters such as **length_of_search_frame** (the length of the search window for finding suitable overlap regions) and **t** (the threshold for the ΔG difference) in **modZN_circular_add_IDT_TM.py** or ***modZN_linear_add_IDT_TM.py** if you need.]

## 📘 User Manual

For detailed usage instructions, parameter explanations, and example workflows, please refer to the **User Manual** provided in:

> 📄 **`modZN_User_Manual.pdf`**



