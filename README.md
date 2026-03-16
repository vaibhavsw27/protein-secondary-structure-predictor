# protein-secondary-structure-predictor
Python implementation of the Chou-Fasman algorithm for predicting protein secondary structure (alpha helices and beta strands) from amino acid sequences.


# Chou-Fasman Secondary Structure Prediction

Python implementation of the Chou-Fasman algorithm for predicting protein secondary structure from amino acid sequences.

This project predicts:
- Alpha helices (H)
- Beta strands (S)

based on amino acid propensities using the classical Chou-Fasman method.

---

## Overview

Proteins fold into secondary structures such as alpha helices and beta sheets.  
The Chou-Fasman algorithm predicts these structures by analyzing the statistical propensities of amino acids to form helices or strands.

This program:

1. Scans the protein sequence
2. Identifies nucleation regions
3. Extends predicted regions
4. Resolves conflicts between helix and strand predictions
5. Outputs predicted secondary structure maps and regions

---

## Features

- Alpha helix prediction
- Beta strand prediction
- Conflict resolution using residue propensities
- Region detection with start/end indices
- Raw and final prediction maps

---

## Requirements

Python 3.x

No external libraries required.

---

## How to Run

Clone the repository:

git clone https://github.com/vaibhavsw27/protein-secondary-structure-predictor.git

Navigate to the folder:

cd protein-secondary-structure-predictor

Run the script:

python chou_fasman.py

You will be prompted to enter a protein sequence.

Example:

Please enter your sequence:
MAQWNQLQQLDTRYLEQLHQLYSDSFPMELR

---

## Output

The program prints:

### Helix Predictions
- Raw helix regions
- Final helix regions
- Helix maps

### Beta Strand Predictions
- Raw strand regions
- Final strand regions
- Strand maps

### Final Secondary Structure

Example:

Final Secondary Structure
---HHHHHHH---SSSSS---HHHH

Where:

H = Alpha helix  
S = Beta strand  
- = No predicted structure

---

## Algorithm

The implementation follows the classical **Chou-Fasman method**, which:

1. Uses amino acid propensities for helix and strand formation
2. Detects nucleation windows
3. Extends regions based on average propensities
4. Resolves helix/strand conflicts using higher residue propensity

---

## Example Input

MAQWNQLQQLDTRYLEQLHQLYSDSFPMELR

---

## File Structure

protein-secondary-structure-predictor
│
├── chou_fasman.py
├── example_sequence.txt
└── README.md
---

## Author

Vaibhav Samrat Waghaye  
B.Tech Computer Engineering
