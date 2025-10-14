# DNA BASICS FOR ALGORITHMIC EXPLORATION

This document provides a **lightweight introducton to DNA and DNA sequencing concepts** that are relevant to the computational modeling.
The document is not a **biology textbook** - as it only covers the parts needed for building and experimenting with DNA sequence algorithms using a randomly generated DNA sequence.

## What is Deoxyribonucleic Acid?

DNA or Deoxyribonucleic acid is a **polymer** composed of two **polynucleotide** chains that coil around each other to form a double helix. The polymer carries genetic instructions for the development, functioning, growth and reproduction of all known organisms and many viruses. DNA and ribonucleic acid (RNA) are **nucleic acids**. Alongside proteins, lipids and complex carbohydrates (polysaccharides), nucleic acids are one of the four major types of macromolecules that are essential for all known forms of life.

## 🧬 How DNA is Built

DNA is made up of two long strands called **polynucleotides**, and each strand is itself a chain of smaller units called **nucleotides**.

### A single nucleotide has three parts

- a nitrogenous base (A, T, C, or G),

- a sugar molecule ***(deoxyribose)***, and

- a phosphate group.

*Nucleotides connect to each other through a **phosphodiester bond**, which links the sugar of one nucleotide to the phosphate of the next. This repeating connection creates a strong **sugar–phosphate backbone**, with the nitrogenous bases projecting out from the chain.*

## More on Nucleotide Bonding:

- Nucleotides are bonded to each other through **condensation reactions**
. Multiple nucleotides bonded together form polynucleotides. Water is released as a by-product.

- Phosphodiester bonds form between nucleotides. A phosphodiester bond forms between the phosphate group of one nucleotide and the pentose sugar of a second nucleotide. Phosphodiester bonds are covalent bonds.

- Phosphodiester bonds create a chain of nucleotides called a polynucleotide chain. These polynucleotide chains also have a sugar-phosphate backbone because of the pentose sugars and the phosphate groups which make up the majority of a nucleotide’s structure.

- Phosphodiester bonds can be broken through hydrolysis reactions.

**DNA becomes double-stranded when two of these polynucleotide chains align in opposite directions (antiparallel) and the bases pair up according to specific rules:**

- A always pairs with T (two hydrogen bonds)

- C always pairs with G (three hydrogen bonds)

The result is the stable, twisted shape known as the DNA double helix.

## What are Monomers?

Monomers are small, single unit molecules that acts as the **building blocks** for more complex structures called polymers.

### What are the four types of monomers?

The four monomers that make up biologically important macromolecules are :-

- Monosaccharides
  - Make up polysaccharides, or carbohydrates
- Fatty Acids
  - Make up triglycerides, or lipids
- Amino Acids
  - Make up polypeptides, or proteins
- Nucleotides
  - Make up polynucleotides, or nucleic acids

## What are Polymers?

Polymers are formed when numerous monomers link together in a repeating chain.

### What are the four main polymers?

The four main biologically important polymers are :-

- Polysaccharides
- Triglycerides
- Polypeptides
- Polynucleotides

### Synthetic Polymers

- Polystyrene
- Polyethylene
- Polybutadiene

## What are Macromolecules?

Macromolecules are large, complex molecules essential for life, formed by joining many smaller repeating units called monomers.

### What are the four main macromolecules?

The four main macromolecules are :-

- Carbohydrates
- Lipids
- Proteins
- Nucleic Acids

### Synthetic Macromolecules

- Styrofoam
- Plastic
- Artificial Rubber

## 🧾 Quick Reference: Monomers → Polymers → Macromolecules

| Macromolecule Class | Monomer (building block) | Polymer (chain)      | Example Macromolecule |
|---------------------|--------------------------|----------------------|------------------------|
| Proteins            | Amino acids              | Polypeptides         | Enzymes, hemoglobin    |
| Carbohydrates       | Monosaccharides (sugars) | Polysaccharides      | Starch, glycogen       |
| Nucleic acids       | Nucleotides (A, C, G, T/U)| Polynucleotides      | DNA, RNA               |
| Lipids*             | (Varied, not always true monomers) | Complex assemblies/Triglycerides | Fats, oils, membranes |

\*Lipids are the odd one out — large, but not always polymers.

### 🔗 Phosphodiester Bond Formation

The **phosphodiester bond** is the covalent linkage that connects nucleotides together to form a continuous **polynucleotide chain** — the backbone of DNA and RNA.

Each nucleotide consists of a sugar (deoxyribose in DNA, ribose in RNA), a phosphate group, and a nitrogenous base. The carbons of the sugar are numbered **1′ to 5′**.

When nucleotides join together:

- The **phosphate group** attached to the **5′ carbon** of one nucleotide’s sugar forms a bond with the **hydroxyl group (-OH)** on the **3′ carbon** of the sugar of the next nucleotide.  
- This **condensation reaction** releases a molecule of **water (H₂O)** as a by-product.  
- The result is a **phosphodiester linkage**, creating a strong **sugar–phosphate backbone** that supports the DNA structure.  
- DNA strands therefore have directionality — a **5′ end** (with a free phosphate group) and a **3′ end** (with a free hydroxyl group).  
- New nucleotides are always added to the **3′ end** during synthesis.

---

#### ⚙️ Energy and Enzymatic Support

Formation of phosphodiester bonds requires **energy**, supplied by the breakdown of **nucleoside triphosphates (NTPs)** or **deoxynucleoside triphosphates (dNTPs)**.  

- When a nucleotide is incorporated into the chain, **two of its three phosphate groups (pyrophosphate)** are cleaved off, releasing the energy needed for bond formation.  
- This process is **enzyme-catalyzed**, primarily by **DNA or RNA polymerases**, which position and link nucleotides using **magnesium ions (Mg²⁺)** as cofactors.

---

#### 🧬 Enzymes Involved

| Enzyme | Role in Phosphodiester Bond Formation |
|---------|----------------------------------------|
| **DNA Polymerase** | Adds new nucleotides to a growing DNA strand during replication. Uses dNTPs and requires Mg²⁺ ions for catalysis. |
| **RNA Polymerase** | Synthesizes RNA using NTPs during transcription. |
| **DNA Ligase** | Seals gaps (“nicks”) in the sugar–phosphate backbone during DNA repair or Okazaki fragment joining on the lagging strand. |
| **Primase** | Synthesizes short RNA primers needed for DNA polymerase to begin synthesis. |

---

#### ⚡ Charge and Structural Implications

- Each **phosphodiester bond** carries a **negative charge** at physiological pH (~7).  
- These negative charges attract **positively charged proteins** (like histones), **metal ions** (e.g., Mg²⁺), and **polyamines** that help stabilize DNA structure.  
- Repulsion between adjacent negative charges also influences DNA’s **conformation** and contributes to its **double helical stability**.

---

#### 🧩 Visual Representation


```bash
5' — Phosphate — Sugar — Phosphate — Sugar — 3'
       |            |            |
      Base         Base         Base
```


The chain grows continuously in the **5′ → 3′ direction**, forming one side of the DNA double helix backbone.

## 🧭 Directionality and Antiparallel Orientation

In molecular biology, directionality refers to the orientation of a nucleic acid strand — defined by the specific numbering of the carbon atoms in the sugar molecule (**deoxyribose** in DNA).

### Each DNA strand has two distinct ends:

**5′ (five-prime) end**:

- The phosphate group is attached to the 5′ carbon of the sugar molecule.

**3′ (three-prime) end**:

- The hydroxyl group (-OH) is attached to the 3′ carbon of the sugar molecule.

When nucleotides join together via phosphodiester bonds, the phosphate on the 5′ carbon of one nucleotide bonds with the –OH on the 3′ carbon of another.
This reaction gives the strand a chemical direction — it runs from the 5′ end → 3′ end.

#### Antiparallel Arrangement in DNA

DNA consists of two complementary strands that run in opposite directions:

- One strand runs 5′ → 3′

- The other runs 3′ → 5′

This antiparallel structure allows the bases of one strand to pair precisely with those of the other, forming the double helix.

### Quick Visualization

```bash
  5′ — A T G C T A G C — 3′
  3′ — T A C G A T C G — 5′
```

Directionality is critical for DNA functions such as **replication** and **transcription**, as enzymes like DNA polymerase can only add new nucleotides to the 3′ end of a growing strand.

## 🧩 Base Pairing and Complementarity

In the DNA double helix, the two strands are held together by **hydrogen bonds** between the **nitrogenous bases** of each strand. This specific pairing of bases is known as **complementary base pairing**.

### Complementary Base Pairing Rules

Each base pairs with only one specific partner:

- **Adenine (A)** always pairs with **Thymine (T)** → forms **two hydrogen bonds**
- **Cytosine (C)** always pairs with **Guanine (G)** → forms **three hydrogen bonds**

These specific interactions are due to the **shape** and **hydrogen bonding sites** of the bases.  
Because of this, **A–T** and **C–G** pairs fit perfectly within the uniform width of the DNA double helix.

### Importance of Complementarity

Complementary base pairing ensures:

- **Accurate DNA replication** → each strand serves as a **template** for building its complementary strand.
- **Stable double helix structure** → consistent base pairing keeps the DNA molecule uniform and strong.
- **Faithful genetic transmission** → when cells divide, genetic information is copied precisely.

\*Each base on one strand has its exact complementary partner on the opposite strand.

### 🧬 Role of Base Pairing in DNA Function

Complementary base pairing is not just a structural property — it is the **core mechanism that allows DNA to replicate faithfully, remain chemically stable, and carry genetic information across generations**.  

---

#### 1. Accurate DNA Replication

During DNA replication, the double helix **unzips**, and each strand serves as a **template** for the synthesis of a new complementary strand.

- Because of the **base pairing rules** (A pairs only with T, and C pairs only with G), each base on the original strand dictates which base will be added to the new strand.  
- This process is carried out by the enzyme **DNA polymerase**, which moves along the template strand, **adding complementary nucleotides** in the 5′ → 3′ direction.  
- As a result, each daughter DNA molecule contains **one original (parental) strand** and **one newly synthesized strand** — a process known as **semi-conservative replication**.  

This mechanism ensures **genetic fidelity**, meaning that the DNA sequence is copied **with minimal errors**, preserving genetic information during cell division.

> 🧩 Example:  
> If one strand has the sequence 5′–ATCGTACG–3′,  
> its complementary strand will automatically be 3′–TAGCATGC–5′.

---

#### 2. Structural Stability of the Double Helix

Base pairing plays a critical role in maintaining the **physical stability** and **three-dimensional structure** of DNA.

- **Hydrogen bonding** between complementary bases holds the two strands together — two bonds between A and T, and three between C and G.  
- The presence of more C–G pairs increases stability because of the **extra hydrogen bond**, which is why **GC-rich regions** are more thermally stable.  
- The pairing creates a **uniform helical structure**, since the distance between A–T and C–G pairs is consistent, keeping the double helix evenly spaced.  
- In addition to hydrogen bonds, **base stacking interactions** — hydrophobic and van der Waals forces between adjacent base pairs — add extra strength and protect the bases from water, stabilizing the molecule even further.

Together, these interactions give DNA both **rigidity and flexibility**, making it stable enough to store information but dynamic enough for replication and transcription.

---

#### 3. Transmission of Genetic Information

Complementary base pairing ensures that genetic information can be **accurately read, copied, and passed on**.

- The **sequence of bases** (A, T, C, G) encodes instructions for building proteins — the foundation of life processes.  
- During **transcription**, one DNA strand (the template strand) is used to synthesize a complementary **messenger RNA (mRNA)** molecule using the same base-pairing principle — except **uracil (U)** replaces **thymine (T)**.  
- During **translation**, the mRNA sequence is read in triplets (codons) to determine the sequence of amino acids in a protein.  
- Because the DNA base pairs are complementary, the information in one strand is **redundantly stored** in the other — a built-in error-correction feature that helps detect and repair damage.

Thus, base pairing provides both the **physical mechanism for genetic replication** and the **logical framework for genetic coding** — enabling the inheritance of traits and continuity of life.

---

#### 🔁 Summary

| Function | How Base Pairing Contributes |
|-----------|------------------------------|
| **Replication** | Ensures each strand serves as a perfect template, allowing precise duplication of genetic material. |
| **Stability** | Hydrogen bonds and base stacking maintain DNA’s uniform double-helical structure and thermal resistance. |
| **Information Transmission** | Preserves and transfers genetic instructions through accurate transcription and translation. |
