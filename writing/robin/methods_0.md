# Faunal references
The faunal reference genomes used for the simulations and further analyses were sourced from the "National Center for Biotechnology Information" *(or is "NCBI" enough?)* **[R1]** (see **table 1**).

*(What should I say about the additions "we" did to the human genome? Like adding phages or something?)*

*(Do we have to justify, why we used these references? For example, I used CanFam3.1, which is outdated, but I used this because I had that ready from the capture design project.)*

*(Should we reorder them and say that the are ordered by distance to human?)*

Num | Common name | Taxon               | Assembly name               | Link
--- | ----------- | ------------------- | --------------------------- | ----
0   | human       | *Homo sapiens*      | modified hg19               | -
1a  | tamarin     | *Saguinus midas*    | ASM2149847v1                | **[R2]**
1b  | tarsier     | *Carlito syrichta*  | Tarsius_syrichta-2.0.1      | **[R3]**
1   | macaque     | *Macaca mulatta*    | Mmul_10                     | **[R4]**
2   | dog         | *Canis lupus*       | CanFam3.1                   | **[R5]** 
3   | bear        | *Ursus americanus*  | gsc_jax_bbear_1.0           | **[R6]**
4   | bison       | *Bison bison bison* | Bison_UMD1.0                | **[R7]**
5   | mouse       | *Mus musculus*      | GRCm39                      | **[R8]**
6   | chicken     | *Gallus gallus*     | bGalGal1.mat.broiler.GRCg7b | **[R9]**
7   | maize       | *Zea mays*          | Zm-B73-REFERENCE-NAM-5.0    | **[R10]**
8   | drosophila  | *Drosophila melanogaster*  | Release 6 plus ISO1 MT | **[R11]** 
9   | yeast       | *Saccharomyces cerevisiae* | R64                    | **[R12]** 
X   | anthrax     | *Bacillus anthracis*       | ASM844v1               | **[R13]**

---

---
# DNA simulation + Mapping
The simulations of ancient Illumina DNA reads were performed using **Gargammel**. 
The nick frequency, length of overhanging ends, probability of deamination of Cs in overhanging double- and single-stranded parts were set to their standard-parameters of 0.03, 0.4, 0.01 and 0.3 *(This came from somewhere, right? Maybe cite this?)* **[R14]**.
The read-length distribution originates from previous observations of real-world Illumina sequencings of ancient *(environmental?)* DNA **[R15]**.
*(Should I say more about what exactly Gargammel does?)*
For each of the 13 reference genomes, 10 Million reads were simulated and 0.2 % of the nucleotides were additionally mutated to a random base.
The simualted reads were mapped to the human reference genomes using **bwa bam2bam**, a modified version of **bwa aln** that takes bams both as input for query sequences and output.
*(Should I even mention the N reference? We decided the third was better, right?)*
Three different versions of the human reference genome were used; unmodified, flipped to "N" and flipped to "third allele".
For the genome flippep to the "third allele", the nucleotides at all diagnostic SNP positions *(how should I reference that?)* were changed to a random base that is neither derived nor ancestral at that position.

---


<!--
# Reference bias
Text
--->

<!--
# Faunal mismapping quantification + program
--->

<!--
# faunalizer program
Text
--->

<!--
# Popgen details
Text
--->

<!--
# Gelabert
Text
--->

<!--
# Competitive mapping
Text
--->

---
# References
Num       | Link
--------- | ----
**[R1]**  | https://www.ncbi.nlm.nih.gov
**[R2]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_021498475.1/
**[R3]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000164805.1/
**[R4]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_003339765.1/
**[R5]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000002285.3/
**[R6]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_020975775.1/
**[R7]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000754665.1/
**[R8]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001635.27/
**[R9]**  | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_016699485.2/
**[R10]** | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_902167145.1/
**[R11]** | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001215.4/
**[R12]** | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000146045.2/
**[R13]** | https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000008445.1/
**[R14]** | Information about gargammel standard damage parameter
**[R15]** | Information about gargammel read-length distribution

---

---
# Tools
Num       | Name      | Version | Link
--------- | --------- | ------- | ----
**[T1]**  | Gargammel | 1.1.2   | https://pubmed.ncbi.nlm.nih.gov/27794556/