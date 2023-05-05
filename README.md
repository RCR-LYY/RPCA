# FAMCSMMA

## Method Description

FAMCSMMA is a fast and efficient method based on RPCA framework to predict potential associations between small molecule drugs and miRNAs . The datasets and codes for the model are released on this page.

## Operating Environment

PyCharm == 2021.1.3       
Python == 3.10  
Windows == 11            
Processor == AMD Ryzen 7 5800H 3.20GHz CPU with 6G RAM       

## Required Packages

numpy == 1.24.1     
matplotlib == 3.6.2    
scipy == 1.10.0     
scikit-learn == 1.2.0   

## Running steps:

**Step 1**: Run the file “GIPK.py”. Its input is the SM--miRNA association matrix (file "SM-miRNA association matrix"), the existing SM similarity matrix (file "SM similarity matrix") and the existing miRNA similarity matrix (file "miRNA similarity matrix"). Its output is the adjacency matrix of the heterogeneous SM-miRNA network involving the GIPK similarity (file "T-GIP.txt"), which is the target matrix for our algorithm.       
**Step 2**: Run the file “Matrix Completion.py”. Its input is the adjacency matrix of the heterogeneous SM-miRNA network (file "T-GIP.txt") and the existing SM similarity matrix (file "SM similarity matrix"). Its output is the SM-miRNA prediction score matrix. The (i,j)-th element of the predictive score matrix, $m'_{ij}$, represents the predictive score between $SM_i$ and $miRNA_j$.

## Contact

If you have any problems or find mistakes, please feel free to contact me: z21070251@s.upc.edu.cn

