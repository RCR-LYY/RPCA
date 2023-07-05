# RPCAΓNR

## Method Description

RPCAΓNR is a fast and efficient method based on RPCA framework to predict potential associations between small molecule drugs and miRNAs . The datasets and codes for the model are released on this page.

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

## File Description

### Dataset: 

This folder contains three datasets (dataset 1, dataset 2, and dataset 3). Taking dataset 1 as an example, it includes the following files:       
**SM number.xlsx**: It contains the CID of 831 SMs and the row number of each SM in the association matrix. For example, the first row of the association matrix represents SM 'CID:137'.   
**miRNA number.xlsx**: It contains the name of 541 human-related miRNAs and the column number of each miRNA in the association matrix. For example, the first column of the association matrix represents miRNA 'hsa-let-7a-1'.   
**SM similarity matrix.txt**: It contains the integrated SM similarity matrix (the dimension is 831 $\times$ 831), where each row and its corresponding column represent a specific SM. The element in the i-th row and j-th column denotes the similarity score between $SM_i$ and $SM_j$.       
**miRNA similarity matrix.txt**: It contains the integrated miRNA similarity (the dimension is 541 $\times$ 541), where each row and its corresponding column represent a specific miRNA. The element in the i-th row and j-th column denotes the similarity score between $miRNA_i$ and $miRNA_j$.     
**SM-miRNA-confirmed associations.txt**: It contains 664 SM-miRNA known associations. For example, the first row in the file, (75,177), denotes that the SMs represented in the 75-th row of the association matrix, 'CID:2578', is associated with the miRNA represented in the 177-th column of the association matrix, 'hsa-mir-200c'.      
**SM-miRNA association matrix.txt**: It is the constructed association matrix (the dimension is 831 $\times$ 541), each row of which represents a specific SM, and each column represents a specific miRNA. The (i,j)-th element of the association matrix, $m_{ij}$, is set to 1 if $SM_i$ is associated with $miRNA_j$, otherwise it is set to 0. The matrix has a total of 664 zeros and 448907 ones.         

### GIPK.py:

This file contains the Python code for calculating the Gaussian Interaction Profile Kernel (GIPK) Similarity of SMs and miRNAs, respectively. 

### Matrix Completion.py:   

This file contains the Python code for our algorithm. 

## Running steps:

**Step 1**: Run the file “GIPK.py”. Its input is the SM--miRNA association matrix (file "SM-miRNA association matrix"), the existing SM similarity matrix (file "SM similarity matrix") and the existing miRNA similarity matrix (file "miRNA similarity matrix"). Its output is the adjacency matrix of the heterogeneous SM-miRNA network involving the GIPK similarity (file "T-GIP.txt"), which is the target matrix for our algorithm.       
**Step 2**: Run the file “Matrix Completion.py”. Its input is the adjacency matrix of the heterogeneous SM-miRNA network (file "T-GIP.txt") and the existing SM similarity matrix (file "SM similarity matrix"). Its output is the SM-miRNA prediction score matrix. The (i,j)-th element of the predictive score matrix, $m'_{ij}$, represents the predictive score between $SM_i$ and $miRNA_j$.

## Contact

If you have any problems or find mistakes, please feel free to contact me: z21070251@s.upc.edu.cn

