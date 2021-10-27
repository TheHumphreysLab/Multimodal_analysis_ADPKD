# **Multimodal single cell analysis on human ADPKD kidneys**
__Yoshiharu Muto, Eryn E. Dixon, Yasuhiro Yoshimura, Haojia Wu, Kohei Omachi, Andrew J. King, Eric Olson, Marvin Gunawan, Jay Kuo, Jennifer Cox, Jeffrey H. Miner, Stephen L. Seliger, Owen M. Woodward, Paul A. Welling, Terry J. Watnick and Benjamin D. Humphreys__  

This work is published on the following [manuscript](https://pubmed.ncbi.nlm.nih.gov/xxxxxxx/)
```
Defining cellular complexity in human autosomal dominant polycystic kidney disease by multimodal single cell analysis

```
Raw fastq files and count tables can be found in GEO: <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185948 <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151302


Welcome to our GitHub repository!  
Here you will find analysis scripts for our manuscript where we integrate paired snRNAseq and snATACseq from 8 ADPKD kidneys and 5 healthy adult kidney cortex samples. Please contact the corresponding author, Dr. Benjamin Humphreys, with questions or comments.  
<br/>
Thanks,  
Yoshi

Visit the Humphrey's lab website:   
www.humphreyslab.com  
<br/>
Check out our interactive datasets with Kidney Interactive mulTiomics (KIT):  
http://humphreyslab.com/SingleCell/
<br/><br/>
Find us on Twitter: 
<br/>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @HumphreysLab</a>
<br/><br/>

**preprocessing workflow**  

1. Preprocessing of individual snRNA-seq data

2. Integration of individual snRNA-seq(control, ADPKD and all data)

3. Preprocessing snATAC-seq data

4. Label transfer from snRNA-seq data to snATAC-seq data and filtering nuclei with low-confident prediciton

5. Differential gene expression analysis on snRNA-seq data

6. Differential accessibility / gene activity / motif enrichment analysis among cell types on snATAC-seq data

7. Differential accessibility / gene activity / motif enrichment analysis between ADPKD and control in each cell type on snATAC-seq data
