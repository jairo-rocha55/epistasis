
This script produces gene pairs related to cancer tumor. 
Download the .vcf files in the annotations/COAD/vcf directories.
Vcf files must be normal tumor mutation pairs.

nohup ./Script

The use of nohup is recommended so that the process do not die even when the console gets disconnected.


Input files are VCF files from TCGA to be processed with Annovar. 
 First with table_annovar procedure
 Then with convert2annovar
 Finally with annotateVariation

If we want to filter:
python SNPFilterAnnovarFiles.py table_annovar_dir files_list.txt 19 -n 20 --maf 0.005 --filtMut --minAD 10 --minRatioAD_DP 0.05

Create SNP x subjects matrix 
python SNPStumorAndnormalFilteredFiles.py directory_of_inputFiles files_list -o outFile --file type_of_file --maf 0.005 --minAD 10 --minRatioAD_DP 0.05

#It also creates a gene matrix from a SNP matrix

After the three special cells in the contingency tables of all gene pairs are analyzed, the output is in 3cell/dataCOAD/GPairsPvaluesTop
