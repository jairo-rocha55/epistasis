#!/bin/sh
#SBATCH --time=1

FILES_DIR=/home/jaume/LUAD38/varscan2/
TABLE_ANNOVAR_DIR=/home/jaume/LUAD38/varscan2/table_annovar
CONVERT_TO_ANNOVAR_DIR=/home/jaume/annovarAvinputs/LUAD38/varscan2/emidio
ANNOTATE_VAR_DIR=/home/jaume/annovarAnnotations/LUAD38/varscan2/table_annovar_emidio
SOFT_DIR=/home/jaume/software/annovar

# For each file in this directory
for FILE in $FILES_DIR/*.vcf; do
        echo 'Table annovar procedure for vcf file: '$FILE
	cd $SOFT_DIR
	NAME_OUT=`basename $FILE`
	perl table_annovar.pl $FILE humandb/ -buildver hg38 -out $TABLE_ANNOVAR_DIR/$NAME_OUT -remove -protocol refGene,knownGene,cytoBand,1000g2015aug_all,exac03,avsnp147,dbnsfp30a -operation g,g,r,f,f,f,f -nastring . -vcfinput -thread 15
	echo 'Convert to Annovar procedure for vcf file: '$NAME_OUT
        perl convert2annovar.pl -format vcf4old -include -withzyg -comment $TABLE_ANNOVAR_DIR/$NAME_OUT.hg19_multianno.vcf -outfile $CONVERT_TO_ANNOVAR_DIR/$NAME_OUT.avinput	
echo 'Annotate Variation Annovar procedure for vcf file: '$NAME_OUT.avinput
        perl annotate_variation.pl -geneanno $CONVERT_TO_ANNOVAR_DIR/$NAME_OUT.avinput -buildver hg38 humandb/ --outfile $ANNOTATE_VAR_DIR/$NAME_OUT
done
echo 'done'
