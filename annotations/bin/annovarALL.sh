#!/bin/sh

SET=$1
MAIN_DIR=/home/jairo/epistasis/annotations
FILES_DIR=$MAIN_DIR/$SET/vcf 
TABLE_ANNOVAR_DIR=$MAIN_DIR/$SET/table_annovar
CONVERT_TO_ANNOVAR_DIR=$MAIN_DIR/$SET/annovar
ANNOTATE_VAR_DIR=$MAIN_DIR/$SET/annovarOut 
SOFT_DIR=. 
HUMAN=hg38
# load in SOFT_DIR/humandb  the   hg19 data 
# load in FILES_DIR the .vcf files
#  vcfList is the list of files in FILES_DIR

FILE=$FILES_DIR/$2
        echo 'Table annovar procedure for vcf file: '$FILE
	#cd $SOFT_DIR
	NAME_OUT=`basename $FILE`
	perl table_annovar.pl $FILE $SOFT_DIR/humandb/ -buildver $HUMAN -out $TABLE_ANNOVAR_DIR/$NAME
_OUT -remove -protocol refGene,knownGene,cytoBand,1000g2015aug_all,exac03,avsnp147,dbnsfp30a -operati
on g,g,r,f,f,f,f -nastring . -vcfinput #-thread 15
	echo 'Convert to Annovar procedure for vcf file: '$NAME_OUT
        perl convert2annovar.pl -format vcf4old -include -withzyg -comment $TABLE_ANNOVAR_DIR/$NAME_O
UT.hg38_multianno.vcf -outfile $CONVERT_TO_ANNOVAR_DIR/$NAME_OUT.anninput
#        rm $TABLE_ANNOVAR_DIR/$NAME_OUT.*
	echo 'Annotate Variation Annovar procedure for vcf file: '$NAME_OUT.anninput
        perl annotate_variation.pl -geneanno $CONVERT_TO_ANNOVAR_DIR/$NAME_OUT.anninput -buildver $HU
MAN $SOFT_DIR/humandb/ --outfile $ANNOTATE_VAR_DIR/$NAME_OUT
#        rm $CONVERT_TO_ANNOVAR_DIR/$NAME_OUT.anninput $ANNOTATE_VAR_DIR/$NAME_OUT.variant_function $
ANNOTATE_VAR_DIR/$NAME_OUT.log $ANNOTATE_VAR_DIR/$NAME_OUT.invalid_input
