


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BS_9WPPT
#7YA	BS_DMKS5N4H	BS_MS2F750M
#chr1	10146	rs375931351	AC	A	225.80	PASS	AC=2;AF=0.333;AN
#=6;DB;DP=112;ExcessHet=1.5490;FS=0.000;MLEAC=3;MLEAF=0.750;MQ=24.10;PG=0,0,0;QD=
#18.82;SOR=1.445;VQSLOD=1.73;culprit=FS;ANN=-|upstream_gene_variant|MODIFIER|DDX1
#1L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogen
#e|||||||||||1863|1||deletion|HGNC|HGNC:37102||||chr1:g.10150del,-|upstream_gene_
#variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_tr
#anscript|||||||||||1722|1||deletion|HGNC|HGNC:37102|YES|||chr1:g.10150del,-|down
#stream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|u
#nprocessed_pseudogene|||||||||||4257|-1||deletion|HGNC|HGNC:38034|YES|||chr1:g.1
#0150del	GT:AD:DP:FT:GQ:PGT:PID:PL:PP	#0/0:26,0:26:lowGQ:0:.:.:0,0,0:0,0,0	
#0/0:64,0:64:lowGQ:0:.:.:0,0,864:0,0,864	1/1:0,12:12:PASS:32:1|1:10144_T_G:237,32
#,0:237,32,0



#!/usr/bin/python
import sys, os, subprocess, string

gzip = 'zcat'



def global_vars():
	# Annovar ouput positions
	global nchr, nid, npos, nref, nalt, npass, ndat, nger, nfam1, nfam2, ninfo
	nchr = 0
	nid = 2
	npos = 1
	nref = 3
	nalt = 4
	npass = 6
	ninfo = 7
	ndat = 8
	nger = 9
	nfam1 = 10
        nfam2 = 11
	global ptype, info_ann, info_start, info_stop
	#ptype = 1                                        
	info_ann = 2
	info_start = 0
	info_stop = 12
	global ann_chr, ann_id, ann_pos, ann_ref, ann_alt, ann_pass, ann_tdat, ann_ger, ann_fam1, ann_fam2, ann_info
	ann_chr = info_start + nchr
	ann_id = info_start + nid
	ann_pos = info_start + npos
	ann_ref = info_start + nref
	ann_alt = info_start + nalt
	ann_pass = info_start + npass
	# Important position in the vcf with information
	ann_tdat = info_start + ndat
	ann_ger = info_start +  nger
	ann_fam1 = info_start + nfam1
        ann_fam2 = info_start + nfam2
	ann_info = info_start + ninfo
	return



def get_mapping_data(data_types, data_values):
	dict = {}
	vtype = data_types.split(':')
	if 'GT' not in vtype:
		print 'ERROR: Incorrect FORMAT position', ann_tdat + 1, 'in VCF'
		exit()
	vdata = data_values.split(':')
	if True or (len(vtype) == len(vdata)):
		for i in range(min(len(vtype),len(vdata))):
			if vtype[i] == 'GT':
				if vdata[i].find('/') >= 0:
					print i, vtype, vdata
					if vdata[i].find('0/0') >= 0: dict[vtype[i]] = 0
					elif vdata[i].find('0') >= 0: dict[vtype[i]] = 1
					else: dict[vtype[i]] = 2
#				elif vdata[i].find('|') >= 0:
#					dict[vtype[i]] = vdata[i].split('|')
				else:
					dict[vtype[i]] = -1 # [vdata[i]]				
			elif vtype[i] == 'AD':
				vnum = []
				v = vdata[i].split(',')
				for num in v:
					try:
						vnum.append(float(num))
					except:
						vnum.append(0.0)
				dict[vtype[i]] = vnum
			elif vtype[i] == 'BQ':   # Not used in this type of files
                                dict[vtype[i]] = vdata[i].split(',')[0]
				
			else:
				try:
					dict[vtype[i]] = float(vdata[i])
				except:
					dict[vtype[i]] = 0.0
	return dict

def get_frequency(data_info):
        info = data_info.split('AF=')
        if (len(info) == 1): return (False)
        info2=info[1].split(";") 
        #info3=info2[0].split(",")
	return(info2[0])     #Allel Frequencies  .  there are several allels


def get_gene(data_info): # for VEP files
        info = data_info.split('ANN=')
        if (len(info) == 1): return(False)
        info2=info[1].split("|")        
        return (info2[1],info2[2],info2[3],info2[7])  #variant_type,consec,gene,genetype


def get_annovar(path, filename0, fileMatrix, fpass=True, maf=1.0, filtMut=True, minAD=0, minRatioAD_DP=0):
	   filename=path+'/'+filename0
	   fileout=filename+'_clean'
	   os.system("gunzip "+filename+".gz")	
	   with open(filename, 'r') as f:
		l = []
		for line in f:
			if line[0] == '#': continue
			v = line.split('\t')
			if len(v) > ann_ger:
				chr = v[ann_chr]
				#rs = v[ann_id]
				rs = '.'
				ipos = v[ann_pos]
                                aref=v[ann_ref]
                                aalt=v[ann_alt]
                                #snp=(chr,rs,ipos,aref,aalt)
				aa = v[ann_ref] + '|' + v[ann_alt]
				#first_data = v[info_ann].split(',')[0]  # to be fix
				#data = first_data.split(':')
				#print data,len(data)
				#gene = data[0]   # for COAD  
				G = get_gene(v[ann_info])
				if G==False:
					print "No ANN info in VEP:",filename,line
				else: (variant_type,consec,gene,genetype) = G
				#if (len(data) >= 5): mut = data[4]
				filterPass = v[ann_pass]
				dic_data = get_mapping_data(v[ann_tdat], v[ann_ger])
				if dic_data['GT'] == -1: print "Error: No '/' en GT ", filename, line, v[ann_tdat], v[ann_ger]

				if len(v) > ann_fam1:
					dic_data_fam1=get_mapping_data(v[ann_tdat], v[ann_fam1])
					if dic_data_fam1['GT'] == -1: print "Error: No '/' en GT", filename, line
				else: 
					dic_data_fam1 = {}
					print "A file with no parents", filename, line
				if len(v) > ann_fam2:
					dic_data_fam2=get_mapping_data(v[ann_tdat], v[ann_fam2])
					if dic_data_fam2['GT'] == -1: print "Error: No '/' en GT", filename, line
				else: 
					dic_data_fam2 = {}
					print "A file with only one parent", filename, line
				if (len(v) > (ann_fam2 + 1)):
					print "A file with more than two parents", filename, line
			else:
				print "A file with no afect", filename, line

				
			#Filter for chrm X
			#if chr == 'X' or chr == 'chrX': continue 
			#print "Filtered by chrm ", chr
			
			#Get frequency data from Annovar file: gnomAD_genome_ALL, ExAC_ALL, 1000g2015aug_all
                        info_freq = get_frequency(v[ann_info])
			if info_freq == False: print "Error: No AF:", filename, line

                        """
			#print dic_info_freq
                        if not bool(dic_info_freq):
                                #print 'Frequency dictionary is empty'
                                freqMin = 0
                                freqMax = 1.0
                        else:
                                freqMin = dic_info_freq
                        	freqMax = dic_info_freq
			#print "Min freq is: ", freq
                        freqMaxstr=str(freqMax)
			"""

			
			if fpass == True and filterPass != 'PASS': continue     #  PASS is needed ??
			#print "Ha pasado los filtros bq, fdp y PASS"

			"""
			#Filter for AF
                        if (float(freqMax) > maf):
                                #print 'Variant ' + str(pos) + ' filtered by freq: ' + str(freqMax) + ' in filename: ' + filename
                                continue
			"""
			#Get the Genotype:
			try :
				iaa = dic_data['GT']
			except:
				iaa = False
				print "No GT :", filename, line
			iaa1 = dic_data_fam1.get('GT',False)
			iaa2 = dic_data_fam2.get('GT',False)
			"""
			#Get the DP and AD for normal and tumor
			fdp = dic_data.get('DP', 0)
                        fdpr = dic_data_ref.get('DP', 0)

                        ad_alt = dic_data.get('AD',[0])[0]
                        ad_ref = dic_data.get('RD',0)
			#print 'RD= ' + str(ad_ref) + ' and AD= ' + str(ad_alt)
			fdpnor = ad_alt + ad_ref
			if (fdpnor==0):
                                #print 'AD,RD '+str(ad_alt) + ','+str(ad_ref)+':'+filename
                                #snp=(chr,rs,ipos,aref,aalt,freqMaxstr)
                                #print snp
                                ratio_ad=0
                        else :
                                ratio_ad = ad_alt / fdpnor

                        adr_alt = dic_data_ref.get('AD',[0])[0]
                        adr_ref = dic_data_ref.get('RD',0)
			fdptum = adr_alt + adr_ref
			if (fdptum==0):
                                #print 'ADt,RDt '+str(adr_alt) + ','+str(adr_ref)+':'+filename
                                #snp=(chr,rs,ipos,aref,aalt,freqMaxstr)
                                #print snp
                                ratio_adr=0
                        else :
                                ratio_adr = adr_alt / fdptum

			#print 'minAD= ' + str(minAD) + ' and minRatioAD_DP= ' + str(minRatioAD_DP)
			if (False and minAD != 0 and minRatioAD_DP != 0):
                        	#Filter if number of reads<10 if GT != 0/0
				if ((fdpnor < minAD and iaa != ['0', '0']) or (fdptum < minAD and iaar != ['0', '0'])):  
					#print 'Variant ' + str(ipos) + ' filtered by DP: ' + str(fdpnor) + '/' + str(fdptum) + ' in filename: ' + filename
					continue 
    			
				#Filter if the ratio of alternate reads < 0.05 if GT != 0/0
				if ((ratio_ad < minRatioAD_DP and iaa != ['0', '0']) or (ratio_adr < minRatioAD_DP and iaar != ['0', '0'])): 
					#print 'Variant ' + str(ipos) + ' filtered by AD: ' + str(ratio_ad) + '/' + str(ratio_adr) + ' in filename: ' + filename
					continue
			"""

			if (gene == "") or (gene.upper()=='UNKNOWN'): continue

			# Include all type of deleterious mutations
			if filtMut:
                                if variant_type == 'synonymous_variant': continue
                                if genetype.find('pseudogene') >= 0: continue

#				if v[ptype] != 'nonsynonymous SNV' and v[ptype] != 'nonframeshift insertion' and v[ptype] != 'nonframeshift deletion' and v[ptype] != 'frameshift insertion' and v[ptype] != 'frameshift deletion' and v[ptype] != 'stopgain' and v[ptype] != 'stoploss': continue 
			mut=variant_type	
			#print chr,ipos,aref,aalt,gene,'_'.join(v[ptype].split()),mut,freqMaxstr,''.join(iaa),ad_ref,ad_alt,''.join(iaar),adr_ref,adr_alt

			"""
			snp=(chr,ipos,aref,aalt,info_freq,gene)
			if dsnps.get(snp,[])==[]: dsnps[snp]=(totsamples+1)*[0] # Dictionary length = number of samples per normals and tumors
			dsnps[snp][j] = iaa
			dsnps[snp][nsamples+j] = iaa1
                        dsnps[snp][2*nsamples+j] = iaa2
			"""
			l.append([chr,ipos,aref,aalt,info_freq,gene,iaa,iaa1,iaa2])
	   try:
              fout=open(fileout, 'w')
           except: 
              print "Can't open ",fileout;
	   for row in l:
		fout.write('\t'.join([str(j) for j in row])+'\n')
	   fout.close()
	   #os.system("gzip " + filename)	
	   """
	print fileMatrix, j
	with open(fileMatrix, 'w') as fout:
	  for snp in dsnps.keys():
	   
           cont=0
           for j in range(nsamples):
             if  dsnps[snp][j] > 0 : cont=cont+1
           if cont > maxmut : continue
           cont=0
           for j in range(nsamples):
             if  dsnps[snp][nsamples + j] > 0 : cont=cont+1
           if cont < minmut  : continue
	   
           fout.write(string.join(snp,':')+'\t')   #+':'+dsnps[snp][totsamples])   # gene name
           for j in range(totsamples):              
            fout.write(" "+str(dsnps[snp][j]))
           fout.write("\n")
	   """
	   return 


def get_options():
	'''
    parse option from command line call
    '''
	import optparse

	desc = 'Script for filtering variants from annovar output'
	parser = optparse.OptionParser('usage: [-h] [-n normal] [-o outfile]', description=desc)
	parser.add_option('-n', dest='normpos', help='Position for checking germline mutation')
	parser.add_option('-o', dest='outfile', help='Output file')
	parser.add_option('--nsubjects', action='store', type='float', dest='nsamples', help='Number of Subjects')
	parser.add_option('--no-germline', action='store_true', dest='rgerm', help='Filtering out germline mutations')
	parser.add_option('--both', action='store_true', dest='both', help='Filtering both positions')
	parser.add_option('-f', '--format', action='store', type='int', dest='pos_format',
					  help='Position of FORMAT information in VCF file')
	parser.add_option('-b', '--bq', action='store', type='float', dest='bq', help='Base Quality Filter')
	parser.add_option('--fa', '--freq', action='store', type='float', dest='fa', help='Allele frequency Filter')
	parser.add_option('-d', '--depth', action='store', type='float', dest='dp', help='Depth Filter')
	parser.add_option('--no-pass', action='store_true', dest='npass', help='Remove PASS filtering')
	parser.add_option('--high', action='store_true', dest='hfvar', help='Select highly frequent variants')
	parser.add_option('--norm', action='store_true', dest='norm', help='Filter for normal positions')
	parser.add_option('--maf', '--maf_freq', action='store', type='float', dest='maf', help='MAF Minimum Allele frequency Filter')
	parser.add_option('--filtMut', action='store_true', dest='filtMut', help='Filter Synonymous variants')
	parser.add_option('--minAD', '--min_AD', action='store', type='float', dest='minAD', help='minAD Minimum Allele Depth')
	parser.add_option('--minRatioAD_DP', '--min_Ratio_AD_DP', action='store', type='float', dest='minRatioAD_DP', help='minRatioAD_DP Minimum Allele Ratio AD_DP')


	(options, args) = parser.parse_args()

	if not os.path.isdir(args[0]):
			print 'ERROR:', args[0], 'input file not found'
			exit()

	rgerm = False
	both = False
	if options.normpos:
		try:
			normpos = int(options.normpos) - 1
			if options.rgerm: rgerm = True
			if options.both: both = True
		except:
			normpos = -1
	else:
		normpos = -1

	if options.norm:
		normFilt = True
	else:
		normFilt = False

	if options.npass:
		npass = False
	else:
		npass = True

	if options.bq:
		try:
			tbq = float(options.bq)
		except:
			tbq = 0
	else:
		tbq = 0

	if options.fa:
		try:
			tfa = float(options.fa)
		except:
			tfa = 0.0
	else:
		tfa = 0.0

	if options.dp:
		try:
			tdp = int(options.dp)
		except:
			tdp = 0
	else:
		tdp = 0

	if options.pos_format:
		global ann_tdat
		ann_tdat = int(options.pos_format) - 1

	if options.hfvar:
		low = False
	else:
		low = True

	if options.filtMut:
                filtMut = True
        else:
                filtMut = False

	if options.maf:
                try:
                        maf = float(options.maf)
                except:
                        maf = 1.0
        else:
                maf = 1.0

	if options.minAD:
                try:
                        minAD = float(options.minAD)
                except:
                        minAD = 0.0
        else:
                minAD = 0.0

	if options.minRatioAD_DP:
                try:
                        minRatioAD_DP = float(options.minRatioAD_DP)
                except:
                        minRatioAD_DP = 0.0
        else:
                minRatioAD_DP = 0.0

	if options.outfile:
		outfile = options.outfile
	else:
		outfile = None

 	if options.nsamples:
                try:
                        nsamples = int(options.nsamples)
                except:
                       nsamples = 0
        else:
                nsamples = 0

	opts = (nsamples, outfile, rgerm, tdp, tfa, tbq, npass, low, normFilt, maf, filtMut, minAD, minRatioAD_DP)

	return args, opts


if __name__ == '__main__':
	global_vars()
	args, opts = get_options()
	if (len(sys.argv) > 1):
           path=sys.argv[1]
	   filename=sys.argv[2]
	else:
		print 'Incorrect input'
		print 'python FilterAnnovarFiles.py path filename  '
		exit()

	nsamples,fileout, rgerm, tdp, tfa, tbq, fpass, low, normFilt, maf, filtMut, minAD, minRatioAD_DP = opts


	get_annovar(path, filename, fileout, fpass, maf, filtMut, minAD, minRatioAD_DP)



