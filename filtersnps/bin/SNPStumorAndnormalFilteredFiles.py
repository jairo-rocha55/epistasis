#!/usr/bin/python
import sys, os, string
import numpy as np



def get_annovar(path, filelist, nsamples, fileout, minmut=0,maxmut=100,  maf=1.0):
	global mafcont, maxmafcont, mingenmutcont
	vmuts = []
	ad=0 
	ratio=0.0
	print minmut, maxmut, maf
	dsnps={}  #dictionary of SNPs
	dsubj={}  #dictionary of subjects
	type_file=""
	totsamples=2*nsamples
	print 'Type_file: ', type_file
	print 'Nsamples: ' , nsamples, 'totsamples', totsamples
	k=0  # number of real subjects , each subject has several samples. 
	mafcont=0  # number of SNPs filtered by maf
	maxmafcont=0 # number of SNPS filtered for far away from maf
	mingenmutcont=0 # number of genes filtered by mingenmut
	with open(filelist, 'r') as list_file:
         for lineList in list_file:
# filelist has three columns: the third column is the subject ID and the first one, the filename
	   NameList=lineList.split()
	   filename=path+'/'+NameList[0]+'.exonic_variant_function'+'_clean'
	   #print filename
	   with open(filename, 'r') as f:
		subjind=dsubj.get(NameList[2],[])
		
		if subjind==[]: 
			dsubj[NameList[2]]=k
			j=k
			k+=1
		else:
			j=subjind
		for line in f:
			v = line.split('\t')
			#print v
			chr = v[0]
			#rs = v[]
			pos = v[1]
                        ref = v[2]
                        alt = v[3]
			gene = v[4].split('/')[0]
			type_mut = v[5]
			filterPASS = v[6]
			freq = v[7] #Max(gnomAD_genome_ALL, ExAC_ALL, 1000g2015aug_all)
			gt_n = v[8]
			rd_n = float(v[9])
			ad_n = float(v[10])
			gt_t = v[11]
                        rd_t = float(v[12])
                        ad_t = float(v[13])
			
			snp=(chr,pos,ref,alt,type_mut) #,freq) #,mut)
	
			#Filter for chrm X
			#if chr == 'X' or chr == 'chrX': continue 
			#print "Filtered by chrm ", chr
			#else: print "Chr is ",chr
			
                        #Filter if number of reads<10 if GT != 0/0
			if (filterPASS!='PASS' and (((rd_n+ad_n) < ad and gt_n != '00') or ((rd_t+ad_t) < ad and gt_t != '00'))):  
				#print 'Variant ' + str(pos) + ' filtered by DP: ' + str(rd_n+ad_n) + '/' + str(rd_t+ad_t) + ' in filename: ' + filename
				continue 
    			
			#Filter if the ratio of alternate reads < 0.05 if GT != 0/0
			if (filterPASS!='PASS' and (( gt_n != '00' and (rd_n+ad_n>0) and(ad_n/(rd_n+ad_n)) < ratio ) or (gt_t != '00' and (rd_t+ad_t>0) and (ad_t/(rd_t+ad_t)) < ratio ))): 
				#print 'Variant ' + str(pos) + ' filtered by AD: ' + str((ad_n/(rd_n+ad_n))) + '/' + str((ad_t/(rd_t+ad_t))) + ' in filename: ' + filename
				continue

			#Filter for AF
			if (float(freq) > maf): 
				if (mafcont==0) : print 'Variant ' + str(pos) + ' filtered by freq: ' + str(freq) + ' in filename: ' + filename
				mafcont=mafcont + 1 
				continue

			# Include all type of deleterious mutations
			if type_mut != 'nonsynonymous_SNV' and type_mut != 'nonframeshift_insertion' and type_mut != 'nonframeshift_deletion' and type_mut != 'frameshift_insertion' and type_mut != 'frameshift_deletion' and type_mut != 'stopgain' and type_mut != 'stoploss': 
				#print 'Variant ' + str(pos) + ' filtered by type of mutation: ' + str(type_mut) + ' in filename: ' + filename
				continue
			
			if dsnps.get(snp,[])==[]: dsnps[snp]=(totsamples+2)*[0] # Dictionary length = number of samples per normals and tumors
                        if gt_n == '00': 
                        	dsnps[snp][j] = 0
                        elif gt_n == '11' or gt_n == '12' or gt_n == '21' or gt_n == '22' or gt_n == '13' or gt_n == '23' or gt_n== '112': 
                        	dsnps[snp][j] = 2 # for normal
                       	elif gt_n == '01' or gt_n == '10' or gt_n == '02' or gt_n == '03' or gt_n == '20' or gt_n == '30'or gt_n == '012':
                        	dsnps[snp][j] = 1
			else:
                        	print 'Case iaa: ', str(gt_n), ' not considered'
                        
			#For normal and tumor 
			if type_file == "" :
                        	if gt_t == '00':
                                	dsnps[snp][nsamples+j] = 0
                        	elif gt_t == '11' or gt_t == '12' or gt_t == '21' or gt_t == '22' or gt_t == '13' or gt_t == '23' or gt_t== '112' :
                                	dsnps[snp][nsamples+j] = 2 # for normal
                       		elif gt_t == '01' or gt_t == '10' or gt_t == '02' or gt_t == '03' or gt_t == '20' or gt_t == '30' or gt_t == '012':
                                	dsnps[snp][nsamples+j] = 1
                        	else:
                                	print 'Case iaa: ', str(gt_t), ' not considered'	
			dsnps[snp][totsamples]=gene
			dsnps[snp][totsamples+1]=float(freq)
	minmut=minmut*nsamples/100
	maxmut=maxmut*nsamples/100
        with open(fileout, 'w') as fout:
	  for snp in dsnps.keys():
           cont=0
           for j in range(k):
             if  dsnps[snp][j] > 0 : cont=cont+1
           freq0=dsnps[snp][totsamples+1]#  MAF	
           if (cont > maxmut) or (cont > ((0.01 + 1.4*freq0)*nsamples)) : continue
           cont=0
           for j in range(k):
             if  dsnps[snp][nsamples + j] > 0 : cont=cont+1
           if cont < minmut  : continue

           fout.write(string.join(snp,':')+':'+str(freq0)+':'+dsnps[snp][totsamples])   # gene name
           for j in range(k):              
            fout.write(" "+str(dsnps[snp][j]))
           for j in range(nsamples,nsamples+k):              
            fout.write(" "+str(dsnps[snp][j]))
           fout.write("\n")
	fout.close()

	#  Group info by genes
	d={}
	maxmafcont=0
	for snp in dsnps.keys():
           cont=0
           for j in range(k):
             if  dsnps[snp][j] > 0 : cont=cont+1
           freq0=dsnps[snp][totsamples+1]	# MAF
           g=dsnps[snp][totsamples]    # gene name	
	   maxmaf = ((0.01 + 1.4*freq0)*nsamples)
           if (cont > maxmut) or (cont > maxmaf)  : 
		if (maxmafcont==0) :
			print string.join(snp,':')+':'+str(freq0)+':'+dsnps[snp][totsamples]
		maxmafcont = maxmafcont + 1 
		continue
           cont=0
           for j in range(k):
             if  dsnps[snp][nsamples + j] > 0 : cont=cont+1
           if cont < minmut  : continue

           d[g]=d.get(g,[np.zeros(nsamples,dtype=int),np.zeros(nsamples,dtype=int)])
           v1=[dsnps[snp][j] for j in range(nsamples)]
           v2=[dsnps[snp][j] for j in range(nsamples,totsamples)]
           for i in range(nsamples):
                    d[g][0][i]=max(d[g][0][i],v1[i])
                    d[g][1][i]=max(d[g][1][i],v2[i])




	print "NumSubjects=",k		

	# print gene matrix
	size=k #nsamples
	maxmutnormal=100
	minmuttumor=5
	with open(fileout+"Gene", 'w') as fout:
	 for k in d.keys():
           if ((sum(d[k][0]) < maxmutnormal*size/100.0)  and (sum(d[k][1]) > minmuttumor*size/100.0)):
		fout.write(str(k) + " " + ' '.join([str(i) for i in d[k][0]]) + " "+ ' '.join([str(i) for i in d[k][1]])+"\n")
           else: 
		if (mingenmutcont == 0) :
			print str(k)
		mingenmutcont = mingenmutcont + 1		
	fout.close()
	return 


def get_options():
	'''
    parse option from command line call
    '''
	import optparse

	desc = 'Script for creating matrix from filtered annovar files'
	parser = optparse.OptionParser('usage: [--maf AF] [--minAD NUM] [--ratio_AD NUM] -o outfile -nsubjects NUM', description=desc)
	parser.add_option('-o', dest='outfile', help='Output file')
	parser.add_option('--nsubjects', action='store', type='float', dest='nsamples', help='Number of Subjects')
        parser.add_option('--minPercMutTumor',  action='store', type='float', dest='MinPercMutTumor', help='Minimun percentage of subjects mutated in tumor')
        parser.add_option('--maxPercMutNormal', action='store', type='float', dest='MaxPercMutNormal', help='Maximum percentage of subjects mutated in normal')
	#parser.add_option('--minAD', '--ad_depth', action='store', type='float', dest='ad', help='Allele Depth Filter')
	#parser.add_option('--file', action='store', dest='type_file', help='Type of VCF file: LUAD, COAD, CHD, 1000g')
	#parser.add_option('--minRatioAD_DP', '--ratio_AD', action='store', type='float', dest='ratio', help='Ratio between AD and AD+RD')
	parser.add_option('--maf', '--maf_freq', action='store', type='float', dest='maf', help='MAF Minimum Allele frequency Filter')

	(options, args) = parser.parse_args()

	if (len(args) < 2):
		print 'Incorrect input'
		print 'python SNPStumorAndnormalFilteredFiles.py path fileList -o outfile [filters]'
		exit()
	else:
		if not os.path.isdir(args[0]):
			print 'ERROR:', args[0], 'input file not found'
			exit()

	if options.outfile:
		outfile = options.outfile
	else:
		outfile = None

 	if options.nsamples:
                try:
                        nsamples = float(options.nsamples)
                except:
                       nsamples = 0.0
        else:
                nsamples = 0.0
	if options.maf:
                try:
                        maf = float(options.maf)
                except:
                        maf = 1.0
        else:
                maf = 1.0
        if options.MinPercMutTumor:
                try:
                        minmut = float(options.MinPercMutTumor)
                except:
                        minmut = 0.0
        else:
                minmut = 0.0

        if options.MaxPercMutNormal:
                try:
                        maxmut = float(options.MaxPercMutNormal)
                except:
                        maxmut = 100
        else:
                maxmut = 100

        ad=0
        ratio=0
	opts = (outfile, minmut, maxmut, int(nsamples), maf)

	return args, opts
'''
        #global type_file
	if options.type_file:
                try:
                	type_file = options.type_file
                except:
			type_file = ""
        else:
                type_file = ""	
	if options.minRatioAD_DP:
                try:
                        ratio = float(options.ratio)
                except:
                        ratio = 0.0
        else:
                ratio = 0.0

	if options.minAD:
                try:
                        ad = float(options.ad)
                except:
                        ad = 0.0
        else:
                ad = 0.0

'''

if __name__ == '__main__':
	args, opts = get_options()
	if (len(sys.argv) > 1):
           path=sys.argv[1]
	   filelist=sys.argv[2]

	outfile, minmut, maxmut, nsamples, maf = opts

	get_annovar(path, filelist, nsamples, outfile, minmut, maxmut, maf)

        print mafcont,': number of SNPs filtered by maf'
        print maxmafcont, ': number of SNPS filtered for far away from maf'
        print mingenmutcont, ': number of genes filtered by mingenmut'

