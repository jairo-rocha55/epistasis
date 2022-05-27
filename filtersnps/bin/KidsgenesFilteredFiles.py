#!/usr/bin/python
import sys, os, string

def get_annovar(path, filelist, nsamples, fileout, minmut=0,maxmut=100,  maf=1.0):
        ad=0 
	ratio=0.0
	dsnps={}
        totsamples=3*nsamples
	print 'Nsamples: ' , nsamples, 'totsamples', totsamples

	j=0
	with open(filelist, 'r') as list_file:
         for lineList in list_file:
	   filename=path+'/'+lineList.rstrip()+'_clean'
	   #print filename
	   with open(filename, 'r') as f:
		for line in f:
			v = line.rstrip().split('\t')
#   0	    1	   2    3       4       5       6       7       8

#chr1	54712	TTTTC	T	0.500	OR4G4P	1	2	0
#chr1	54714	TTCTTTCTTTC	*,T	0.667,0.167	OR4G4P	2	2	1
#chr1	54720	CTTTCT	C	0.167	OR4G4P	0	0	1
			#print v
			chr = v[0]
			#rs = v[]
			pos = v[1]
                        ref = v[2]
                        alt = v[3]
			gene = v[5]  #.split('/')[0]
			type_mut = v[5]
			mut = v[5]
			freq = v[4] #Max(gnomAD_genome_ALL, ExAC_ALL, 1000g2015aug_all)
			iaa = v[6]
			#rd_n = float(v[9])
			#ad_n = float(v[10])
			iaa1 = v[7]
			if (len(v) == 9):
				iaa2 = v[8]
			else:   iaa2 = "-1"
			
			snp=gene #(chr,pos,ref,alt)
	
			#Filter for chrm X
			#if chr == 'X' or chr == 'chrX': continue 
			#print "Filtered by chrm ", chr
			#else: print "Chr is ",chr
			
 			if dsnps.get(snp,[])==[]: dsnps[snp]=(totsamples+1)*[0] # Dictionary length = number of samples per afecto and parents 
			dsnps[snp][j] = max(iaa,dsnps[snp][j])
			dsnps[snp][nsamples+j] = max(iaa1,dsnps[snp][nsamples+j])
                        dsnps[snp][2*nsamples+j] = max(iaa2,dsnps[snp][2*nsamples+j])
			#dsnps[snp][totsamples]=gene
	   j+=1
#	minmut=minmut*nsamples/100
#       maxmut=maxmut*nsamples/100
        with open(fileout, 'w') as fout:
	  for snp in dsnps.keys():
#           cont=0
#           for j in range(nsamples):
#             if  dsnps[snp][j] > 0 : cont=cont+1
#           if cont > maxmut : continue
#           cont=0
#           for j in range(nsamples):
#             if  dsnps[snp][nsamples + j] > 0 : cont=cont+1
#           if cont < minmut  : continue

           fout.write(snp) ###(string.join(snp,':')+':'+dsnps[snp][totsamples])   # gene name
           for j in range(totsamples):              
            fout.write(" "+str(dsnps[snp][j]))
           fout.write("\n")
		
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
