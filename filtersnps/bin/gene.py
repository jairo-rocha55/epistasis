import sys
import numpy as np


def get_dic(filename,s1,s2,l):
	d={}
	f=open(filename)
	for line in f:
		v=line.rstrip().split()
		g=v[0].split(':')[-1]
		d[g]=d.get(g,[np.zeros(l,dtype=int),np.zeros(l,dtype=int)])
		v1=[int(v[i]) for i in range(s1,s1+l)]
		v2=[int(v[i]) for i in range(s2,s2+l)]
		for i in range(l):
			d[g][0][i]=max(d[g][0][i],v1[i])
			d[g][1][i]=max(d[g][1][i],v2[i])
	return d


if __name__ == '__main__':
	filename=sys.argv[1]
	size=int(sys.argv[2])
	minmuttumor=int(sys.argv[3])
	maxmutnormal=int(sys.argv[4])
	d=get_dic(filename,1,size+1,size)
	for k in d.keys():
           if ((sum(d[k][0]) < maxmutnormal*size/100.0)  and (sum(d[k][1]) > minmuttumor*size/100.0)):
		print k,' '.join([str(i) for i in d[k][0]]),' '.join([str(i) for i in d[k][1]])
	
