#!/usr/bin/python
import sys
import pylab as plt


def read_survdata(filename,ptime,pcens,pgs=None):
	dgs={}
	lines=open(filename).readlines()
	data=[i.rstrip().split() for i in lines]
	if pgs:
		for i in data:
			dgs[i[pgs]]=dgs.get(i[pgs],[])
			dgs[i[pgs]].append([int(i[ptime]),int(i[pcens])])
	else:
		dgs['ALL']=[]
		for i in data:	
			dgs['ALL'].append([int(i[ptime]),int(i[pcens])])
	for k in dgs.keys():
		dgs[k].sort(key=lambda x: (x[0],x[1]))
	return dgs


def kaplan(survdata):
	h_coords=[]
	v_coords=[]
	lost=[]
	x=0
	y=1
	for i in survdata:
		if i[1]!=1:
			lost.append([i[0],y])
		else:
			h_coords.append([i[0],y])
			y=y*len(survdata[survdata.index(i)+1:])/float(len(survdata[survdata.index(i):
]))
			v_coords.append([i[0],h_coords[-1][-1],y])
			break
	newsurvdata=survdata[survdata.index(i)+1:]
	l=0
	while len(newsurvdata)>0:
		l=l+1
		newsurvdata,y,h_coords,v_coords,lost=survloop(newsurvdata,y,h_coords,v_coords,lost)
	return h_coords,v_coords,lost


def survloop(newsurvdata,y,h_coords,v_coords,lost):
	for j in newsurvdata:
		if j[1]!=1:
			lost.append([j[0],y])
		else:
			h_coords.append([j[0],y])
			y=y*len(newsurvdata[newsurvdata.index(j)+1:])/float(len(newsurvdata[newsurvdata.index(j):]))
			v_coords.append([j[0],h_coords[-1][-1],y])
			break
	newsurvdata=newsurvdata[newsurvdata.index(j)+1:]
	return newsurvdata,y,h_coords,v_coords,lost


def plot_kaplan(dks,d=0.01,width=3):
	colors=['b','r','b','c','m','g','y']
	c=0
	fig=plt.figure()
	ax=fig.add_subplot(111)
	fig.subplots_adjust(bottom=0.1)
	fig.subplots_adjust(top=0.95)
	fig.subplots_adjust(left=0.1)
	fig.subplots_adjust(right=0.75)
	for k in dks.keys():
		start=0
		kaplan=dks[k]
		for i in kaplan[0]:
			ax.hlines(i[1],start,i[0],linewidths=width,color=colors[c])
			start=i[0]
		if kaplan[-1][-1][0]>kaplan[0][-1][0]:
			ax.hlines(kaplan[-1][-1][1],kaplan[0][-1][0],kaplan[-1][-1][0],linewidths=width,color=colors[c])
		for i in kaplan[1]:
			ax.vlines(i[0],i[2]-(width*1.7/1000),i[1]+(width*1.7/1000),linewidths=width,color=colors[c])
		for i in kaplan[2]:
			ax.vlines(i[0],i[1]-d,i[1]+d,color=colors[c])
		c=c+1

	#ax.tick_params(axis='x',length=15,width=3,direction='out',labelsize=30)
	#ax.tick_params(axis='y',length=15,width=3,direction='out',labelsize=30)
	#ax.spines['bottom'].set_position(['outward',10])
	#ax.spines['left'].set_position(['outward',10])
	ax.spines['left'].set_linewidth(width)
	ax.spines['right'].set_linewidth(width)
	ax.spines['top'].set_linewidth(width)
	ax.spines['bottom'].set_linewidth(width)
	#ax.spines['left'].set_bounds(0,1)
	plt.xlim(0,)
	plt.ylim(0,1)
	plt.savefig('testfif.pdf')
	plt.show()


if __name__ == '__main__':
	filename=sys.argv[1]
	ptime=int(sys.argv[2])-1
	pcens=int(sys.argv[3])-1
	pgs=None
	if len(sys.argv)>4:
		pgs=int(sys.argv[4])-1
	dgs=read_survdata(filename,ptime,pcens,pgs)
	dks={}
	for k in dgs.keys():
		dks[k]=kaplan(dgs[k])
		#print dks[k]
	plot_kaplan(dks)
