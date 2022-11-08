import sys
import numpy as np
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
#from lifelines.utils import median_survival_times

PVALUE=0.0

def read_file(file):
	t1=[]
	t2=[]
	c1=[]
	c2=[]
	with open(file, 'r') as f:
		for line in f:
			v = line.split(' ') 
			#print v
			if v[1] == 'YES': #Mutation=group1
				c1.append(int(v[2]))
				t1.append(int(v[3].rstrip()))
			elif v[1] == 'NO': #Not mutation=group 2
				c2.append(int(v[2]))
				t2.append(int(v[3].rstrip()))		
			else:
				print ('Error: not mutation found')
	#print t1,t2,c1,c2
	return t1,t2,c1,c2

def plotKaplan(T, lab):
	kmf=KaplanMeierFitter()
	ax=kmf.fit(T, label=lab).plot(ax=ax)

def kaplanMeier(T,E):
	kmf = KaplanMeierFitter()
	if len(T)!=0 and len(E)!=0:
		kmf.fit(T, event_observed=E)
		return kmf
	else:
		print ('Empty set for Kaplan Meier')
		sys.exit()

if __name__ == '__main__':
	f='TmpgroupsSurv.txt' #sys.argv[1]
	T1,T2,C1,C2 = read_file(f)
	results = logrank_test(T1, T2, event_observed_A=C1, event_observed_B=C2)
	#print(results.__dict__)
	#results.print_summary()
	#print(results.p_value)  
	#print(results.test_statistic)

	kmf1 = kaplanMeier(T1,C1)
	kmf2 = kaplanMeier(T2,C2)
	PVALUE = results.p_value
	if (kmf1.median_survival_time_ > kmf2.median_survival_time_):
		MORESURV="G1"
	else:
		MORESURV="G2" 
	print(results.p_value, results.test_statistic, kmf1.median_survival_time_, kmf2.median_survival_time_, len(T1), len(T2))
