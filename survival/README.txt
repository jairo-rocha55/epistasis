Here there are the routines for survival prob calculations. 
Execute (cd bin; ./ScriptCOAD) to calculate the survival p-values for the 
top pairs of 3cell/. Two groups are calculated: the ones with both genes mutated and the other subjects.  The pvalue is the significance of the survival difference between the two groups.

For each pair, a temporal file TmpgroupsSurv.txt is generated with the list of subjects.
Then the routine  survival.py is called for the Tmp file. The pvalue is copied and the next pair is processed.

The command python3 survivalplot.py TmpgroupSurv.txt is used for plotting the survival distributions using as input the Tmp file. 

 
