for FILE in `cat dataCOAD/subjectsJAG2RASA4B1`; do 
     fgrep -f dataCOAD/tmp  $FILE >> dataCOAD/SNPsSubjects
done

