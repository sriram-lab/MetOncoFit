for i in `cat csv.lst`
do
    #grep -v "GENE" $i.csv > temp
    #cat header temp > $i.csv
    python $1\_ensemble.py $i.train.csv > $i.pred\_$1
done
