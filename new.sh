#!/bin/bash -ex
file=`ls ToRun/`
echo $file
if [ -f ToRun/$file ]
then
echo ToRun/$file
    year=`grep "year" ToRun/$file | awk '{print $2}'`
    week=`grep "week" ToRun/$file | awk '{print $2}'`
    sqlite1=`grep "run1" ToRun/$file | awk '{print $2}'`
    sqlite2=`grep "run2" ToRun/$file | awk '{print $2}'`
cp ToRun/$file RunFiles/.
rm ToRun/$file
echo "./L1RateValidation_2022.sh $sqlite1 $sqlite2 $week $year"
./L1RateValidation_2022.sh $sqlite1 $sqlite2 $week $year 
git commit -a -m "clean ToRun files"
git push
else
echo "No new files"
fi

