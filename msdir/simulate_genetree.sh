#!/bin/bash 
echo Hello World
#not working tho
./ms 50 1 -t 30 | tail -n +7 > ../data/repo2/temp.txt
uniq -c ../data/repo2/temp.txt | sed s/" "/:/g | sed s/:::::://g | sed 's/^/1/' | sed -e 's/\(.\)/\1 /g' > ../data/repo2/genetree_50_t30.txt

# for genetree
./genetree ../data/repo2/genetree_50_t30.txt 30 1 123123 -f surf_50t30.txt -g 20 40 11 -2


uniq -c gt100_t10.txt | sed s/" "/:/g | sed s/:::::://g | sed 's/^/1/' | sed -e 's/\(.\)/\1 /g' > gt100_t10.csv

sed s/\"//g temp.csv | sed -e 's/\(.\)/\1 /g'