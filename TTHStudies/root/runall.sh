#! /bin/sh


#combineCards.py Name1=SL_cat1_New.txt Name2=SL_cat2_New.txt Name3=SL_cat3_New.txt > SL_New.txt
#combineCards.py Name1=DL_cat4_New.txt > DL_New.txt 
#combineCards.py Name1=SL_New.txt Name2=DL_New.txt > COMB_New.txt

#combineCards.py Name1=SL_cat1a.txt Name2=SL_cat1b.txt Name3=SL_cat2.txt Name4=SL_cat3.txt Name5=SL_cat4.txt Name6=SL_cat5.txt Name7=DL_cat6.txt Name8=DL_cat7.txt > COMB.txt

#combine -M Asymptotic -t -1 -n COMB     -d COMB_New.txt 
#combine -M Asymptotic -t -1 -n SL       -d SL_New.txt 
#combine -M Asymptotic -t -1 -n DL       -d DL_New.txt 

combine -M Asymptotic -t -1 -n SL_cat1  -d SL_cat1_New.txt 
#combine -M Asymptotic -t -1 -n SL_cat2  -d SL_cat2_New.txt 
#combine -M Asymptotic -t -1 -n SL_cat3  -d SL_cat3_New.txt 
#combine -M Asymptotic -t -1 -n SL_cat4  -d SL_cat4_New.txt 
#combine -M Asymptotic -t -1 -n SL_cat5  -d SL_cat5_New.txt 
#combine -M Asymptotic -t -1 -n DL_cat6  -d DL_cat4_New.txt 
