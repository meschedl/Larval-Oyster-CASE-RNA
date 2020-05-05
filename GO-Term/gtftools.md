Trying to use [gtftools](http://www.genemine.org/gtftools.php) to calculate the length of the transcripts. Ultimately did not work.

```
wget http://www.genemine.org/codes/GTFtools_0.6.9.zip
unzip  GTFtools_0.6.9.zip
cd  GTFtools_0.6.9
chmod u+x gtftools.py
```
Needs this other python program to run so have to download that one
```
wget https://github.com/python/cpython/blob/3.8/Lib/argparse.py
chmod u+x argparse.py

mv GTFtools_0.6.9/gtftools.py .


ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/restring/C_Vir_ST_merged_NO_A_test.gtf .
```


```
python gtftools.py -l gene_length.gtf C_Vir_ST_merged_NO_A_test.gtf
```
ok keeps giving me this error  
File "gtftools.py", line 283  
  exon={}  
        ^  
TabError: inconsistent use of tabs and spaces in indentation  

so I guess there's something wrong with the code. I will try downloading a different version I guess

```
rm gtftools.py
rm -r GTFtools_0.6.9
rm GTFtools_0.6.9.zip

wget http://www.genemine.org/codes/GTFtools_0.6.5.zip
unzip GTFtools_0.6.5.zip
cd GTFtools_0.6.5/
chmod u+x gtftools.py
cd ..
mv GTFtools_0.6.5/gtftools.py .

python gtftools.py -l gene_length.gtf C_Vir_ST_merged_NO_A_test.gtf
```

Another error.

File "gtftools.py", line 785  
  print 'Analyzing chromosomes:'  
                               ^  
SyntaxError: Missing parentheses in call to 'print'. Did you mean print('Analyzing chromosomes:')?

maybe I am not in the right version of python?

```
python --version
Python 3.6.8 :: Anaconda, Inc.

conda activate python27

python gtftools.py -l gene_length.gtf C_Vir_ST_merged_NO_A_test.gtf
```
File "gtftools.py", line 4, in <module>  
   import argparse  
 File "/home/mschedl/Working-CASE-RNA/histat/stringtie/restring/GOGO/argparse.py", line 7  
   `<!DOCTYPE html>`  
   ^  
SyntaxError: invalid syntax

ok why does it say that
```
nano argparse.py   
```
this is WROOOOOOOOOOONG it's the html of the github page not the code
```
rm argparse.py
```
what about if I use the raw code link?
```
wget https://raw.githubusercontent.com/python/cpython/3.8/Lib/argparse.py

chmod u+x argparse.py
nano argparse.py
```
looks right!
```
python gtftools.py -l gene_length.gtf C_Vir_ST_merged_NO_A_test.gtf
```
Nope  
Traceback (most recent call last):  
  File "gtftools.py", line 4, in <module>  
    import argparse  
  File "/home/mschedl/Working-CASE-RNA/histat/stringtie/restring/GOGO/argparse.py", line 606  
    raise ValueError("invalid nargs value") from None  
                                               ^  
SyntaxError: invalid syntax

ok... now is this the wrong version for this one?

I will try going back to the other version of gtf tools...
```
wget http://www.genemine.org/codes/GTFtools_0.6.9.zip
unzip  GTFtools_0.6.9.zip
cd  GTFtools_0.6.9
chmod u+x gtftools.py

cd ..
mv GTFtools_0.6.9/gtftools.py .

python gtftools.py -l gene_length.gtf C_Vir_ST_merged_NO_A_test.gtf
```
nope, same error.


Could be that the two programs need different versions of python. or something else is wrong. But I don't know what to do if it's telling me their code is wrong.

Guess I'll try to do something different.




might need to have chromosomes labeled as 1,2,3 etc not NC_007175.2


search and replace
