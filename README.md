# MinimalResolution

The algorithm is explained in the pdf file MinimalResolution.pdf

The codes can be compiled with GCC. The GNU Multiple Precision Arithmetic Library should be installed.

To compile, run the following batch files:
sh st_compiling
sh BPtable_complile
sh BP_compile

To get the minimal resolution for BP/I, for t<=50, s<=21 (say), run
./mr_st 25 21
To get the structure maps of the BP Hopf algebroid for t<=50, run
./BPtab 25
To get the minimal resolution for BP, for t<=50, s<=20, run
./mr_BP 25 20

Warning:
The first input parameter is half of t.
The three executalbes are dependent, and should be run in the above order. 
The s for the minimal resolution for BP/I should be at least one larger than the s for that of BP.

Any mistake of the input could result in unpredictible behaviour, usually a break-down of the program such as a segmentation error.
