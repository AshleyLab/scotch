# scotch

```
Usage ./scotch.sh [workingDir] [bedsDir] [bam] [id] [fastaRef] [gatkJAR] [rfsDir]
        workingDir      absolute path to directory where Scotch should put intermediate files and results
        bedsDir         directory with .bed files listing the regions Scotch should examine
        bam             BAM file for which Scotch should call indels
        id              a name for this sample
        fastaRef        FASTA reference for genome build
        gatkJAR         Genome Analysis Toolkit JAR file
        rfsDir          directory with computed region features
More at https://github.com/AshleyLab/scotch.
```

Running ./scotch.sh creates a directory, ${workingDir}/do/ and populates it with executables to run to execute Scotch's pipeline. 

```
			./scotch.sh			X
--------------------------------------------------------
1. BAM preparation	   
                           rmdup.sh
		      split.${c}.sh
		     unclip.${c}.sh
--------------------------------------------------------
2. Calculate features
              gd.${c}.sh gn.${c}.sh gr.${c}.sh
--------------------------------------------------------
3. Compile features 
                    compile.${c}.sh
--------------------------------------------------------
4. Predict
                    predict.${c}.sh 
