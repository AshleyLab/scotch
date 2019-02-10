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

*Note*: `$workingDir` should be an absolute path. 
