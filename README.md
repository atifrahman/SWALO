# SWALO
Scaffolding with assembly likelihood optimization

##Installation

To install SWALO, downlaod the latest version from [http://atifrahman.github.io/SWALO/](http://atifrahman.github.io/SWALO/) and run the following commands.
```
tar xzf swalo-[version]-beta.tar.gz
cd swalo
make
```

##Running

#### Using shell script

SWALO can be run by copying the script `runSwalo` for Bowtie 2 (recommended) to the desired directory, modifying the filenames and parameters, and running 
```
.\runSwalo
```
This will run all the steps for mapping and scaffolding. The scaffolds will be in the file *scaffolds.fa* in fasta format.

#### Running commands separately

The directories containing Bowtie 2 and SWALO need to be added to the $PATH variable.

###### Mapping reads using Bowtie 2

To map mate-pair reads (reads from jumping library) using Bowtie 2 run
```
bowtie2-build <contigFile_Fasta> contigs
bowtie2 -k 5 -x contigs -X <maxInsertSize> --rf -p <numberThreads> -1 <readFile1_Fastq> -2 <readFile2_Fastq> -S <mapFile_SAM>
```
For paired-end reads (regular or fosmid libraries) delete the `--rf`. For non repeat rich genomes `-a` can be used instead of `-k 5`. The `<maxInsertSize>` can be set something large. This is unlikely to effect results as cut-off points will be learned from data during scaffolding but setting a large value will increase the time needed to map reads.

###### Scaffolding using SWALO

For scaffolding run the following commands
```
bowtie2convert <mapFile_SAM> <contigFile_Fasta> <maxInsertSize>
align <contigFile_Fasta>
swalo <contigFile_Fasta> <minContigLength> --jump [options]
```
For paired-end reads (regular or fosmid libraries) delete the `--jump`. The `<minContigLength>` is the minimum length of contigs that will be used for learning insert size distribution. Other options are

`--dist <mean> <standardDev> or -d <mean> <standardDev>` : Use Normal distribution with `<mean>` and `<standardDev>` instead of learning insert size distribution. Recommended if number of inserts to learn the distribution is less than 100,000.

`--conservative or -c` : Use a conservative mode (please see paper). Recommended if standard deviation of insert library is high (>1000). 
