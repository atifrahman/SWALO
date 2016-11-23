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

To create index of the contig file and map mate-pair reads (reads from jumping library) using Bowtie 2 run
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
*scaffolds.fa* will contain the results in fasta format. For paired-end reads (regular or fosmid libraries) delete the `--jump`. The `<minContigLength>` is the minimum length of contigs that will be used for learning insert size distribution. Other options are

`--dist <mean> <standardDev> or -d <mean> <standardDev>` : Use Normal distribution with `<mean>` and `<standardDev>` instead of learning insert size distribution. Recommended if number of inserts to learn the distribution is less than 100,000.

`--conservative or -c` : Use a conservative mode (please see paper). Recommended if standard deviation of insert library is high (>1000). 

###### Using SWALO with Bowtie

Reads can also be mapped with Bowtie although number of reads mapped can be quite low for poor quality libraries. This can be done by running the script `.\runSwalo_bowtie` as discussed above or run
```
bowtie-build <contigFile_Fasta> contigs_bowtie
bowtie -k 5 -v <numberMismatches> contigs_bowtie <readFile1_Fastq> -S <mapFile1_SAM>
bowtie -k 5 -v <numberMismatches> contigs_bowtie <readFile2_Fastq> -S <mapFile2_SAM>
bowtieconvert <mapFile1_SAM> <mapFile2_SAM> <contigFile_Fasta> <maxInsertSize> <insertSizeLimit(>=maxInsertSize)> --jump
```
As before `-a` can be used instead of `-k 5` and `--jump` need to be deleted for regular or fosmid libraries. The last two steps `align` and `swalo` remain the same.

###### Important

If Bowtie is used either use a single thread for mapping or sort the results before running `bowtieconvert` as reads can be in different order if multiple threads are used.


#### Combining multiple insert libraries
There are two approaches for multiple insert size libraries. The first is to build the scaffold graph separately and then combine the graphs before actually merging contigs. The other approach is to do the scaffolding hierarchically. We recommend the first approach unless some insert libraries have large standard deviations in which case we recommend scaffolding using all other libraries using the first approach and then scaffold using libraries with high standard deviation together using the conservative mode.

For this first run SWALO on each library separately in different directories as discussed above. This will generate three files with names staring with the same `<prefix>` and a text file `prefixes.txt` containing the `<prefix>`. Now create a new directory, copy into the new directory the three files starting with `<prefix>` from each directory and create a new text file `prefixes.txt` in the new directory containing one prefix in each line. Then in the new directory run
```
swaloFile <contigFile_Fasta> [options]
```
This will generate the scaffolds. Add `--conservative or -c` to run in conservative mode. 

An example script to do this is provided with SWALO which will generate results in `combined_bowtie2` with `bowtie2` and `long_bowtie2` being the source directories. 

#### SWALO Help

Please email your questions, comments, suggestions, and bug reports to atif DOT bd AT gmail DOT com

