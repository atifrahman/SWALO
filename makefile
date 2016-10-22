all:
	g++ bowtie2convert.cpp -lm -o bowtie2convert
	g++ bowtieconvert.cpp -lm -o bowtieconvert
	g++ pAlign.cpp -lm -lpthread -o align
	g++ swalo.cpp -lm -lpthread -o swalo
	g++ swalo_file.cpp -lm -lpthread -o swaloFile
