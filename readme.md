Installation

install python modules
pip install -r requirements.txt

mkdir files/
mkdir files/data/
mkdir files/libs/
mkdir files/output/


-> shortstack
git clone https://github.com/MikeAxtell/ShortStack.git files/libs/

- samtools 1.2  
    wget  "https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2"  && tar xfj samtools-1.2.tar.bz2  && (cd samtools-1.2 && make)
    export PATH=/home/trigo/runs3/miRNA-MITE/samtools-1.2:${PATH}
- ViennaRNA
    wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_17_04/viennarna_2.4.3-1_i386.deb
    sudo dpkg -i viennarna_2.4.3-1_i386.deb
OR
    wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.2.tar.gz
    tar -zxvf ViennaRNA-2.4.3.tar.gz
    cd ViennaRNA-2.4.3
    ./configure
    make
    sudo make install



nohup ./files/libs/ShortStack/ShortStack --readfile files/data/T --outdir files/output/sstack_genesseqs_mirna21 --genomefile files/output/genes_sequences.fasta --bowtie_cores 3  & 