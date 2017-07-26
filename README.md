# INTEGRATE-Vis

INTEGRATE-Vis is a gene fusion visualization tool. It is written in Python.

### Prerequisites

Please make sure you have installed the following tools:

  - [Python](https://www.python.org/)
  - [CMake](https://cmake.org/)
  - [GCC](https://gcc.gnu.org/)
  - [Matplotlib](http://matplotlib.org/)
  - [gtfToGenePred](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred)

If not, please install these languages or tools.
Note: Matplotlib can also be installed through using [EPDFree](https://www.enthought.com/products/epd/).

### Installation

Download INTEGRATE-Vis at https://github.com/ChrisMaherLab/INTEGRATE-Vis.

Run the installation script:

```sh
$ cd INTEGRATE-Vis.1.0.0
$ chmod +x install.sh
$ ./install.sh -o /opt/bin/
```

Note that you can choose wherever you like to install the software. It can be different from "/opt/bin/".

### Input

If you type the following:

```sh
$ ./Integrate-Vis.py --help
```
or

```sh
python ./Integrate-vis.py --help
```

you will see the usages for 4 sub utils and explanations, which are for the 4 types of figures that INTEGRATE-Vis plots. The utils includes: structure, domain, exon expression, and gene expression in a cohort.

You can run the following commands:

    Integrate-vis structure <parameters>
    Integrate-vis domain    <parameters>
    Integrate-vis exon-exp  <parameters>
    Integrate-vis gene-exp  <parameters>

For example, you can run
```sh
Integrate-vis structure --help
```
to see what input values or files are needed for the structure util.

Input file formats include BEDPE, BAM, GTF, and TSV.

The BEDPE format for gene fusions follows the standardized format provided by The ICGC-TCGA DREAM Somatic Mutation Calling - RNA Challenge ([SMC-RNA](http://dreamchallenges.org/)).

The GTF files can be downloaded from Ensembl:

GRCh37, e.g., v75: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

GRCh38, e.g., v86: ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz

Ideogram.tsv can be created by the following R commands:

```sh
library(IdeoViz)
ideo <- getIdeo("hg38")
write.table(ideo,"Desktop/ideogram.hg38.tsv", sep="\t", quote=F, row.names = F)
```

Domain_table.tsv can be created by the python script domain_table.prep.py under the src/ dir.

Both ideogram.hg38.tsv and domain.table.tsv have also been included under the data/ dir.

### output

Output files are figures in PDF format.

### Important

The chromosome names in the reference genome, the gene models, and the fusions should be consistent.

### Example command lines:

Plotting structure, domain, and exon expression for single sample:

```sh
Integrate-vis.py structure -b SAMPLE_NAME.bedpe -s SAMPLE_NAME -d ideogram.hg38.txt -r Homo_sapiens.GRCh38.85.fa -g Homo_sapiens.GRCh38.85.gtf (-m SAMPLE_NAME.sort.bam) -o OUTPUT_DIR -k
Integrate-vis.py domain -b SAMPLE_NAME.bedpe -s SAMPLE_NAME -d domain.table.txt -r Homo_sapiens.GRCh38.85.fa -g Homo_sapiens.GRCh38.85.gtf -o OUTPUT_DIR -k
Integrate-vis.py exon-exp -b SAMPLE_NAME.bedpe -s SAMPLE_NAME -m SAMPLE_NAME.sort.bam -r Homo_sapiens.GRCh38.85.fa -g Homo_sapiens.GRCh38.85.gtf -o ./OUTPUT_DIR -k
```

Plotting gene partner expression in cohort:

a. preperation steps to generate COHORT_NAME.fusions.tsv and COHORT_NAME.gene_expression.tsv (Only need to run once).

```sh
pd_fusion_converter.py -r Homo_sapiens.GRCh38.85.fa -g Homo_sapiens.GRCh38.85.gtf -o OUTPUT_DIR -k -a COHORT_NAME.fusion.bedpe.dir.tsv -c COHORT_NAME
pd_expression_converter.py -a COHORT_NAME.expression.files.dir.tsv -g 1 -e 7 -o OUTPUT_DIR -c COHORT_NAME
```

The above commands generate COHORT_NAME.fusions.tsv and COHORT_NAME.gene_expression.tsv under OUTPUT_DIR.


b. running the gene-exp util:
```sh
Integrate-vis.py gene-exp -f OUTPUT_DIR/COHORT_NAME.fusions.tsv -e OUTPUT_DIR/COHORT_NAME.gene_expression.tsv -t COHORT_NAME.type.tsv -g Homo_sapiens.GRCh38.85.gtf -m "Read count" -c COHORT_NAME -o OUTPUT_DIR -k
```

The input files COHORT_NAME.fusion.bedpe.dir.tsv, COHORT_NAME.expression.files.dir.tsv, COHORT_NAME.type.tsv are for the directories of the BEDPE files, the directories of the gene expression files (Can be from any tool), and the sample types. Examples of these files can be found under data/ dir.

### Enjoy!
