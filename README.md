# metaClean_long-reads
[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) workflow to eliminate low-quality reads and human contamination from raw long-read metagenomes. 

The user provides a path to raw (long-read) metagenomes to clean. The workflow runs `chopper` for adaptor trimming and removal of low-quality reads (>20QC) and `minimap2/samtools` to remove human contamination. This workflow eliminates raw reads and intermediate files to improve space efficiency. The output of this workflow is clean compressed metagenomes. 

## Installation

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
   
3. Create a conda environment
   
<pre><code>conda create -n snakemake-7.3.7 python=3.10 -y 
conda activate snakemake-7.3.7</code></pre> 

3. Install snakemake v7.3.7 (using conda or mamba, a faster conda alternative)

<pre><code>conda create -n snakemake-7.3.7 -c conda-forge -c bioconda snakemake=7.3.7</code></pre>

If you don't have access to bioconda, you can also try installing this using `pip`

<pre><code>pip install snakemake==7.3.7</code></pre>

4. Clone the repository

<pre><code>git clone https://github.com/luisgv467/metaClean_long-reads.git</code></pre>

## How to Run

### Prepare your files before running

1. Download the latest version of the human genome (GRCh38)

<pre><code>wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz</code></pre>

2. Build a `Minimap2` index

<pre><code>#Create the directory where you will store your index
mkdir index
cd index

#Build the index
minimap2 -t 24 -I 100G -d human.mmi ../GCF_000001405.40_GRCh38.p14_genomic.fna </code></pre>

## Run snakemake workflow

1. We need the basename for the raw metagenomic assemblies that we plan to process. This workflow assumes that the raw-metagenomes are compressed and have a file extension `.fastq.gz`.

<pre><code> ls path_to_raw_reads/*_2.fastq.gz | sed -e 's|path_to_raw_reads/||' -e 's|_2.fastq.gz||' > metagenome_paths.txt </pre></code> 

2. Modify the `config/config.yml` with your input and output locations

<pre><code># Required arguments

#File with the basenames of the samples you plan to process
input:
    "/data/pam/lg21g/scratch/Programs/snakemake/MetaClean_long-reads/metagenome_paths.txt"

#Path to raw long-reads fastqs
input_path:
    "/data/pam/lg21g/scratch/chapters/chapter_1/Resp_microbiome_db/raw_long-reads/"

#Output directory to store results
output:
    "/data/pam/lg21g/scratch/chapters/chapter_1/Resp_microbiome_db/clean_long_reads/"</pre></code>

3. Give execution permissions

This workflow counts the number of clean long reads as the last step. However, the script needs execution permissions to work:

<pre><code>chmod +x scripts/count_reads.sh</pre></code>

5. Run Snakemake

If using a server:
<pre><code>snakemake --use-conda -k -j 100 --profile config/lfs --rerun-incomplete --latency-wait 120 --scheduler greedy </pre></code>

If running locally:
<pre><code>snakemake --use-conda -k -j 10 --rerun-incomplete --latency-wait 60 </pre></code>

## Output

The clean fastqs will be located in the directory `clean_long_reads/`. You will also find a file called `Number_reads.txt` with the read count of all your clean metagenomes.  









