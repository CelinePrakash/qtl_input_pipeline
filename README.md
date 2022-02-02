# qtl_input_pipeline
A Pipeline to prepare Restriction site-associated DNA sequencing (RAD-seq) data for quantitative trait locus (QTL) analysis

Please refer to the included PDF for full documentation

## Requirements
* python > 2.7
* pandas
* numpy
* matplotlib
* seaborn

## Usage
Example command:
`python QTL_input_pipeline.py -i INPUT_TABLE.tsv -g ADAPTERS.fa -r REFERENCE.fa -j 0.9 -s 0.9 -m THREADS -c CROSSTYPE -z PARENTS_LIST -p GRANDPARENTS_LIST -o OUTPUT_DIR -v OUTPUT.vcf`

Please refer to `python QTL_input_pipeline.py --help` for a complete list of options.
