######### Get protein domain info for all coding genes in python using Biomart API ###########
from biomart import BiomartServer
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Get ensembl id')
parser.add_argument('gtf_dir', type=str,
                    help='Directory to gene model file')
parser.add_argument('-out', type=str, default="domain_table.txt",
                    help='Output domain table file directory')
parser.add_argument('-version', type=str, default="GRCh38",
                    help='Annotation version, GRCh37 or GRCh38')

if __name__ == '__main__':
    args = parser.parse_args()
    # set server
    server = BiomartServer("http://grch37.ensembl.org/biomart") if args.version == "GRCh37" else BiomartServer("http://useast.ensembl.org/biomart")
    print("Using biomart server for %s" % args.version)
    human_genes = server.datasets['hsapiens_gene_ensembl']

    all_ID = subprocess.check_output( "awk '$3 == \"gene\"' %s | grep protein_coding | cut -f9 | cut -d ' ' -f2 | cut -c 2-16 | sort -u" % args.gtf_dir, shell = True).rstrip('\n').split('\n')
    print("Got %s protein coding gene IDs" % len(all_ID))
    batches = [all_ID[i:i + 300] for i in xrange(0, len(all_ID), 300)]
    print("Divided into %s batches of size 300" % len(batches))
    counter = 1
    # do the querying in batch of 300
    with open(args.out, 'w') as output_file:
        for IDs in batches:
            print("Querying batch %s out of %s" % (counter, len(batches)))
            counter += 1;
            response_interpro = human_genes.search({
              'filters': {
                  'ensembl_gene_id': IDs
              },
              'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'interpro_start', 'interpro_end', 'interpro', 'interpro_description']
            })
            response_interpro = response_interpro.text.rstrip('\n').encode('UTF-8').split('\n')
            response_interpro[:] = [line.split('\t') for line in response_interpro]
            response_pfam = human_genes.search({
              'filters': {
                  'ensembl_gene_id': IDs
              },
              'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'pfam_start', 'pfam_end', 'pfam']
            })
            response_pfam = response_pfam.text.rstrip('\n').encode('latin-1').split('\n')
            response_pfam[:] = [line.split('\t') for line in response_pfam]
            response_combined = [pfam[:4] + interpro[4:] for pfam in response_pfam for interpro in response_interpro if pfam[0] == interpro[0] and pfam[1] == interpro[1] and pfam[2] == interpro[2] and pfam[3] == interpro[3] and (pfam[2] != '') and (interpro[2] != '')]
            for line in response_combined:
                output_file.write("\t".join(line) + '\n')

