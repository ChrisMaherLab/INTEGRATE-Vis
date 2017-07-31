##################################################
## Please run the following two commands first: ##
##################################################

######### Get the IDs of all coding genes from gtf ###################################################
# cat Homo_sapiens.GRCh38.85.gtf | awk '$3 == "gene"' | grep protein_coding > all_coding_genes.txt
# cut -f9 all_coding_genes.txt | cut -d ' ' -f2 | cut -c 2-16 | sort -u > all_coding_genes_ID.txt


#########################################
## Next run the following python code  ##
#########################################

######### Get protein domain info for all these genes in python using Biomart API ####################
from biomart import BiomartServer
server = BiomartServer("http://useast.ensembl.org/biomart")
human_genes = server.datasets['hsapiens_gene_ensembl']

f = open('all_coding_genes_ID.txt')
all_ID = f.read()
all_ID = all_ID.split('\n')
f.close()

# do the querying in batch of 300
for i in range(len(all_ID)/300):
    output_file = open('domain_table.txt', 'a+')
    print("querying: " + str(i * 300) + " to " + str((i+1) * 300))
    response_interpro = human_genes.search({
      'filters': {
          'ensembl_gene_id': all_ID[i*300:(i+1)*300]
      },
      'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'interpro_start', 'interpro_end', 'interpro', 'interpro_description']
    })
    response_interpro = response_interpro.text.rstrip('\n').encode('latin-1').split('\n')
    response_interpro[:] = [line.split('\t') for line in response_interpro]
    response_pfam = human_genes.search({
      'filters': {
          'ensembl_gene_id': all_ID[i*300:(i+1)*300]
      },
      'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'pfam_start', 'pfam_end', 'pfam']
    })
    response_pfam = response_pfam.text.rstrip('\n').encode('latin-1').split('\n')
    response_pfam[:] = [line.split('\t') for line in response_pfam]
    response_combined = [pfam[:4] + interpro[4:] for pfam in response_pfam for interpro in response_interpro if pfam[0] == interpro[0] and pfam[1] == interpro[1] and pfam[2] == interpro[2] and pfam[3] == interpro[3] and (pfam[2] != '') and (interpro[2] != '')]
    for line in response_combined:
        output_file.write("\t".join(line) + '\n')
    output_file.close()

# do the remaining few ones
output_file = open('domain_table.txt', 'a+')
response_interpro = human_genes.search({
  'filters': {
      'ensembl_gene_id': all_ID[19800:]
  },
  'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'interpro_start', 'interpro_end', 'interpro', 'interpro_description']
})
response_interpro = response_interpro.text.rstrip('\n').encode('latin-1').split('\n')
response_interpro[:] = [line.split('\t') for line in response_interpro]
response_pfam = human_genes.search({
  'filters': {
      'ensembl_gene_id': all_ID[19800:]
  },
  'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'pfam_start', 'pfam_end', 'pfam']
})
response_pfam = response_pfam.text.rstrip('\n').encode('latin-1').split('\n')
response_pfam[:] = [line.split('\t') for line in response_pfam]
response_combined = [pfam[:4] + interpro[4:] for pfam in response_pfam for interpro in response_interpro if pfam[0] == interpro[0] and pfam[1] == interpro[1] and pfam[2] == interpro[2] and pfam[3] == interpro[3] and (pfam[2] != '') and (interpro[2] != '')]
for line in response_combined:
    output_file.write("\t".join(line) + '\n')
output_file.close()
