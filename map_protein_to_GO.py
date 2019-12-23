import re


def Transcript_and_SwissProt_ID_to_dict(blast_filename):
    '''
    Get the transcript ID's (qseqid) and SwissProt ID's (sseqid) without
    a version number and add them to the to the transcript_to_protein 
    dictionary.
    Args:
        blast_filename: a tab separated BLAST output file
    Returns:
        A dictionary (transcript_to_protein)
    '''
    transcript_to_protein = {}
    with open(blast_filename, "r") as blast_file:
      for line in blast_file:
        qseqid, sseqid, pident, length, mismatch, gapopen, \
        qstart, qend, sstart, send, evalue, bitscore = line.rstrip().split("\t")
        transcript, isoform = qseqid.split("|")
        gi_type, gi, sp_type, sp, sp_name = sseqid.split("|")
        sp_id, sp_version = sp.split(".")
        if float(pident) > 99:
            transcript_to_protein[transcript] = sp_id
    return transcript_to_protein


def ProtID_and_GOTerm_to_dict(gene_to_go_filename):
    '''
    Get the protein ID's and their corresponding GO terms and add 
    them to the gene_to_go dictionary
    Args:
        gene_to_go_filename: a GAF file of GO annotations
    Returns:
        A dictionary (gene_to_go)
    '''
    gene_to_go = {}
    with open(gene_to_go_filename, "r") as gene_to_go_file:
        for line in gene_to_go_file:
            db, object_id, object_symbol, qualifier, go_id, *others = line.split("\t")

            '''
            Check if both protein and GO IDs have a value before adding
            '''
            if object_id and go_id:
                if object_id in gene_to_go:
                    gene_to_go[object_id].add(go_id) # Append to Set Data Structure
                else:
                    gene_to_go[object_id] = {go_id} # Assign a Set Data Structure
    # Convert the set structure to List
    gene_new = {key:sorted(list(item)) for key, item in gene_to_go.items()}
    
    return gene_new

def GOID_to_dict(go_terms_filename):
    '''
    Load the GO ID's with their names and add them to the previous dictionary
    Args:
        go_terms_filename: a gene ontology file
    Returns:
        A dictionary (go_to_desc)
    '''
    go_to_desc = {}
    with open(go_terms_filename, "r") as go_terms_file:
        terms = go_terms_file.read()
        terms = re.findall(r"\[Term]\n(.*?)\n\n", terms, re.DOTALL)

    for term in terms:
        go_id = re.search(r"^id:\s+(GO:\d+?)\n", term)
        go_name = re.search(r"^name:\s+(.+?)\n", term, re.M)

        '''
        Check if both ID and name have a value before adding
        '''
        if go_id and go_name:
            go_to_desc[go_id.group(1)] = go_name.group(1)
    
    return go_to_desc

'''
Loop through the file for differential expression and get the protein ID, 
GO term, and GO name and stitch them together.
Print the final results to REPORT output: report file
Args:
    diff_exp_filename: a Trinity output consisting of a matrix of values differentially expressed genes
Returns:
    Output results in TSV file
'''

def Results_to_Report(diff_exp_filename, transcript_to_protein, gene_to_go, go_to_desc):
    report_file = open(report_filename, "w")
    with open(diff_exp_filename, "r") as diff_exp_file:
        diff_exp_file.readline()  # skip header
        for line in diff_exp_file:
            transcript, sp_ds, sp_hs, sp_log, sp_plat = line.rstrip().split("\t")
            protein = transcript_to_protein.get(transcript, "NA")
            go_id = gene_to_go.get(protein, "NA")
            go_desc = ''
            for item in go_id:                
                go_desc += '\t' + item + '\t' + go_to_desc.get(item, "NA") +'\n'+'\t'*5

            report_file.write("\t".join([transcript, protein, sp_ds, sp_hs,
                                         sp_log, sp_plat, go_desc])+'\n')
    report_file.close()



if __name__ == "__main__":
    gene_to_go_filename = <'GAF file'>
    blast_filename = <'BLAST output file'>
    diff_exp_filename = <'Trinity Assembly Matrix'>
    go_terms_filename = <'GO file'>
    report_filename = <'Output.TSV'>

    transcript_to_protein = Transcript_and_SwissProt_ID_to_dict(blast_filename)
    gene_to_go = ProtID_and_GOTerm_to_dict(gene_to_go_filename)
    go_to_desc = GOID_to_dict(go_terms_filename)
    go_to_desc = Results_to_Report(diff_exp_filename, transcript_to_protein, 
                                   gene_to_go, go_to_desc)
