% read the Genbank genome into Matlab
% or handle the fasta sequence text file if you prefer.
genome=genbankread('sequence.gb');

% genome.CDS contains all the coding sequences and their information
% e.g., genome.CDS(1).protein_id (e.g., WP_005685692.1) is the ID that is in the carveme model.genes
% And genome.CDS(1).text contains the locus tag (e.g., LGG_RS00005) that is used in the KBase
% model. Write something here to unify the IDs. Can use one ID format in
% model.genes and then add the other format in model.(something). I think
% using the locus tag is the more usual convention and at the same time noting the gene ID
% and protein ID in other fields in the model. The KBase model already have
% the gene ID (which is in model.proteinisncbigeneID and geneEntrezID) but not the protein ID 

% Step 1:
% for each gene in the CarveMe model, add the protein ID in e.g., model.proteinIDs,
% find the corresponding locus tag and replace the protein IDs by the locus tag

% Step 2:
% rewrite the carveme model as SBML file using writeCbModel