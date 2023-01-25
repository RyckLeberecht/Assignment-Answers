require "./BLASTforOrthologues.rb"
# The fasta files have be put into the same folder this file is in!
BLASTforOrthologues.find_orthologues("Spombe_peptides.fa","Arabidopsis_nucleotides.fa")
