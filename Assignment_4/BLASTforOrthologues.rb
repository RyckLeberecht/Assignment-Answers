#!/usr/bin/env ruby
require 'fileutils'
require 'bio'
require 'stringio'

# @author Ryck Leberecht
# A class to make use of BLAST to discover putative Orthologues 
class BLASTforOrthologues
       
    # Guesses the type of sequence present inside a fasta file
    # @param [String] file a file with FASTA format
    # @return [String] the type of sequence
    def self.sequence_type(fasta_file)
      # Open the file and read the fasta entry
      fasta_entry = Bio::Sequence.new(Bio::FlatFile.auto(fasta_file).next_entry.seq)
  
      # Use the `guess` method to guess the type of the sequence
      seq_type = fasta_entry.guess
      
      # Return the appropriate sequence type
      if seq_type == Bio::Sequence::NA
        return 'nucl'
      elsif seq_type == Bio::Sequence::AA
        return 'prot'
      else
        abort("Could not determine sequence type")
      end
    end

    # Creates a blast data base
    # @param [String] fastafile a file with fasta sequences
    # @return [String] a path to the created data base 
    def self.create_blast_db(file)
        # Check if the "Database" folder already exists
        unless File.directory?("Database")
            # Create the "Database" folder if it doesn't exist
            FileUtils.mkdir_p("Database")
        end
      
        # Copy the input file to the "Database" folder
        FileUtils.cp(file, "Database")

        # Determine the type of input
        input_type = BLASTforOrthologues.sequence_type(file)
        # Move to the "Database" folder
        Dir.chdir("Database") do
            # Run makeblastdb command
            system("makeblastdb -in #{file} -dbtype #{input_type} -out #{File.basename(file)}")
        end
        # Return the path to the created database
        path = "./Database/#{File.basename(file)}"
        #return path.gsub(".fa", "")
        return path
    end
    
    # Find the blast top hit
    # @param [String] fasta_sequenze a fasta sequenze
    # @param [Integer] bit_treshold a Interger to evaluate significance
    # @param [Integer] evalue_threshold a Integer to evaluate e
    # @return [String] top_hit a String with the name of the top hit
    # @note Both the bit treshold and e-value treshold are based on (Pearson, 2013; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/)
    # @note This paper claims, that a bit treshold of 50 almost garantees the significance of hits. While the bit_treshold is based on the evalue, the paper also claims that a good cutoff is at e < 1e-10
    # @note It also states that a smaller database usually needs a lower bit_treshold to make hits more likely.
    def self.find_tophit(query,factory, bit_treshold = 0, evalue_threshold = 0)
        # Do the homology search 
        report = factory.query(query)
        # Filter for the top hit
        top_hit = report.hits.select.first { |hit| hit.bit_score > bit_treshold && hit.evalue < evalue_threshold }
        if top_hit
          # Extract the name of the homologue gene/protein
          top_hit = top_hit.definition.match(/^[^\s]+/)[0]
          return top_hit
        else
          return top_hit
        end      
    end
  
    # Find all orthologues
    # @param [File] sequenze_1 a fasta (.fa) file
    # @param [File] sequenze_1 a fasta (.fa) file
    # @return [Hash] a hash filled with all orthologues
    # @note This code takes ages and therefore I break it after a defined amount of checked genes
    def self.find_orthologues(sequenze_1, sequenze_2)
        # Create an empty hash for the orthologues
        orthologues = {}
        # Create the db paths
        path_1 = BLASTforOrthologues.create_blast_db(sequenze_1)
        path_2 = BLASTforOrthologues.create_blast_db(sequenze_2)
        # Create the "Blast" factories
        factory_1 = Bio::Blast.local('tblastn', "#{path_2}")
        factory_2 = Bio::Blast.local('blastx', "#{path_1}")
        
        k = 0 # used to break the loop earlier
        key = 0 # identifiers for orthologue pairs
        Bio::FlatFile.auto(sequenze_1).each_entry do |query|
          # Extract the name of the tested gene
          entry_id_stripped = query.entry_id.strip
          gene_name_1 = entry_id_stripped.split("|")[0]
          # Find the best homologue
          top_hit = BLASTforOrthologues.find_tophit(query, factory_1)
          #puts gene_name_1
          #puts query
          #puts top_hit
          if top_hit && !top_hit.empty?
            Bio::FlatFile.auto(sequenze_2).each_entry do |query_2|
              # Extract the gene name to check if it is identical to the top_hit
              entry_id_stripped = query_2.entry_id.strip
              gene_name_2 = entry_id_stripped.split("|")[0]
              if gene_name_2 == top_hit
                # Find the best homologue
                hit_back = BLASTforOrthologues.find_tophit(query_2, factory_2)
                hit_back = hit_back.split("|")[0]
                #puts hit_back
                if hit_back == gene_name_1
                  # Save the orthologue genes
                  key += 1
                  puts "Found orthologue: #{gene_name_1} + #{top_hit}"
                  orthologues[key.to_s] = [gene_name_1, top_hit]
                end
              end
            end
          end
          k += 1
          if k > 5 # break the loop after a specific amount of runs 
            # Write the orthologues to a txt file
            File.open("orthologues.txt", "w") do |file|
              orthologues.each do |key, value|
                file.write "Orthologue pair number #{key}: #{value}\n"
              end
            end
            break
          end
        end
        return orthologues
    end
end
