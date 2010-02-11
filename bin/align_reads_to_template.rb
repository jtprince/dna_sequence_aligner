#!/usr/bin/ruby

require 'bio'
require 'optparse'


CLUSTALW_OPTS = %w(gapopen gapext dnamatrix type)

=begin
# http://align.genome.jp/clustalw/clustalw_help.html
 >>HELP 8 <<      Help for command line parameters

                DATA (sequences)

-INFILE=file.ext                             :input sequences.



                VERBS (do things)

-OPTIONS	    :list the command line parameters
-HELP  or -CHECK    :outline the command line params.
-ALIGN              :do full multiple alignment 
-TREE               :calculate NJ tree.
-BOOTSTRAP(=n)      :bootstrap a NJ tree (n= number of bootstraps; def. = 1000).
-CONVERT            :output the input sequences in a different file format.


                PARAMETERS (set things)

***General settings:****
-INTERACTIVE :read command line, then enter normal interactive menus
-QUICKTREE   :use FAST algorithm for the alignment guide tree
-TYPE=       :PROTEIN or DNA sequences
-NEGATIVE    :protein alignment with negative values in matrix
-OUTFILE=    :sequence alignment file name
-OUTPUT=     :GCG, GDE, PHYLIP, PIR or NEXUS
-OUTORDER=   :INPUT or ALIGNED
-CASE        :LOWER or UPPER (for GDE output only)
-SEQNOS=     :OFF or ON (for Clustal output only)
-SEQNO_RANGE=:OFF or ON (NEW: for all output formats) 
-RANGE=m,n   :sequence range to write starting m to m+n. 

***Fast Pairwise Alignments:***
-KTUPLE=n    :word size
-TOPDIAGS=n  :number of best diags.
-WINDOW=n    :window around best diags.
-PAIRGAP=n   :gap penalty
-SCORE       :PERCENT or ABSOLUTE


***Slow Pairwise Alignments:***
-PWMATRIX=    :Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename
-PWDNAMATRIX= :DNA weight matrix=IUB, CLUSTALW or filename
-PWGAPOPEN=f  :gap opening penalty        
-PWGAPEXT=f   :gap opening penalty


***Multiple Alignments:***
-NEWTREE=      :file for new guide tree
-USETREE=      :file for old guide tree
-MATRIX=       :Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename
-DNAMATRIX=    :DNA weight matrix=IUB, CLUSTALW or filename
-GAPOPEN=f     :gap opening penalty        
-GAPEXT=f      :gap extension penalty
-ENDGAPS       :no end gap separation pen. 
-GAPDIST=n     :gap separation pen. range
-NOPGAP        :residue-specific gaps off  
-NOHGAP        :hydrophilic gaps off
-HGAPRESIDUES= :list hydrophilic res.    
-MAXDIV=n      :% ident. for delay
-TYPE=         :PROTEIN or DNA
-TRANSWEIGHT=f :transitions weighting


***Trees:***
-OUTPUTTREE=nj OR phylip OR dist OR nexus
-SEED=n        :seed number for bootstraps.
-KIMURA        :use Kimura's correction.   
-TOSSGAPS      :ignore positions with gaps.
-BOOTLABELS=node OR branch :position of bootstrap values in tree display

=end

class Fasta
  class Entry
    attr_accessor :header, :biosequence
    def initialize(header, biosequence)
      @header = header
      @biosequence = biosequence
    end
  end
end

# returns an Array of Fasta::Entry objects
def fasta_sequences(file)
  # remove any commented lines
  clean_string = IO.read(file).split("\n").reject {|line| line =~ /^\#/ }.join("\n")
    io = StringIO.new clean_string
  objects = []
  Bio::FlatFile.auto(io) do |ff|
    ff.each_entry do |entry|
      objects << Fasta::Entry.new( entry.definition, Bio::Sequence::NA.new(entry.seq) )
    end
  end
  objects
end

# returns the index of the starting run of good chars
def find_bad_end(iupac_concensus_string, min_length)
  good_char_count = 0
  char_index = 0
  iupac_concensus_string.each_char do |char|
    if char =~ /[^\?\-Nn]/
      good_char_count += 1
      if good_char_count >= min_length
        break
      end
    else
      good_char_count = 0
    end
    char_index += 1
  end
  char_index - (good_char_count - 1)
end

# returns (start, length) where min_length reads are correct
def find_good_section(iupac_concensus_string, min_length)
  start = find_bad_end(iupac_concensus_string, min_length)
  from_end = find_bad_end(iupac_concensus_string.reverse, min_length)
  length = iupac_concensus_string.length - start - from_end
  [start, length]
end

def hash_opts_to_clustalopts(hash)
  array = []
  hash.each do |k,v|
    if CLUSTALW_OPTS.include?(k.to_s)
      array << "-#{k}=#{v}"
    end
  end
  array
end

def lstrip_dash(string)
  chr = first_non_dash_char(string)
  string[chr..-1]
end

def strip_dash(string)
  ls = lstrip_dash(string)
  lstrip_dash(ls.reverse).reverse
end

def first_non_dash_char(string)
  char_cnt = 0
  string.each_char do |char|
    if char != '-'
      break
    end
    char_cnt += 1
  end
  char_cnt
end

def clustal_align(bioseqs, factory)
  al = Bio::Alignment.new(bioseqs)
  al.do_align(factory)
end

def consensus(files, opt={})
  factory = Bio::ClustalW.new
  clustal_opts = hash_opts_to_clustalopts(opt)
  factory.options = clustal_opts
  files = [files] unless files.is_a?(Array) 
  fasta_entries = files.inject([]) {|ar, file| ar.push( *fasta_sequences(file) ) }
  bioseqs = fasta_entries.map(&:biosequence)
  template = bioseqs.shift
  start_length = []
  trimmed = bioseqs.map do |bseq|
    clust_al = clustal_align([template, bseq], factory)
    cl_cons = clust_al.consensus
    aligned_string = clust_al[1].to_s
    (st, len) = find_good_section(aligned_string, opt[:fidelity_length])
    pristine_read = aligned_string[st, len].gsub('-','')
    better_al = clustal_align([template, Bio::Sequence::NA.new(pristine_read)], factory)
    p better_al[1].to_s
  end
end


opt = {
  :type => 'DNA',
  :gapopen => 20,
  :gapext => 20,
  :dnamatrix => 'IUB',  # "IUB" || "CLUSTALW"
  :fidelity_length => 10,
}

op = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <template>.fasta <read>.txt ..."
  op.separator " " 
  op.on("--fidelity-length <#{opt[:fidelity_length]}>", Integer, "min length of correct reads on the ends") {|v| opt[:fidelity_length] = v }
  op.separator " " 
  op.on("--type <DNA|PROTEIN>", "type of bio sequences") {|v| opt[:type] = v }
  op.on("--gapopen <float>", Float, "gap opening penalty") {|v| opt[:gapopen] = v }
  op.on("--gapext <float>", Float, "gap extension penalty") {|v| opt[:gapext] = v }
  op.on("--dnamatrix <IUB|CLUSTALW>", "DNA weight matrix IUB|CLUSTALW") {|v| opt[:dnamatrix] = v }
  op.separator " " 
  op.separator "the first sequence is assumed the template sequence"

end

if $0 == __FILE__

  if ARGV.size == 0
    puts op
    exit
  end

  consensus(ARGV.map, opt)

end




