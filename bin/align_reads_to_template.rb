#!/usr/bin/ruby

require 'bio'
require 'optparse'

require 'bio/alignment/dna_sequence'

# returns an Array of Bio::FastaFormat objects.
# The method will remove any commented lines first
def fasta_entries(file)
  clean_string = IO.read(file).split("\n").reject {|line| line =~ /^\#/ }.join("\n")
    io = StringIO.new clean_string
  objects = []
  Bio::FlatFile.auto(io) do |ff|
    ff.each_entry do |entry|
      objects << entry
    end
  end
  objects
end

opt = Bio::Alignment::DNASequenceReads::ALIGN_OPTS.dup

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

if ARGV.size == 0
  puts op
  exit
end

dumpfile = "mydump.rbdump"

palign = 
  if File.exist?(dumpfile)
    Marshal.load(IO.read(dumpfile))
  else
    files = ARGV.map
    fasta_entries = files.inject([]) {|ar, file| ar.push( *fasta_entries(file) ) }
    bioseqs = fasta_entries.map {|entry| Bio::Sequence::NA.new(entry.seq) }
    labels = fasta_entries.map {|entry| entry.definition }

    pairwise = Bio::Alignment::DNASequenceReads.align_pairwise(bioseqs, opt)
    File.open(dumpfile, 'w') {|out| out.print Marshal.dump(pairwise) }
  end

Bio::Alignment::DNASequenceReads.merge_pairwise(palign)


#Bio::Alignment::DNASequenceReads.print_align(STDOUT, align_s, labels)




