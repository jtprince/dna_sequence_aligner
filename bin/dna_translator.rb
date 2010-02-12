#!/usr/bin/ruby

require 'bio'
require 'optparse'

opt = {
  :frameshift => 0
}
op = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <dnaseq>.fasta"
  op.on("-s", "--shift <int>", Integer, "frameshift") {|v| opt[:frameshift] = v }
end
op.parse!

frameshift = opt[:frameshift]
p frameshift

if ARGV.size == 0
  puts op
  exit
end

file = ARGV.shift
string = IO.read(file).split("\n").reject {|line| line =~ /^\#/ }.join("\n")
st = StringIO.new(string)
ff = Bio::FlatFile.auto(st)

seqs = []
ff.each_entry do |entry|
  seq = entry.seq
  seqs << seq[frameshift..-1]
end

length = 70

seqs.each do |seq|
  bsq = Bio::Sequence::NA.new(seq)
  protseq = bsq.translate
  start = 0
  loop do
    break if start >= seq.length
    frag = seq[start, length]
    puts frag
    prot_line = (start...(start+length)).to_a.map do |x|
      if x % 3 == 0
        prot_i = x / 3
        protseq[prot_i,1]
      else
        " "
      end
    end.join
    puts prot_line
    start += length
  end
  print "NUM START/STOP CODONS: "
  puts protseq.to_s.split("").select {|v| v == '*' }.size
end


