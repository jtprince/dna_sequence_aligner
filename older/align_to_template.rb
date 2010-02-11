#!/usr/bin/ruby

require 'bio'
require 'optparse'

def printv(*args)
  if $VERBOSE
    print(*args) ; $stdout.flush
  end
end

def putsv(*args)
  if $VERBOSE
    puts(*args) ; $stdout.flush
  end
end


def print_align(io, sequences, labels, opts={})
  opts = {:cutoff => 100, :start => 0, :chars => 20}.merge(opts)
  (start, length, chars) = opts.values_at(:start, :cutoff, :chars)

  loop do
    fin = false
    sequences.zip(labels) do |string, label|
      fin = (start >= string.length )
      break if fin
      io.puts "#{label.exactly_chars(chars)} : #{string[start,length]}"
    end
    io.puts " "
    break if fin
    start += length
  end
end

class String

  # returns [% same chars, % same letters (template), % same letters self]
  def percent_similar_to(template)
    num_same = 0
    num_same_letters = 0
    num_letters_in_template = 0
    num_letters_in_self = 0
    (0...(template.size)).each do |i|
      if letters = (self[i,1] =~ /[A-Za-z]/)
        num_letters_in_self += 1 
      end
      if template[i,1] =~ /[A-Za-z]/
        num_letters_in_template += 1
      end
      if self[i] == template[i]
        num_same += 1 
        if letters
          num_same_letters += 1
        end
      end
    end
    [[num_same, template.size], [num_same_letters, num_letters_in_template], [num_same_letters, num_letters_in_self]].map {|a,b| (a.to_f / b) * 100 }
  end

  def exactly_chars(n)
    at_least = "%#{n}s" % self
    at_least[0,n]
  end

end

def seqs_and_defs(file)
  ff = Bio::FlatFile.auto(file)
  na_seq_objs = []
  definitions = []
  ff.each_entry do |entry|
    definitions << entry.definition
    na_seq_objs << Bio::Sequence::NA.new(entry.seq.to_s)
  end
  [na_seq_objs, definitions]
end

outfile = "aligned.txt"

$VERBOSE = 3
opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} template.fasta"
  op.separator "output: aligned.txt"
  op.separator "if template.ANNOTATED.fasta, then strips leading '#' lines and writes template.fasta"
  op.separator " "
end
opts.parse!


if ARGV.size != 2
  puts opts
  exit
end

template = ARGV.shift

#seqs = seqs[184,10]
#definitions = definitions[184,10]

#PSIM = %w(%chars_to_template %letters_to_template %letters_to_self)

factory = Bio::ClustalW.new

pass_threshold = []

printv "performing pairwise alignments [on #{seqs.size} seqs]: "
seqs.zip(definitions) do |seq, df|
  if seq.to_s !~ /[^N]/i
    printv '- '
    next
  end
  align = Bio::Alignment.new([template_seq, seq])
  result = align.do_align(factory)
  (template_s, seq_s) = result.map do |seq|
    seq.to_s
  end
  psimilar = seq_s.percent_similar_to(template_s)
  printv( ("%.0f" % psimilar.last) + ' ')
  if psimilar.last > opt[:min]
    pass_threshold << [df, seq]
    #printv '*'
  else
    #printv '.'
  end
end
putsv "Done!"

abort "none found above threshold! #{opt[:min]}" if pass_threshold.size == 0
putsv "Found #{pass_threshold.size} sequence(s) above #{opt[:min]}% identical"

pass_threshold << [template_def, template_seq]

multi_align = Bio::Alignment.new( pass_threshold.map {|pair| pair.last } )
m_result = multi_align.do_align(factory).strip

labels = pass_threshold.map {|pair| pair.first }
aligned_seqs = m_result.map {|seq| seq.to_s }

File.open(outfile, 'w') do |out|
  print_align(out, aligned_seqs, labels)
end

