#!/usr/bin/ruby

outfile = "ANALYZE.FASTA"

if ARGV.size == 0
  puts "usage: #{File.basename(__FILE__)} <file>.fasta ..."
  puts "comments (starting with '#') are ok"
  puts "outputs: #{outfile}"
  exit
end

all_text = ARGV.map do |file|
  IO.read(file).split("\n").reject {|line| line =~ /^\#/ }.join("\n")
end.join("\n")

File.open(outfile, 'w') do |out|
  out.print all_text
end

