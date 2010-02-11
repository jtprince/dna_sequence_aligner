
require 'spec/more'
require 'align_reads_to_template'

describe 'aligning' do

  before do
    @string = 'AAAATTTTGGGGGCCCCCC'
    @conc = '--A?A-AT?TTGGGGGCCCAAC?C---'
    @testcase = "testcase.fasta"
  end

  it 'removes bad ends' do
    (start, len) = find_bad_ends(@conc, 4)
    p @conc[start, len]
    puts "OUTPUT: "
    p start
    p len
    puts "EO"
    1.is 1
  end

  it 'aligns' do
    concensus(@testcase)
    1.is 1
  end
end
