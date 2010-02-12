
require File.dirname(__FILE__) + '/spec_helper'
require 'bio/alignment/dna_sequence'

DNAReads = Bio::Alignment::DNASequenceReads

describe 'aligning' do

  before do
    @string = 'AAAATTTTGGGGGCCCCCC'
    @conc = '--A?A-AT?TTGGGGGCCCAAC?C---'
    @testcase = "testcase.fasta"

    @pa = [ ["--ABCDEFGHIJKLMNOP", 
             "-----DEFGHIJK-MN--"], 
            ["--ABCDEFGHIJKLM-NOP", 
             "--ABCDE---IJKLMZNOP"],
            ["--ABCDEFGHIJKLMNOP", 
             "-------------LMNOP"],
            ["--ABCDEFGHIJKLMNOP", 
             "--ABCDEFGHIJKLMN--"],
            ["--ABCDEFGHIJKLMNOP", 
             "--ABC------JKLM--P"],
            ["--ABC--DEFGHIJKLMNOP", 
              "--ABCZZDEFGHIJKLMNOP"],
    ]
    @template = "--ABC--DEFGHIJKLM-NO"
    @aligned = ["-------DEFGHIJK-M-N-", 
                "--ABC--DE---IJKLMZNO", 
                "---------------LM-NO", 
                "--ABC--DEFGHIJKLM-N-",
                "--ABC--------JKLM---", 
                "--ABCZZDEFGHIJKLM-NO"
                ] 

    @labels = %w(one two three four five six)

  end

  it 'removes bad ends' do
    (start, len) = DNAReads.find_good_section(@conc, 4)
    @conc[start, len].is "TTGGGGGCCCAAC"
  end

  it 'aligns pairwise' do
    (template, others) = DNAReads.merge_pairwise(@pa)
    template.is @template
    @aligned.enums others
  end

  it 'can create a good consensus string' do
    (string, stats) = DNAReads.consensus_string_and_stats([@template, *@aligned])
    string.is "  ===^^==========^=="
    stats.enums [2, 3, 15, 0, 0, 0]
    (string, stats) = DNAReads.consensus_string_and_stats([@template, "-------DEFGHIJK-M-N-"])
    string.is "  ...  ========.= =."
    stats.enums [5, 0, 10, 5, 0, 0]
  end

  xit 'prints useful printout' do
    st = StringIO.new
    DNAReads.print_align(st, @aligned, @labels, :template => @template, :template_label => "template", :chars => 8)
    puts " "
    puts st.string
    1.is 1
  end
end
