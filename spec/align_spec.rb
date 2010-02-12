
require File.dirname(__FILE__) + '/spec_helper'
require 'bio/alignment/dna_sequence'

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
  end

  it 'removes bad ends' do
    (start, len) = Bio::Alignment::DNASequenceReads.find_good_section(@conc, 4)
    @conc[start, len].is "TTGGGGGCCCAAC"
  end

  it 'aligns pairwise' do
    (template, others) = Bio::Alignment::DNASequenceReads.merge_pairwise(@pa)
    template.is @template
    @aligned.enums others
  end
end
