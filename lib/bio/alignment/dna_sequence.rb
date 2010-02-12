

module Bio
  module Alignment
    module DNASequenceReads

      module_function
      CLUSTALW_OPTS = %w(gapopen gapext dnamatrix type)

      # returns the index of the starting run of good chars
      def find_start_good_section(iupac_concensus_string, min_length)
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
        start = find_start_good_section(iupac_concensus_string, min_length)
        from_end = find_start_good_section(iupac_concensus_string.reverse, min_length)
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

      ALIGN_OPTS = {
        :type => 'DNA',
        :gapopen => 20,
        :gapext => 20,
        :dnamatrix => 'IUB',  # "IUB" || "CLUSTALW"
        :fidelity_length => 10,
      }

      # returns high quality pairwise alignments
      # based on the fidelity_length option
      def align_pairwise(bioseqs, opt={})
        factory = Bio::ClustalW.new
        clustal_opts = hash_opts_to_clustalopts(opt)
        factory.options = clustal_opts
        template = bioseqs.shift
        start_length = []
        pairwise_aligns = bioseqs.map do |bseq|
          clust_al = clustal_align([template, bseq], factory)
          cl_cons = clust_al.consensus
          aligned_string = clust_al[1].to_s
          (st, len) = find_good_section(aligned_string, opt[:fidelity_length])
          pristine = aligned_string[st, len].gsub('-','')  # pristine read (ends removed)
          clustal_align([template.to_s, Bio::Sequence::NA.new(pristine)], factory)
        end
      end

      # assumes all were aligned to the same template (the first of a pair)
      def merge_pairwise(aligns)
        ps = aligns.map do |align| 
          seqs = []
          align.each do |bioseq|
            seqs << bioseq.to_s
          end
          seqs
        end
        template = []
        #m,x,n
        x = 2
        ftemp = ps.first.first
        nmax = ps.map {|pair| pair.first.size }.max
        mmax = ps.size
        mar = (0...mmax).to_a
        others = mar.map { [] }
        ns = mar.map { 0 }
        tn = 0
        on = 0
        (0...nmax).each do |n|
          (t_dsh, t_no_dsh) = mar.partition do |m| 
            # this is RUBY 1.8 ONLY!!
            ps[m][0][ns[m]] == 45  # '-' is ascii 45
          end

          # if a template has a dash, all other off-templates need a dash
          if t_dsh.size > 0
            template[tn] = 45
            t_no_dsh.each do |m|
              # don't update these guys counter
              others[m][tn] = 45
            end
            t_dsh.each do |m|
              others[m][tn] = ps[m][1][ns[m]]
              ns[m] += 1
            end
          else # no dashes in the template
            t_no_dsh.each do |m|
              others[m][tn] = ps[m][1][ns[m]]
            end
            template[tn] = ps[0][0][ns[0]]
            ns.map!{|v| v+1 } 
          end
          tn += 1
        end
        [cs_to_s(template), others.map! {|ar| cs_to_s(ar) } ]
      end

      def cs_to_s(ar)
        ar.map {|v| v.nil? ? '-' : v.chr }.join
      end

      # adjust all pairwise alignments to fit each other

      #consensus_template = []
      #max_length = pairs_of_strings.map {|pair| pair.first.size }.max
      #(0...max_length).each do |n|
      #  pairs_of_strings.map {|pair| pair.map {|st| st[n] } }
      #end


      def exactly_chars(string, n)
        at_least = "%#{n}s" % self
        at_least[0,n]
      end

      def print_align(io, sequences, labels, opts={})
        opts = {:cutoff => 100, :start => 0, :chars => 20}.merge(opts)
        (start, length, chars) = opts.values_at(:start, :cutoff, :chars)
        spacer = "  "

        loop do
          fin = false

          ## the counters at the top of the line
          start_s = start.to_s
          finish_s = (start + length).to_s
          count_line_gap = length - (start_s.size + finish_s.size)

          count_line = [start_s, " " * count_line_gap, finish_s].join
          io.puts [exactly_chars("", chars), spacer, count_line].join

          sequences.zip(labels) do |string, label|
            fin = (start >= string.length )
            break if fin
            io.puts "#{exactly_chars(label, chars)}#{spacer}#{string[start,length]}"
          end
          io.puts " "
          break if fin
          start += length
        end
      end

    end
  end
end
