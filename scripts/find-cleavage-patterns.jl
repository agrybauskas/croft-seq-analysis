#!/usr/bin/env julia

using ArgParse
using BioAlignments
using BioSequences
using CSVFiles
using DataFrames
using Formatting
using Gumbo
using Printf
using Statistics
using PyCall

arg_parser = ArgParseSettings()

@add_arg_table arg_parser begin
    "--bam"
        arg_type = String
        default = ""
        help = "BAM alignment file as input"
    "--reference"
        arg_type = String
        help = "FASTA reference file"
    "--seq"
        arg_type = String
        default = ""
        help = "Possible cleavage sequence, but only for visualization"
    "--depth-change"
        arg_type = Float64
        default = 7.0
        help = "Coverage depth from background that has to be reached in order "*
               "the interval to be included"
    "--bg-include"
        action = :store_true
        help = "Includes background into calculations. Removes those targets " *
               "that exceed background error values."
    "--cleavage-ratio"
        arg_type = Float64
        default = 0.95
        help = "Ratio of reads that belong to cleaved reads. Works when " *
               "--bg-include flag is on."
    "--nt-around-cleavage"
        arg_type = Int64
        default = 20
        help = "Nucleotide count forward and backwards from cleavage position."
    "--reading-frame"
        arg_type = Int64
        default = 2
        help = "Reading frame for cleavage detection."
    "--left-bp"
        arg_type = Int64
        default = 19
        help = "Nucleotide count forwards from cleavage position."
    "--right-bp"
        arg_type = Int64
        default = 5
        help = "Nucleotide count backwards from cleavage position."
    "--ignore-R1"
        action = :store_true
        help = "Ignore R1 reads."
    "--ignore-R2"
        action = :store_true
        help = "Ignore R2 reads."
    "--only-cleavage-position"
        action = :store_true
        help = "Only cleavage positions are printed out in first three BED columns."
    "--top"
        arg_type = Int64
        default = 1000000
        help = "Regions with highest depth"
    "--include"
        arg_type = String
        default = ""
        help = "Chromosome to be included."
    "--exclude"
        arg_type = String
        default = ""
        help = "Chromosomes to be excluded."
    "--region"
        arg_type = String
        default = ""
        help = "Analyses only specified chromosome region. Used for debugging."
    "--gap-open-penalty"
        arg_type = Int64
        help = "Alignment gap-open penalty score."
        default = -20
    "--gap-extend-penalty"
        arg_type = Int64
        help = "Alignment gap extension penalty score."
        default = -2
    "--mapq-min"
        arg_type = Int64
        help = "Minimal MAPQ score for reads to be included in the analysis."
        default = 1
    "--emboss"
        action = :store_true
        help = "Use EMBOSS aligner."
    "--log-file"
        arg_type = String
        default = ""
        help = "Prints out warnings and errors to the file for debugging."
end

parsed_args = parse_args(ARGS, arg_parser)

Base.@kwdef mutable struct CleavageCandidate
    strand::String = ""
    pos::Int64 = -1
    centered_pos::Int64 = -1
    coverage::Complex{Int64} = -1
    centered_coverage::Complex{Int64} = -1
    score = nothing
    target_seq_aln::String = ""
    cleavage_seq_aln::String = ""
    cleavage_seq_start::Int64 = -1
end

#
# Reads BAM format file.
#

function read_bam(bam::String;
                  ignore_r1_reads::Bool=false,
                  ignore_r2_reads::Bool=false,
                  mapq_min::Int64=1,
                  exclude_white_reads::Bool=true,
                  region::Dict{String,Union{String,Int64}}=Dict(),
                  bg_include::Bool=false)::Tuple{Dict{String,Array{Complex{Int64},1}},Dict{String,Dict{Int64,Bool}},Dict{String,Array{Int64,1}}, Int64, Int64}
    # Reading bam.
    bam_reader::BioAlignments.BAM.Reader{IOStream} = open(BAM.Reader, bam)

    read_start_quantity::Dict{String,Array{Complex{Int64},1}} = Dict()
    bg_coverage::Dict{String,Array{Int64,1}} = Dict()
    read_start_pos::Dict{String,Dict{Int64,Bool}} = Dict()
    ref_length::Dict{String,Int64} = Dict()
    total_ref_length::Int64 = 0
    total_read_start_quantity::Int64 = 0

    bam_record::BioAlignments.BAM.Record = BAM.Record()
    while !eof(bam_reader)
        read!(bam_reader, bam_record)

        if BAM.ismapped(bam_record)
            refname::String = BAM.refname(bam_record)

            # Flag explantion:
            # 1. read paired
            # 2. read mapped in proper pair
            # 3. read unmapped
            # 4. mate unmapped
            # 5. read reverse strand
            # 6. mate reverse strand
            # 7. first in pair
            # 8. second in pair
            # 9. not primary alignment
            # 10. read fails platform/vendor quality checks
            # 11. read is PCR or optical duplicate
            # 12. supplementary alignment
            flags::Array{Int64,1} = digits(BAM.flag(bam_record), base=2, pad=12)

            # Removes R1 or R2 accordingly.
            if ignore_r1_reads && flags[7] == 1
                continue
            end

            if ignore_r2_reads && flags[7] != 1
                continue
            end

            # Including only certain region. For debugging purposes usually.
            # TODO: refactor nesting.
            if length(region) > 0
                if(refname != region["chr"])
                    continue
                end

                if BAM.position(bam_record) < region["start"] || BAM.position(bam_record) > region["end"]
                    continue
                end
            end

            if ! haskey(read_start_quantity, refname)
                ref_length[refname] = BAM.reflen(bam_record)
                read_start_quantity[refname] = zeros(ref_length[refname])
                total_ref_length += ref_length[refname]
            end
            if ! haskey(read_start_pos, refname)
                read_start_pos[refname] = Dict()
            end
            if bg_include && ! haskey(bg_coverage, refname)
                bg_coverage[refname] = zeros(ref_length[refname])
            end

            # Determines how to modify coverage from record position, direction
            # and length.
            left_pos::Int64 = BAM.position(bam_record)
            right_pos::Int64 = BAM.rightposition(bam_record)
            mapq::UInt8 = BAM.mappingquality(bam_record)

            # Filtering out white reads.
            if (mapq == 0 && exclude_white_reads) || (mapq > 0 && mapq < mapq_min)
                continue
            end

            is_reverse::Bool = false
            if flags[5] == 1
                is_reverse = true
            end

            # TODO: a bit too deep nesting.
            # Using complex numbers to differentiate forward read start from the
            # reverse one. Real part is for forward read start and complex -- for
            # reverse.
            if is_reverse
                # Normalizing reverse and forward positions.
                if right_pos < ref_length[refname]
                    read_start_quantity[refname][right_pos+1] += 1im
                    read_start_pos[refname][right_pos+1] = true
                    total_read_start_quantity += 1
                    if bg_include
                        for base_pos = [left_pos+1:right_pos+1;]
                            bg_coverage[refname][base_pos] += 1
                        end
                    end
                end
            else
                read_start_quantity[refname][left_pos] += 1
                read_start_pos[refname][left_pos] = true
                total_read_start_quantity += 1
                if bg_include
                    for base_pos = [left_pos:right_pos;]
                        bg_coverage[refname][base_pos] += 1
                    end
                end
            end
        end
    end

    close(bam_reader)

    return read_start_quantity, read_start_pos, bg_coverage,
           total_ref_length, total_read_start_quantity
end

#
# Detects regions that satisfies the accumulation of read start quantities
# (possible cleavage sites).
#

function find_cleavage_regions(read_start_pos::Dict{String,Dict{Int64,Bool}},
                               read_start_quantity::Dict{String,Array{Complex{Int64},1}};
                               reading_frame::Int64=2,
                               depth_change::Float64=4.0)::Dict{String,Array{Tuple{Int64,Int64},1}}
    # Copy has to be done here, because delete!() acts on already allocated
    # dictiory values and there is no delete().
    read_start_pos_cp::Dict{String,Dict{Int64,Bool}} = deepcopy(read_start_pos)
    cleavage_regions::Dict{String,Array{Tuple{Int64,Int64},1}} = Dict()
    for refname::String in sort(collect(keys(read_start_pos_cp)))
        ref_len::Int64 = length(read_start_quantity[refname])

        while length(read_start_pos_cp[refname]) > 0
            # For early exit from while loop.
            skip_left_pos::Bool = true
            skip_right_pos::Bool = true

            pos::Int64 = collect(keys(read_start_pos_cp[refname]))[1]

            # Checks left boundaries.
            left_boundary::Int64, right_boundary::Int64 =
                initiate_boundary_pos(pos, reading_frame, max=ref_len)

            while left_boundary >= 1 &&
                  real(sum(read_start_quantity[refname][left_boundary:right_boundary])) + imag(sum(read_start_quantity[refname][left_boundary:right_boundary])) >= depth_change
                left_boundary -= 1
                right_boundary -= 1
                skip_left_pos = false
            end

            final_left_pos::Int64 = left_boundary + 1

            # Checks right boundaries.
            # Reset boundary positions.
            left_boundary, right_boundary =
                initiate_boundary_pos(pos, reading_frame, max=ref_len)

            while right_boundary <= ref_len &&
                  real(sum(read_start_quantity[refname][left_boundary:right_boundary])) + imag(sum(read_start_quantity[refname][left_boundary:right_boundary])) >= depth_change
                left_boundary += 1
                right_boundary += 1
                skip_right_pos = false
            end

            final_right_pos::Int64 = left_boundary - 1

            if skip_left_pos && skip_right_pos
                delete!(read_start_pos_cp[refname], pos)
                continue
            end

            if !haskey(cleavage_regions, refname)
                cleavage_regions[refname] = []
            end

            push!(cleavage_regions[refname], (final_left_pos, final_right_pos))

            # Remove every position that is in between.
            for i in [final_left_pos:right_boundary-1;]
                delete!(read_start_pos_cp[refname], i)
            end
        end
    end

    return cleavage_regions
end

#
# Depending on the reading frame length (even or odd), places the boundary
# around the cleavage position.
#

function initiate_boundary_pos(pos::Int64,
                               reading_frame::Int64;
                               max::Int64=-1)::Tuple{Int64,Int64}
    # FIXME: include situations where reading frame is larger than genome.
    is_reading_frame_even::Bool = reading_frame % 2 == 0

    left_pos::Int64 = 0
    if is_reading_frame_even
        left_pos = pos - Int((reading_frame / 2))
    else
        left_pos = pos - Int((reading_frame - 1) / 2)
    end

    if left_pos < 1
        left_pos = 1
    end

    right_pos::Int64 = left_pos + reading_frame - 1

    if max != -1 && right_pos > max
        right_pos = max
    end

    return left_pos, right_pos
end

#
# Defines the center of the reading frame.
#

function center_reading_frame_pos(pos::Int64, reading_frame::Int64)::Int64
    is_reading_frame_even::Bool = reading_frame % 2 == 0
    if is_reading_frame_even
        return pos + Int((reading_frame / 2)) - 1
    else
        return pos + Int((reading_frame - 1) / 2)
    end
end

#
# Aligns supplied by CLI target sequence to the potential cleavage region.
#

function align_seq(seq1::String,
                   seq2::String,
                   scoremodel::AffineGapScoreModel{Int64};
                   log_file::String="",
                   do_emboss::Bool=false)::Tuple{String,String,Int64}
    # TODO: add types when it will be used in the future.
    if do_emboss
        # Use get(o, i-1)
        emboss_alignment = emboss_scoremodel.align(string(seq1), string(seq2))
        emboss_seq, _, emboss_ref_seq, _ = split(string(emboss_alignment[1].__str__()), "\n")
        emboss_score = floor(Int, emboss_alignment[1].score)
        # NOTE: debugging.
        if log_file != ""
            open(log_file, "a+") do f
                println(f, emboss_alignment[1].__str__())
            end
        end
        return String(emboss_seq), String(emboss_ref_seq), emboss_score
    else
        # TODO: types should be simplified.
        pair_alignment::PairwiseAlignmentResult{Int64,BioSequence{DNAAlphabet{4}},BioSequence{DNAAlphabet{4}}} =
            pairalign(GlobalAlignment(), DNASequence(seq1),
                      DNASequence(seq2), scoremodel, banded=true)
        aln::PairwiseAlignment{BioSequence{DNAAlphabet{4}},BioSequence{DNAAlphabet{4}}} =
            alignment(pair_alignment)
        alignment_score::Int64 = score(pair_alignment)

        if log_file != ""
            open(log_file, "a+") do f
                println(f, pair_alignment)
            end
        end

        # NOTE: should the type be declared to x and y?
        return String(DNASequence([x for (x, _) in aln])), String(DNASequence([y for (_, y) in aln])), alignment_score
    end
end

function update_if_better!(original_cleavage_candidate::CleavageCandidate,
                           cleavage_candidate::CleavageCandidate;
                           log_file::String="")::Nothing
    if log_file != ""
        open(log_file, "a+") do f
            println(f, "Cleavage candidate")
            println(f, "Pos: " * string(cleavage_candidate.pos))
            println(f, "Centered RF pos: " * string(cleavage_candidate.centered_pos))
            println(f, "Strand: " * cleavage_candidate.strand)
            println(f, "Score: " * string(cleavage_candidate.score))
            println(f, "Centered coverage: " * string(cleavage_candidate.centered_coverage))
            println(f, "Coverage: " * string(cleavage_candidate.coverage))
            println(f, "."^80);
        end
    end

    # Checks strandedness. Double stranded has highest priority.
    is_lower_stranded::Bool =
        original_cleavage_candidate.centered_coverage != -1 &&
        Int(real(cleavage_candidate.centered_coverage) != 0) + Int(imag(cleavage_candidate.centered_coverage) != 0) < Int(real(original_cleavage_candidate.centered_coverage) != 0) + Int(imag(original_cleavage_candidate.centered_coverage) != 0)

    if is_lower_stranded
        return nothing
    end

    # If strandedness is equal, then centered coverage is checked.
    is_equal_stranded::Bool =
        original_cleavage_candidate.centered_coverage != -1 &&
        Int(real(cleavage_candidate.centered_coverage) != 0) + Int(imag(cleavage_candidate.centered_coverage) != 0) == Int(real(original_cleavage_candidate.centered_coverage) != 0) + Int(imag(original_cleavage_candidate.centered_coverage) != 0)
    is_centered_coverage_lower::Bool =
        original_cleavage_candidate.centered_coverage != -1 &&
        real(cleavage_candidate.centered_coverage) + imag(cleavage_candidate.centered_coverage) < real(original_cleavage_candidate.centered_coverage) + imag(original_cleavage_candidate.centered_coverage)

    if is_equal_stranded && is_centered_coverage_lower
        return nothing
    end

    # If centered coverage is equal, then overall coverage is checked.
    is_centered_coverage_equal::Bool =
        original_cleavage_candidate.centered_coverage != -1 &&
        real(cleavage_candidate.centered_coverage) + imag(cleavage_candidate.centered_coverage) == real(original_cleavage_candidate.centered_coverage) + imag(original_cleavage_candidate.centered_coverage)
    is_coverage_lower::Bool =
        original_cleavage_candidate.coverage != -1 &&
        real(cleavage_candidate.coverage) + imag(cleavage_candidate.coverage) < real(original_cleavage_candidate.coverage) + imag(original_cleavage_candidate.coverage)

    if is_equal_stranded && is_centered_coverage_equal && is_coverage_lower
        return nothing
    end

    # If coverage is equal, then the alignment score is compared.
    is_coverage_equal::Bool =
        original_cleavage_candidate.coverage != -1 &&
        real(cleavage_candidate.coverage) + imag(cleavage_candidate.coverage) == real(original_cleavage_candidate.coverage) + imag(original_cleavage_candidate.coverage)
    is_score_lower_or_equal::Bool =
        original_cleavage_candidate.score != nothing &&
        cleavage_candidate.score <= original_cleavage_candidate.score

    if is_equal_stranded && is_centered_coverage_equal && is_coverage_equal && is_score_lower_or_equal
        return nothing
    end

    original_cleavage_candidate.strand = cleavage_candidate.strand
    original_cleavage_candidate.pos = cleavage_candidate.pos
    original_cleavage_candidate.centered_pos = cleavage_candidate.centered_pos
    original_cleavage_candidate.coverage = cleavage_candidate.coverage
    original_cleavage_candidate.centered_coverage = cleavage_candidate.centered_coverage
    original_cleavage_candidate.score = cleavage_candidate.score
    original_cleavage_candidate.target_seq_aln = cleavage_candidate.target_seq_aln
    original_cleavage_candidate.cleavage_seq_aln = cleavage_candidate.cleavage_seq_aln
    original_cleavage_candidate.cleavage_seq_start = cleavage_candidate.cleavage_seq_start

    return nothing
end

#
# Identifies sequence position from potential cleavage regions in the reference
# sequence.
#

function choose_best_cleavage_positions(reference_cleavage_regions::Dict{String,Array{Tuple{Int64,Int64},1}},
                                        read_start_quantity::Dict{String,Array{Complex{Int64},1}},
                                        refname::String,
                                        reading_frame::Int64,
                                        reference_seq::String,
                                        target_seq::String,
                                        scoremodel::AffineGapScoreModel{Int64};
                                        bg_coverage::Dict{String,Array{Int64,1}}=Dict(),
                                        cleavage_ratio::Float64=0.95,
                                        left_bp::Int64=length(target_seq),
                                        right_bp::Int64=length(target_seq),
                                        log_file::String="",
                                        do_emboss::Bool=false)::Dict{Int64,Dict{String,Union{Int64,String}}}
    ref_len::Int64 = length(reference_seq)
    bg_include::Bool = length(bg_coverage) > 0 ? true : false

    target_seq_rc::String =
        string(reverse_complement(DNASequence(target_seq)))

    best_cleavage_positions::Dict{Int64,Dict{String,Union{Int64,String}}} = Dict()
    for cleavage_regions::Tuple{Int64,Int64} in reference_cleavage_regions[refname]
        best_cleavage_candidate::CleavageCandidate = CleavageCandidate()

        for region_pos::Int64 in [cleavage_regions[1]:cleavage_regions[2];]
            centered_reading_frame_pos::Int64 =
                center_reading_frame_pos(region_pos, reading_frame)

            # The check on cleavage ratio against background coverage is
            # performed if the flag is on.
            if bg_include
                current_cleavage_ratio::Float64 =
                    cleavage_ratio_rf(read_start_quantity[refname],
                                      bg_coverage[refname],
                                      region_pos,
                                      reading_frame)

                if current_cleavage_ratio > 1.0
                    println("chrom: ", refname)
                    println("pos: ", region_pos)
                end

                if log_file != ""
                    open(log_file, "a+") do f
                        println(f, "Background coverage")
                        println(f, "Chromosome: " * refname)
                        println(f, "Region pos start: " * string(region_pos))
                        println(f, "Region pos end: " * string(region_pos+reading_frame-1))
                        println(f, "Region pos centered: " * string(centered_reading_frame_pos))
                        println(f, "Cleavage ratio: " * string(cleavage_ratio))
                        println(f, "Actual ratio: " * string(current_cleavage_ratio))
                        println(f, "Conditional: " * string(current_cleavage_ratio < cleavage_ratio))
                        println(f, "."^80);
                    end
                end

                if current_cleavage_ratio < cleavage_ratio
                    continue
                end
            end

            # Forward sequence alignment.
            if region_pos - left_bp + 1 > 0 && region_pos + right_bp - 1 <= ref_len
                cleavage_seq::String =
                    string(reverse_complement(DNASequence(reference_seq[region_pos-left_bp+1:region_pos+right_bp-1])))
                target_seq_aln::String, cleavage_seq_aln::String, aln_score::Int64 =
                    align_seq(target_seq_rc, cleavage_seq, scoremodel;
                              log_file=log_file, do_emboss=do_emboss)

                # Reverse complementing after alignment.
                target_seq_aln =
                    string(reverse_complement(DNASequence(target_seq_aln)))
                cleavage_seq_aln =
                    string(reverse_complement(DNASequence(cleavage_seq_aln)))

                cleavage_candidate::CleavageCandidate =
                    CleavageCandidate(strand="+",
                                      pos=region_pos,
                                      centered_pos=centered_reading_frame_pos,
                                      coverage=sum(read_start_quantity[refname][region_pos:region_pos+reading_frame-1]),
                                      centered_coverage=read_start_quantity[refname][centered_reading_frame_pos],
                                      score=aln_score,
                                      target_seq_aln=target_seq_aln,
                                      cleavage_seq_aln=cleavage_seq_aln,
                                      cleavage_seq_start=region_pos-left_bp+1)

                update_if_better!(best_cleavage_candidate, cleavage_candidate;
                                  log_file=log_file)
            end

            # Reverse sequence alignment.
            if region_pos - right_bp + 1 > 0 && region_pos + left_bp - 1 <= ref_len
                cleavage_seq_rc::String =
                    reference_seq[region_pos-right_bp+1:region_pos+left_bp-1]

                target_seq_aln_rc::String, cleavage_seq_aln_rc::String,  aln_score_rc::Int64 =
                    align_seq(target_seq_rc, cleavage_seq_rc, scoremodel;
                              log_file=log_file, do_emboss=do_emboss)

                cleavage_candidate_rc::CleavageCandidate =
                    CleavageCandidate(strand="-",
                                      pos=region_pos,
                                      centered_pos=centered_reading_frame_pos,
                                      coverage=sum(read_start_quantity[refname][region_pos:region_pos+reading_frame-1]),
                                      centered_coverage=read_start_quantity[refname][centered_reading_frame_pos],
                                      score=aln_score_rc,
                                      target_seq_aln=target_seq_aln_rc,
                                      cleavage_seq_aln=cleavage_seq_aln_rc,
                                      cleavage_seq_start=region_pos-right_bp+1)

                update_if_better!(best_cleavage_candidate, cleavage_candidate_rc;
                                  log_file=log_file)
            end
        end

        if best_cleavage_candidate.cleavage_seq_aln == "" ||
            best_cleavage_candidate.target_seq_aln == "" ||
            best_cleavage_candidate.strand == "" ||
            best_cleavage_candidate.pos == -1 ||
            best_cleavage_candidate.cleavage_seq_start == -1
            continue
        end

        best_cleavage_positions[best_cleavage_candidate.pos] =
            Dict("strand" => best_cleavage_candidate.strand,
                 "cleavage_seq_aln" => best_cleavage_candidate.cleavage_seq_aln,
                 "target_seq_aln" => best_cleavage_candidate.target_seq_aln,
                 "cleavage_seq_start" => best_cleavage_candidate.cleavage_seq_start)

        if log_file != ""
            open(log_file, "a+") do f
                println(f, "Best cleavage candidate")
                println(f, "Strand: " * best_cleavage_candidate.strand)
                println(f, "Best pos: " * string(best_cleavage_candidate.pos))
                println(f, "Best centered pos: " * string(best_cleavage_candidate.centered_pos))
                println(f, "Best score: " * string(best_cleavage_candidate.score))
                println(f, "Best centered coverage: " * string(best_cleavage_candidate.centered_coverage))
                println(f, "Best coverage: " * string(best_cleavage_candidate.coverage))
                println(f, "-"^80)
            end
        end
    end

    return best_cleavage_positions
end

function cleavage_ratio_rf(read_start_quantity::Array{Complex{Int64},1},
                           bg_coverage::Array{Int64,1},
                           pos::Int64,
                           reading_frame::Int64)::Float64
    cleavage_nucl_count::Int64 = 0
    total_nucl_count::Int64 = sum(bg_coverage[pos:pos+reading_frame-1])
    for i::Int64 in [pos:pos+reading_frame-1;]
        cleavage_nucl_count +=
            (reading_frame - (i - pos + 1) + 1) * real(read_start_quantity[i])+
            (i - pos + 1) * imag(read_start_quantity[i]);
    end
    return cleavage_nucl_count / total_nucl_count
end

function find_region_boundaries(target_seq_aln::String,
                                cleavage_seq_aln::String,
                                cleavage_seq_start::Int64)::Tuple{Int64,Int64}
    start_pos::Int64 = -1
    end_pos::Int64 = -1
    do_find_aln_start_pos::Bool = true

    for aln_pos::Int64 = [1:length(target_seq_aln);]
        # Determines gaps in target sequence alignment.
        if target_seq_aln[aln_pos] != '-'
            if do_find_aln_start_pos
                start_pos = cleavage_seq_start + aln_pos - 2
                do_find_aln_start_pos = false
            end

            end_pos = cleavage_seq_start + aln_pos - 1 - count(r"-", replace(replace(cleavage_seq_aln, r"^[-]*" => ""), r"[-]*$" => ""))
        end
    end

    return start_pos, end_pos
end

#
# Parses the alignment of the target sequence and the cleavage region that is
# suitable for visualization.
#

function parse_alignment(seq::String,
                         ref_seq::String;
                         strand::String="+",
                         log_file::String="")::Tuple{Dict{Int64,Bool},Array{Int64,1},Array{Int64,1},Int64,Int64,Dict{Int64,Int64},Int64}
    # Searches for trailing gaps and gaps in the middle of the sequences.
    seq_gaps::Array{Dict{String,Int64},1} =
        [Dict("gap_pos" => match.offset,
              "gap_len" => length(match.match))
         for match::RegexMatch in eachmatch(r"[-]+", seq)]
    ref_seq_gaps::Array{Dict{String,Int64},1} =
        [Dict("gap_pos" => match.offset,
              "gap_len" => length(match.match))
         for match::RegexMatch in eachmatch(r"[-]+", ref_seq)]

    # Calculates gaps.
    trailing_gap_start::Int64 =
        length(seq_gaps) > 0 &&
        seq_gaps[1]["gap_pos"] == 1 ? seq_gaps[1]["gap_len"] : 0
    trailing_gap_end::Int64 =
        length(seq_gaps) > 0 &&
        seq_gaps[end]["gap_pos"] + seq_gaps[end]["gap_len"] - 1 == length(seq) ?
        seq_gaps[end]["gap_len"] : 0

    # Removes trailing gaps (trimming from both sides).
    seq_trimmed =
        seq[1+trailing_gap_start:length(seq)-trailing_gap_end]
    ref_seq_trimmed =
        ref_seq[1+trailing_gap_start:length(ref_seq)-trailing_gap_end]

    # Both calculates mismatches and removes
    seq_no_gaps::String = ""
    ref_seq_no_gaps::String = ""

    mismatch_quantity::Int64 = 0
    for i::Int64 = [1:length(seq_trimmed);]
        if seq_trimmed[i] == '-' && ref_seq_trimmed[i] == '-'
            continue
        end

        if seq_trimmed[i] != '-'
            ref_seq_no_gaps *= ref_seq_trimmed[i]
        end

        if ref_seq_trimmed[i] != '-'
            seq_no_gaps *= seq_trimmed[i]
        end

        if seq_trimmed[i] != ref_seq_trimmed[i]  &&
            seq_trimmed[i] != '-' && ref_seq_trimmed[i] != '-' &&
            seq_trimmed[i] != 'N' && ref_seq_trimmed[i] != 'N'
            mismatch_quantity += 1
        end
    end

    # Finding bars.
    seq_bars::Array{Dict{String,Int64},1} =
        [Dict("bar_pos" => match.offset,
              "bar_len" => length(match.match))
         for match::RegexMatch in eachmatch(r"[ACGT]+", seq_no_gaps)]

    bar_pos::Array{Int64,1} = []
    bar_len::Array{Int64,1} = []
    for seq_bar::Dict{String,Int64} in seq_bars
        push!(bar_pos, seq_bar["bar_pos"] - 1)
        push!(bar_len, seq_bar["bar_len"])
    end

    # Redefining gap positions.
    seq_gaps =
        [Dict("gap_pos" => match.offset,
              "gap_len" => length(match.match))
         for match::RegexMatch in eachmatch(r"[-]+", seq_no_gaps)]

    # Finding gap positions.
    gap_pos::Dict{Int64,Bool} = Dict()
    for seq_gap in seq_gaps
        gap_pos[seq_gap["gap_pos"]] = true
    end

    # Redefining ref gap positions.
    ref_seq_gaps =
        [Dict("gap_pos" => match.offset,
              "gap_len" => length(match.match))
         for match::RegexMatch in eachmatch(r"[-]+", ref_seq_no_gaps)]

    ref_gap_pos::Dict{Int64,Int64} = Dict()
    for ref_seq_gap::Dict{String,Int64} in ref_seq_gaps
        ref_gap_pos[ref_seq_gap["gap_pos"]] = ref_seq_gap["gap_len"]
    end

    ref_gap_len::Array{Int64,1} =
        map((x) -> ref_gap_pos[x], sort(collect(keys(ref_gap_pos))))

    gap_quantity::Int64 = length(seq_gaps)
    overall_gap_quantity::Int64 =
        length(collect(eachmatch(r"-", seq_no_gaps))) +
        length(collect(eachmatch(r"-", ref_seq_no_gaps)))

    if log_file != ""
        final_seq::String = replace(ref_seq_no_gaps, r"-" => "")

        open(log_file, "a+") do f
            println(f, "seq:     ", seq)
            println(f, "ref_seq: ", ref_seq)
            println(f, "seq_trimmed:     ", " "^trailing_gap_start * seq_trimmed)
            println(f, "ref_seq_trimmed: ", " "^trailing_gap_start * ref_seq_trimmed)
            println(f, "final_seq: ", " "^trailing_gap_start * final_seq)
            println(f, "mismatch_quantity: ", mismatch_quantity)
            println(f, "bar_pos: ", bar_pos)
            println(f, "bar_len: ", bar_len)
            println(f, "gap_pos: ", gap_pos)
            println(f, "seq_gaps: ", seq_gaps)
            println(f, "ref_gap_pos: ", ref_gap_pos)
            println(f, "gap_quantity: ", gap_quantity)
            println(f, "overall_gap_quantity:", overall_gap_quantity)
            println(f, "trailing_gap_start: ", trailing_gap_start)
            println(f, "trailing_gap_end: ", trailing_gap_end)
            println(f, "-"^80)
        end
    end

    return gap_pos, bar_pos, bar_len, gap_quantity, mismatch_quantity,
           ref_gap_pos, overall_gap_quantity
end

function main()::Nothing
    # Argument parsing.
    options::Dict{String,Union{String,Float64,Int64,Bool}} = Dict()
    for (key::String,val::Union{String,Float64,Int64,Bool}) in parsed_args
        options[key] = val
    end

    if options["emboss"]
        emboss_align = pyimport("Bio.Align")
        emboss_scoremodel = emboss_align.PairwiseAligner()
        emboss_scoremodel.mode = "global"
        emboss_scoremodel.open_gap_score = -10.0
        emboss_scoremodel.extend_gap_score = -0.5
        emboss_scoremodel.target_end_gap_score = -10.0
        emboss_scoremodel.target_extend_gap_score = -0.5
    end

    bam::String = options["bam"]
    reference_file::String = options["reference"]
    reading_frame::Int64 = options["reading-frame"]
    target_seq::String = options["seq"]
    depth_change::Float64 = options["depth-change"]
    bg_include::Bool = options["bg-include"]
    cleavage_ratio::Float64 = options["cleavage-ratio"]
    ignore_r1::Bool = options["ignore-R1"]
    ignore_r2::Bool = options["ignore-R2"]
    top::Int64 = options["top"] # TODO: implement!
    only_cleavage_pos::Bool = options["only-cleavage-position"]
    chr_include_str::String = options["include"]
    chr_exclude_str::String = options["exclude"]
    region_str::String = options["region"]
    gap_open_penalty::Int64 = options["gap-open-penalty"]
    gap_extend_penalty::Int64 = options["gap-extend-penalty"]
    mapq_min::Int64 = options["mapq-min"]
    nt_around_cleavage::Int64 = options["nt-around-cleavage"]
    left_bp::Int64 = options["left-bp"]
    right_bp::Int64 = options["right-bp"]
    log_file::String = options["log-file"]
    do_emboss_alignment::Bool = options["emboss"]

    # Cleans up a log file if it is defined.
    if log_file != ""
        open(log_file, "w") do f
            println(f,"")
        end
    end

    # Searching around cleavage is important, because it has to be taken into
    # account that both '+' and '-' target sequence alignments are possible.
    if left_bp > 0 && right_bp > 0
        if left_bp >= right_bp
            nt_around_cleavage = left_bp
        else
            nt_around_cleavage = right_bp
        end
    end

    chr_include::Dict{String,Bool} = Dict()
    for chr::String in split(chr_include_str, ",")
        if chr != ""
            chr_include[chr] = true
        end
    end
    chr_exclude::Dict{String,Bool} = Dict()
    for chr::String in split(chr_exclude_str, ",")
        if chr != ""
            chr_exclude[chr] = true
        end
    end

    region::Dict{String,Union{String,Int64}} = Dict()
    if region_str != ""
        region_chr::String = split(region_str, ":")[1]
        region_start::Int64 = parse(Int64, split(split(region_str, ":")[2], "-")[1])
        region_end::Int64 = parse(Int64, split(split(region_str, ":")[2], "-")[2])
        region = Dict("chr" => region_chr,
                      "start" => region_start,
                      "end" => region_end)
    end

    fmt::FormatExpr =
        FormatExpr("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n")

    if bam != ""
        read_start_quantity::Dict{String,Array{Complex{Int64},1}},
        read_start_pos::Dict{String,Dict{Int64,Bool}},
        bg_coverage::Dict{String,Array{Int64,1}},
        total_ref_length::Int64,
        total_read_start_quantity::Int64 =
            read_bam(bam, ignore_r1_reads=ignore_r1, ignore_r2_reads=ignore_r2,
                     mapq_min=mapq_min, region=region, bg_include=bg_include)

        reference::BioSequences.FASTA.Reader =
            FASTA.Reader(Base.open(reference_file, "r"),
                         index=index="$reference_file.fai")
        scoremodel::AffineGapScoreModel{Int64} =
            AffineGapScoreModel(EDNAFULL, match=5, mismatch=-4,
                                gap_open=gap_open_penalty,
                                gap_extend=gap_extend_penalty)

        cleavage_regions::Dict{String,Array{Tuple{Int64,Int64},1}} =
            find_cleavage_regions(read_start_pos, read_start_quantity,
                                  reading_frame=reading_frame,
                                  depth_change=depth_change)

        # Clearing memory.
        read_start_pos = Dict()

        # Calculates ratio of read starts per whole reference length. It will
        # be useful when determining depth change automaticaly.
        read_start_ratio_per_ref::Float64 =
            total_read_start_quantity / total_ref_length;

        output_data::Array{Tuple{Int64,Any},1} = []

        # Cleavage regions data parsing.
        region_counter::Int64 = 1
        for refname::String in sort(collect(keys(cleavage_regions)))
            reference_seq::String = string(sequence(reference[refname]))
            ref_len::Int64 = length(reference_seq)

            # Including/excluding chromosomes.
            if length(chr_include) > 0 && !haskey(chr_include, refname)
                continue
            end
            if length(chr_exclude) > 0 && haskey(chr_exclude, refname)
                continue
            end

            # Best cleavage positions are defined as these positions that have
            # the potential target sequences that are analogous to target with
            # the best alignment score within the cleavage region.
            # TODO: refactor function, because it has too many arguments.
            cleavage_positions::Dict{Int64,Dict{String,Union{Int64,String}}} =
                choose_best_cleavage_positions(cleavage_regions,
                                               read_start_quantity,
                                               refname,
                                               reading_frame,
                                               reference_seq,
                                               target_seq,
                                               scoremodel,
                                               bg_coverage=bg_coverage,
                                               cleavage_ratio=cleavage_ratio,
                                               left_bp=left_bp,
                                               right_bp=right_bp,
                                               log_file=log_file)

            # Extracts additional information on the alignment that will be used
            # later.
            for pos::Int64 in sort(collect(keys(cleavage_positions)))
                gap_pos::Dict{Int64,Bool},
                bar_pos::Array{Int64,1},
                bar_len::Array{Int64,1},
                gap_quantity::Int64,
                mismatch_quantity::Int64,
                ref_gap_pos::Dict{Int64,Int64},
                overall_gap_quantity::Int64 =
                    parse_alignment(
                        cleavage_positions[pos]["target_seq_aln"],
                        cleavage_positions[pos]["cleavage_seq_aln"],
                        strand=cleavage_positions[pos]["strand"],
                        log_file=log_file
                )

                start_pos::Int64, end_pos::Int64 =
                    find_region_boundaries(
                        cleavage_positions[pos]["target_seq_aln"],
                        cleavage_positions[pos]["cleavage_seq_aln"],
                        cleavage_positions[pos]["cleavage_seq_start"]
                )

                coverage::Int64 =
                    real(sum(read_start_quantity[refname][pos:pos+reading_frame-1])) +
                    imag(sum(read_start_quantity[refname][pos:pos+reading_frame-1]))
                strand::String = cleavage_positions[pos]["strand"]

                if ! only_cleavage_pos
                    push!(output_data, (
                        coverage,
                        format(fmt,
                               refname,
                               string(start_pos),
                               string(end_pos),
                               "region_" * string(region_counter),
                               "COV=" * string(coverage) *
                               ";MM=" * string(mismatch_quantity) *
                               ";GAPS=" * string(gap_quantity) *
                               ";GAP_POS=" * join(sort(collect(keys(gap_pos))), ",") *
                               ";REF_GAP_POS=" * join(sort(collect(keys(ref_gap_pos))), ",") *
                               ";REF_GAP_LEN=" * join(map(x->ref_gap_pos[x], sort(collect(keys(ref_gap_pos)))), ",") *
                               ";OVERALL_GAP_QUANTITY=" * string(overall_gap_quantity),
                               strand,
                               string(start_pos),
                               string(end_pos),
                               "120,120,120",
                               string(length(bar_pos)),
                               join(bar_len, ","),
                               join(bar_pos, ",")))
                    )
                elseif only_cleavage_pos
                    push!(output_data, (
                        coverage,
                        format(fmt,
                               refname,
                               string(start_pos),
                               string(start_pos+1),
                               "region_" * string(region_counter),
                               "COV=" * string(coverage) *
                               ";MM=" * string(mismatch_quantity) *
                               ";GAPS=" * string(gap_quantity) *
                               ";GAP_POS=" * join(sort(collect(keys(gap_pos))), ",") *
                               ";REF_GAP_POS=" * join(sort(collect(keys(ref_gap_pos))), ",") *
                               ";REF_GAP_LEN=" * join(map(x->ref_gap_pos[x], sort(collect(keys(ref_gap_pos)))), ",") *
                               ";OVERALL_GAP_QUANTITY=" * string(overall_gap_quantity),
                               strand,
                               string(start_pos),
                               string(end_pos),
                               "120,120,120",
                               string(length(bar_pos)),
                               join(bar_len, ","),
                               join(bar_pos, ",")))
                    )
                end

                region_counter += 1
            end
        end
    end

    sort!(output_data, by=x->x[1], rev=true)

    for output_row::Tuple{Int64,String} in output_data
        print(output_row[2])
    end

    return nothing
end

main()
