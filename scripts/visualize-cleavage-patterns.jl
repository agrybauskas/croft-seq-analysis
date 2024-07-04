#!/usr/bin/env julia

using ArgParse
using BioAlignments
using BioSequences
using CSVFiles, DataFrames
using Formatting
using Gumbo
using Printf
using Statistics
using PyCall

arg_parser = ArgParseSettings()

@add_arg_table arg_parser begin
    "--bed"
        help = "BED alignment file as input. Used only for the conversion from BED to HTML"
        required = true
    "--reference"
        help = "FASTA reference file"
        required = true
    "--seq"
        default = ""
        help = "Possible cleavage sequence, but only for visualization"
    "--show-gaps"
        help = "Show gaps using '|' symbol."
        action = :store_true
    "--show-next-seq"
        help = "Shows sequence next to the cleavage site."
        default = ""
    "--move-to-back"
        help = "Moves selected chromosome to the back of the list."
        action = :store_true
    "--top"
        arg_type = Int64
        help = "Regions with highest depth"
        default = 1000000
    "--no-alt-sequence"
        help = "Searches for alternative target positions and counts them."
        action = :store_true
    "--include"
        default = ""
        help = "Chromosome to be included."
    "--exclude"
        default = ""
        help = "Chromosomes to be excluded."
    "--output-format"
        default = "aln"
        help = "Output format: HTML, BED."
end

parsed_args = parse_args(ARGS, arg_parser)

function read_bed(bed::String)::Array{Dict{String,Any},1}
    output_data::Array{Dict{String,Any},1} = []
    bed_data::DataFrame =
        DataFrame(load(File{format"CSV"}(bed), delim='\t', header_exists=false))
    for bed_row::DataFrameRow{DataFrame,DataFrames.Index} in eachrow(bed_data)
        chrom::String = bed_row["Column1"]
        chrom_start::Int64 = bed_row["Column2"]
        chrom_end::Int64 = bed_row["Column3"]
        strand::String = bed_row["Column6"]
        depth_gaps_and_mms::Array{String,1} = split(bed_row["Column5"], ";")
        depth::Float64 = parse(Float64, split(depth_gaps_and_mms[1], "COV=")[2])
        gaps::String = split(depth_gaps_and_mms[3], "GAPS=")[2]
        mismatches::String = split(depth_gaps_and_mms[2], "MM=")[2]
        overall_gap_quantity::String = split(depth_gaps_and_mms[7], "OVERALL_GAP_QUANTITY=")[2]

        gap_pos::Dict{Int64,Bool} = Dict()
        for pos::String in split(split(depth_gaps_and_mms[4], "GAP_POS=")[2], ",")
            if pos != ""
                gap_pos[parse(Int64, pos)] = true
            end
        end

        ref_gap_pos::Dict{Int64,Int64} = Dict()
        ref_gap_pos_list::Array{String,1} =
            split(split(depth_gaps_and_mms[5], "REF_GAP_POS=")[2], ",")
        ref_gap_len_list::Array{String,1} =
            split(split(depth_gaps_and_mms[6], "REF_GAP_LEN=")[2], ",");
        for i::Int64 in [1:length(ref_gap_pos_list);]
            if ref_gap_pos_list[i] != ""
                ref_gap_pos[parse(Int64, ref_gap_pos_list[i])] =
                    parse(Int64, ref_gap_len_list[i])
            end
        end

        bar_len::Array{Int64,1} =
            [Base.parse(Int64, len) for len::String in split(string(bed_row["Column11"]), ",")]
        bar_pos::Array{Int64,1} =
            [Base.parse(Int64, len) for len::String in split(string(bed_row["Column12"]), ",")]

        push!(output_data,
              Dict("chrom" => chrom,
                   "chrom_start" => chrom_start,
                   "chrom_end" => chrom_end,
                   "strand" => strand,
                   "depth" => depth,
                   "gaps" => gaps,
                   "mismatches" => mismatches,
                   "gap_pos" => gap_pos,
                   "ref_gap_pos" => ref_gap_pos,
                   "bar_len" => bar_len,
                   "bar_pos" => bar_pos,
                   "overall_gap_quantity" => overall_gap_quantity))
    end

    return sort!(output_data, by=x->x["depth"], rev=true)
end

function ruller(seq::String, breaks::Int64=5, padding::String=".")::String
    seq_length::Int64 = length(seq)
    ruller::Array{Char,1} =
        ['.' for _ in [1:seq_length+length(string(seq_length))-1;]]
    for i::Int64 in [1:seq_length;]
        if i == 1 || i % breaks == 0
            for j::Int64 in [0:length(string(i))-1;]
                ruller[i+j] = string(i)[j+1]
            end
        end
    end
    return join(ruller)
end

function generate_seq_dots(seq::String, ref_seq::String)::String
    seq_dots::String = ""
    for i::Int64 in [1:length(seq);]
        if seq[i] == ref_seq[i]
            seq_dots = seq_dots * '.'
        else
            seq_dots = seq_dots * seq[i]
        end
    end
    return seq_dots
end

function generate_seq_dots_html(seq::String,
                                ref_seq::String;
                                gap_pos::Dict{Int64,Bool}=Dict(),
                                ref_gap_pos::Dict{Int64,Int64}=Dict())::String
    seq_dots_html::String = ""

    # Calculates the number of gaps that should be added.
    ref_gap_length::Int64 = 0
    if length(ref_gap_pos) > 0
        ref_gap_length = sum(map(x->ref_gap_pos[x], collect(keys(ref_gap_pos))))
    end

    # HACK: have to find more robust method than adding and substracting nt.
    if length(ref_gap_pos) > 0
        final_seq::String = ""
        for i::Int64 in [1:length(seq);]
            if haskey(ref_gap_pos, i)
                final_seq = final_seq * '-'^ref_gap_pos[i] * seq[i]
            else
                final_seq = final_seq * seq[i]
            end
        end
    else
        final_seq = seq
    end

    # Adding reference gaps.
    for i::Int64 in [1:length(ref_seq);]
        if final_seq[i] == ref_seq[i] || final_seq[i] == 'N' || ref_seq[i] == 'N'
            if length(gap_pos) > 0 && haskey(gap_pos, i) && i != 1
                seq_dots_html =
                    seq_dots_html * "<a style=\"border-left-style: solid;\">" *
                    final_seq[i] * "</a>"
            else
                seq_dots_html =
                    seq_dots_html * final_seq[i]
            end
        else
            gap_css_property::String = ""
            if length(gap_pos) > 0 && haskey(gap_pos, i) && i != 1
                gap_css_property = gap_css_property * "border-left-style: solid;"
            end

            if final_seq[i] == 'A'
                seq_dots_html =
                    seq_dots_html * "<a style=\"" * gap_css_property *
                    "background: ForestGreen;\">" * final_seq[i] * "</a>"
            elseif final_seq[i] == 'C'
                seq_dots_html =
                    seq_dots_html * "<a style=\"" * gap_css_property *
                    "background: GoldenRod;\">" * final_seq[i] * "</a>"
            elseif final_seq[i] == 'G'
                seq_dots_html =
                    seq_dots_html * "<a style=\"" * gap_css_property *
                    "background: FireBrick;\">" * final_seq[i] * "</a>"
            elseif final_seq[i] == 'T'
                seq_dots_html =
                    seq_dots_html * "<a style=\"" * gap_css_property *
                    "background: SteelBlue;\">" * final_seq[i] * "</a>"
            elseif final_seq[i] == '-'
                seq_dots_html =
                    seq_dots_html * "<a style=\"" * gap_css_property *
                    "font-weight: bold;\">" * final_seq[i] * "</a>"
            end
        end
    end

    return seq_dots_html
end

function mismatch_count(target_seq::String, ref_seq::String)::Int64
    mismatch_quantity::Int64 = 0
    for i::Int64 = [1:length(target_seq);]
        if target_seq[i] != ref_seq[i]  &&
            target_seq[i] != '-' && ref_seq[i] != '-' &&
            target_seq[i] != 'N' && ref_seq[i] != 'N'
            mismatch_quantity += 1
        end
    end
    return mismatch_quantity
end

function main()
    # Argument parsing.
    options::Dict{String,Union{String,Float64,Int64,Bool}} = Dict()
    for (key::String,val::Union{String,Float64,Int64,Bool}) in parsed_args
        options[key] = val
    end

    bed::String = options["bed"]

    reference_file::String = options["reference"]

    target_seq::String = options["seq"]
    do_show_gaps::Bool = options["show-gaps"]
    show_next_seq::String = options["show-next-seq"]
    next_seq_len::Int64 = length(show_next_seq)
    move_to_back::Bool = options["move-to-back"]
    top::Int64 = options["top"]
    no_alt_sequence::Bool = options["no-alt-sequence"]
    chr_include_str::String = options["include"]
    chr_exclude_str::String = options["exclude"]
    output_format::String = lowercase(options["output-format"])

    chr_include::Dict{String,Bool} = Dict()
    chr_exclude::Dict{String,Bool} = Dict()
    for chr in split(chr_include_str, ",")
        if chr != ""
            chr_include[chr] = true
        end
    end
    for chr in split(chr_exclude_str, ",")
        if chr != ""
            chr_exclude[chr] = true
        end
    end

    if output_format == "html"
        if no_alt_sequence
            fmt = FormatExpr("{:<6} {:<40} {} {:<18} {:<18} {:<18} {:<18} {:<18}\n")
        else
            fmt = FormatExpr("{:<6} {:<40} {} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18}\n")
        end
    elseif output_format == "bed"
        fmt = FormatExpr("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n")
    end

    output_data = []

    if output_format == "html"
        bed_data::Array{Dict{String,Any},1} = read_bed(bed)
        reference::BioSequences.FASTA.Reader =
            FASTA.Reader(Base.open(reference_file, "r"),
                         index=index="$reference_file.fai")

        ref_dic::Dict{String,String} = Dict()
        ref_len::Dict{String,Int64} = Dict()
        for ref in reference
            ref_dic[seqname(ref)] = string(DNASequence(sequence(ref)))
            ref_len[seqname(ref)] = length(ref_dic[seqname(ref)])
        end

        rank::Int64 = 1
        for bed_row::Dict{String,Any} in bed_data
            # Exits early if reaches top limit.
            if rank > top
                break
            end

            chrom_name = bed_row["chrom"]
            chrom_length::Int64 = ref_len[chrom_name]
            start_pos::Int64 = bed_row["chrom_start"]
            end_pos::Int64 = bed_row["chrom_end"]
            strand::String = bed_row["strand"]
            depth::Float64 = bed_row["depth"]
            gap_quantity::String = bed_row["gaps"]
            mismatch_quantity::Int64 = parse(Int64, bed_row["mismatches"])
            gap_pos::Dict{Int64,Bool} = bed_row["gap_pos"]
            ref_gap_pos::Dict{Int64,Int64} = bed_row["ref_gap_pos"]
            bar_pos::Array{Int64,1} = bed_row["bar_pos"]
            bar_len::Array{Int64,1} = bed_row["bar_len"]
            overall_gap_quantity::String = bed_row["overall_gap_quantity"]

            final_seq::String = ""

            # Searches for exact sequences in the reference genome.
            # TODO: should be refactored.
            alt_seq_count::Int64 = 0
            if ! no_alt_sequence
                for ref_name::String in sort(collect(keys(ref_dic)))
                    alt_seq_count +=
                        length(collect(eachmatch(Regex(string(ref_dic[chrom_name][start_pos+1:end_pos]) * "|" *
                        string(reverse_complement(DNASequence(ref_dic[chrom_name][start_pos+1:end_pos])))), ref_dic[ref_name]; overlap=true)))
                end
            end

            for i::Int64 = [1:length(bar_pos);]
                bar_start::Int64 = bar_pos[i] + 1
                bar_end::Int64 = bar_start + bar_len[i] - 1
                final_seq =
                    final_seq *
                    ref_dic[chrom_name][start_pos+bar_start:start_pos+bar_end]
            end

            # For reverse strands gap positions changes with reference to
            # a visualization.
            if strand == "-"
                final_seq = string(reverse_complement(DNASequence(final_seq)))

                # Gap positions has to be reset and recalculated for "-" strand.
                gap_pos_old::Dict{Int64,Bool} = gap_pos
                gap_pos = Dict()
                for pos::Int64 in sort(collect(keys(gap_pos_old)))
                    gap_pos[length(final_seq)-pos+2] = gap_pos_old[pos]
                end

                # Reference gap positions has to be reset and recalculated for
                # "-" strand.
                ref_gap_pos_old::Dict{Int64,Int64} = ref_gap_pos
                ref_gap_pos = Dict()
                for pos::Int64 in sort(collect(keys(ref_gap_pos_old)))
                    ref_gap_pos[length(final_seq)-pos+2] = ref_gap_pos_old[pos]
                end
            end

            # Adds sequences if --show-next-seq is ON.
            if next_seq_len > 0
                final_seq_add::String = "";
                final_seq_start::Int64 = 0;
                final_seq_end::Int64 = 0;

                if strand == "+"
                    final_seq_start = end_pos + 1
                    final_seq_end =
                        end_pos + next_seq_len > chrom_length ?
                        end_pos + (chrom_length - end_pos) :
                        end_pos + next_seq_len

                    final_seq_add =
                        string(ref_dic[chrom_name][final_seq_start:final_seq_end])
                else
                    final_seq_start =
                        start_pos - next_seq_len < 1 ?
                        1 :
                        start_pos - next_seq_len + 1
                    final_seq_end = start_pos

                    final_seq_add =
                        string(reverse_complement(DNASequence(ref_dic[chrom_name][final_seq_start:final_seq_end])))
                end

                final_seq *= final_seq_add
                mismatch_quantity += mismatch_count(final_seq_add, show_next_seq)
            end

            if final_seq == ""
                continue
            end

            seq_dots_html =
                generate_seq_dots_html(final_seq,
                                       length(show_next_seq) > 0 ?
                                           target_seq * show_next_seq :
                                           target_seq;
                                       gap_pos=gap_pos,
                                       ref_gap_pos=ref_gap_pos)

            if no_alt_sequence
                push!(output_data,
                      ("<span>" *
                       format(fmt,
                              string(rank),
                              chrom_name * ':' * string(start_pos) * '-' * string(end_pos),
                              seq_dots_html,
                              strand,
                              string(round(depth; digits=3)),
                              string(gap_quantity),
                              overall_gap_quantity,
                              string(mismatch_quantity)) *
                       "</span>"))
            else
                push!(output_data,
                      ("<span>" *
                       format(fmt, string(rank),
                              chrom_name * ':' * string(start_pos) * '-' * string(end_pos),
                              seq_dots_html,
                              strand,
                              string(round(depth; digits=3)),
                              string(gap_quantity),
                              overall_gap_quantity,
                              string(mismatch_quantity),
                              string(alt_seq_count)) *
                       "</span>"))
            end

            rank += 1
        end
    end

    if output_format == "html"
        if no_alt_sequence
            output_html =
                format(fmt, "Rank", "Position", ruller(target_seq * show_next_seq),
                       "Strand", "Depth", "Gap count", "Gap length sum", "Mismatches") *
                format(fmt, "", "Target-sequence", target_seq * show_next_seq,
                       "", "", "", "", "")
        else
            output_html =
                format(fmt, "Rank", "Position", ruller(target_seq * show_next_seq),
                       "Strand", "Depth", "Gap count", "Gap length sum",
                       "Mismatches", "Alt. seq. count") *
                format(fmt, "", "Target-sequence", target_seq * show_next_seq,
                       "", "", "", "", "", "")
        end

        for output_row in output_data
            output_html = output_html * output_row
        end

        print(parsehtml("<pre>" * output_html * "</pre>"))
    end
end

main()
