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
        help = "BED alignment files as input."
        required = true
    "--top"
        arg_type = Int
        help = "Regions with highest depth."
        default = 1000000
    "--include"
        default = ""
        help = "Chromosome to be included."
    "--exclude"
        default = ""
        help = "Chromosomes to be excluded."
    "--operator"
        default = "AND"
        help = "Operator used to combine BED files. OR, AND or NOT can be chosen."
end

parsed_args = parse_args(ARGS, arg_parser)

# Argument parsing.
options = Dict()
for (key,val) in parsed_args
    options[key] = val
end

function read_bed(bed)
    output_data = []
    bed_data = DataFrame(load(File{format"CSV"}(bed), delim='\t', header_exists=false))
    for bed_row in eachrow(bed_data)
        chrom = bed_row["Column1"]
        chrom_start = bed_row["Column2"]
        chrom_end = bed_row["Column3"]
        strand = bed_row["Column6"]
        depth_gaps_and_mms = split(bed_row["Column5"], ";")
        depth = parse(Float64, split(depth_gaps_and_mms[1], "COV=")[2])
        gaps = split(depth_gaps_and_mms[3], "GAPS=")[2]
        mismatches = split(depth_gaps_and_mms[2], "MM=")[2]
        overall_gap_quantity = split(depth_gaps_and_mms[7], "OVERALL_GAP_QUANTITY=")[2]

        gap_pos = Dict()
        for pos in split(split(depth_gaps_and_mms[4], "GAP_POS=")[2], ",")
            if pos != ""
                gap_pos[parse(Int, pos)] = true
            end
        end

        ref_gap_pos = Dict()
        ref_gap_pos_list = split(split(depth_gaps_and_mms[5], "REF_GAP_POS=")[2], ",")
        ref_gap_len_list = split(split(depth_gaps_and_mms[6], "REF_GAP_LEN=")[2], ",");
        for i in [1:length(ref_gap_pos_list);]
            if ref_gap_pos_list[i] != ""
                ref_gap_pos[parse(Int, ref_gap_pos_list[i])] = parse(Int, ref_gap_len_list[i])
            end
        end

        bar_len = [Base.parse(Int, len) for len in split(string(bed_row["Column11"]), ",")]
        bar_pos = [Base.parse(Int, len) for len in split(string(bed_row["Column12"]), ",")]
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

    return output_data
end

function join_bed(bed_datas, fmt; operator="and")
    # Key is generated that is composed of chromosome, strand, start and end
    # positions in order to quicly find visited cleavage sites.
    combined_visited_bed_data = Dict()
    output_data = []

    bed_counter = 1
    for bed_data in bed_datas
        for bed_row in bed_data
            chr = bed_row["chrom"]
            strand = bed_row["strand"]
            start_pos = bed_row["chrom_start"]
            end_pos = bed_row["chrom_end"]
            bar_pos = bed_row["bar_pos"]
            bar_len = bed_row["bar_len"]

            key = join([chr, strand, start_pos, end_pos, join(keys(bar_pos), ","), join(bar_len, ",")], ";" )

            if haskey(combined_visited_bed_data, key)
                combined_visited_bed_data[key] = Dict("mean_coverage_sum" => combined_visited_bed_data[key]["mean_coverage_sum"] + bed_row["depth"],
                                                      "bed_counter" => combined_visited_bed_data[key]["bed_counter"] + 1,
                                                      "bed_row" => bed_row,
                                                      "bed_first_idx" => combined_visited_bed_data[key]["bed_first_idx"])
            else
                combined_visited_bed_data[key] = Dict("mean_coverage_sum" => bed_row["depth"],
                                                      "bed_counter" => 1,
                                                      "bed_row" => bed_row,
                                                      "bed_first_idx" => bed_counter)
            end
        end

        bed_counter += 1
    end

    region_counter = 1
    for bed_key in collect(keys(combined_visited_bed_data))
        bed_counter = combined_visited_bed_data[bed_key]["bed_counter"]
        # A bit weird logic, but it should work. By identifying where the first
        # region occurred and by tracking the overall quantity of occurences, it
        # is possible to remove regions.
        bed_first_idx = combined_visited_bed_data[bed_key]["bed_first_idx"]

        if operator == "and"
            if bed_counter < length(bed_datas)
                continue
            end
        end

        if operator == "not" && ! (bed_first_idx == 1 && bed_counter == 1)
            continue
        end

        bed_row = combined_visited_bed_data[bed_key]["bed_row"]
        mean_coverage_sum = combined_visited_bed_data[bed_key]["mean_coverage_sum"]
        chr = bed_row["chrom"]
        strand = bed_row["strand"]
        start_pos = bed_row["chrom_start"]
        end_pos = bed_row["chrom_end"]
        gap_quantity = bed_row["gaps"]
        mismatch_quantity = bed_row["mismatches"]
        gap_pos = bed_row["gap_pos"]
        ref_gap_pos = bed_row["ref_gap_pos"]
        overall_gap_quantity = bed_row["overall_gap_quantity"]
        bar_pos = bed_row["bar_pos"]
        bar_len = bed_row["bar_len"]
        mean_coverage = mean_coverage_sum / bed_counter

        push!(output_data, (mean_coverage, format(fmt, chr, string(start_pos), string(end_pos), "region_" * string(region_counter), "COV=" * string(mean_coverage) * ";MM=" * string(mismatch_quantity) * ";GAPS=" * string(gap_quantity) * ";GAP_POS=" * join(sort(collect(keys(gap_pos))), ",") * ";REF_GAP_POS=" * join(sort(collect(keys(ref_gap_pos))), ",") * ";REF_GAP_LEN=" * join(map(x->ref_gap_pos[x], sort(collect(keys(ref_gap_pos)))), ",") * ";OVERALL_GAP_QUANTITY=" * string(overall_gap_quantity), strand, string(start_pos), string(end_pos), "120,120,120", string(length(bar_pos)), join(bar_len, ","), join(bar_pos, ","))))

        region_counter += 1
    end

    return output_data
end

function main()
    # Argument parsing.
    options = Dict()
    for (key,val) in parsed_args
        options[key] = val
    end

    beds = options["bed"]
    top = options["top"]
    chr_include_str = options["include"]
    chr_exclude_str = options["exclude"]
    operator = lowercase(options["operator"])

    chr_include = Dict()
    chr_exclude = Dict()
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

    fmt = FormatExpr("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n")

    bed_datas = []

    for bed in split(beds, ",")
        push!(bed_datas, read_bed(bed))
    end

    output_data = join_bed(bed_datas, fmt; operator=operator)

    sort!(output_data, by=x->x[1], rev=true) # Probably do not have to sort.

    for output_row in output_data
        print(output_row[2])
    end
end

main()
