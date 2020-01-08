#include <algorithm>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/dream_index/binning_directory.hpp>
#include <seqan3/search/dream_index/concept.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include "minimizer3.h"

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;               // instead of dna5
};

struct cmd_arguments
{
    std::filesystem::path gene_file;
    std::filesystem::path exp_file;
    std::filesystem::path path_in{"./"};
    float expression{1.0};
    bool compressed = false;
    uint8_t k{20};
    uint16_t window_size{60};
    uint64_t shape;
    uint64_t seed{0x8F3F73B5CF1C9ADE};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Search through an IBF.";
    parser.add_positional_option(args.gene_file, "Please provide a fasta file.");
    parser.add_option(args.exp_file, 'x', "expression", "A tab seperated file containing expression information per "
                      "transcript given. By Default: Search for Existence in experiments",
                      seqan3::option_spec::DEFAULT,seqan3::input_file_validator{{"tsv"}});
    parser.add_option(args.path_in, 'i', "in", "Directory where input files can be found.");
    parser.add_option(args.expression, 'e', "expression_value", "Which expression level should be considered during a "
                      "search.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
    parser.add_flag(args.compressed, 'c', "compressed", "If set ibf is compressed. Default: Not compressed.");
}

template <typename number_type, typename range_type>
number_type to_number(range_type && range)
{
    std::string str;
    number_type num;
    std::ranges::copy(range, std::back_inserter(str));
    auto res = std::from_chars(&str[0], &str[0] + str.size(), num);
    if (res.ec != std::errc{})
    {
        seqan3::debug_stream << "Could not cast '" << range << "' to a valid number\n";
        throw std::invalid_argument{"CAST ERROR"};
    }
    return num;
}

template <class BD>
int do_search(BD bd, cmd_arguments & args)
{
    std::vector<uint32_t> counter;
    std::vector<uint32_t> results;
    std::vector<float> expression;
    std::vector<seqan3::dna4_vector> seqs;

    seqan3::sequence_file_input<my_traits> input_file{args.gene_file};
    for (auto & rec : input_file)
    {
        seqs.push_back(seqan3::get<seqan3::field::SEQ>(rec));
    }

    if (args.exp_file != "")
    {
        std::ifstream file{args.exp_file.string()};
        if (file.is_open())
        {
            std::string line;
            while (std::getline(file, line))
            {
                auto splitted_line = line | std::view::split('\t');
                auto it = splitted_line.begin(); // move to 1rst column
                expression.push_back(to_number<double>(*std::next(it, 1)));
            }

            if (expression.size() != seqs.size())
            {
                seqan3::debug_stream << "Error! Number of given expression levels do not match number of sequences.\n";
                return -1;
            }
        }
    }
    else
    {
        expression.assign(seqs.size(),args.expression);
    }

    std::ifstream is{args.path_in.string() + "IBF_" + std::to_string(expression[0]), std::ios::binary};
    seqan3::debug_stream << "IBF_" + std::to_string(expression[0])<< "\n";
    cereal::BinaryInputArchive iarchive{is};
    iarchive(bd);

    uint32_t minimizer_length;
    counter.resize(bd.get_bins(), 0);
    results.resize(bd.get_bins(), 0);
    for (unsigned i = 0; i < expression.size(); i++)
    {
        if ((i > 0) && (expression[i] != expression[i-1]))
        {
            std::ifstream is{"IBF_" + std::to_string(expression[i]), std::ios::binary};
            seqan3::debug_stream << "IBF_" + std::to_string(expression[i])<< "\n";
            cereal::BinaryInputArchive iarchive{is};
            iarchive(bd);
        }

        minimizer_length = 0;
        for (auto & minHash : compute_minimizer(seqs[i], args.k, args.window_size, args.shape, args.seed))
        {
            std::transform (counter.begin(), counter.end(), bd.get(minHash).begin(), counter.begin(), std::plus<int>());
            ++minimizer_length;
        }

        for(unsigned j = 0; j < counter.size(); j++)
        {
            if ( (counter[j] >= minimizer_length/2))
                results[j] = results[j] + 1;
        }
        counter.clear();
        counter.assign(bd.get_bins(), 0);
    }

    seqan3::debug_stream << "Results: " << results << "\n";
    return 0;
}

int main(int const argc, char const ** argv)
{

    seqan3::argument_parser parser("needle-search", argc, argv);
    cmd_arguments args{};
    initialize_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }

    if (args.compressed)
    {
        seqan3::binning_directory_compressed bd;
        return do_search(bd, args);
    }
    else
    {
        seqan3::binning_directory bd;
        return do_search(bd, args);
    }
}
