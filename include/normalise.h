#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/trim.hpp>
#include <seqan3/std/ranges>

using seqan3::operator""_dna4;
seqan3::phred42 average{};
uint16_t min_length{};

//!\brief Use dna4 instead of default dna5
struct traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};


// std::views::filter takes a function object (a lambda in this case) as input that returns a boolean
auto minimum_quality_filter = std::views::filter([] (auto const & rec)
{
    auto qual = seqan3::get<seqan3::field::qual>(rec) | std::views::transform([] (auto q) { return q.to_phred(); });
    double sum = std::accumulate(qual.begin(), qual.end(), 0);
    return ((sum / std::ranges::size(qual) >= average) &
            (std::ranges::size(seqan3::get<seqan3::field::seq>(rec)) >= min_length));
});

void normalise(std::vector<std::filesystem::path> sequence_files, std::filesystem::path path_out, seqan3::phred42 avg,
               uint16_t min_len, seqan3::phred42 quality)
{
    average = avg;
    min_length = min_len;
    for (auto & file : sequence_files )
    {
        seqan3::sequence_file_output fout{path_out.string() + "Normalised_" + file.string()};
        seqan3::sequence_file_input<traits> fin{file};
        using record_type = decltype(fin)::record_type;
        std::vector<record_type> records{};
        for (auto & rec : fin | minimum_quality_filter)
        {
            auto vec = rec | seqan3::views::trim(quality);
            records.push_back(std::move(rec));
        }

        for (auto & rec : records)
        {
            fout.push_back(rec);
        }
    }
}
// | seqan3::views::trim(quality)
