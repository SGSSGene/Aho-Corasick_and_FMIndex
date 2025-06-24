#include "SearchHammingAAA.h"

#include <channel/channel.h>
#include <clice/clice.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/fmindex/diskStorage.h>
#include <fmindex-collection/search/SearchHammingSM.h>
#include <fmt/format.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <ranges>
#include <thread>

namespace {
auto cliHelp = clice::Argument { .args     = {"-h", "--help"},
                                 .desc     = "prints the help page"
};

auto cliReferenceFile = clice::Argument { .args = {"-r", "--ref", "--reference"},
                                          .desc = "reference file (must be fasta file)",
                                          .value = std::string{},
                                          .tags = {"required"},
};
auto cliInputFile = clice::Argument { .args = {"-i", "--input"},
                                      .desc = "input file, must be a fasta file",
                                      .value = std::string{},
                                      .tags = {"required"},
};
auto cliErrors = clice::Argument { .args = {"-k", "--errors"},
                                   .desc = "number of errors that are allowed during the search",
                                   .value = size_t{0},
};

auto cliAAA = clice::Argument { .args = {"--aaa"},
                                .desc = "number of allowed ambiguous amino acids",
                                .value = size_t{0},
};

auto cliThreads = clice::Argument { .args = {"-t", "--threads"},
                                    .desc = "number of threads",
                                    .value = size_t{1},
};

auto cliCountOnly = clice::Argument { .args = {"--count-only"},
                                      .desc = "counts the number of results, but doesn't write anything to a file"
};
}

struct Timer {
  std::chrono::high_resolution_clock::time_point start_;

    Timer() : start_(std::chrono::high_resolution_clock::now()) {}
    void reset() { start_ = std::chrono::high_resolution_clock::now(); }
    double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start_).count();
    }
};


//!TODO quick hack, to get aminoacid, scoring scheme in here
template <typename index_t, fmc::Sequence query_t, typename callback_t>
void search(index_t const& _index, query_t const& _query, size_t _errorsMM, size_t _errorsAAA, callback_t _callback) {
    using cursor_t = fmc::select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    static thread_local auto cache = std::tuple<size_t, size_t, fmc::search_scheme::Scheme, fmc::search_scheme::Scheme>{std::numeric_limits<size_t>::max(), 0, {}, {}};
    // check if last scheme has correct errors and length, other wise generate it
    auto& [error, length, ss, search_scheme] = cache;
    if (error != _errorsMM + _errorsAAA) { // regenerate everything
        ss            = fmc::search_scheme::generator::h2(_errorsMM+_errorsAAA+2, 0, _errorsMM + _errorsAAA);
        length        = _query.size();
        search_scheme = limitToHamming(fmc::search_scheme::expand(ss, length));
    } else if (length != _query.size()) {
        length        = _query.size();
        search_scheme = limitToHamming(fmc::search_scheme::expand(ss, length));
    }

    static thread_local auto scoringMatrix = [&]() {
        // by default the diagonal is set to 0 and the rest to 1
        using Alphabet = ivs::delimited_alphabet<ivs::aa27>;
        static_assert(Alphabet::size() == index_t::Sigma, "This is a hack, Alphabet should always be equal to the one used in the index");

        auto matrix = fmc::search_hamming_aaa::ScoringMatrix<index_t::Sigma>{};

        for (auto base : Alphabet::ambiguous_bases()) {
            for (auto alt : Alphabet::base_alternatives(base)) {
                matrix.setCost(base, alt, 0);
                matrix.setCost(alt, base, 0);
            }
        }

        return matrix;
    }();


    fmc::search_hamming_aaa::search(_index, _query, _errorsMM, _errorsAAA, search_scheme, scoringMatrix, _callback);
}


template <typename Alphabet>
auto loadFastaFile(std::filesystem::path const& path) -> std::vector<std::vector<uint8_t>> {
    auto res = std::vector<std::vector<uint8_t>>{};
    for (auto record : ivio::fasta::reader {{.input = path}}) {
        res.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));

        if (auto pos = ivs::verify_rank(res.back()); pos) {
            throw std::runtime_error{fmt::format("input file has unexpected letter in record {}: {} at position {}", res.size(), record.id, *pos)};
        }
    }
    return res;

}

template <typename Alphabet>
auto loadOrConstructIndex(std::filesystem::path const& _inputFile, std::vector<std::vector<uint8_t>> const& ref, size_t _nbrThreads) {
    using Index = fmc::BiFMIndex<Alphabet::size()>;

    auto indexPath = std::filesystem::path{_inputFile.string() + ".idx"};
    // try to load index from file
    if (std::filesystem::exists(indexPath)) {
        fmt::print("loading index from file\n");
        return fmc::loadIndex<Index>(indexPath);
    }

    // create index
    auto index = Index{ref, /*.samplingRate=*/16, /*.threadNbr=*/_nbrThreads};
    fmt::print("created index\n");

    // store index
    fmt::print("saving index to disk\n");
    saveIndex(index, indexPath);
    return index;
}

static void app() {
    // define Alphabet and Index type
    using Alphabet = ivs::delimited_alphabet<ivs::aa27>;
    auto totalTimeTimer = Timer{};
    auto timer = Timer{};

    // load index from fasta file
    auto ref = loadFastaFile<Alphabet>(*cliReferenceFile);
    fmt::print("loaded reference with {} entries in {} seconds\n", ref.size(), timer.elapsed());
    timer.reset();

    // load index, either load from file, or create from fasta file
    auto index = loadOrConstructIndex<Alphabet>(*cliReferenceFile, ref, *cliThreads);
    fmt::print("loaded or constructed index in {} seconds\n", timer.elapsed());
    timer.reset();

    auto queries = loadFastaFile<Alphabet>(*cliInputFile);
    fmt::print("loaded {} queries in {} seconds\n", queries.size(), timer.elapsed());
    timer.reset();

    // search each entry
    fmt::print("executing search\n");

    std::atomic_size_t totalHits{};
    std::mutex mutex;
    size_t nextQuery{};

    // run searches in parallel
    auto workers = channel::workers{*cliThreads, [&]() {
        size_t hitsPerThread{};
        auto output = std::stringstream{};
        do {
            size_t queryId{};
            {
                auto g = std::lock_guard{mutex};
                if (nextQuery == queries.size()) break;
                queryId = nextQuery;
                ++nextQuery;
            }
            auto& query = queries[queryId];

            // run the actual search, results are report by a call to the given lambda
            // uses a scoring matrix
            search(index, query, *cliErrors, *cliAAA, [&](auto cursor, size_t error) {
                for (auto sa_entry : cursor) {
                    hitsPerThread += 1;
                    auto [entry, offset] = index.locate(sa_entry);
                    auto [refId, refPos] = entry;
                    if (!cliCountOnly) {
                        auto g = std::lock_guard{mutex};
                        fmt::print("Found needle #{}({}) at pos {} in protein {}\n", nextQuery, ivs::view_rank_to_char<Alphabet>(query), refPos+offset, ivs::view_rank_to_char<Alphabet>(ref[refId]));
                    }
                }
            });
        } while(true);
        totalHits += hitsPerThread;
    }};
    workers.join();

    fmt::print("found {} hits in {} seconds (total time: {} seconds)\n", totalHits.load(), timer.elapsed(), totalTimeTimer.elapsed());
}

int main(int argc, char** argv) {
    clice::parse({
        .args = {argc, argv},
        .allowDashCombi  = true, // default false, -a -b -> -ab
        .helpOpt         = true, // default false, registers a --help option and generates help page
        .catchExceptions = true, // default false, catches exceptions and prints them to the command line and exists with code 1
        .run = app,
    });
}
