#include "SearchHammingAAA.h"

#include <channel/value_mutex.h>
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

namespace fmc = fmindex_collection;


namespace {
auto cliHelp = clice::Argument { .args     = {"-h", "--help"},
                                 .desc     = "prints the help page"
};

auto cliReferenceFile = clice::Argument { .args = {"-r", "--ref", "--reference"},
                                          .desc = "reference file (must be fasta file)",
                                          .value = std::string{},
                                          .tags = {"required"},
};
auto cliNoIndex = clice::Argument { .args = {"--noindex"},
                                    .desc = "neither read nor write the reference index file"
};
auto cliInputFile = clice::Argument { .args = {"-i", "--input"},
                                      .desc = "input file, must be a fasta file",
                                      .value = std::string{},
                                      .tags = {"required"},
};
auto cliOutputFile = clice::Argument { .args = {"-o", "--output"},
                                      .desc = "output file",
                                      .value = std::string{}
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
auto loadOrConstructIndex(std::filesystem::path const& _inputFile, bool _skipLoadingAndSaving, size_t _nbrThreads) {
    using String = fmc::string::FlattenedBitvectors_512_64k<Alphabet::size()>;
    using Index = fmc::BiFMIndex<String>;

    auto indexPath = std::filesystem::path{_inputFile.string() + ".idx"};
    // try to load index from file
    if (!_skipLoadingAndSaving && std::filesystem::exists(indexPath)) {
        fmt::print("loading index from file\n");
        return fmc::loadIndex<Index>(indexPath);
    }

    // load index from fasta file
    auto ref = std::vector<std::vector<uint8_t>>{};
    for (auto record : ivio::fasta::reader {{.input = _inputFile}}) {
        ref.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));

        if (auto pos = ivs::verify_rank(ref.back()); pos) {
            throw std::runtime_error{fmt::format("input file has unexpected letter in record {}: {} at position {}", ref.size(), record.id, *pos)};
        }
    }
    fmt::print("loaded reference with {} entries\n", ref.size());

    // create index
    auto index = Index{ref, /*.samplingRate=*/16, /*.threadNbr=*/_nbrThreads};
    fmt::print("created index\n");

    // store index
    if (!_skipLoadingAndSaving) {
        fmt::print("saving index to disk\n");
        saveIndex(index, indexPath);
    }
    return index;
}

int main(int argc, char** argv) {
    try {
        // parse and run clice commands
        if (auto failed = clice::parse(argc, argv); failed) {
            fmt::print(stderr, "parsing failed {}\n", *failed);
            return 1;
        }
        if (auto ptr = std::getenv("CLICE_COMPLETION"); ptr) {
            return 0;
        }

        // print help if requested
        if (cliHelp) {
            fmt::print("{}", clice::generateHelp());
            return 0;
        }

        // check number of threads are valid
        if (*cliThreads < 1) {
            fmt::print("invalid number of threads given: {}\n", *cliThreads);
            return 1;
        }

        // define Alphabet and Index type
        using Alphabet = ivs::delimited_alphabet<ivs::aa27>;

        // load index, either load from file, or create from fasta file
        auto index = loadOrConstructIndex<Alphabet>(*cliReferenceFile, cliNoIndex, *cliThreads);


        // load index from fasta file
        auto ref = std::vector<std::vector<uint8_t>>{};
        auto ref_as_str = std::vector<std::string>{};
        for (auto record : ivio::fasta::reader {{.input = *cliReferenceFile}}) {
            ref_as_str.emplace_back(record.seq);
            ref.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));

            if (auto pos = ivs::verify_rank(ref.back()); pos) {
                throw std::runtime_error{fmt::format("input file has unexpected letter in record {}: {} at position {}", ref.size(), record.id, *pos)};
            }
        }
        fmt::print("loaded reference with {} entries\n", ref.size());


        fmt::print("loading queries\n");

        //!TODO not working with gcc12
        //auto queries = std::ranges::to<std::vector>(ivio::fasta::reader{{.input = *cliInputFile}});
        auto queries = std::vector<ivio::fasta::record>{};
        for (auto record_view : ivio::fasta::reader{{.input = *cliInputFile}}) {
            queries.emplace_back(record_view);
        }

        // run searches in parallel
        auto outputBuffers = std::vector<std::stringstream>{};
        outputBuffers.resize(*cliThreads);

        // search each entry
        fmt::print("executing search\n");


        // Number of queries occupied by a thread
        size_t chunk_size = 10;

        auto lastProcessedQuery = channel::value_mutex<size_t>{0};

        std::atomic_size_t totalHits{};

        auto workers = std::vector<std::jthread>{};
        for (size_t j{0}; j < *cliThreads; ++j) {
            workers.emplace_back([&, j]() {
                size_t hitsPerThread{};
                do {
                    // lock query and fetch processing range
                    auto [g, v] = *lastProcessedQuery;
                    if (*v == queries.size()) {
                        break;
                    }
                    size_t startId = *v;
                    size_t endId = startId + chunk_size;
                    if (endId > queries.size()) {
                        endId = queries.size();
                    }
                    *v = endId;
//                    fmt::print("starting {}-{}/{} ({}%)\n", startId, endId, queries.size(), startId*100./queries.size());
                    g.unlock();

                    // process search
                    for (size_t i{startId}; i < endId; ++i) {
                        auto const& record = queries[i];


                        auto query = ivs::convert_char_to_rank<Alphabet>(record.seq);
                        if (auto pos = ivs::verify_rank(query); pos) {
                            fmt::print("skipping record {}, it has an invalid character at position {}\n", record.id, *pos);
                            continue;
                        }

                        // run the actual search, results are report by a call to the given lambda
                        // uses a scoring matrix
                        search(index, query, *cliErrors, *cliAAA, [&](auto cursor, size_t error) {
                            for (auto sa_entry : cursor) {
                                hitsPerThread += 1;
                                auto [refId, refPos] = index.locate(sa_entry);
                                if (!cliCountOnly) {
//                                        outputBuffers[j] << std::cout << "Found needle #" << h.needle_index << "(" << peptides[h.needle_index] << ") at pos " << h.query_pos << " in protein " << prot << "\n";
                                    outputBuffers[j] << fmt::format("Found needle #{}({}) at pos {} in protein {}\n", i, queries[i].seq, refPos, ref_as_str[refId]);
//                                        outputBuffers[j] << fmt::format("{} {} {} {}\n", i, refId, refPos, error);
                                }
                            }
                        });
                    }
                } while(true);
                totalHits += hitsPerThread;
            });
        }
        workers.clear(); // join all threads

        fmt::print("totalHits: {}\n", totalHits.load());

        if (!cliCountOnly && cliOutputFile) {
            auto ofs = std::ofstream{*cliOutputFile};
            for (auto const& ob : outputBuffers) {
                ofs << ob.str();
            }
        }
    } catch (std::exception const& e) {
        fmt::print(stderr, "error {}\n", e.what());
    }
}
