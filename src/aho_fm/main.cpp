#include <channel/value_mutex.h>
#include <clice/clice.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/fmindex/diskStorage.h>
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
                                      .value = std::string{},
                                      .tags = {"required"},
};
auto cliErrors = clice::Argument { .args = {"-k", "--errors"},
                                   .desc = "number of errors that are allowed during the search",
                                   .value = size_t{0},
};

auto cliThreads = clice::Argument { .args = {"-t", "--threads"},
                                    .desc = "number of threads",
                                    .value = size_t{1},
};

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
        size_t chunk_size = 100;

        auto lastProcessedQuery = channel::value_mutex<size_t>{0};

        auto workers = std::vector<std::jthread>{};
        for (size_t j{0}; j < *cliThreads; ++j) {
            workers.emplace_back([&, j]() {
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
                    g.unlock();

                    // process search
                    for (size_t i{startId}; i < endId; ++i) {
                        auto const& record = queries[i];

                        fmt::print("starting {}/{} ({}%)\n", i, queries.size(), i*100./queries.size());

                        auto query = ivs::convert_char_to_rank<Alphabet>(record.seq);
                        if (auto pos = ivs::verify_rank(query); pos) {
                            fmt::print("skipping record {}, it has an invalid character at position {}\n", record.id, *pos);
                            continue;
                        }

                        // run the actual search, results are report by a call to the given lambda
                        fmc::search</*.editdistance=*/false>(index, query, *cliErrors, [&](auto cursor, size_t error) {
                            for (auto i : cursor) {
                                auto [refId, refPos] = index.locate(i);
                                outputBuffers[j] << fmt::format("{} {} {} {}\n", i, refId, refPos, error);
                            }
                        });
                    }
                } while(true);
            });
        }
        workers.clear(); // join all threads

        auto ofs = std::ofstream{*cliOutputFile};
        for (auto const& ob : outputBuffers) {
            ofs << ob.str();
        }
    } catch (std::exception const& e) {
        fmt::print(stderr, "error {}\n", e.what());
    }
}
