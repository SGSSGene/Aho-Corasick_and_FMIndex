// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "fmindex-collection/concepts.h"
#include "fmindex-collection/fmindex/BiFMIndexCursor.h"
#include "fmindex-collection/search_scheme/Scheme.h"
#include "fmindex-collection/search/SelectCursor.h"

// search with hamming distance and using a scoring matrix

namespace fmc::search_hamming_aaa {

template <size_t QuerySigma, size_t RefSigma = QuerySigma>
struct ScoringMatrix {
    std::array<std::array<uint8_t, RefSigma>, QuerySigma> matrix;
    std::array<std::vector<size_t>, QuerySigma> noCostList;
    std::array<std::vector<size_t>, QuerySigma> costList;

    ScoringMatrix() {
        for (size_t y{1}; y < RefSigma; ++y) {
            for (size_t x{1}; x < QuerySigma; ++x) {
                if (x == y) {
                    setCost(x, y, 0);
                } else {
                    setCost(x, y, 1);
                }
            }
        }
    }

    void setCost(size_t queryRank, size_t refRank, size_t cost) {
        std::erase(noCostList[queryRank], refRank);
        std::erase(costList[queryRank], refRank);

        matrix[queryRank][refRank] = cost;
        if (cost == 1) {
            costList[queryRank].push_back(refRank);
        } else {
            if (queryRank != refRank) {
                noCostList[queryRank].push_back(refRank);
            }
        }
    }
};

template <typename index_t, typename search_scheme_t, Sequence query_t, typename delegate_t, typename SM>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = select_cursor_t<index_t>;

    index_t const& index;

    decltype(search_scheme_t::pi) const& pi;
    decltype(search_scheme_t::l) const& l;
    decltype(search_scheme_t::u) const& u;

    query_t const& query;
    SM const& sm;
    delegate_t const& delegate;

    mutable size_t remainingMM{};
    mutable size_t remainingAAA{};

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, SM const& _sm, delegate_t const& _delegate, size_t _remainingMM, size_t _remainingAAA) noexcept
        : index {_index}
        , pi{_search.pi}
        , l{_search.l}
        , u{_search.u}
        , query{_query}
        , sm{_sm}
        , delegate  {_delegate}
        , remainingMM{_remainingMM}
        , remainingAAA{_remainingAAA}
    {
        auto cur       = cursor_t{index};
        searchPart(cur, 0, 0);
    }

    auto extendCursor(cursor_t const& cur, size_t part) const {
        if (part == 0 or pi[part-1] < pi[part]) {
            return cur.extendRight();
        }
        return cur.extendLeft();
    }

    void searchPart(cursor_t const& cur, size_t e, size_t part) const noexcept {
        if (cur.count() == 0) {
            return;
        }
        if (part == query.size()) {
            if (l[part-1] <= e and e <= u[part-1]) {
                delegate(cur, e);
            }
            return;
        }
        if (e > u[part]) {
            return;
        }

        auto rank = query[pi[part]];
        auto cursors = extendCursor(cur, part);

        // search matches
        if (l[part] <= e) {
            searchPart(cursors[rank], e, part+1);
        }

        // search substitutes
        if (e+1 <= u[part]) {
            if (remainingAAA > 0) {
                remainingAAA -= 1;
                for (auto i : sm.noCostList[rank]) {
                    searchPart(cursors[i], e+1, part+1);
                }
                remainingAAA += 1;
            }

            if (remainingMM > 0) {
                remainingMM -= 1;
                if (remainingAAA == 0) {
                    for (auto i : sm.noCostList[rank]) {
                        searchPart(cursors[i], e+1, part+1);
                    }
                }
                for (auto i : sm.costList[rank]) {
                    searchPart(cursors[i], e+1, part+1);
                }
                remainingMM += 1;
            }
        }
    }

};

template <typename index_t, Sequence query_t, typename search_schemes_t, typename delegate_t, typename SM>
void search(index_t const & index, query_t && query, size_t errorsMM, size_t errorsAAA, search_schemes_t const & search_scheme, SM const& sm, delegate_t && delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    for (size_t j{0}; j < search_scheme.size(); ++j) {
        Search{index, search_scheme[j], query, sm, delegate, errorsMM, errorsAAA};
    }
}
}
