#pragma once

#include <tuple>
#include <vector>
#include <iostream>
#include "needleman_wunsch.hpp"

namespace alignment {

template <typename AlignType, typename ScoreSystem>
class hirschberg
{
  typedef typename ScoreSystem::align_type align_type;
  typedef typename ScoreSystem::score_type score_type;

  ScoreSystem score_of;

public:

  template <typename RandomAccessIterator1, typename RandomAccessIterator2>
  auto align_score(RandomAccessIterator1 itx, RandomAccessIterator1 etx,
                   RandomAccessIterator2 ity, RandomAccessIterator2 ety) const
  {
    const auto N = std::distance(itx, etx), M = std::distance(ity, ety);
    std::vector<score_type> row0(M+1);
    std::vector<score_type> row1(M+1);

    for (auto col = 1; col <= M; ++col) {
      row0[col] = row0[col-1] + score_of.ins();
    }

    for (auto row = 1; row <= N; ++row) {
      row1[0] = row0[0] + score_of.del();

      for (auto col = 1; col <= M; ++col) {
        score_type up   = row0[col-0] + score_of.del();
        score_type left = row1[col-1] + score_of.ins();
        score_type diag = row0[col-1] + score_of(itx[row-1], ity[col-1]);

        row1[col] = std::max(up, std::max(left, diag));
      }

      row0.swap(row1);
    }

    return row0;
  }

public:

  template <typename RandomAccessIterator1, typename RandomAccessIterator2>
  auto operator()(RandomAccessIterator1 itx, RandomAccessIterator1 etx,
                  RandomAccessIterator2 ity, RandomAccessIterator2 ety,
                  const align_type &GAP_SYMBOL) const
  {
    const auto N = std::distance(itx, etx);
    const auto M = std::distance(ity, ety);

    if (N == 0) {
      return std::make_tuple(
          M * score_of.ins(),
          std::vector<AlignType>(M, GAP_SYMBOL),
          std::vector<AlignType>(ity, ety));
    }
    else if (M == 0) {
      return std::make_tuple(
          N * score_of.del(),
          std::vector<AlignType>(itx, etx),
          std::vector<AlignType>(N, GAP_SYMBOL));
    }
    else if (N == 1 || M == 1) {
      needleman_wunsch<AlignType, ScoreSystem> nw;
      return nw(itx, etx, ity, ety, GAP_SYMBOL);
    }

    typedef std::reverse_iterator<RandomAccessIterator1> RIter1;
    typedef std::reverse_iterator<RandomAccessIterator2> RIter2;

    auto xmid = N/2;
    auto l_score = align_score(itx, itx + xmid, ity, ety);
    auto r_score = align_score(RIter1(etx), RIter1(itx + xmid), RIter2(ety), RIter2(ity));

    // Y partition
    auto ymid = 0;
    auto ymax = l_score[ymid] + r_score[M - ymid];
    for (auto idx = 1; idx <= M; ++idx) {
      auto v = l_score[idx] + r_score[M - idx];
      if (v > ymax) {
        ymax = v;
        ymid = idx;
      }
    }

    score_type sc0, sc1;
    std::vector<AlignType> Z0, Z1, W0, W1;
    std::tie(sc0, Z0, W0) = operator()(itx, itx + xmid, ity, ity + ymid, GAP_SYMBOL);
    std::tie(sc1, Z1, W1) = operator()(itx + xmid, etx, ity + ymid, ety, GAP_SYMBOL);

    Z0.insert(Z0.end(), Z1.begin(), Z1.end());
    W0.insert(W0.end(), W1.begin(), W1.end());

    return std::make_tuple(sc0 + sc1, Z0, W0);
  }

  template <typename Container1, typename Container2>
  auto operator()(const Container1 &X, const Container2 &Y, const AlignType &GAP_SYMBOL) const {
    return operator()(std::begin(X), std::end(X), std::begin(Y), std::end(Y), GAP_SYMBOL);
  }
};

}
