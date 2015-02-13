#include "basic_score.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

namespace alignment
{

template <typename AlignType, typename ScoreSystem>
struct needleman_wunsch
{
  typedef typename ScoreSystem::align_t align_t;
  typedef typename ScoreSystem::score_t score_t;

private:
  ScoreSystem score_of;

  struct slot
  {
    enum direction_t
    {
      UP, LEFT, UP_LEFT, UNDEFINED
    };

    slot() : slot(direction_t::UNDEFINED, 0) { }
    slot(direction_t direction, score_t score) : score_(score), direction_(direction) { }

    direction_t &direction() { return direction_; }
    const direction_t &direction() const { return direction_; }

    score_t &score() { return score_; }
    const score_t &score() const { return score_; }

  private:
    score_t score_;
    direction_t direction_;
  };

public:

  template<typename RandomAccessIterator>
  auto operator()(RandomAccessIterator itx, RandomAccessIterator etx,
                  RandomAccessIterator ity, RandomAccessIterator ety,
                  const align_t &GAP_SYMBOL) const
  {
    auto N = std::distance(itx, etx), M = std::distance(ity, ety);
    std::vector<std::vector<slot>> score(N + 1, std::vector<slot>(M + 1));

    // Fill top row with GAP multiples
    for (auto col = 1; col <= M; ++col) {
      score[0][col].score() = score[0][col - 1].score() + score_of.ins();
      score[0][col].direction() = slot::LEFT;
    }

    for (auto row = 1; row <= N; ++row) {
      // Fill left column with GAP multiples
      score[row][0].score() = score[row - 1][0].score() + score_of.del();
      score[row][0].direction() = slot::UP;

      for (auto col = 1; col <= M; ++col) {
        auto up   = score[row - 1][col - 0].score() + score_of.del();
        auto left = score[row - 0][col - 1].score() + score_of.ins();
        auto diag = score[row - 1][col - 1].score() + score_of(itx[row - 1], ity[col - 1]);
        auto max = std::max(up, std::max(left, diag));

        score[row][col].score() = max;

        // TODO: return all best alignments?
        // It may be convenient to know all best alignments. Collecting all directions which have
        // maximum score is one (the?) way to generate them.
        if      (max == diag) { score[row][col].direction() = slot::UP_LEFT; }
        else if (max == up)   { score[row][col].direction() = slot::UP; }
        else if (max == left) { score[row][col].direction() = slot::LEFT; }
        else                  { assert(!"max is neither up, left or diagonal!"); }
      }
    }

    std::vector<std::pair<AlignType, AlignType>> result;
    auto row = N, col = M;
    while (row > 0 || col > 0) {
      switch (score[row][col].direction()) {
      case slot::LEFT:
        --col;
        result.emplace_back(GAP_SYMBOL, ity[col]);
        break;
      case slot::UP:
        --row;
        result.emplace_back(itx[row], GAP_SYMBOL);
        break;
      case slot::UP_LEFT:
        --row, --col;
        result.emplace_back(itx[row], ity[col]);
        break;
      default:
        assert(!"score has undefined direction (should never happen!)");
        break;
      }
    }

    std::reverse(result.begin(), result.end());

    return std::make_pair(score[N][M].score(), result);
  }

  template <typename Container>
  auto operator()(const Container &X, const Container &Y, const AlignType &GAP_SYMBOL) {
    return operator()(std::begin(X), std::end(X), std::begin(Y), std::end(Y), GAP_SYMBOL);
  }
};

}
