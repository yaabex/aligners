#pragma once

#include <algorithm>

namespace alignment {

template <typename AlignType, typename Matcher = std::equal_to<AlignType>, typename ScoreType = std::ptrdiff_t>
struct basic_score {
  /** Typedef for the alignment type. */
  typedef AlignType align_t;

  /** Typedef for the matcher type.
  *
  * Must implement `bool operator()(const align_t &, const align_t &)`.
  */
  typedef Matcher   matcher_t;

  /** Typedef for the score type.
  *
  * Must implement addition and be comparable.
  */
  typedef ScoreType score_t;

  /** A simple scoring system.
  *
  * \param match The score for a match.
  * \param mismatch The score for a mismatch.
  * \param insertion The score for an insertion (gap in first sequence).
  * \param deletion The score for a deletion (gap in second sequence).
  */
  basic_score(score_t match, score_t mismatch, score_t insertion, score_t deletion)
      : match_(match)
      , mismatch_(mismatch)
      , insertion_(insertion)
      , deletion_(deletion) {}

  /** Scores the given element combination.
  *
  * If the element matches, returns the value of `match()`. Otherwise, returns `mismatch()`.
  *
  * \param src The element of the first sequence.
  * \param dst The element of the second sequence.
  */
  score_t operator()(const align_t &src, const align_t &dst) const {
    if (matcher_(src, dst)) {
      return match();
    }
    return mismatch();
  }

  /** The score of a match in the sequences. */
  score_t match() const { return match_; }

  /** The score of a mismatch in the sequences. */
  score_t mismatch() const { return mismatch_; }

  /** Gives the score of a gap in the first sequence. */
  score_t ins() const { return insertion_; }

  /** Gives the score of a gap in the second sequence. */
  score_t del() const { return deletion_; }

private:
  score_t match_;
  score_t mismatch_;
  score_t insertion_;
  score_t deletion_;
  matcher_t matcher_;
};

}