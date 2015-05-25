#pragma once

#include <algorithm>

namespace alignment {

template <typename AlignType, typename Matcher = std::equal_to<AlignType>, typename ScoreType = std::ptrdiff_t>
struct basic_score {
  /** Typedef for the alignment type. */
  typedef AlignType align_type;

  /** Typedef for the matcher type.
  *
  * Must implement `bool operator()(const align_type &, const align_type &)`.
  */
  typedef Matcher   matcher_type;

  /** Typedef for the score type.
  *
  * Must implement addition and be comparable.
  */
  typedef ScoreType score_type;

  /** A simple scoring system.
  *
  * \param match The score for a match.
  * \param mismatch The score for a mismatch.
  * \param insertion The score for an insertion (gap in first sequence).
  * \param deletion The score for a deletion (gap in second sequence).
  */
  basic_score(score_type match, score_type mismatch, score_type insertion, score_type deletion)
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
  score_type operator()(const align_type &src, const align_type &dst) const {
    if (matcher_(src, dst)) {
      return match();
    }
    return mismatch();
  }

  /** The score of a match in the sequences. */
  score_type match() const { return match_; }

  /** The score of a mismatch in the sequences. */
  score_type mismatch() const { return mismatch_; }

  /** Gives the score of a gap in the first sequence. */
  score_type ins() const { return insertion_; }

  /** Gives the score of a gap in the second sequence. */
  score_type del() const { return deletion_; }

private:
  score_type match_;
  score_type mismatch_;
  score_type insertion_;
  score_type deletion_;
  matcher_type matcher_;
};

}