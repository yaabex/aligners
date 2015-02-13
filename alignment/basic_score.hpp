#include <algorithm>

namespace alignment {

using std::ptrdiff_t;

template <typename AlignType, typename Matcher = std::equal_to<AlignType>, typename ScoreType = ptrdiff_t>
struct basic_score {
  typedef AlignType align_t;
  typedef Matcher   matcher_t;
  typedef ScoreType score_t;

  basic_score(score_t match, score_t mismatch, score_t insertion, score_t deletion)
      : match_(match)
      , mismatch_(mismatch)
      , insertion_(insertion)
      , deletion_(deletion) {}

  score_t operator()(const align_t &src, const align_t &dst) const {
    if (matcher_(src, dst)) {
      return match_;
    }
    return mismatch_;
  }

  score_t ins() const { return insertion_; }
  score_t del() const { return deletion_; }

private:
  score_t match_;
  score_t mismatch_;
  score_t insertion_;
  score_t deletion_;
  matcher_t matcher_;
};

}