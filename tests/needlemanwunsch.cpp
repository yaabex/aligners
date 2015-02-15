#include "alignment/basic_score.hpp"
#include "alignment/needleman_wunsch.hpp"

#include <gtest/gtest.h>

using alignment::basic_score;
using alignment::needleman_wunsch;
using std::string;
using std::vector;

template <int Match, int Mismatch, int Insertion, int Deletion>
struct simple_score : basic_score<char> {
  simple_score() : basic_score(Match, Mismatch, Insertion, Deletion) {}
};

using score_matrix_1 = simple_score<1, -1, -1, -1>;
using score_matrix_2 = simple_score<2, -1, -2, -2>;

template <typename ContainerType>
auto to_string(const ContainerType &C) {
  return string(std::begin(C), std::end(C));
}

TEST(NeedlemanWunsch, EmptyLeft) {
  needleman_wunsch<char, score_matrix_1> nw;

  auto res = nw(string(""), string("ACGT"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, score_matrix_1().ins() * 4);
  ASSERT_EQ(Z, "----");
  ASSERT_EQ(W, "ACGT");
}

TEST(NeedlemanWunsch, EmptyRight) {
  needleman_wunsch<char, score_matrix_1> nw;

  auto res = nw(string("ACGT"), string(""), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, score_matrix_1().del() * 4);
  ASSERT_EQ(Z, "ACGT");
  ASSERT_EQ(W, "----");
}

TEST(NeedlemanWunsch, EmptyBoth) {
  needleman_wunsch<char, score_matrix_1> nw;

  auto res = nw(string(""), string(""), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, 0);
  ASSERT_EQ(Z, "");
  ASSERT_EQ(W, "");
}

TEST(NeedlemanWunsch, WikipediaExample1) {
  needleman_wunsch<char, score_matrix_1> nw;

  auto res = nw(string("GCATGCU"), string("GATTACA"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  EXPECT_EQ(score, 0);
  EXPECT_TRUE(Z == "GCATG-CU" || Z == "GCA-TGCU" || Z == "GCAT-GCU");
  EXPECT_EQ(W, "G-ATTACA");
}

TEST(NeedlemanWunsch, WikipediaExample2) {
  needleman_wunsch<char, score_matrix_2> nw;

  auto res = nw(string("AGTACGCA"), string("TATGC"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  EXPECT_EQ(score, 1);
  EXPECT_EQ(Z, "AGTACGCA");
  EXPECT_EQ(W, "--TATGC-");
}
