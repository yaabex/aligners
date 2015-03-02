#include "alignment/basic_score.hpp"
#include "alignment/hirschberg.hpp"
#include "alignment/needleman_wunsch.hpp"

#include <gtest/gtest.h>

using std::string;
using std::vector;

template <int Match, int Mismatch, int Insertion, int Deletion>
struct simple_score : alignment::basic_score<char> {
  simple_score() : basic_score(Match, Mismatch, Insertion, Deletion) {}
};

using score_matrix_1 = simple_score<1, -1, -1, -1>;
using score_matrix_2 = simple_score<2, -1, -2, -2>;
alignment::needleman_wunsch<char, score_matrix_1> nw1;
alignment::needleman_wunsch<char, score_matrix_2> nw2;
alignment::hirschberg<char, score_matrix_1> hb1;
alignment::hirschberg<char, score_matrix_2> hb2;


template <typename ContainerType>
auto to_string(const ContainerType &C) {
  return string(std::begin(C), std::end(C));
}

TEST(NeedlemanWunsch, EmptyLeft) {
  auto res = nw1(string(""), string("ACGT"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, score_matrix_1().ins() * 4);
  ASSERT_EQ(Z, "----");
  ASSERT_EQ(W, "ACGT");
}

TEST(NeedlemanWunsch, EmptyRight) {
  auto res = nw1(string("ACGT"), string(""), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, score_matrix_1().del() * 4);
  ASSERT_EQ(Z, "ACGT");
  ASSERT_EQ(W, "----");
}

TEST(NeedlemanWunsch, EmptyBoth) {
  auto res = nw1(string(""), string(""), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, 0);
  ASSERT_EQ(Z, "");
  ASSERT_EQ(W, "");
}

TEST(NeedlemanWunsch, WikipediaExample1) {
  auto res = nw1(string("GCATGCU"), string("GATTACA"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  EXPECT_EQ(score, 0);
  EXPECT_TRUE(Z == "GCATG-CU" || Z == "GCA-TGCU" || Z == "GCAT-GCU");
  EXPECT_EQ(W, "G-ATTACA");
}

TEST(NeedlemanWunsch, WikipediaExample2) {
  auto res = nw2(string("AGTACGCA"), string("TATGC"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  EXPECT_EQ(score, 1);
  EXPECT_EQ(Z, "AGTACGCA");
  EXPECT_EQ(W, "--TATGC-");
}

TEST(Hirschberg, EmptyLeft) {
  auto res = hb1(string(""), string("ACGT"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, score_matrix_1().ins() * 4);
  ASSERT_EQ(Z, "----");
  ASSERT_EQ(W, "ACGT");
}

TEST(Hirschberg, EmptyRight) {
  auto res = hb1(string("ACGT"), string(""), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, score_matrix_1().del() * 4);
  ASSERT_EQ(Z, "ACGT");
  ASSERT_EQ(W, "----");
}

TEST(Hirschberg, EmptyBoth) {
  auto res = hb1(string(""), string(""), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  ASSERT_EQ(score, 0);
  ASSERT_EQ(Z, "");
  ASSERT_EQ(W, "");
}

TEST(Hirschberg, WikipediaExample1) {
  auto res = hb1(string("GCATGCU"), string("GATTACA"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  EXPECT_EQ(score, 0);
  EXPECT_TRUE(Z == "GCATG-CU" || Z == "GCA-TGCU" || Z == "GCAT-GCU");
  EXPECT_EQ(W, "G-ATTACA");
}

TEST(Hirschberg, WikipediaExample2) {
  auto res = hb2(string("AGTACGCA"), string("TATGC"), '-');
  auto score = std::get<0>(res);
  auto Z = to_string(std::get<1>(res));
  auto W = to_string(std::get<2>(res));

  EXPECT_EQ(score, 1);
  EXPECT_EQ(Z, "AGTACGCA");
  EXPECT_EQ(W, "--TATGC-");
}
