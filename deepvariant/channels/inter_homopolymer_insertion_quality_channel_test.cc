/*
 * Copyright 2025 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "deepvariant/channels/inter_homopolymer_insertion_quality_channel.h"

#include <cstdint>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "tensorflow/core/platform/test.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace {

using ::nucleus::genomics::v1::Read;

// Creates a DeepVariantCall proto for testing.
DeepVariantCall CreateDvCall() {
  DeepVariantCall dv_call;
  return dv_call;
}

// Creates a Read proto for testing with aligned_sequence and aligned_quality.
Read CreateReadWithT0(absl::string_view aligned_sequence,
                      absl::string_view t0_value,
                      const std::vector<int>& aligned_quality = {}) {
  Read read;
  read.set_aligned_sequence(aligned_sequence.data());
  std::string quality_string;
  quality_string.reserve(aligned_sequence.size());
  if (aligned_quality.empty()) {
    for (int i = 0; i < aligned_sequence.size(); ++i) {
      quality_string.push_back(static_cast<char>(20));  // Default quality
    }
  } else {
    for (int qual : aligned_quality) {
      quality_string.push_back(static_cast<char>(qual));
    }
  }
  read.set_aligned_quality(quality_string);

  auto& info = (*read.mutable_info())["t0"];
  info.add_values()->set_string_value(t0_value.data());
  return read;
}

class InterHomopolymerInsertionQualityChannelTest : public ::testing::Test {
 protected:
  PileupImageOptions options_;
  InterHomopolymerInsertionQualityChannel channel_{10, options_};
};

TEST_F(InterHomopolymerInsertionQualityChannelTest, GetT0QualityValues) {
  // Sequence: AAAA, t0: "5555" (phred score 20 for all bases)
  Read read = CreateReadWithT0("AAAA", "5555");
  // Scaled color for base quality 20 should be 54
  std::vector<std::uint8_t> expected_qualities = {54, 54, 54, 54};

  std::vector<std::uint8_t> result = channel_.GetT0QualityValues(read);
  ASSERT_EQ(result.size(), 4);
  EXPECT_EQ(result, expected_qualities);
}

TEST_F(InterHomopolymerInsertionQualityChannelTest, GetT0QualityValuesEmptyT0) {
  // Sequence: AAAA, empty t0 tag
  Read read = CreateReadWithT0("AAAA", "");
  std::vector<std::uint8_t> expected_qualities = {0, 0, 0,
                                                  0};  // Should return 0s

  std::vector<std::uint8_t> result = channel_.GetT0QualityValues(read);
  ASSERT_EQ(result.size(), 4);
  EXPECT_EQ(result, expected_qualities);
}

TEST_F(InterHomopolymerInsertionQualityChannelTest, FillReadBase) {
  // Sequence: ATGC, t0: "5432" (phred scores 20, 19, 18, 17)
  Read read = CreateReadWithT0("ATGC", "5432");
  DeepVariantCall dv_call = CreateDvCall();
  std::vector<std::string> alt_alleles = {};
  std::vector<unsigned char> data(10, 0);  // Simulate a data row of size 10

  // Expected scaled colors:
  // 20 -> 54  (floor(255 * 20 / 93) = 54)
  // 19 -> 52  (floor(255 * 19 / 93) = 52)
  // 18 -> 49  (floor(255 * 18 / 93) = 49)
  // 17 -> 46  (floor(255 * 17 / 93) = 46)
  std::vector<unsigned char> expected_colors = {54, 52, 49, 46};

  // Test fill at various read_index and col combinations
  channel_.FillReadBase(data, 0, 'A', 'A', 20, read, 0, dv_call, alt_alleles);
  EXPECT_EQ(data[0], expected_colors[0]);
  channel_.FillReadBase(data, 1, 'T', 'T', 20, read, 1, dv_call, alt_alleles);
  EXPECT_EQ(data[1], expected_colors[1]);
  channel_.FillReadBase(data, 5, 'G', 'G', 20, read, 2, dv_call, alt_alleles);
  EXPECT_EQ(data[5], expected_colors[2]);
  channel_.FillReadBase(data, 9, 'C', 'C', 20, read, 3, dv_call, alt_alleles);
  EXPECT_EQ(data[9], expected_colors[3]);

  // Test out-of-bounds read_index (should fill with 0)
  channel_.FillReadBase(data, 3, 'A', 'A', 20, read, 4, dv_call, alt_alleles);
  EXPECT_EQ(data[3], 0);
}

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
