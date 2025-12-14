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

#include "deepvariant/channels/homopolymer_deletion_quality_channel.h"

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
Read CreateReadWithQualities(absl::string_view aligned_sequence,
                             const std::vector<int>& aligned_quality,
                             const std::vector<int>& tp_values = {}) {
  Read read;
  read.set_aligned_sequence(aligned_sequence.data());
  std::string quality_string;
  quality_string.reserve(aligned_quality.size());
  for (int qual : aligned_quality) {
    quality_string.push_back(static_cast<char>(qual));
  }
  read.set_aligned_quality(quality_string);
  if (!tp_values.empty()) {
    auto& info = (*read.mutable_info())["tp"];
    for (int tp : tp_values) {
      info.add_values()->set_int_value(tp);
    }
  }
  return read;
}

class HomopolymerDeletionQualityChannelTest : public ::testing::Test {
 protected:
  PileupImageOptions options_;
  HomopolymerDeletionQualityChannel channel_{10, options_};
};

TEST_F(HomopolymerDeletionQualityChannelTest, HomoPolymerInDelQualityBasic) {
  // Sequence: GGG, qualities: 20, 20, 20. TP: -1, -1, -1 (deletion errors)
  Read read = CreateReadWithQualities("GGG", {20, 20, 20}, {-1, -1, -1});
  // The homopolymer is GGG. It has three deletion errors.
  // Deletion error probability: 3 * 10^(-20/10) = 0.03
  // Q = -10 * log10(0.03) = 15.22. Scaled color: 255 * (15.22/93) = 41
  std::vector<unsigned char> expected_qualities = {41, 41, 41};

  std::vector<std::uint8_t> result =
      channel_.HomoPolymerInDelQuality(read, true);  // true for deletion
  ASSERT_EQ(result.size(), 3);
  EXPECT_EQ(result, expected_qualities);
}

TEST_F(HomopolymerDeletionQualityChannelTest, HomoPolymerInDelQualityMixedTP) {
  // Sequence: GGGT, qualities: 20, 20, 20, 20. TP: -1, 1, -1, 0
  // Only first G and third G should contribute to deletion quality
  Read read = CreateReadWithQualities("GGGT", {20, 20, 20, 20}, {-1, 1, -1, 0});
  // The first homopolymer is GGG. It has two deletion errors at index 0 and 2.
  // The deletion error probability is 2 * 10^(-20/10) = 0.02.
  // Q = -10 * log10(0.02) = 16.98. This is truncated to 16.
  // Scaled color is 255 * (16/93) = 43.
  std::vector<unsigned char> expected_qualities = {
      43, 43, 43, 255};  // Last T has no TP, so max quality

  std::vector<std::uint8_t> result =
      channel_.HomoPolymerInDelQuality(read, true);  // true for deletion
  EXPECT_EQ(result, expected_qualities);
}

TEST_F(HomopolymerDeletionQualityChannelTest, FillReadBase) {
  // Sequence: GGG, qualities: 20, 20, 20. TP: -1, -1, -1
  Read read = CreateReadWithQualities("GGG", {20, 20, 20}, {-1, -1, -1});
  DeepVariantCall dv_call = CreateDvCall();
  std::vector<std::string> alt_alleles = {};
  std::vector<unsigned char> data(10, 0);  // Simulate a data row of size 10

  // Expected color based on 3 Gs with deletion errors. Calculation is similar
  // to the insertion case. Deletion error probability: 3 * 10^(-20/10) = 0.03
  // Q = -10 * log10(0.03) = 15.22. Scaled color: 255 * (15.22/93) = 41
  unsigned char expected_color = 41;

  // Test fill at various read_index and col combinations
  channel_.FillReadBase(data, 0, 'G', 'G', 20, read, 0, dv_call, alt_alleles);
  EXPECT_EQ(data[0], expected_color);
  channel_.FillReadBase(data, 1, 'G', 'G', 20, read, 1, dv_call, alt_alleles);
  EXPECT_EQ(data[1], expected_color);
  channel_.FillReadBase(data, 5, 'G', 'G', 20, read, 2, dv_call, alt_alleles);
  EXPECT_EQ(data[5], expected_color);

  // Test out-of-bounds read_index (should fill with 0)
  channel_.FillReadBase(data, 3, 'G', 'G', 20, read, 3, dv_call, alt_alleles);
  EXPECT_EQ(data[3], 0);
}

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
