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

#include "deepvariant/channels/homopolymer_indel_quality_channel.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/channels/channel_utils.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/log/log.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

namespace {

bool ShouldLogHmerRead(const Read& read) {
  const char* target_substr = std::getenv("DV_DEBUG_READ_SUBSTR");
  if (target_substr == nullptr || target_substr[0] == '\0') {
    return false;
  }
  return read.fragment_name().find(target_substr) != std::string::npos;
}

}  // namespace

HomopolymerInDelQualityChannel::HomopolymerInDelQualityChannel(
    int width, const PileupImageOptions& options)
    : Channel(width, options) {}

// Explanation of what TP means in Ultima Genomics data:
// Given read R of size N,
// Cram will have a QUAL vector of size N, and a TP vector of size N.
//
// TP[i] refers to the direction and magnitude of the error which is encoded by
// QUAL[i]
//
// For example:
// If read has sequence ...TAAAAAG...,
// decided_homopolymer_size for the A homopolymer is len(AAAAA) = 5.
// Let i be an index in the read in the range of this homopolymer: [start, end].
// When TP[i] = 1, it means QUAL[i] will contain the PHRED score for the
// homopolymer to be 5 + 1 = 6.
// When TP[i] = -1, it means QUAL[i] will contain the PHRED score for the
// homopolymer to be 5 - 1 = 4.
// When TP[i] = 2, it means QUAL[i] will contain the PHRED score for the
// homopolymer to be 5 + 2 = 7.
std::vector<int8_t> HomopolymerInDelQualityChannel::GetTPValues(
    const Read& read) {
  const bool inspect_tp_read = ShouldLogHmerRead(read);
  std::vector<int8_t> int_tps(read.aligned_sequence().size());
  if (!read.info().contains("tp")) {
    if (inspect_tp_read) {
      LOG(WARNING) << "[TP_DEBUG] read=" << read.fragment_name()
                   << " | has_tp=0"
                   << " | seq_size=" << read.aligned_sequence().size();
    }
    // Return a vector of zeros if tp tag is not present.
    return int_tps;
  }
  const auto& tps = read.info().at("tp").values();

  if (tps.empty()) {
    if (inspect_tp_read) {
      LOG(WARNING) << "[TP_DEBUG] read=" << read.fragment_name()
                   << " | has_tp=1"
                   << " | tp_values_size=0"
                   << " | seq_size=" << read.aligned_sequence().size();
    }
    return int_tps;
  }

  if (inspect_tp_read) {
    LOG(WARNING) << "[TP_DEBUG] read=" << read.fragment_name()
                 << " | has_tp=1"
                 << " | tp_values_size=" << tps.size()
                 << " | seq_size=" << read.aligned_sequence().size();
    for (int i = 0; i < tps.size() && i < 8; ++i) {
      const auto& tpv = tps[i];
      std::string kind = "UNKNOWN";
      if (tpv.kind_case() == nucleus::genomics::v1::Value::kIntValue) {
        kind = "INT";
      } else if (tpv.kind_case() == nucleus::genomics::v1::Value::kNumberValue) {
        kind = "NUMBER";
      } else if (tpv.kind_case() == nucleus::genomics::v1::Value::kStringValue) {
        kind = "STRING";
      }
      LOG(WARNING) << "[TP_DEBUG] read=" << read.fragment_name()
                   << " | idx=" << i
                   << " | kind=" << kind
                   << " | int_value=" << tpv.int_value()
                   << " | number_value=" << tpv.number_value()
                   << " | string_size=" << tpv.string_value().size();
    }
  }

  // Decode TP values into a full vector first. If the read sequence is trimmed
  // by region query, we may need to apply an offset before slicing to the
  // current aligned_sequence() length.
  std::vector<int8_t> decoded_tps;

  // TP may be encoded as a single string where each char stores a signed
  // offset relative to 'A' (matching haplotype's SRead decoding path).
  if (tps.size() == 1 &&
      tps[0].kind_case() == nucleus::genomics::v1::Value::kStringValue &&
      !tps[0].string_value().empty()) {
    const std::string& encoded_tps = tps[0].string_value();
    decoded_tps.resize(encoded_tps.size());
    for (int i = 0; i < encoded_tps.size(); i++) {
      decoded_tps[i] = static_cast<int8_t>(
          static_cast<int>(static_cast<unsigned char>(encoded_tps[i])) -
          static_cast<int>('A'));
    }
  } else {
    decoded_tps.resize(tps.size());
    for (int i = 0; i < tps.size(); i++) {
      const nucleus::genomics::v1::Value& tp_value = tps[i];
      if (tp_value.kind_case() == nucleus::genomics::v1::Value::kIntValue) {
        decoded_tps[i] = static_cast<int8_t>(tp_value.int_value());
      } else if (tp_value.kind_case() ==
                 nucleus::genomics::v1::Value::kNumberValue) {
        decoded_tps[i] = static_cast<int8_t>(tp_value.number_value());
      } else if (tp_value.kind_case() ==
                     nucleus::genomics::v1::Value::kStringValue &&
                 !tp_value.string_value().empty()) {
        decoded_tps[i] = static_cast<int8_t>(
            static_cast<int>(
                static_cast<unsigned char>(tp_value.string_value()[0])) -
            static_cast<int>('A'));
      }
    }
  }

  // If the read is trimmed, aligned_sequence can be a suffix of the original
  // read represented by tp. In that case, use leading hard-clip length as
  // offset into decoded_tps.
  int tp_offset = 0;
  if (read.has_alignment()) {
    for (const auto& cigar_elt : read.alignment().cigar()) {
      if (cigar_elt.operation() == nucleus::genomics::v1::CigarUnit::CLIP_HARD) {
        tp_offset += cigar_elt.operation_length();
      } else {
        break;
      }
    }
  }

  // In some region-query code paths, tp is left untrimmed while
  // aligned_sequence/aligned_quality are trimmed to the overlapping window,
  // and CIGAR may not carry hard-clips for that trimming. In this case,
  // fallback to aligning tp to the rightmost (same-length) window.
  if (tp_offset == 0 && decoded_tps.size() > int_tps.size()) {
    tp_offset = decoded_tps.size() - int_tps.size();
  }

  if (tp_offset < 0 || tp_offset >= decoded_tps.size()) {
    tp_offset = 0;
  }

  for (int i = 0; i < int_tps.size(); i++) {
    int decoded_index = tp_offset + i;
    if (decoded_index >= 0 && decoded_index < decoded_tps.size()) {
      int_tps[i] = decoded_tps[decoded_index];
    }
  }

  if (inspect_tp_read) {
    int non_zero_decoded = 0;
    for (int i = 0; i < decoded_tps.size(); ++i) {
      if (decoded_tps[i] != 0) ++non_zero_decoded;
    }
    int non_zero_final = 0;
    for (int i = 0; i < int_tps.size(); ++i) {
      if (int_tps[i] != 0) ++non_zero_final;
    }
    LOG(WARNING) << "[TP_DEBUG] read=" << read.fragment_name()
                 << " | decoded_size=" << decoded_tps.size()
                 << " | tp_offset=" << tp_offset
                 << " | decoded_non_zero=" << non_zero_decoded
                 << " | final_non_zero=" << non_zero_final;
    for (int i = 0; i < int_tps.size() && i < 20; ++i) {
      LOG(WARNING) << "[TP_DEBUG] read=" << read.fragment_name()
                   << " | final_tp[" << i << "]=" << static_cast<int>(int_tps[i]);
    }
  }

  return int_tps;
}

std::vector<std::uint8_t> HomopolymerInDelQualityChannel::HomoPolymerWeighted(
    const Read& read) {
  // Generates a vector reflecting the number of repeats observed
  std::vector<std::uint8_t> homopolymer_weighted(read.aligned_sequence().size(),
                                                 1);
  const absl::string_view seq = read.aligned_sequence();

  if (seq.empty()) {
    return homopolymer_weighted;
  }

  int i = 0;
  while (i < seq.size()) {
    int hmer_length = 1;
    char current_base = seq[i];
    int j = i + 1;

    // Count consecutive identical bases
    while (j < seq.size() && seq[j] == current_base) {
      hmer_length++;
      j++;
    }

    // Fill all positions in this homopolymer with its length
    for (int k = i; k < j; k++) {
      // Maximum value of uint8_t is 255.
      homopolymer_weighted[k] = std::min(
          hmer_length, static_cast<int>(std::numeric_limits<uint8_t>::max()));
    }

    i = j;
  }

  return homopolymer_weighted;
}

std::vector<std::uint8_t>
HomopolymerInDelQualityChannel::HomoPolymerInDelQuality(const Read& read,
                                                        bool is_deletion) {
  // Decode base quality and tp tags to generate a vector indicating the
  // probability of not having an hmer deletion/insertion for each homopolymer
  // in the read. The probabilities are encoded in phred-scores.
  // The parameter is_deletion determines the direction of the computed error
  // rate per homopolymer.
  const int quality_cap = options_.base_quality_cap();
  std::vector<std::uint8_t> hmer_directed_qualities(
      read.aligned_sequence().size(),
      channels::internal::BaseQualityColor(quality_cap, quality_cap));

  std::string seq(read.aligned_sequence());
  auto hmer_lengths = HomoPolymerWeighted(read);
  auto tps = GetTPValues(read);
  const bool log_this_read = ShouldLogHmerRead(read);

  // If tp tag is not present, return default quality values
  if (tps.empty() || tps.size() != seq.size()) {
    if (log_this_read) {
      LOG(WARNING) << "[HMER_INDEL_READ] read=" << read.fragment_name()
                   << "/" << read.read_number()
                   << " | is_deletion=" << is_deletion
                   << " | status=missing_or_mismatched_tp"
                   << " | seq_size=" << seq.size()
                   << " | tp_size=" << tps.size();
    }
    return hmer_directed_qualities;
  }

  if (log_this_read) {
    LOG(WARNING) << "[HMER_INDEL_READ] read=" << read.fragment_name()
                 << "/" << read.read_number()
                 << " | is_deletion=" << is_deletion
                 << " | status=begin"
                 << " | seq_size=" << seq.size()
                 << " | qual_size=" << read.aligned_quality().size()
                 << " | quality_cap=" << quality_cap;
  }

  int i = 0;
  while (i < hmer_lengths.size()) {
    int hmer_length = hmer_lengths[i];
    float hmer_directed_error_prob = 0;

    // Iterate quality encodings for hmer [i], and sum them up to upward
    // direction or downward direction quality PHRED scores.
    for (int j = 0; j < hmer_length; j++) {
      const int pos = i + j;
      const int tp_value = tps[pos];
      const int encoded_hmer_qual = read.aligned_quality()[pos];
      float error_prob_component = 0;
      bool contributes = false;

      if (tp_value == 0) {
        if (log_this_read) {
          LOG(WARNING) << "[HMER_INDEL_BASE] read=" << read.fragment_name()
                       << "/" << read.read_number()
                       << " | is_deletion=" << is_deletion
                       << " | hmer_start=" << i
                       << " | hmer_length=" << hmer_length
                       << " | pos=" << pos
                       << " | tp=0"
                       << " | encoded_qual=" << encoded_hmer_qual
                       << " | contributes=0"
                       << " | error_prob_component=0";
        }
        continue;
      }

      bool is_deletion_err = tp_value < 0;
      if (is_deletion_err == is_deletion) {
        error_prob_component = std::pow(10, (encoded_hmer_qual / -10.0));
        hmer_directed_error_prob += error_prob_component;
        contributes = true;
      }

      if (log_this_read) {
        LOG(WARNING) << "[HMER_INDEL_BASE] read=" << read.fragment_name()
                     << "/" << read.read_number()
                     << " | is_deletion=" << is_deletion
                     << " | hmer_start=" << i
                     << " | hmer_length=" << hmer_length
                     << " | pos=" << pos
                     << " | tp=" << tp_value
                     << " | encoded_qual=" << encoded_hmer_qual
                     << " | contributes=" << contributes
                     << " | error_prob_component=" << error_prob_component;
      }
    }

    int hmer_directed_quality =
        hmer_directed_error_prob == 0
            ? quality_cap
            : static_cast<int>(-10 * std::log10(hmer_directed_error_prob));
    // Clamp to valid range
    if (hmer_directed_quality > quality_cap) {
      hmer_directed_quality = quality_cap;
    }

    if (log_this_read) {
      LOG(WARNING) << "[HMER_INDEL_HMER] read=" << read.fragment_name()
                   << "/" << read.read_number()
                   << " | is_deletion=" << is_deletion
                   << " | hmer_start=" << i
                   << " | hmer_length=" << hmer_length
                   << " | total_error_prob=" << hmer_directed_error_prob
                   << " | computed_quality=" << hmer_directed_quality
                   << " | color_value="
                   << static_cast<int>(channels::internal::BaseQualityColor(
                          hmer_directed_quality, quality_cap));
    }

    for (int j = 0; j < hmer_length; j++) {
      hmer_directed_qualities[i + j] =
          channels::internal::BaseQualityColor(hmer_directed_quality,
                                               quality_cap);
    }
    i += hmer_length;
  }
  if (log_this_read) {
    LOG(WARNING) << "[HMER_INDEL_READ] read=" << read.fragment_name()
                 << "/" << read.read_number()
                 << " | is_deletion=" << is_deletion
                 << " | status=end"
                 << " | output_size=" << hmer_directed_qualities.size();
  }
  return hmer_directed_qualities;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
