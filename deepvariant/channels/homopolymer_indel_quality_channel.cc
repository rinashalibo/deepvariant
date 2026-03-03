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
  std::vector<int8_t> int_tps(read.aligned_sequence().size());
  bool log_this_read = (read.fragment_name().find("2749846422") != std::string::npos);
  
  if (log_this_read) {
    LOG(WARNING) << "[GET_TP_VALUES] START | read=" << read.fragment_name()
               << " | has_tp_tag=" << read.info().contains("tp")
               << " | seq_size=" << read.aligned_sequence().size();
  }
  
  if (!read.info().contains("tp")) {
    // Return a vector of zeros if tp tag is not present.
    if (log_this_read) {
      LOG(WARNING) << "[GET_TP_VALUES] NO_TP_TAG | read=" << read.fragment_name();
    }
    return int_tps;
  }
  const auto& tps = read.info().at("tp").values();

  if (log_this_read) {
    LOG(WARNING) << "[GET_TP_VALUES] TP_TAG_FOUND | read=" << read.fragment_name()
               << " | tps.empty()=" << tps.empty()
               << " | tps.size()=" << tps.size()
               << " | expected_size=" << int_tps.size();
  }

  if (tps.empty()) {
    if (log_this_read) {
      LOG(WARNING) << "[GET_TP_VALUES] EMPTY_TP_VALUES | read=" << read.fragment_name();
    }
    return int_tps;
  }

  for (int i = 0; i < tps.size() && i < int_tps.size(); i++) {
    int_tps[i] = tps[i].int_value();
    if (log_this_read && i < 20) {  // Log first 20 values only
      LOG(WARNING) << "[GET_TP_VALUES] VALUE | read=" << read.fragment_name()
                 << " | index=" << i
                 << " | tp_value=" << static_cast<int>(int_tps[i]);
    }
  }
  
  if (log_this_read) {
    LOG(WARNING) << "[GET_TP_VALUES] COMPLETE | read=" << read.fragment_name()
               << " | populated_count=" << (tps.size() < int_tps.size() ? tps.size() : int_tps.size());
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
  std::vector<std::uint8_t> hmer_directed_qualities(
      read.aligned_sequence().size(),
      channels::internal::BaseQualityColor(kMaxQScore));
  
  bool log_this_read = (read.fragment_name().find("2749846422") != std::string::npos);
  if (log_this_read) {
    LOG(WARNING) << "[HMER_INDEL_QUALITY] START | read=" << read.fragment_name() 
             << " | is_deletion=" << is_deletion 
             << " | seq_size=" << read.aligned_sequence().size()
             << " | vector_size=" << hmer_directed_qualities.size()
             << " | kMaxQScore=" << kMaxQScore;
  }

  std::string seq(read.aligned_sequence());
  auto hmer_lengths = HomoPolymerWeighted(read);
  auto tps = GetTPValues(read);

  if (log_this_read) {
    LOG(WARNING) << "[HMER_INDEL_QUALITY] TP_TAG_CHECK | read=" << read.fragment_name()
             << " | tps.empty()=" << tps.empty()
             << " | tps.size()=" << tps.size()
             << " | seq.size()=" << seq.size()
             << " | tps=";
    for (size_t k = 0; k < tps.size() && k < 20; k++) {
      LOG(WARNING) << "tp[" << k << "]="  << static_cast<int>(tps[k]) << " ";
    }
  }

  // If tp tag is not present, return default quality values
  if (tps.empty() || tps.size() != seq.size()) {
    if (log_this_read) {
      LOG(WARNING) << "[HMER_INDEL_QUALITY] RETURNING_DEFAULT | read=" << read.fragment_name()
                 << " | reason=no_tp_tag";
    }
    return hmer_directed_qualities;
  }

  if (log_this_read) {
    LOG(WARNING) << "[HMER_INDEL_QUALITY] READ_QUALITIES | read=" << read.fragment_name()
               << " | aligned_quality_size=" << read.aligned_quality().size()
               << " | aligned_quality_values=";
    for (size_t k = 0; k < read.aligned_quality().size() && k < 20; k++) {
      LOG(WARNING) << "q[" << k << "]="  << static_cast<int>(read.aligned_quality()[k]) << " ";
    }
  }

  int i = 0;
  while (i < hmer_lengths.size()) {
    int hmer_length = hmer_lengths[i];
    float hmer_directed_error_prob = 0;

    if (log_this_read) {
      LOG(WARNING) << "[HMER_INDEL_QUALITY] HMER_GROUP | read=" << read.fragment_name()
                 << " | hmer_start=" << i
                 << " | hmer_length=" << hmer_length
                 << " | is_deletion=" << is_deletion;
    }

    // Iterate quality encodings for hmer [i], and sum them up to upward
    // direction or downward direction quality PHRED scores.
    for (int j = 0; j < hmer_length; j++) {
      int pos = i + j;
      std::uint8_t encoded_hmer_qual = read.aligned_quality()[pos];
      int tp_value = tps[pos];
      bool is_deletion_err = tp_value < 0;
      
      if (log_this_read) {
        LOG(WARNING) << "[HMER_INDEL_QUALITY] BASE | read=" << read.fragment_name()
                   << " | pos=" << pos
                   << " | encoded_qual=" << static_cast<int>(encoded_hmer_qual)
                   << " | tp=" << tp_value
                   << " | is_deletion_err=" << is_deletion_err
                   << " | matches_direction=" << (is_deletion_err == is_deletion);
      }

      if (tp_value == 0) {
        if (log_this_read) {
          LOG(WARNING) << "[HMER_INDEL_QUALITY] SKIPPING_ZERO_TP | pos=" << pos;
        }
        continue;
      }
      
      if (is_deletion_err == is_deletion) {
        float error_prob = std::pow(10, (encoded_hmer_qual / -10.0));
        if (log_this_read) {
          LOG(WARNING) << "[HMER_INDEL_QUALITY] ADDING_ERROR_PROB | pos=" << pos
                     << " | encoded_qual=" << static_cast<int>(encoded_hmer_qual)
                     << " | error_prob_component=" << error_prob;
        }
        hmer_directed_error_prob += error_prob;
      }
    }

    int hmer_directed_quality =
        hmer_directed_error_prob == 0
            ? kMaxQScore
            : static_cast<int>(-10 * std::log10(hmer_directed_error_prob));
    
    if (log_this_read) {
      LOG(WARNING) << "[HMER_INDEL_QUALITY] FINAL_QUALITY | read=" << read.fragment_name()
                 << " | hmer_start=" << i
                 << " | hmer_length=" << hmer_length
                 << " | total_error_prob=" << hmer_directed_error_prob
                 << " | computed_quality=" << hmer_directed_quality
                 << " | color_value=" << static_cast<int>(channels::internal::BaseQualityColor(hmer_directed_quality, 40));
    }

    // Clamp to valid range
    if (hmer_directed_quality > kMaxQScore) {
      if (log_this_read) {
        LOG(WARNING) << "[HMER_INDEL_QUALITY] CLAMPING | original=" << hmer_directed_quality
                   << " | to_max=" << kMaxQScore;
      }
      hmer_directed_quality = kMaxQScore;
    }

    for (int j = 0; j < hmer_length; j++) {
      int pos = i + j;
      // Use max_quality of 40 for homopolymer quality color scaling (Ultima standard)
      std::uint8_t color_value = channels::internal::BaseQualityColor(hmer_directed_quality, 40);
      hmer_directed_qualities[pos] = color_value;
      if (log_this_read) {
        LOG(WARNING) << "[HMER_INDEL_QUALITY_FULL] read=" << read.fragment_name()
                   << " | pos=" << pos
                   << " | quality=" << hmer_directed_quality
                   << " | color_value=" << static_cast<int>(color_value);
      }
    }
    i += hmer_length;
  }
  
  if (log_this_read) {
    LOG(WARNING) << "[HMER_INDEL_QUALITY] COMPLETE | read=" << read.fragment_name()
               << " | is_deletion=" << is_deletion
               << " | final_vector_size=" << hmer_directed_qualities.size();
  }
  
  return hmer_directed_qualities;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
