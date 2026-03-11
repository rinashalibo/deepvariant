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

#include "deepvariant/channels/homopolymer_insertion_quality_channel.h"

#include <cstdlib>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/channels/homopolymer_indel_quality_channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/log/log.h"

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

HomopolymerInsertionQualityChannel::HomopolymerInsertionQualityChannel(
    int width, const PileupImageOptions& options)
    : HomopolymerInDelQualityChannel(width, options) {}

void HomopolymerInsertionQualityChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  const bool log_this_read = ShouldLogHmerRead(read);
  
  if (!homopolymer_insertion_quality_vector_.has_value()) {
    if (log_this_read) {
      LOG(WARNING) << "[HMER_INS_QUALITY] BEFORE_COMPUTE | read=" << read.fragment_name()
                  << "/" << read.read_number()
                  << " | read.aligned_sequence().size()=" << read.aligned_sequence().size()
                  << " | read.aligned_quality().size()=" << read.aligned_quality().size();
    }
    homopolymer_insertion_quality_vector_ =
        HomoPolymerInDelQuality(read, false);  // false = insertion
     if (log_this_read) {
       LOG(WARNING) << "[HMER_INS_QUALITY] COMPUTED_VECTOR | read=" << read.fragment_name()
                   << "/" << read.read_number()
                   << " | seq_length=" << read.aligned_sequence().size()
                   << " | vector_size=" << homopolymer_insertion_quality_vector_->size();
       if (read.aligned_sequence().size() != homopolymer_insertion_quality_vector_->size()) {
         LOG(WARNING) << "[HMER_INS_QUALITY] SIZE_MISMATCH! read=" << read.fragment_name()
                     << "/" << read.read_number()
                     << " | seq_length=" << read.aligned_sequence().size()
                     << " | vector_size=" << homopolymer_insertion_quality_vector_->size();
       }
       for (size_t i = 0; i < homopolymer_insertion_quality_vector_->size(); i++) {
         LOG(WARNING) << "[HMER_INS_QUALITY_FULL] read=" << read.fragment_name()
                     << "/" << read.read_number()
                     << " | index=" << i
                     << " | value=" << static_cast<int>((*homopolymer_insertion_quality_vector_)[i]);
       }
     }
  }
  
  bool is_valid_index = (read_index >= 0 && 
                        read_index < homopolymer_insertion_quality_vector_->size());
  if (is_valid_index) {
    data[col] = (*homopolymer_insertion_quality_vector_)[read_index];
    if (log_this_read) {
      LOG(WARNING) << "[HMER_INS_QUALITY] read=" << read.fragment_name()
                  << "/" << read.read_number()
                  << " | read_index=" << read_index 
                  << " | col=" << col 
                  << " | value=" << static_cast<int>(data[col])
                  << " | vector_size=" << homopolymer_insertion_quality_vector_->size();
    }
  } else {
    data[col] = 0;
    if (log_this_read) {
      LOG(WARNING) << "[HMER_INS_QUALITY] OUT_OF_BOUNDS! read=" << read.fragment_name()
                  << "/" << read.read_number()
                  << " | read_index=" << read_index 
                  << " | col=" << col 
                  << " | vector_size=" << homopolymer_insertion_quality_vector_->size();
    }
  }
}

void HomopolymerInsertionQualityChannel::FillRefBase(
    std::vector<unsigned char>& ref_data, int col, char ref_base,
    const std::string& ref_bases) {
  ref_data.push_back(0);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
