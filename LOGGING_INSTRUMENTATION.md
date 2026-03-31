# DeepVariant Pileup Processing Instrumentation

## Overview
Added detailed logging throughout the pileup image processing pipeline to track data flow and identify where divergence occurs between haplotype and deepvariant repos.

## Instrumentation Points

### 1. **PILEUP_START** - BuildPileupForOneSample Entry
**File**: `deepvariant/pileup_image_native.cc` (lines ~335)

Logs at the beginning of pileup assembly after read downsampling:
- Variant location (reference_name:start)
- Total reads available vs. sampled
- Details of first 10 sampled reads:
  - Read name and read_number
  - Mapping quality
  - Alignment position

**Purpose**: Verify which reads are selected for the pileup

### 2. **ENCODE_READ_START** - Per-Read Processing Entry
**File**: `deepvariant/pileup_image_native.cc` (lines ~375)

Logs when EncodeRead is called for each read:
- Read fragment_name/read_number
- Sequence and quality length
- Alignment position

**Purpose**: Track which reads enter the encoding pipeline

### 3. **ENCODE_READ Result** - Post-EncodeRead Status
**File**: `deepvariant/pileup_image_native.cc` (lines ~382-383)

Logs:
- If EncodeRead returned nullptr (likely due to low quality at call site)
- If EncodeRead succeeded, number of channels created

**Purpose**: Identify which reads fail to encode

### 4. **ENCODE_READ_DETAIL** - Detailed Encoding Diagnostics
**File**: `deepvariant/pileup_image_native.cc` (lines ~510)

Detailed logs within EncodeRead function:
- Mapping quality check (actual vs. minimum required)
- image_start_pos
- ref_bases size and content
- Number of channels
- Result of CalculateChannels call

**Purpose**: Understand reference bases and channel setup for each read

### 5. **Channel Value Changes** - Per-CIGAR-Unit Processing
**File**: `deepvariant/pileup_channel_lib.cc` (lines ~150-193)

Logs during channel calculation:
- CIGAR operation type (MATCH, MISMATCH, INSERT, DELETE)
- Reference position and read position
- Column index in image
- Read base, ref base, base quality
- For each channel:
  - Channel index and enum value
  - Value before and after FillReadBase

**Purpose**: Trace exactly which bases are written to which channels

### 6. **PILEUP_END** - Final Assembly Summary
**File**: `deepvariant/pileup_image_native.cc` (lines ~465-478)

Logs at BuildPileupForOneSample return:
- Total rows (references + reads + empty)
- Breakdown of rows:
  - Reference band rows
  - Read rows
  - Empty rows
- Image width
- Details of first 5 reads in final pileup:
  - Haplotype index
  - Allele support group
  - Alignment position

**Purpose**: Verify final pileup composition and read ordering

## Log Format Convention

All logs follow this pattern:
```
=== STAGE_NAME | key_info ===
  sub_info1: value1
  sub_info2: value2
```

## Comparison Strategy

To find divergence between repos:
1. Run make_examples with both repos on same CRAM
2. Enable logging: `--alsologtostderr --v=2`
3. Capture logs to files
4. Compare PILEUP_START sections:
   - Are same reads sampled?
   - Are they in same order?
5. If yes, compare ENCODE_READ sections:
   - Do all reads successfully encode in both?
   - Are ref_bases identical?
6. If yes, compare ENCODE_READ_DETAIL:
   - Are mapping qualities checked?
   - Are image_start_pos values same?
7. If yes, compare Channel calculations:
   - Are channel values being written identically?
8. If yes, compare PILEUP_END:
   - Are final row compositions identical?

## Running with Logging

```bash
docker run \
  -v /data:/data \
  deepvariant-local \
  /opt/deepvariant/bin/make_examples \
  --ref <reference.fasta> \
  --reads <input.cram> \
  --examples <output.tfrecord> \
  --regions "chr1:1000-2000" \
  --mode calling \
  --alsologtostderr \
  --v=2 \
  2>&1 | tee dv_detailed.log
```

Then filter logs by stage:
```bash
grep "=== PILEUP" dv_detailed.log
grep "=== ENCODE_READ" dv_detailed.log
grep "=== Channel" dv_detailed.log | head -50
```
