// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef COMMON_DATAMODEL_MFTMCHMATCHTRUTH_H_
#define COMMON_DATAMODEL_MFTMCHMATCHTRUTH_H_
#include <Framework/AnalysisDataModel.h>
#include <cstdint>
#include <vector>
namespace o2::aod {
namespace mftmchmatchtruth {
DECLARE_SOA_COLUMN(ResidualRowId, residualRowId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(MchMcParticleId, mchMcParticleId, int64_t);
DECLARE_SOA_COLUMN(MftMcParticleId, mftMcParticleId, int64_t);
DECLARE_SOA_COLUMN(HasMchMcParticle, hasMchMcParticle, bool);
DECLARE_SOA_COLUMN(HasMftMcParticle, hasMftMcParticle, bool);
DECLARE_SOA_COLUMN(IsCorrectMatch, isCorrectMatch, bool);
// -1 = missing/invalid particle label, 0 = distinct particles, 1 = same particle.
DECLARE_SOA_COLUMN(IsFalseMatch, isFalseMatch, bool);
DECLARE_SOA_COLUMN(MatchTruthStatus, matchTruthStatus, int8_t);
DECLARE_SOA_COLUMN(MchMcMask, mchMcMask, uint8_t);
DECLARE_SOA_COLUMN(MftMcMask, mftMcMask, uint8_t);
DECLARE_SOA_COLUMN(MchMcCollisionId, mchMcCollisionId, int64_t);
DECLARE_SOA_COLUMN(MftMcCollisionId, mftMcCollisionId, int64_t);
DECLARE_SOA_COLUMN(MchPdgCode, mchPdgCode, int);
DECLARE_SOA_COLUMN(MftPdgCode, mftPdgCode, int);
DECLARE_SOA_COLUMN(MchMcPx, mchMcPx, float);
DECLARE_SOA_COLUMN(MchMcPy, mchMcPy, float);
DECLARE_SOA_COLUMN(MchMcPz, mchMcPz, float);
DECLARE_SOA_COLUMN(MftMcPx, mftMcPx, float);
DECLARE_SOA_COLUMN(MftMcPy, mftMcPy, float);
DECLARE_SOA_COLUMN(MftMcPz, mftMcPz, float);
DECLARE_SOA_COLUMN(MchPhysicalPrimary, mchPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(MftPhysicalPrimary, mftPhysicalPrimary, bool);
}
DECLARE_SOA_TABLE(MftMchMatchTruths, "AOD", "FWDMPMCLBL",
 mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
 mftmchmatchtruth::MchMcParticleId, mftmchmatchtruth::MftMcParticleId,
 mftmchmatchtruth::HasMchMcParticle, mftmchmatchtruth::HasMftMcParticle, mftmchmatchtruth::IsCorrectMatch,
 mftmchmatchtruth::MchMcMask, mftmchmatchtruth::MftMcMask,
 mftmchmatchtruth::MchMcCollisionId, mftmchmatchtruth::MftMcCollisionId,
 mftmchmatchtruth::MchPdgCode, mftmchmatchtruth::MftPdgCode,
 mftmchmatchtruth::MchMcPx, mftmchmatchtruth::MchMcPy, mftmchmatchtruth::MchMcPz,
 mftmchmatchtruth::MftMcPx, mftmchmatchtruth::MftMcPy, mftmchmatchtruth::MftMcPz,
 mftmchmatchtruth::MchPhysicalPrimary, mftmchmatchtruth::MftPhysicalPrimary,
 mftmchmatchtruth::IsFalseMatch, mftmchmatchtruth::MatchTruthStatus);
DECLARE_SOA_TABLE(MixedMftMchMatchTruths, "AOD", "FWDMPMIXLBL",
 mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
 mftmchmatchtruth::MchMcParticleId, mftmchmatchtruth::MftMcParticleId,
 mftmchmatchtruth::HasMchMcParticle, mftmchmatchtruth::HasMftMcParticle, mftmchmatchtruth::IsCorrectMatch,
 mftmchmatchtruth::MchMcMask, mftmchmatchtruth::MftMcMask,
 mftmchmatchtruth::MchMcCollisionId, mftmchmatchtruth::MftMcCollisionId,
 mftmchmatchtruth::MchPdgCode, mftmchmatchtruth::MftPdgCode,
 mftmchmatchtruth::MchMcPx, mftmchmatchtruth::MchMcPy, mftmchmatchtruth::MchMcPz,
 mftmchmatchtruth::MftMcPx, mftmchmatchtruth::MftMcPy, mftmchmatchtruth::MftMcPz,
 mftmchmatchtruth::MchPhysicalPrimary, mftmchmatchtruth::MftPhysicalPrimary,
 mftmchmatchtruth::IsFalseMatch, mftmchmatchtruth::MatchTruthStatus);
namespace mftmchmatchrelation {
DECLARE_SOA_COLUMN(PairOrigin, pairOrigin, int8_t);
DECLARE_SOA_COLUMN(AncestryStatus, ancestryStatus, int);
DECLARE_SOA_COLUMN(AncestryComplete, ancestryComplete, bool);
DECLARE_SOA_COLUMN(IsMchMuon, isMchMuon, bool);
DECLARE_SOA_COLUMN(IsMftAncestorOfMch, isMftAncestorOfMch, bool);
DECLARE_SOA_COLUMN(IsMftToMuonDecay, isMftToMuonDecay, bool);
DECLARE_SOA_COLUMN(IsDirectMftToMuonDecay, isDirectMftToMuonDecay, bool);
DECLARE_SOA_COLUMN(AncestryDepth, ancestryDepth, int);
DECLARE_SOA_COLUMN(MchMcProcess, mchMcProcess, int);
DECLARE_SOA_COLUMN(MftMcProcess, mftMcProcess, int);
DECLARE_SOA_COLUMN(MchProducedByGenerator, mchProducedByGenerator, bool);
DECLARE_SOA_COLUMN(MftProducedByGenerator, mftProducedByGenerator, bool);
DECLARE_SOA_COLUMN(MchMcVertexX, mchMcVertexX, float);
DECLARE_SOA_COLUMN(MchMcVertexY, mchMcVertexY, float);
DECLARE_SOA_COLUMN(MchMcVertexZ, mchMcVertexZ, float);
DECLARE_SOA_COLUMN(MftMcVertexX, mftMcVertexX, float);
DECLARE_SOA_COLUMN(MftMcVertexY, mftMcVertexY, float);
DECLARE_SOA_COLUMN(MftMcVertexZ, mftMcVertexZ, float);
DECLARE_SOA_COLUMN(MchMotherIds, mchMotherIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MchMotherPdgs, mchMotherPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(MftMotherIds, mftMotherIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MftMotherPdgs, mftMotherPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(AncestryIds, ancestryIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(AncestryPdgs, ancestryPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(AncestryProcesses, ancestryProcesses, std::vector<int>);
DECLARE_SOA_COLUMN(AncestryGeneratorFlags, ancestryGeneratorFlags, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(AncestryVertexX, ancestryVertexX, std::vector<float>);
DECLARE_SOA_COLUMN(AncestryVertexY, ancestryVertexY, std::vector<float>);
DECLARE_SOA_COLUMN(AncestryVertexZ, ancestryVertexZ, std::vector<float>);
}
// Row-aligned with MatchPlaneResiduals, DownstreamMatchPlaneResiduals and MftMchMatchTruths.
DECLARE_SOA_TABLE(MftMchMatchRelations, "AOD", "FWDMPMCREL",
 mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
 mftmchmatchrelation::PairOrigin,
 mftmchmatchrelation::AncestryStatus,
 mftmchmatchrelation::AncestryComplete,
 mftmchmatchrelation::IsMchMuon,
 mftmchmatchrelation::IsMftAncestorOfMch,
 mftmchmatchrelation::IsMftToMuonDecay,
 mftmchmatchrelation::IsDirectMftToMuonDecay,
 mftmchmatchrelation::AncestryDepth,
 mftmchmatchrelation::MchMcProcess,
 mftmchmatchrelation::MftMcProcess,
 mftmchmatchrelation::MchProducedByGenerator,
 mftmchmatchrelation::MftProducedByGenerator,
 mftmchmatchrelation::MchMcVertexX,
 mftmchmatchrelation::MchMcVertexY,
 mftmchmatchrelation::MchMcVertexZ,
 mftmchmatchrelation::MftMcVertexX,
 mftmchmatchrelation::MftMcVertexY,
 mftmchmatchrelation::MftMcVertexZ,
 mftmchmatchrelation::MchMotherIds,
 mftmchmatchrelation::MchMotherPdgs,
 mftmchmatchrelation::MftMotherIds,
 mftmchmatchrelation::MftMotherPdgs,
 mftmchmatchrelation::AncestryIds,
 mftmchmatchrelation::AncestryPdgs,
 mftmchmatchrelation::AncestryProcesses,
 mftmchmatchrelation::AncestryGeneratorFlags,
 mftmchmatchrelation::AncestryVertexX,
 mftmchmatchrelation::AncestryVertexY,
 mftmchmatchrelation::AncestryVertexZ);
} // namespace o2::aod
#endif
