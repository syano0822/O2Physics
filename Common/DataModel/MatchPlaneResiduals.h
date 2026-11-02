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

// Propagation, conversions, timing and chi2 definitions based on
// PWGDQ/Tasks/global-muon-matcher.cxx (0dfb35775). No realignment is applied here.
// Phi residuals in angular chi2 are wrapped to [-pi, pi) for periodicity.
/// \file MatchPlaneResiduals.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_DATAMODEL_MATCHPLANERESIDUALS_H_
#define COMMON_DATAMODEL_MATCHPLANERESIDUALS_H_
#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>
namespace o2::aod
{
namespace matchplaneresidual
{
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(MchCollisionId, mchCollisionId, int64_t);
DECLARE_SOA_COLUMN(MftCollisionId, mftCollisionId, int64_t);
DECLARE_SOA_COLUMN(MatchingPlaneZ, matchingPlaneZ, float);
DECLARE_SOA_COLUMN(BzAtMftCenter, bzAtMftCenter, float);
DECLARE_SOA_COLUMN(MchX, mchX, double);
DECLARE_SOA_COLUMN(MchY, mchY, double);
DECLARE_SOA_COLUMN(MchPhi, mchPhi, double);
DECLARE_SOA_COLUMN(MchTanl, mchTanl, double);
DECLARE_SOA_COLUMN(MchInvQPt, mchInvQPt, double);
DECLARE_SOA_COLUMN(MftX, mftX, double);
DECLARE_SOA_COLUMN(MftY, mftY, double);
DECLARE_SOA_COLUMN(MftPhi, mftPhi, double);
DECLARE_SOA_COLUMN(MftTanl, mftTanl, double);
DECLARE_SOA_COLUMN(MftInvQPt, mftInvQPt, double);
DECLARE_SOA_COLUMN(DeltaX, deltaX, double);
DECLARE_SOA_COLUMN(DeltaY, deltaY, double);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, double);
DECLARE_SOA_COLUMN(DeltaPhiWrapped, deltaPhiWrapped, double);
DECLARE_SOA_COLUMN(DeltaTanl, deltaTanl, double);
DECLARE_SOA_COLUMN(DeltaInvQPt, deltaInvQPt, double);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, double);
DECLARE_SOA_COLUMN(DeltaPt, deltaPt, double);
DECLARE_SOA_COLUMN(Chi2All, chi2All, double);
DECLARE_SOA_COLUMN(Chi2XYPhiTanl, chi2XYPhiTanl, double);
DECLARE_SOA_COLUMN(Chi2XY, chi2XY, double);
DECLARE_SOA_COLUMN(IsFinite, isFinite, bool);
DECLARE_SOA_COLUMN(MchPairMassPV, mchPairMassPV, double);
DECLARE_SOA_COLUMN(HasMchPairMassPV, hasMchPairMassPV, bool);
DECLARE_SOA_COLUMN(IsMchPairJpsi, isMchPairJpsi, bool);
DECLARE_SOA_COLUMN(MchPairPartnerId, mchPairPartnerId, int64_t);
DECLARE_SOA_COLUMN(NumMchJpsiPartners, numMchJpsiPartners, int);
DECLARE_SOA_COLUMN(AbsorberStepCm, absorberStepCm, double);
// 0=success; 1..5=joint transport/covariance status; 101..106=standalone MFT transport failure.
DECLARE_SOA_COLUMN(TransportStatus, transportStatus, int);
} // namespace matchplaneresidual
DECLARE_SOA_TABLE(DownstreamMatchPlaneResiduals, "AOD", "FWDMPDOWN",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::AbsorberStepCm,
                  matchplaneresidual::TransportStatus);

namespace matchplanemixing
{
DECLARE_SOA_COLUMN(MixingStreamId, mixingStreamId, uint64_t);
DECLARE_SOA_COLUMN(TargetFrameId, targetFrameId, uint64_t);
DECLARE_SOA_COLUMN(DonorFrameId, donorFrameId, uint64_t);
} // namespace matchplanemixing
// Optional diagnostic matrices; row-major 5x5. No change to matching.
namespace matchplanemixdiag
{
DECLARE_SOA_COLUMN(IsMixed, isMixed, bool);
DECLARE_SOA_COLUMN(MchCov, mchCov, std::vector<double>);
DECLARE_SOA_COLUMN(MftCov, mftCov, std::vector<double>);
DECLARE_SOA_COLUMN(ResidualCov, residualCov, std::vector<double>);
} // namespace matchplanemixdiag
DECLARE_SOA_TABLE(MftMchMixingDiagnostics, "AOD", "FWDMPMIXDG",
                  matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  matchplaneresidual::MatchingPlaneZ, matchplanemixdiag::IsMixed,
                  matchplanemixdiag::MchCov, matchplanemixdiag::MftCov, matchplanemixdiag::ResidualCov,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
// Mixed rows use independent tables; both planes carry source IDs and pool diagnostics.
namespace matchplanemixing
{
DECLARE_SOA_COLUMN(DonorMchTrackId, donorMchTrackId, int64_t);
DECLARE_SOA_COLUMN(DonorCollisionId, donorCollisionId, int64_t);
DECLARE_SOA_COLUMN(MixedGroupIndex, mixedGroupIndex, int);
DECLARE_SOA_COLUMN(NumMixedGroups, numMixedGroups, int);
DECLARE_SOA_COLUMN(MixDeltaBc, mixDeltaBc, int64_t);
DECLARE_SOA_COLUMN(MixDeltaPhi, mixDeltaPhi, double);
DECLARE_SOA_COLUMN(MixDeltaR, mixDeltaR, double);
DECLARE_SOA_COLUMN(MixDeltaAttemptsRel, mixDeltaAttemptsRel, double);
DECLARE_SOA_COLUMN(MixDeltaZ, mixDeltaZ, double);
DECLARE_SOA_COLUMN(TargetMatchAttempts, targetMatchAttempts, int);
DECLARE_SOA_COLUMN(DonorMatchAttempts, donorMatchAttempts, int);
DECLARE_SOA_COLUMN(TargetMftMultiplicity, targetMftMultiplicity, int);
DECLARE_SOA_COLUMN(DonorMftMultiplicity, donorMftMultiplicity, int);
DECLARE_SOA_COLUMN(TargetMchMultiplicity, targetMchMultiplicity, int);
DECLARE_SOA_COLUMN(DonorMchMultiplicity, donorMchMultiplicity, int);
DECLARE_SOA_COLUMN(TargetCentrality, targetCentrality, double);
DECLARE_SOA_COLUMN(DonorCentrality, donorCentrality, double);
DECLARE_SOA_COLUMN(TargetCbMultiplicity, targetCbMultiplicity, int);
DECLARE_SOA_COLUMN(DonorCbMultiplicity, donorCbMultiplicity, int);
DECLARE_SOA_COLUMN(MixingWeight, mixingWeight, double);
// Combine frame ID with the output job/stream identity when merging jobs.

DECLARE_SOA_COLUMN(RetainedFrames, retainedFrames, std::vector<uint64_t>);
DECLARE_SOA_COLUMN(RetainedEvents, retainedEvents, uint64_t);
DECLARE_SOA_COLUMN(RetainedMftStates, retainedMftStates, uint64_t);
DECLARE_SOA_COLUMN(EvictedEvents, evictedEvents, uint64_t);
DECLARE_SOA_COLUMN(OversizedEvents, oversizedEvents, uint64_t);
DECLARE_SOA_COLUMN(PoolWasReset, poolWasReset, bool);
} // namespace matchplanemixing
// One state row per input DF, including empty/rejected DFs. Also synchronizes
// MC label-cache lifetime with the exact donor-history frames owned upstream.
DECLARE_SOA_TABLE(MftMchMixingPoolStates, "AOD", "FWDMPMIXPS",
                  matchplanemixing::TargetFrameId, matchplanemixing::RetainedFrames,
                  matchplanemixing::RetainedEvents, matchplanemixing::RetainedMftStates,
                  matchplanemixing::EvictedEvents, matchplanemixing::OversizedEvents,
                  matchplanemixing::PoolWasReset, matchplanemixing::MixingStreamId);
DECLARE_SOA_TABLE(RawMixedMatchPlaneResiduals, "AOD", "FWDMPMIXRAW",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::AbsorberStepCm,
                  matchplaneresidual::TransportStatus,
                  matchplanemixing::DonorMchTrackId,
                  matchplanemixing::DonorCollisionId,
                  matchplanemixing::MixedGroupIndex,
                  matchplanemixing::NumMixedGroups,
                  matchplanemixing::MixDeltaBc,
                  matchplanemixing::MixDeltaPhi,
                  matchplanemixing::MixDeltaR,
                  matchplanemixing::MixDeltaAttemptsRel,
                  matchplanemixing::MixDeltaZ,
                  matchplanemixing::TargetMatchAttempts,
                  matchplanemixing::DonorMatchAttempts,
                  matchplanemixing::TargetMftMultiplicity,
                  matchplanemixing::DonorMftMultiplicity,
                  matchplanemixing::TargetMchMultiplicity,
                  matchplanemixing::DonorMchMultiplicity,
                  matchplanemixing::TargetCentrality,
                  matchplanemixing::DonorCentrality,
                  matchplanemixing::TargetCbMultiplicity,
                  matchplanemixing::DonorCbMultiplicity,
                  matchplanemixing::MixingWeight,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
DECLARE_SOA_TABLE(MixedMatchPlaneResiduals, "AOD", "FWDMPMIX",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::AbsorberStepCm,
                  matchplaneresidual::TransportStatus,
                  matchplanemixing::DonorMchTrackId,
                  matchplanemixing::DonorCollisionId,
                  matchplanemixing::MixedGroupIndex,
                  matchplanemixing::NumMixedGroups,
                  matchplanemixing::MixDeltaBc,
                  matchplanemixing::MixDeltaPhi,
                  matchplanemixing::MixDeltaR,
                  matchplanemixing::MixDeltaAttemptsRel,
                  matchplanemixing::MixDeltaZ,
                  matchplanemixing::TargetMatchAttempts,
                  matchplanemixing::DonorMatchAttempts,
                  matchplanemixing::TargetMftMultiplicity,
                  matchplanemixing::DonorMftMultiplicity,
                  matchplanemixing::TargetMchMultiplicity,
                  matchplanemixing::DonorMchMultiplicity,
                  matchplanemixing::TargetCentrality,
                  matchplanemixing::DonorCentrality,
                  matchplanemixing::TargetCbMultiplicity,
                  matchplanemixing::DonorCbMultiplicity,
                  matchplanemixing::MixingWeight,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
DECLARE_SOA_TABLE(RawMatchPlaneResiduals, "AOD", "FWDMPRERAW",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite);
DECLARE_SOA_TABLE(MatchPlaneResiduals, "AOD", "FWDMPRESID",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::MchPairMassPV,
                  matchplaneresidual::HasMchPairMassPV,
                  matchplaneresidual::IsMchPairJpsi,
                  matchplaneresidual::MchPairPartnerId,
                  matchplaneresidual::NumMchJpsiPartners);
DECLARE_SOA_TABLE(JpsiMatchPlaneResiduals, "AOD", "FWDMPJPSI",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::MchPairMassPV,
                  matchplaneresidual::HasMchPairMassPV,
                  matchplaneresidual::IsMchPairJpsi,
                  matchplaneresidual::MchPairPartnerId,
                  matchplaneresidual::NumMchJpsiPartners);
// One row per input MCH(-MID) track, including rejected tracks for selection QA.
// IDs must be combined with the input file and DF identifier.
namespace matchplaneeventsel
{
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MchCollisionId, mchCollisionId, int64_t);
DECLARE_SOA_COLUMN(SelectionMode, selectionMode, int8_t); // 0=all, 1=int8/sel8, 2=UPC candidate
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);
DECLARE_SOA_COLUMN(HasMchBC, hasMchBC, bool);
DECLARE_SOA_COLUMN(MchGlobalBC, mchGlobalBC, uint64_t);
DECLARE_SOA_COLUMN(NumMchMidTracks, numMchMidTracks, int);
DECLARE_SOA_COLUMN(MchSel8, mchSel8, bool);
DECLARE_SOA_COLUMN(HasFV0, hasFV0, bool);
DECLARE_SOA_COLUMN(MaxFV0Amplitude, maxFV0Amplitude, float);
DECLARE_SOA_COLUMN(NumMchTracks, numMchTracks, int); // MCH-MID + MCH-only at original collisionId
DECLARE_SOA_COLUMN(NumMftTracks, numMftTracks, int); // All MFT at same collisionId, before covariance selection
DECLARE_SOA_COLUMN(PassTrackMultiplicity, passTrackMultiplicity, bool);
DECLARE_SOA_COLUMN(PassMchUnlikeSign, passMchUnlikeSign, bool);
DECLARE_SOA_COLUMN(PassFIT, passFIT, bool); // UDHelpers cleanFIT in configured BC window; true outside UPC mode
DECLARE_SOA_COLUMN(PassMuonQuality, passMuonQuality, bool);
DECLARE_SOA_COLUMN(IsMuonPropagated, isMuonPropagated, bool);
DECLARE_SOA_COLUMN(MuonCutEta, muonCutEta, float);
DECLARE_SOA_COLUMN(MuonCutPt, muonCutPt, float);
DECLARE_SOA_COLUMN(MftAbsorberStepCm, mftAbsorberStepCm, double); // effective material delta-z; 0=whole absorber/inactive
} // namespace matchplaneeventsel
DECLARE_SOA_TABLE(MatchPlaneEventSels, "AOD", "FWDMPEVSEL",
                  matchplaneeventsel::MchTrackId,
                  matchplaneeventsel::MchCollisionId,
                  matchplaneeventsel::SelectionMode,
                  matchplaneeventsel::IsSelected,
                  matchplaneeventsel::HasMchBC,
                  matchplaneeventsel::MchGlobalBC,
                  matchplaneeventsel::NumMchMidTracks,
                  matchplaneeventsel::MchSel8,
                  matchplaneeventsel::HasFV0,
                  matchplaneeventsel::MaxFV0Amplitude,
                  matchplaneeventsel::NumMchTracks,
                  matchplaneeventsel::NumMftTracks,
                  matchplaneeventsel::PassTrackMultiplicity,
                  matchplaneeventsel::PassMchUnlikeSign,
                  matchplaneeventsel::PassFIT,
                  matchplaneeventsel::PassMuonQuality,
                  matchplaneeventsel::IsMuonPropagated,
                  matchplaneeventsel::MuonCutEta,
                  matchplaneeventsel::MuonCutPt,
                  matchplaneeventsel::MftAbsorberStepCm);
// One row per selected same-collision unordered MCH pair, independent of MFT candidates.
namespace mchdimuonpv
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId1, mchTrackId1, int64_t);
DECLARE_SOA_COLUMN(MchTrackId2, mchTrackId2, int64_t);
DECLARE_SOA_COLUMN(PairCharge, pairCharge, int8_t);
DECLARE_SOA_COLUMN(VertexX, vertexX, float);
DECLARE_SOA_COLUMN(VertexY, vertexY, float);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
DECLARE_SOA_COLUMN(Px1, px1, double);
DECLARE_SOA_COLUMN(Py1, py1, double);
DECLARE_SOA_COLUMN(Pz1, pz1, double);
DECLARE_SOA_COLUMN(Px2, px2, double);
DECLARE_SOA_COLUMN(Py2, py2, double);
DECLARE_SOA_COLUMN(Pz2, pz2, double);
DECLARE_SOA_COLUMN(Mass, mass, double);
DECLARE_SOA_COLUMN(PairPt, pairPt, double);
DECLARE_SOA_COLUMN(PairRapidity, pairRapidity, double);
} // namespace mchdimuonpv
DECLARE_SOA_TABLE(MchDimuonsPV, "AOD", "FWDMCHPAIR",
                  mchdimuonpv::CollisionId, mchdimuonpv::MchTrackId1, mchdimuonpv::MchTrackId2,
                  mchdimuonpv::PairCharge, mchdimuonpv::VertexX, mchdimuonpv::VertexY, mchdimuonpv::VertexZ,
                  mchdimuonpv::Px1, mchdimuonpv::Py1, mchdimuonpv::Pz1,
                  mchdimuonpv::Px2, mchdimuonpv::Py2, mchdimuonpv::Pz2,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity);
// ML selections are reconstructed from real-data observables only.
namespace mlmatchedpv
{
DECLARE_SOA_COLUMN(FrameId, frameId, uint64_t);
DECLARE_SOA_COLUMN(MlGoodScore, mlGoodScore, double);
DECLARE_SOA_COLUMN(HasMftCovariance, hasMftCovariance, bool);
DECLARE_SOA_COLUMN(Status, status, int); // 0=success, 1=MCH PV failure, 2=invalid momentum, 3=nonfinite matched PV state
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(Px, px, double);
DECLARE_SOA_COLUMN(Py, py, double);
DECLARE_SOA_COLUMN(Pz, pz, double);
DECLARE_SOA_COLUMN(MftTrackId1, mftTrackId1, int64_t);
DECLARE_SOA_COLUMN(MftTrackId2, mftTrackId2, int64_t);
DECLARE_SOA_COLUMN(MlGoodScore1, mlGoodScore1, double);
DECLARE_SOA_COLUMN(MlGoodScore2, mlGoodScore2, double);
} // namespace mlmatchedpv
DECLARE_SOA_TABLE(MlMatchedMuonsPV, "AOD", "FWDMLMUONPV",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  mlmatchedpv::MlGoodScore, mlmatchedpv::HasMftCovariance, mlmatchedpv::Status,
                  mlmatchedpv::Charge, mlmatchedpv::Px, mlmatchedpv::Py, mlmatchedpv::Pz);
DECLARE_SOA_TABLE(MlMatchedDimuonsPV, "AOD", "FWDMLPAIRPV",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId,
                  mchdimuonpv::MchTrackId1, mlmatchedpv::MftTrackId1, mchdimuonpv::MchTrackId2, mlmatchedpv::MftTrackId2,
                  mlmatchedpv::MlGoodScore1, mlmatchedpv::MlGoodScore2, mchdimuonpv::PairCharge,
                  mchdimuonpv::Px1, mchdimuonpv::Py1, mchdimuonpv::Pz1, mchdimuonpv::Px2, mchdimuonpv::Py2, mchdimuonpv::Pz2,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity);

// Raw-pair alignment QA; MC can reject a raw mixed pair later. Join the final
// mixed table by (MchTrackId,MftTrackId,MatchingPlaneZ,DonorMchTrackId).
namespace mixedalignment
{
DECLARE_SOA_COLUMN(Applied, applied, bool);
DECLARE_SOA_COLUMN(Status, status, int);
DECLARE_SOA_COLUMN(DonorTransportStatus, donorTransportStatus, int);
DECLARE_SOA_COLUMN(DonorParameters, donorParameters, std::vector<double>);
DECLARE_SOA_COLUMN(TargetParameters, targetParameters, std::vector<double>);
DECLARE_SOA_COLUMN(MftBeforeParameters, mftBeforeParameters, std::vector<double>);
DECLARE_SOA_COLUMN(AnchorDelta, anchorDelta, std::vector<double>);
DECLARE_SOA_COLUMN(ResidualBefore, residualBefore, std::vector<double>);
DECLARE_SOA_COLUMN(Rotation, rotation, std::vector<double>);
DECLARE_SOA_COLUMN(DonorSeedMomentum, donorSeedMomentum, double);
} // namespace mixedalignment
DECLARE_SOA_TABLE(MftMchMixedAlignments, "AOD", "FWDMPMIXAL",
                  matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  matchplanemixing::DonorMchTrackId, matchplaneresidual::MatchingPlaneZ,
                  mixedalignment::Applied, mixedalignment::Status, mixedalignment::DonorTransportStatus,
                  mixedalignment::DonorParameters, mixedalignment::TargetParameters,
                  mixedalignment::MftBeforeParameters, mixedalignment::AnchorDelta, mixedalignment::ResidualBefore,
                  mixedalignment::Rotation, mixedalignment::DonorSeedMomentum,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
// Pair-specific covariance diagnostics. Matrices are row-major. Joint order:
// [MFT(x,y,phi,tanl,q/pt), MCH(x,y,phi,tanl,q/pt)] at matching plane.
namespace mftmchpaircov
{
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(PairCovValid, pairCovValid, bool);
DECLARE_SOA_COLUMN(PairCovStatus, pairCovStatus, int);
DECLARE_SOA_COLUMN(SeedMomentum, seedMomentum, double);
DECLARE_SOA_COLUMN(OriginalMftInvQPt, originalMftInvQPt, double);
DECLARE_SOA_COLUMN(UsesMchCharge, usesMchCharge, bool);
DECLARE_SOA_COLUMN(RankAll, rankAll, int);
DECLARE_SOA_COLUMN(RankFour, rankFour, int);
DECLARE_SOA_COLUMN(RankXY, rankXY, int);
DECLARE_SOA_COLUMN(JointCov, jointCov, std::vector<double>);
DECLARE_SOA_COLUMN(ResidualCov, residualCov, std::vector<double>);
DECLARE_SOA_COLUMN(SeedCov, seedCov, std::vector<double>);
DECLARE_SOA_COLUMN(SeedMchCrossCov, seedMchCrossCov, std::vector<double>);
DECLARE_SOA_COLUMN(PairJacobian, pairJacobian, std::vector<double>);
DECLARE_SOA_COLUMN(UpstreamSigmaELoss, upstreamSigmaELoss, double);
DECLARE_SOA_COLUMN(JacobianScale, jacobianScale, double);
} // namespace mftmchpaircov
DECLARE_SOA_TABLE(MftMchPairCovariances, "AOD", "FWDMPPAIRCV",
                  mftmchpaircov::MchTrackId,
                  mftmchpaircov::MftTrackId,
                  mftmchpaircov::PairCovValid,
                  mftmchpaircov::PairCovStatus,
                  mftmchpaircov::SeedMomentum,
                  mftmchpaircov::OriginalMftInvQPt,
                  mftmchpaircov::UsesMchCharge,
                  mftmchpaircov::RankAll,
                  mftmchpaircov::RankFour,
                  mftmchpaircov::RankXY,
                  mftmchpaircov::JointCov,
                  mftmchpaircov::ResidualCov,
                  mftmchpaircov::SeedCov,
                  mftmchpaircov::SeedMchCrossCov,
                  mftmchpaircov::PairJacobian,
                  mftmchpaircov::UpstreamSigmaELoss,
                  mftmchpaircov::JacobianScale);
namespace selectedmatchingmuon
{
DECLARE_SOA_COLUMN(Method, method, int8_t);     // 0=MCHPV,1=chi2,2=ML
DECLARE_SOA_COLUMN(Score, score, double);       // chi2 or good-source score; MCHPV NaN
DECLARE_SOA_COLUMN(RawScore, rawScore, double); // chi2 or good logit; MCHPV NaN
DECLARE_SOA_COLUMN(DcaX, dcaX, double);
DECLARE_SOA_COLUMN(DcaY, dcaY, double);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, double);
} // namespace selectedmatchingmuon
DECLARE_SOA_TABLE(SelectedMatchingMuons, "AOD", "FWDEVALMUON",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  selectedmatchingmuon::Method, selectedmatchingmuon::Score, selectedmatchingmuon::RawScore,
                  mlmatchedpv::Charge, mlmatchedpv::Px, mlmatchedpv::Py, mlmatchedpv::Pz,
                  selectedmatchingmuon::DcaX, selectedmatchingmuon::DcaY, selectedmatchingmuon::DcaZ);

} // namespace o2::aod
#endif // COMMON_DATAMODEL_MATCHPLANERESIDUALS_H_
