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

#ifndef COMMON_CORE_MFTMCHMATCHING_H_
#define COMMON_CORE_MFTMCHMATCHING_H_

// Shared data/MC matching implementation. Keep MC annotations independent of
// reconstructed-track selection, propagation and inference.

// -----------------------------------------------------------------------------
// MftMchMatchingWorkflow
// -----------------------------------------------------------------------------
// Single implementation shared by the data and MC entry points.
// Keep truth labels outside this workflow: they must not affect matching.
#include "PWGUD/Core/UDHelpers.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/VarManager.h"

// -----------------------------------------------------------------------------
// GlobalMuonMatchingPlane
// -----------------------------------------------------------------------------
// Propagation, conversions, timing and chi2 definitions based on
// PWGDQ/Tasks/global-muon-matcher.cxx (0dfb35775). No realignment is applied here.
// Phi residuals in angular chi2 are wrapped to [-pi, pi) for periodicity.
#include "Common/Core/fwdtrackUtilities.h"
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MFTTracking/Constants.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackParam.h>
#include <ReconstructionDataFormats/TrackFwd.h>
#include <Math/MatrixFunctions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TVector2.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>
namespace o2::aod::fwdtrackutils
{
struct GlobalMuonMatchingPlane {
  using SMatrix55Sym = o2::track::SMatrix55Sym;
  using SMatrix55Std = o2::track::SMatrix55Std;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;
  static constexpr int MchMidTrackType = 3;
  static constexpr int ThetaAbsBoundaryDeg = 3;
  static constexpr double SlopeResolutionZ = 535.;
  double mBzAtMftCenter{0.};
  std::map<std::string, std::function<std::tuple<double, int>(const o2::track::TrackParCovFwd&, const o2::track::TrackParCovFwd&)>> mMatchingFunctionMap;
  GlobalMuonMatchingPlane() { initMatchingFunctions(); }
  template <class T>
  o2::mch::TrackParam fwdToMch(const T& fwdtrack)
  {
    // Convert Forward Track parameters and covariances matrix to the
    // MCH track format.

    // Parameter conversion
    const double x2 = fwdtrack.getPhi();
    const double x3 = fwdtrack.getTanl();
    const double x4 = fwdtrack.getInvQPt();

    const auto sinX2 = std::sin(x2);
    const auto cosX2 = std::cos(x2);

    const double alpha1 = cosX2 / x3;
    const double alpha3 = sinX2 / x3;
    const double alpha4 = x4 / std::sqrt(x3 * x3 + sinX2 * sinX2);

    const auto kNorm = std::sqrt(x3 * x3 + sinX2 * sinX2);
    const auto kNorm3 = kNorm * kNorm * kNorm;

    // Covariances matrix conversion
    SMatrix55Std jacobian;
    SMatrix55Sym covariances;

    covariances(0, 0) = fwdtrack.getCovariances()(0, 0);
    covariances(0, 1) = fwdtrack.getCovariances()(0, 1);
    covariances(0, 2) = fwdtrack.getCovariances()(0, 2);
    covariances(0, 3) = fwdtrack.getCovariances()(0, 3);
    covariances(0, 4) = fwdtrack.getCovariances()(0, 4);

    covariances(1, 1) = fwdtrack.getCovariances()(1, 1);
    covariances(1, 2) = fwdtrack.getCovariances()(1, 2);
    covariances(1, 3) = fwdtrack.getCovariances()(1, 3);
    covariances(1, 4) = fwdtrack.getCovariances()(1, 4);

    covariances(2, 2) = fwdtrack.getCovariances()(2, 2);
    covariances(2, 3) = fwdtrack.getCovariances()(2, 3);
    covariances(2, 4) = fwdtrack.getCovariances()(2, 4);

    covariances(3, 3) = fwdtrack.getCovariances()(3, 3);
    covariances(3, 4) = fwdtrack.getCovariances()(3, 4);

    covariances(4, 4) = fwdtrack.getCovariances()(4, 4);

    jacobian(0, 0) = 1;

    jacobian(1, 2) = -sinX2 / x3;
    jacobian(1, 3) = -cosX2 / (x3 * x3);

    jacobian(2, 1) = 1;

    jacobian(3, 2) = cosX2 / x3;
    jacobian(3, 3) = -sinX2 / (x3 * x3);

    jacobian(4, 2) = -x4 * sinX2 * cosX2 / kNorm3;
    jacobian(4, 3) = -x3 * x4 / kNorm3;
    jacobian(4, 4) = 1 / kNorm;
    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    std::array<double, 15> cov = {covariances(0, 0), covariances(1, 0), covariances(1, 1), covariances(2, 0), covariances(2, 1), covariances(2, 2), covariances(3, 0), covariances(3, 1), covariances(3, 2), covariances(3, 3), covariances(4, 0), covariances(4, 1), covariances(4, 2), covariances(4, 3), covariances(4, 4)};
    std::array<double, 5> param = {fwdtrack.getX(), alpha1, fwdtrack.getY(), alpha3, alpha4};

    o2::mch::TrackParam convertedTrack(fwdtrack.getZ(), param.data(), cov.data());
    return {convertedTrack};
  }

  o2::track::TrackParCovFwd mchToFwd(const o2::mch::TrackParam& mchParam)
  {
    // Convert a MCH Track parameters and covariances matrix to the
    // Forward track format. Must be called after propagation though the absorber

    o2::track::TrackParCovFwd convertedTrack;

    // Parameter conversion
    const double alpha1 = mchParam.getNonBendingSlope();
    const double alpha3 = mchParam.getBendingSlope();
    const double alpha4 = mchParam.getInverseBendingMomentum();

    const double x2 = std::atan2(-alpha3, -alpha1);
    const double x3 = -1. / std::sqrt(alpha3 * alpha3 + alpha1 * alpha1);
    const double x4 = alpha4 * -x3 * std::sqrt(1 + alpha3 * alpha3);

    const auto kNorm = alpha1 * alpha1 + alpha3 * alpha3;
    const auto kNorm32 = kNorm * std::sqrt(kNorm);
    const auto slopeLen = std::sqrt(alpha3 * alpha3 + 1);

    // Covariances matrix conversion
    SMatrix55Std jacobian;
    SMatrix55Sym covariances;

    covariances(0, 0) = mchParam.getCovariances()(0, 0);
    covariances(0, 1) = mchParam.getCovariances()(0, 1);
    covariances(0, 2) = mchParam.getCovariances()(0, 2);
    covariances(0, 3) = mchParam.getCovariances()(0, 3);
    covariances(0, 4) = mchParam.getCovariances()(0, 4);

    covariances(1, 1) = mchParam.getCovariances()(1, 1);
    covariances(1, 2) = mchParam.getCovariances()(1, 2);
    covariances(1, 3) = mchParam.getCovariances()(1, 3);
    covariances(1, 4) = mchParam.getCovariances()(1, 4);

    covariances(2, 2) = mchParam.getCovariances()(2, 2);
    covariances(2, 3) = mchParam.getCovariances()(2, 3);
    covariances(2, 4) = mchParam.getCovariances()(2, 4);

    covariances(3, 3) = mchParam.getCovariances()(3, 3);
    covariances(3, 4) = mchParam.getCovariances()(3, 4);

    covariances(4, 4) = mchParam.getCovariances()(4, 4);

    jacobian(0, 0) = 1;

    jacobian(1, 2) = 1;

    jacobian(2, 1) = -alpha3 / kNorm;
    jacobian(2, 3) = alpha1 / kNorm;

    jacobian(3, 1) = alpha1 / kNorm32;
    jacobian(3, 3) = alpha3 / kNorm32;

    jacobian(4, 1) = -alpha1 * alpha4 * slopeLen / kNorm32;
    jacobian(4, 3) = alpha3 * alpha4 * (1 / (std::sqrt(kNorm) * slopeLen) - slopeLen / kNorm32);
    jacobian(4, 4) = slopeLen / std::sqrt(kNorm);

    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    // Set output
    convertedTrack.setX(mchParam.getNonBendingCoor());
    convertedTrack.setY(mchParam.getBendingCoor());
    convertedTrack.setZ(mchParam.getZ());
    convertedTrack.setPhi(x2);
    convertedTrack.setTanl(x3);
    convertedTrack.setInvQPt(x4);
    convertedTrack.setCharge(mchParam.getCharge());
    convertedTrack.setCovariances(covariances);

    return convertedTrack;
  }

  void initMatchingFunctions()
  {
    using SVector2 = ROOT::Math::SVector<double, 2>;
    using SVector4 = ROOT::Math::SVector<double, 4>;
    using SVector5 = ROOT::Math::SVector<double, 5>;

    using SMatrix44 = ROOT::Math::SMatrix<double, 4>;
    using SMatrix45 = ROOT::Math::SMatrix<double, 4, 5>;
    using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
    using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;

    // Define built-in matching functions
    //________________________________________________________________________________
    mMatchingFunctionMap["matchALL"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating all parameters: X,Y, phi, tanl & q/pt

      SMatrix55Sym hK, vK;
      SVector5 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl(), mftTrack.getInvQPt()),
        rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      vK(2, 2) = mftTrack.getCovariances()(2, 2);
      vK(3, 3) = mftTrack.getCovariances()(3, 3);
      vK(4, 4) = mftTrack.getCovariances()(4, 4);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;
      hK(2, 2) = 1.0;
      hK(3, 3) = 1.0;
      hK(4, 4) = 1.0;

      // Covariance of residuals
      SMatrix55Std invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Update Parameters
      rKKminus1 = mK - hK * globalMuonTrackParameters; // Residuals of prediction
      rKKminus1[2] = TVector2::Phi_mpi_pi(rKKminus1[2]);

      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 5};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating positions & angles

      SMatrix45 hK;
      SMatrix44 vK;
      SVector4 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl()),
        rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      vK(2, 2) = mftTrack.getCovariances()(2, 2);
      vK(3, 3) = mftTrack.getCovariances()(3, 3);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;
      hK(2, 2) = 1.0;
      hK(3, 3) = 1.0;

      // Covariance of residuals
      SMatrix44 invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      rKKminus1 = mK - hK * globalMuonTrackParameters;
      rKKminus1[2] = TVector2::Phi_mpi_pi(rKKminus1[2]);

      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 4};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Calculate Matching Chi2 - X and Y positions

      SMatrix25 hK;
      SMatrix22 vK;
      SVector2 mK(mftTrack.getX(), mftTrack.getY()), rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;

      // Covariance of residuals
      SMatrix22 invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      rKKminus1 = mK - hK * globalMuonTrackParameters;
      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return reduced chi2
      return {matchChi2Track, 2};
    };
  }

  o2::track::TrackParCovFwd propagateToZMch(const o2::track::TrackParCovFwd& muon, const double z)
  {
    auto mchTrack = fwdToMch(muon);

    float absFront = -90.f;
    float absBack = -505.f;

    if (muon.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else {
      // all other cases
      o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
    }

    return mchToFwd(mchTrack);
  }

  o2::track::TrackParCovFwd propagateToZMft(const o2::track::TrackParCovFwd& mftTrack, const double z)
  {
    o2::track::TrackParCovFwd trackExtrap{mftTrack};
    trackExtrap.propagateToZ(z, mBzAtMftCenter);
    return trackExtrap;
  }

  template <class EVT, class BC, class TMUON, class TMFT>
  bool isMftMchTimeCompatible(EVT const& collisions,
                              BC const& bcs,
                              TMUON const& mchTrack,
                              TMFT const& mftTrack)
  {
    if (!mchTrack.has_collision() || !mftTrack.has_collision()) {
      return false;
    }

    const auto& collMch = collisions.rawIteratorAt(mchTrack.collisionId());
    const auto& bcMch = bcs.rawIteratorAt(collMch.bcId());
    const auto& collMft = collisions.rawIteratorAt(mftTrack.collisionId());
    const auto& bcMft = bcs.rawIteratorAt(collMft.bcId());

    int64_t deltaBc = static_cast<int64_t>(bcMft.globalBC()) - static_cast<int64_t>(bcMch.globalBC());
    double deltaBcNS = o2::constants::lhc::LHCBunchSpacingNS * deltaBc;
    double deltaTrackTime = mftTrack.trackTime() - mchTrack.trackTime() + deltaBcNS;
    double trackTimeResTot = mftTrack.trackTimeRes() + mchTrack.trackTimeRes();

    return std::fabs(deltaTrackTime) <= trackTimeResTot;
  }

  // Equivalent to VarManager::PropagateMuon(..., kToVertex) with its default
  // zero coordinate shifts. Kept here to avoid a Common -> PWGDQ dependency.
  template <class T, class C>
  o2::dataformats::GlobalFwdTrack propagateToVertexMchNative(const T& track, const C& collision)
  {
    o2::globaltracking::MatchGlobalFwd matching;
    auto mch = matching.FwdtoMCH(getTrackParCovFwd3DShift(track, 0.f, 0.f, 0.f, track));
    o2::mch::TrackExtrap::extrapToVertex(mch, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    return matching.MCHtoFwd(mch);
  }

  template <class T, class C>
  bool pDcaCut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    constexpr double AbsorberEndZ = 505.;
    constexpr double RadToDeg = 180. / o2::constants::math::PI;
    double thetaAbs = std::atan(mchTrack.rAtAbsorberEnd() / AbsorberEndZ) * RadToDeg;

    // propagate muon track to vertex
    auto mchTrackAtVertex = propagateToVertexMchNative(mchTrack, collision);

    // double pUncorr = mchTrack.p();
    double p = mchTrackAtVertex.getP();

    double pDCA = mchTrack.pDca();
    double sigmaPDCA = (thetaAbs < ThetaAbsBoundaryDeg) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = SlopeResolutionZ * slopeRes * p;
    double sigmaPDCAWithRes = std::sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
    return pDCA <= nSigmaPDCA * sigmaPDCAWithRes;
  }
  template <class T, class C>
  bool isGoodMuon(const T& mchTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    // chi2 cut
    if (mchTrack.chi2() > chi2Cut) {
      return false;
    }

    // momentum cut
    if (mchTrack.p() < pCut) {
      return false; // skip low-momentum tracks
    }

    // transverse momentum cut
    if (mchTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mchTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    // RAbs cut
    double rAbs = mchTrack.rAtAbsorberEnd();
    if ((rAbs < rAbsCut[0] || rAbs > rAbsCut[1])) {
      return false;
    }

    // pDCA cut
    return pDcaCut(mchTrack, collision, nSigmaPdcaCut);
  }
  template <class TMCH, class C>
  o2::track::TrackParCovFwd propagateToVertexMch(const TMCH& muon,
                                                 const C& collision)
  {
    auto mchTrack = fwdToMch(fwdtrackutils::getTrackParCovFwd(muon, muon));
    o2::mch::TrackExtrap::extrapToVertex(mchTrack,
                                         collision.posX(),
                                         collision.posY(),
                                         collision.posZ(),
                                         collision.covXX(),
                                         collision.covYY());
    return mchToFwd(mchTrack);
  }
  template <class TMUON, class C, class TagConfig>
  void getTaggedMuons(C const& collisions,
                      TMUON const& muonTracks,
                      std::vector<int64_t>& taggedMuons, TagConfig& configMuonTagging)
  {
    taggedMuons.clear();
    for (const auto& muonTrack : muonTracks) {

      // only consider MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) != MchMidTrackType) {
        continue;
      }

      // only select MCH-MID tracks associated to a collision
      if (!muonTrack.has_collision()) {
        continue;
      }

      const auto& collision = collisions.rawIteratorAt(muonTrack.collisionId());

      // select MCH tracks with strict quality cuts
      if (!isGoodMuon(muonTrack, collision,
                      configMuonTagging.cfgMuonTaggingTrackChi2MchUp,
                      configMuonTagging.cfgMuonTaggingPMchLow,
                      configMuonTagging.cfgMuonTaggingPtMchLow,
                      {configMuonTagging.cfgMuonTaggingEtaMchLow, configMuonTagging.cfgMuonTaggingEtaMchUp},
                      {configMuonTagging.cfgMuonTaggingRabsLow, configMuonTagging.cfgMuonTaggingRabsUp},
                      configMuonTagging.cfgMuonTaggingPdcaUp)) {
        continue;
      }

      // propagate MCH track to the vertex
      auto mchTrackAtVertex = propagateToVertexMch(muonTrack, collision);

      // propagate the track from the vertex to the first MFT plane
      const auto& extrapToMFTfirst = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[0]);
      double rFront = std::sqrt(extrapToMFTfirst.getX() * extrapToMFTfirst.getX() + extrapToMFTfirst.getY() * extrapToMFTfirst.getY());
      if (rFront < configMuonTagging.cfgMuonTaggingRadiusAtMftFrontLow.value || rFront > configMuonTagging.cfgMuonTaggingRadiusAtMftFrontUp.value) {
        continue;
      }

      // propagate the track from the vertex to the last MFT plane
      const auto& extrapToMFTlast = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[9]);
      double rBack = std::sqrt(extrapToMFTlast.getX() * extrapToMFTlast.getX() + extrapToMFTlast.getY() * extrapToMFTlast.getY());
      if (rBack < configMuonTagging.cfgMuonTaggingRadiusAtMftBackLow.value || rBack > configMuonTagging.cfgMuonTaggingRadiusAtMftBackUp.value) {
        continue;
      }

      int64_t muonTrackIndex = muonTrack.globalIndex();
      taggedMuons.emplace_back(muonTrackIndex);
    }
  }
};
} // namespace o2::aod::fwdtrackutils

// -----------------------------------------------------------------------------
// MftMchEventMixing
// -----------------------------------------------------------------------------
#include <cmath>
#include <cstdint>
#include <limits>
namespace o2::analysis::mftmch::mixing
{
struct Config {
  int depth = -1;
  int64_t minDeltaBc = 3564;
  double maxDeltaPhi = 3.14159265358979323846 / 10.; // negative disables position-phi selection.
  double maxDeltaR = 10.;
  double maxDeltaAttempts = .1;
  double maxDeltaZ = 1.;
  bool symmetricCuts = true;              // false preserves signed cuts in global-muon-matcher.
  double maxDeltaMftMultiplicityRel = .1; // -1 disables; symmetric relative to target.
  int maxDeltaMchMultiplicity = 0;        // -1 disables; 0 requires equal counts.
  double maxDeltaCentrality = -1.;        // percentile points; -1 disables.
  double maxDeltaMchPhiPV = -1.;          // wrapped momentum-azimuth difference (rad); -1 disables.
  double maxDeltaMchMomentumRel = -1.;    // |donor pPV - target pPV| / target pPV; negative disables.
  double maxDeltaMchEtaPV = -1.;          // momentum-pseudorapidity difference; -1 disables.
  double maxDeltaCbMultiplicityRel = -1.; // PVMults::multNTracksPV, relative to target.
  bool mftMultiplicityStrict = false;     // true excludes the configured boundary.
  bool mchRadiusStrict = false;
  bool mchEtaStrict = false;
  double maxDeltaMchPtRel = -1.; // |donor pT(PV) - target pT(PV)| / target pT(PV); negative disables.
  bool cbMultiplicityStrict = false;
  bool vertexZStrict = false;
};
struct Event {
  int64_t collision = -1, mch = -1, bc = 0;
  int run = -1, attempts = 0;
  double phi = 0., radius = 0., z = 0.;
  int mftMultiplicity = -1, mchMultiplicity = -1;
  int cbMultiplicity = -1;
  double centrality = std::numeric_limits<double>::quiet_NaN();
  double mchPhiPV = std::numeric_limits<double>::quiet_NaN();
  double mchEtaPV = std::numeric_limits<double>::quiet_NaN();
  double mchMomentumPV = std::numeric_limits<double>::quiet_NaN();
  double mchPtPV = std::numeric_limits<double>::quiet_NaN();
  uint64_t frame = 0; // DPL timeslice; row IDs are scoped to this frame and job.
};
struct Deltas {
  int64_t bc = 0;
  double phi = 0., radius = 0., attempts = 0., z = 0.;
};
inline bool compatible(const Event& target, const Event& donor, const Config& c, Deltas& d)
{
  if (target.collision < 0 || donor.collision < 0 || (target.frame == donor.frame && target.collision == donor.collision) ||
      target.run != donor.run || donor.attempts <= 0) {
    return false;
  }
  if (c.maxDeltaMftMultiplicityRel >= 0.) {
    if (target.mftMultiplicity < 0 || donor.mftMultiplicity < 0) {
      return false;
    }
    if (target.mftMultiplicity == 0) {
      if (donor.mftMultiplicity != 0) {
        return false;
      }
    } else {
      const double relative = std::abs(static_cast<double>(donor.mftMultiplicity) - target.mftMultiplicity) / target.mftMultiplicity;
      if (c.mftMultiplicityStrict ? relative >= c.maxDeltaMftMultiplicityRel : relative > c.maxDeltaMftMultiplicityRel) {
        return false;
      }
    }
  }
  if (c.maxDeltaCbMultiplicityRel >= 0.) {
    if (target.cbMultiplicity < 0 || donor.cbMultiplicity < 0 ||
        (target.cbMultiplicity == 0 ? donor.cbMultiplicity != 0 : (c.cbMultiplicityStrict ? std::abs(static_cast<double>(donor.cbMultiplicity) - target.cbMultiplicity) >= c.maxDeltaCbMultiplicityRel * target.cbMultiplicity : std::abs(static_cast<double>(donor.cbMultiplicity) - target.cbMultiplicity) > c.maxDeltaCbMultiplicityRel * target.cbMultiplicity))) {
      return false;
    }
  }
  if (c.maxDeltaMchMultiplicity >= 0 &&
      (target.mchMultiplicity < 0 || donor.mchMultiplicity < 0 ||
       std::abs(donor.mchMultiplicity - target.mchMultiplicity) > c.maxDeltaMchMultiplicity)) {
    return false;
  }
  if (c.maxDeltaCentrality >= 0. &&
      (!std::isfinite(target.centrality) || !std::isfinite(donor.centrality) ||
       target.centrality < 0. || target.centrality > 100. || donor.centrality < 0. || donor.centrality > 100. ||
       std::abs(donor.centrality - target.centrality) > c.maxDeltaCentrality)) {
    return false;
  }
  if (c.maxDeltaMchMomentumRel >= 0. &&
      (!std::isfinite(target.mchMomentumPV) || !std::isfinite(donor.mchMomentumPV) ||
       target.mchMomentumPV <= 0. || donor.mchMomentumPV <= 0. ||
       std::abs(donor.mchMomentumPV - target.mchMomentumPV) > c.maxDeltaMchMomentumRel * target.mchMomentumPV)) {
    return false;
  }
  if (c.maxDeltaMchPtRel >= 0. &&
      (!std::isfinite(target.mchPtPV) || !std::isfinite(donor.mchPtPV) ||
       target.mchPtPV <= 0. || donor.mchPtPV <= 0. ||
       std::abs(donor.mchPtPV - target.mchPtPV) > c.maxDeltaMchPtRel * target.mchPtPV)) {
    return false;
  }
  // Select donor MCH directions, never the borrowed MFT or the mixed-pair residual.
  if (c.maxDeltaMchPhiPV >= 0. &&
      (!std::isfinite(target.mchPhiPV) || !std::isfinite(donor.mchPhiPV) ||
       std::abs(std::remainder(donor.mchPhiPV - target.mchPhiPV, 2. * 3.14159265358979323846)) > c.maxDeltaMchPhiPV)) {
    return false;
  }
  if (c.maxDeltaMchEtaPV >= 0. &&
      (!std::isfinite(target.mchEtaPV) || !std::isfinite(donor.mchEtaPV) ||
       (c.mchEtaStrict ? std::abs(donor.mchEtaPV - target.mchEtaPV) >= c.maxDeltaMchEtaPV : std::abs(donor.mchEtaPV - target.mchEtaPV) > c.maxDeltaMchEtaPV))) {
    return false;
  }
  d.bc = target.bc >= donor.bc ? target.bc - donor.bc : donor.bc - target.bc;
  d.phi = std::abs(std::remainder(donor.phi - target.phi, 2. * 3.14159265358979323846));
  d.radius = donor.radius - target.radius;
  d.attempts = target.attempts > 0 ? static_cast<double>(donor.attempts - target.attempts) / target.attempts : 0.;
  d.z = std::abs(donor.z - target.z);
  return d.bc >= c.minDeltaBc &&
         (c.maxDeltaPhi < 0. || (std::isfinite(d.phi) && d.phi <= c.maxDeltaPhi)) &&
         (c.maxDeltaZ < 0. || (std::isfinite(d.z) && (c.vertexZStrict ? d.z < c.maxDeltaZ : d.z <= c.maxDeltaZ))) &&
         (c.maxDeltaR < 0. || (std::isfinite(d.radius) && (c.mchRadiusStrict ? (c.symmetricCuts ? std::abs(d.radius) : d.radius) < c.maxDeltaR : (c.symmetricCuts ? std::abs(d.radius) : d.radius) <= c.maxDeltaR))) &&
         (c.maxDeltaAttempts < 0. || (std::isfinite(d.attempts) && (c.symmetricCuts ? std::abs(d.attempts) : d.attempts) <= c.maxDeltaAttempts));
}
} // namespace o2::analysis::mftmch::mixing

// -----------------------------------------------------------------------------
// MftMchMixingHistory
// -----------------------------------------------------------------------------
#include <deque>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <utility>
namespace o2::analysis::mftmch::mixing
{
// Item must own its state and expose stateCount. Never store AO2D iterators.
template <class Item>
struct History {
  std::deque<Item> events;
  size_t states = 0, evicted = 0, oversized = 0;
  int run = -1;
  double solenoid = NAN, dipole = NAN;
  void clear()
  {
    events.clear();
    states = 0;
  }
  bool setRun(int value)
  {
    if (run == value)
      return false;
    clear();
    run = value;
    solenoid = dipole = NAN;
    return true;
  }
  bool setField(double sol, double dip)
  {
    if (solenoid == sol && dipole == dip)
      return false;
    const bool reset = std::isfinite(solenoid) || !events.empty();
    clear();
    solenoid = sol;
    dipole = dip;
    return reset;
  }
  bool push(Item item, size_t maxEvents, size_t maxStates)
  {
    if (item.stateCount > maxStates || maxEvents == 0) {
      ++oversized;
      return false;
    }
    while (!events.empty() && (events.size() >= maxEvents || states + item.stateCount > maxStates)) {
      states -= events.front().stateCount;
      events.pop_front();
      ++evicted;
    }
    states += item.stateCount;
    events.push_back(std::move(item));
    return true;
  }
};
} // namespace o2::analysis::mftmch::mixing

#include <Framework/TimingInfo.h>
#include <deque>
#include <fstream>
#include <sstream>
#include <random>
// -----------------------------------------------------------------------------
// MftMchMixingAlignment
// -----------------------------------------------------------------------------
#include <array>
#include <algorithm>
#include <cmath>
namespace o2::analysis::mftmch::alignment
{
using Vec = std::array<double, 3>;
using Rotation = std::array<double, 9>;
inline double dot(const Vec& a, const Vec& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec cross(const Vec& a, const Vec& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
inline Vec rotate(const Rotation& r, const Vec& v)
{
  Vec out{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out[i] += r[3 * i + j] * v[j];
  return out;
}
// Shortest proper rotation (no extra roll) from the donor to the target direction.
inline bool directionRotation(Vec donor, Vec target, Rotation& r)
{
  double a = std::sqrt(dot(donor, donor)), b = std::sqrt(dot(target, target));
  if (!std::isfinite(a) || !std::isfinite(b) || a <= 0. || b <= 0.)
    return false;
  for (int i = 0; i < 3; ++i) {
    donor[i] /= a;
    target[i] /= b;
  }
  double c = std::clamp(dot(donor, target), -1., 1.);
  r = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  if (c < -1. + 1.e-12) {
    Vec axis = cross(donor, std::abs(donor[0]) < .8 ? Vec{1., 0., 0.} : Vec{0., 1., 0.});
    double n = std::sqrt(dot(axis, axis));
    for (double& v : axis)
      v /= n;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r[3 * i + j] = 2. * axis[i] * axis[j] - (i == j ? 1. : 0.);
    return true;
  }
  Vec v = cross(donor, target);
  Rotation k{0., -v[2], v[1], v[2], 0., -v[0], -v[1], v[0], 0.};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      r[3 * i + j] += k[3 * i + j];
      for (int l = 0; l < 3; ++l)
        r[3 * i + j] += k[3 * i + l] * k[3 * l + j] / (1. + c);
    }
  return true;
}
inline Vec transformPosition(const Rotation& r, const Vec& position, const Vec& donor, const Vec& target)
{
  Vec v{};
  for (int i = 0; i < 3; ++i)
    v[i] = position[i] - donor[i];
  v = rotate(r, v);
  for (int i = 0; i < 3; ++i)
    v[i] += target[i];
  return v;
}
} // namespace o2::analysis::mftmch::alignment

#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/MatchPlaneResiduals.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/MatchMFTFT0.h"
#include "Common/DataModel/MatchMFTMuonData.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DataFormatsGlobalTracking/RecoContainer.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/InitContext.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TDatabasePDG.h>
#include <TGeoGlobalMagField.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSymEigen.h>
#include <TVector2.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <vector>

namespace o2::analysis::mftmch
{
using namespace std;

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = aod::Collisions;
using MyBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedToFT0>;
using MyMUONs = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMFTs = aod::MFTTracks;

using MyCollision = MyCollisions::iterator;
using MyBC = MyBCs::iterator;
using MyMUON = MyMUONs::iterator;
using MyMFT = MyMFTs::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

inline float mMu = TDatabasePDG::Instance()->GetParticle(13)->Mass();
inline int mRunNumber;
/*
  TLorentzVector muon1LV;
  TLorentzVector muon2LV;
  TLorentzVector dimuonLV;
*/
inline unordered_map<int, vector<int64_t>> map_mfttracks;
inline unordered_map<int, vector<int64_t>> map_muontracks;
inline unordered_map<int, bool> map_collisions;
inline unordered_map<int, bool> map_has_mfttracks_collisions;
inline unordered_map<int, bool> map_has_muontracks_collisions;
inline unordered_map<int, float> map_vtxz;
inline unordered_map<int, int> map_nmfttrack;

// Adapt raw FIT records to the found-detector interface used by UDHelpers.
// Missing detectors pass, exactly as in udhelpers::cleanFIT.
struct MatchingFITRecord {
  aod::FV0As const& fv0s;
  aod::FT0s const& ft0s;
  aod::FDDs const& fdds;
  int64_t fv0 = -1, ft0 = -1, fdd = -1;
  bool has_foundFV0() const { return fv0 >= 0; }
  bool has_foundFT0() const { return ft0 >= 0; }
  bool has_foundFDD() const { return fdd >= 0; }
  auto foundFV0() const { return fv0s.rawIteratorAt(fv0); }
  auto foundFT0() const { return ft0s.rawIteratorAt(ft0); }
  auto foundFDD() const { return fdds.rawIteratorAt(fdd); }
};

// Downstream counterpart of O2 TrackExtrap::extrapToVertexWithoutBranson.
// Material/ELoss/MCS helpers below are adapted from O2 TrackExtrap.cxx,
// revision 024f05624a, retaining its linear absorber-path / fixed-p budget
// approximation, two half-loss corrections, and covariance conventions.
// Changes: upstream -> downstream, E + loss -> E - loss; add MCS at the EXIT.
// This is a study extension, not an official O2 downstream transport API.
struct MftAbsorberExtrap {
  using TrackParam = o2::mch::TrackParam;
  static constexpr double SMuMass = 0.105658; // Same mass as O2 TrackExtrap
  static constexpr double AbsorberBegin = -90.;
  static constexpr double AbsorberEnd = -505.;
  enum class Status { Success,
                      InvalidInput,
                      FieldPropagation,
                      Geometry,
                      RangeOut,
                      Nonfinite,
                      StepLimit,
                      Count };
  struct NoiseOptions {
    double lossPull = 0.;
    bool includeLossVariance = true;
    bool includeMCS = true;
  };
  struct TransportInfo {
    double pathLength = 0., energyLoss = 0., momentumIn = 0., momentumOut = 0.;
    double energyLossSigma = 0.;
    int materialSteps = 0;
    Status status = Status::Success;
  };

  static bool finite(const TrackParam& track)
  {
    if (!std::isfinite(track.getZ()) || !track.hasCovariances() || !std::isfinite(track.p()) || track.p() <= 0.) {
      return false;
    }
    for (int i = 0; i < 5; ++i) {
      if (!std::isfinite(track.getParameters()(i, 0))) {
        return false;
      }
      for (int j = 0; j < 5; ++j) {
        if (!std::isfinite(track.getCovariances()(i, j))) {
          return false;
        }
      }
    }
    return true;
  }

  // One material interval. With the full -90 -> -505 interval, this retains
  // the previous whole-absorber operation order and representative momentum.
  static bool transportMaterialSlice(TrackParam& track, double zTarget, TransportInfo& result, const NoiseOptions* noise = nullptr)
  {
    TrackParam transported(track);
    // 2. Field-only pilot determines the exit of the straight material chord,
    //    analogous to the unconstrained upstream pilot in official TrackExtrap.
    const double xyzIn[3] = {transported.getNonBendingCoor(), transported.getBendingCoor(), transported.getZ()};
    TrackParam pilot(transported);
    if (!o2::mch::TrackExtrap::extrapToZ(pilot, zTarget)) {
      result.status = Status::FieldPropagation;
      return false;
    }
    double in[3] = {xyzIn[0], xyzIn[1], xyzIn[2]};
    double out[3] = {pilot.getNonBendingCoor(), pilot.getBendingCoor(), zTarget};
    double pathLength = 0., f0 = 0., f1 = 0., f2 = 0., meanRho = 0., loss = 0., lossVariance = 0.;
    const double momentumIn = transported.p();
    if (!getAbsorberCorrectionParam(in, out, momentumIn, pathLength, f0, f1, f2, meanRho, loss, lossVariance) ||
        !std::isfinite(pathLength) || !std::isfinite(loss) || loss < 0. ||
        !std::isfinite(lossVariance) || lossVariance < 0. ||
        !std::isfinite(f0) || !std::isfinite(f1) || !std::isfinite(f2) || f0 < 0.) {
      result.status = Status::Geometry;
      return false;
    }
    // O2's Landau-width approximation adds sigma (not variance) along the
    // material path. Preserve that convention across user intervals: adding
    // independent sigma_step^2 would make straggling vanish as step -> 0.
    const double sigmaStep = std::sqrt(lossVariance);
    if (result.materialSteps > 0) {
      lossVariance = sigmaStep * (2. * result.energyLossSigma + sigmaStep);
    }
    if (noise) {
      loss += noise->lossPull * sigmaStep;
      if (!noise->includeLossVariance)
        lossVariance = 0.;
    }
    // 3. Range-out guard. Unlike upstream reconstruction, E must DECREASE.
    if (std::hypot(momentumIn, SMuMass) - loss <= SMuMass) {
      result.status = Status::RangeOut;
      return false;
    }
    // 4. Half loss -> field-map/covariance propagation -> MCS -> half loss.
    //    correctELossEffectInAbsorber accepts signed ELoss: negative subtracts E.
    //    The variance is POSITIVE in both propagation directions.
    correctELossEffectInAbsorber(transported, -0.5 * loss, 0.5 * lossVariance);
    if (!finite(transported) || !o2::mch::TrackExtrap::extrapToZCov(transported, zTarget)) {
      result.status = Status::FieldPropagation;
      return false;
    }
    // The O2 MCS helper assumes parameters at zOut. zOut-zIn is negative
    // for a forward muon: position/slope cross terms must therefore be negative.
    if (!noise || noise->includeMCS)
      addMCSEffectInAbsorber(transported, -pathLength, f0, f1, f2);
    correctELossEffectInAbsorber(transported, -0.5 * loss, 0.5 * lossVariance);
    if (!finite(transported)) {
      result.status = Status::Nonfinite;
      return false;
    }
    result.pathLength += pathLength;
    result.energyLoss += loss;
    result.energyLossSigma += sigmaStep;
    ++result.materialSteps;
    track = transported;
    return true;
  }

  // Shorten a positive-step interval at its first material boundary. A small
  // navigation margin crosses the boundary, avoiding repeated zero steps.
  // The material budget still integrates every medium crossed by the chord.
  static bool limitAtMaterialBoundary(const TrackParam& track, double& zTarget, TransportInfo& result)
  {
    TrackParam pilot(track);
    if (!o2::mch::TrackExtrap::extrapToZ(pilot, zTarget)) {
      result.status = Status::FieldPropagation;
      return false;
    }
    double point[3] = {track.getNonBendingCoor(), track.getBendingCoor(), track.getZ()};
    double direction[3] = {pilot.getNonBendingCoor() - point[0], pilot.getBendingCoor() - point[1], zTarget - point[2]};
    double length = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
    if (!std::isfinite(length) || length <= 0.) {
      result.status = Status::Nonfinite;
      return false;
    }
    for (double& component : direction) {
      component /= length;
    }
    if (!gGeoManager->InitTrack(point, direction)) {
      result.status = Status::Geometry;
      return false;
    }
    gGeoManager->FindNextBoundary(length);
    double boundary = gGeoManager->GetStep();
    if (!std::isfinite(boundary) || boundary < 0.) {
      result.status = Status::Geometry;
      return false;
    }
    constexpr double navigationMargin = 1.e-6; // cm, same scale as official material integration
    if (boundary + navigationMargin < length) {
      double fraction = std::min(1., (boundary + navigationMargin) / length);
      zTarget = point[2] + direction[2] * length * fraction;
    }
    return true;
  }

  // Input: original MFT converted to TrackParam upstream of absorber.
  // stepCm==0: previous whole-absorber model, unchanged at zEnd=-505.
  // stepCm>0: maximum |delta z| per material update, shortened at boundaries.
  // RK magnetic integration remains adaptive inside every interval.
  // Additional material downstream of -505 is intentionally NOT modelled here.
  // Transactional: leave input unchanged on any failure.
  static bool extrapToAbsorberEndWithoutBranson(TrackParam& track, double zEnd = AbsorberEnd,
                                                TransportInfo* info = nullptr, double stepCm = 0., const NoiseOptions* noise = nullptr)
  {
    TransportInfo result;
    auto fail = [&](Status status) {
      result.status = status;
      if (info) {
        *info = result;
      }
      return false;
    };
    if (!finite(track) || track.getZ() < AbsorberBegin || !std::isfinite(zEnd) || zEnd > AbsorberEnd ||
        !std::isfinite(stepCm) || stepCm < 0.) {
      return fail(Status::InvalidInput);
    }
    if (!gGeoManager) {
      return fail(Status::Geometry);
    }
    TrackParam transported(track);
    if (!o2::mch::TrackExtrap::extrapToZCov(transported, AbsorberBegin) || !finite(transported)) {
      return fail(Status::FieldPropagation);
    }
    result.momentumIn = transported.p();
    constexpr int maxMaterialSteps = 1000000;
    do {
      if (result.materialSteps >= maxMaterialSteps) {
        return fail(Status::StepLimit);
      }
      double nextZ = stepCm == 0. ? AbsorberEnd : std::max(AbsorberEnd, transported.getZ() - stepCm);
      if (stepCm > 0. && !limitAtMaterialBoundary(transported, nextZ, result)) {
        return fail(result.status);
      }
      if (!(nextZ < transported.getZ())) {
        return fail(Status::StepLimit);
      }
      if (!transportMaterialSlice(transported, nextZ, result, noise)) {
        return fail(result.status);
      }
    } while (transported.getZ() > AbsorberEnd);
    if (zEnd < AbsorberEnd && !o2::mch::TrackExtrap::extrapToZCov(transported, zEnd)) {
      return fail(Status::FieldPropagation);
    }
    if (!finite(transported)) {
      return fail(Status::Nonfinite);
    }
    result.momentumOut = transported.p();
    track = transported;
    if (info) {
      *info = result;
    }
    return true;
  }

  // Supporting calculations copied from the official private helpers.
  static bool getAbsorberCorrectionParam(double trackXYZIn[3], double trackXYZOut[3], double pTotal,
                                         double& pathLength, double& f0, double& f1, double& f2,
                                         double& meanRho, double& totalELoss, double& sigmaELoss2)
  {
    /// Parameters used to correct for Multiple Coulomb Scattering and energy loss in absorber
    /// Calculated assuming a linear propagation from trackXYZIn to trackXYZOut (order is important)
    /// pathLength:  path length between trackXYZIn and trackXYZOut (cm)
    /// f0:          0th moment of z calculated with the inverse radiation-length distribution
    /// f1:          1st moment of z calculated with the inverse radiation-length distribution
    /// f2:          2nd moment of z calculated with the inverse radiation-length distribution
    /// meanRho:     average density of crossed material (g/cm3)
    /// totalELoss:  total energy loss in absorber
    /// sigmaELoss2: square of energy loss fluctuation in absorber

    // Check whether the geometry is available
    if (!gGeoManager) {
      LOG(warning) << "geometry is missing";
      return false;
    }

    // Initialize starting point and direction
    pathLength = TMath::Sqrt((trackXYZOut[0] - trackXYZIn[0]) * (trackXYZOut[0] - trackXYZIn[0]) +
                             (trackXYZOut[1] - trackXYZIn[1]) * (trackXYZOut[1] - trackXYZIn[1]) +
                             (trackXYZOut[2] - trackXYZIn[2]) * (trackXYZOut[2] - trackXYZIn[2]));
    if (pathLength < TGeoShape::Tolerance()) {
      LOG(warning) << "path length is too small";
      return false;
    }
    double b[3] = {(trackXYZOut[0] - trackXYZIn[0]) / pathLength, (trackXYZOut[1] - trackXYZIn[1]) / pathLength, (trackXYZOut[2] - trackXYZIn[2]) / pathLength};
    TGeoNode* currentnode = gGeoManager->InitTrack(trackXYZIn, b);
    if (!currentnode) {
      LOG(warning) << "starting point out of geometry";
      return false;
    }

    // loop over absorber slices and calculate absorber's parameters
    f0 = f1 = f2 = meanRho = totalELoss = 0.;
    double sigmaELoss(0.);
    double zB = trackXYZIn[2];
    double remainingPathLength = pathLength;
    do {

      // Get material properties
      TGeoMaterial* material = currentnode->GetVolume()->GetMedium()->GetMaterial();
      double rho = material->GetDensity(); // material density (g/cm3)
      double x0 = material->GetRadLen();   // radiation-length (cm-1)
      double atomicZ = material->GetZ();   // Z of material
      double atomicZoverA(0.);             // Z/A of material
      if (material->IsMixture()) {
        TGeoMixture* mixture = static_cast<TGeoMixture*>(material);
        double sum(0.);
        for (int iel = 0; iel < mixture->GetNelements(); ++iel) {
          sum += mixture->GetWmixt()[iel];
          atomicZoverA += mixture->GetWmixt()[iel] * mixture->GetZmixt()[iel] / mixture->GetAmixt()[iel];
        }
        atomicZoverA /= sum;
      } else {
        atomicZoverA = atomicZ / material->GetA();
      }

      // Get path length within this material
      gGeoManager->FindNextBoundary(remainingPathLength);
      double localPathLength = gGeoManager->GetStep() + 1.e-6;
      // Check if boundary within remaining path length. If so, make sure to cross the boundary to prepare the next step
      if (localPathLength >= remainingPathLength) {
        localPathLength = remainingPathLength;
      } else {
        currentnode = gGeoManager->Step();
        if (!currentnode) {
          LOG(warning) << "navigation failed";
          return false;
        }
        if (!gGeoManager->IsEntering()) {
          // make another small step to try to enter in new absorber slice
          gGeoManager->SetStep(0.001);
          currentnode = gGeoManager->Step();
          if (!gGeoManager->IsEntering() || !currentnode) {
            LOG(warning) << "navigation failed";
            return false;
          }
          localPathLength += 0.001;
        }
      }

      // calculate absorber's parameters
      double zE = b[2] * localPathLength + zB;
      double dzB = zB - trackXYZIn[2];
      double dzE = zE - trackXYZIn[2];
      f0 += localPathLength / x0;
      f1 += (dzE * dzE - dzB * dzB) / b[2] / b[2] / x0 / 2.;
      f2 += (dzE * dzE * dzE - dzB * dzB * dzB) / b[2] / b[2] / b[2] / x0 / 3.;
      meanRho += localPathLength * rho;
      // Vacuum has no energy loss; avoid log(0) in the official formula.
      if (rho > 0.) {
        totalELoss += betheBloch(pTotal, localPathLength, rho, atomicZ, atomicZoverA);
        sigmaELoss += energyLossFluctuation(pTotal, localPathLength, rho, atomicZoverA);
      }

      // prepare next step
      zB = zE;
      remainingPathLength -= localPathLength;
    } while (remainingPathLength > TGeoShape::Tolerance());

    meanRho /= pathLength;
    sigmaELoss2 = sigmaELoss * sigmaELoss;

    return true;
  }

  static void addMCSEffectInAbsorber(TrackParam& param, double signedPathLength, double f0, double f1, double f2)
  {
    /// Add to the track parameter covariances the effects of multiple Coulomb scattering
    /// signedPathLength must have the sign of (zOut - zIn) where all other parameters are assumed to be given at zOut.

    // absorber related covariance parameters
    double bendingSlope = param.getBendingSlope();
    double nonBendingSlope = param.getNonBendingSlope();
    double inverseBendingMomentum = param.getInverseBendingMomentum();
    double alpha2 = 0.0136 * 0.0136 * inverseBendingMomentum * inverseBendingMomentum * (1.0 + bendingSlope * bendingSlope) /
                    (1.0 + bendingSlope * bendingSlope + nonBendingSlope * nonBendingSlope); // velocity = 1
    double pathLength = TMath::Abs(signedPathLength);
    double varCoor = alpha2 * (pathLength * pathLength * f0 - 2. * pathLength * f1 + f2);
    double covCorrSlope = TMath::Sign(1., signedPathLength) * alpha2 * (pathLength * f0 - f1);
    double varSlop = alpha2 * f0;

    // Set MCS covariance matrix
    TMatrixD newParamCov(param.getCovariances());
    // Non bending plane
    newParamCov(0, 0) += varCoor;
    newParamCov(0, 1) += covCorrSlope;
    newParamCov(1, 0) += covCorrSlope;
    newParamCov(1, 1) += varSlop;
    // Bending plane
    newParamCov(2, 2) += varCoor;
    newParamCov(2, 3) += covCorrSlope;
    newParamCov(3, 2) += covCorrSlope;
    newParamCov(3, 3) += varSlop;

    // Set momentum related covariances if B!=0
    if (o2::mch::TrackExtrap::isFieldON()) {
      // compute derivative d(q/Pxy) / dSlopeX and d(q/Pxy) / dSlopeY
      double dqPxydSlopeX = inverseBendingMomentum * nonBendingSlope / (1. + nonBendingSlope * nonBendingSlope + bendingSlope * bendingSlope);
      double dqPxydSlopeY = -inverseBendingMomentum * nonBendingSlope * nonBendingSlope * bendingSlope /
                            (1. + bendingSlope * bendingSlope) / (1. + nonBendingSlope * nonBendingSlope + bendingSlope * bendingSlope);
      // Inverse bending momentum (due to dependences with bending and non bending slopes)
      newParamCov(4, 0) += dqPxydSlopeX * covCorrSlope;
      newParamCov(0, 4) += dqPxydSlopeX * covCorrSlope;
      newParamCov(4, 1) += dqPxydSlopeX * varSlop;
      newParamCov(1, 4) += dqPxydSlopeX * varSlop;
      newParamCov(4, 2) += dqPxydSlopeY * covCorrSlope;
      newParamCov(2, 4) += dqPxydSlopeY * covCorrSlope;
      newParamCov(4, 3) += dqPxydSlopeY * varSlop;
      newParamCov(3, 4) += dqPxydSlopeY * varSlop;
      newParamCov(4, 4) += (dqPxydSlopeX * dqPxydSlopeX + dqPxydSlopeY * dqPxydSlopeY) * varSlop;
    }

    // Set new covariances
    param.setCovariances(newParamCov);
  }

  static double betheBloch(double pTotal, double pathLength, double rho, double atomicZ, double atomicZoverA)
  {
    /// Returns the mean total momentum energy loss of muon with total momentum='pTotal'
    /// in the absorber layer of lenght='pathLength', density='rho', Z='atomicZ' and mean Z/A='atomicZoverA'
    /// This is the parameterization of the Bethe-Bloch formula inspired by Geant.

    static constexpr double mK = 0.307075e-3; // [GeV*cm^2/g]
    static constexpr double me = 0.511e-3;    // [GeV/c^2]
    static constexpr double x0 = 0.2 * 2.303; // density effect first junction point * 2.303
    static constexpr double x1 = 3. * 2.303;  // density effect second junction point * 2.303

    double bg = pTotal / SMuMass; // beta*gamma
    double bg2 = bg * bg;
    double maxT = 2 * me * bg2; // neglecting the electron mass

    // mean exitation energy (GeV)
    double mI = (atomicZ < 13) ? (12. * atomicZ + 7.) * 1.e-9 : (9.76 * atomicZ + 58.8 * TMath::Power(atomicZ, -0.19)) * 1.e-9;

    // density effect
    double x = TMath::Log(bg);
    double lhwI = TMath::Log(28.816 * 1e-9 * TMath::Sqrt(rho * atomicZoverA) / mI);
    double d2(0.);
    if (x > x1) {
      d2 = lhwI + x - 0.5;
    } else if (x > x0) {
      double r = (x1 - x) / (x1 - x0);
      d2 = lhwI + x - 0.5 + (0.5 - lhwI - x0) * r * r * r;
    }

    return pathLength * rho * (mK * atomicZoverA * (1 + bg2) / bg2 * (0.5 * TMath::Log(2 * me * bg2 * maxT / (mI * mI)) - bg2 / (1 + bg2) - d2));
  }

  static double energyLossFluctuation(double pTotal, double pathLength, double rho, double atomicZoverA)
  {
    /// Returns the total momentum energy loss fluctuation of muon with total momentum='pTotal'
    /// in the absorber layer of lenght='pathLength', density='rho', A='atomicA' and Z='atomicZ'

    static constexpr double mK = 0.307075e-3; // GeV.g^-1.cm^2

    double p2 = pTotal * pTotal;
    double beta2 = p2 / (p2 + SMuMass * SMuMass);

    double fwhm = 2. * mK * rho * pathLength * atomicZoverA / beta2; // FWHM of the energy loss Landau distribution

    return fwhm / TMath::Sqrt(8. * log(2.)); // gaussian: fwmh = 2 * srqt(2*ln(2)) * sigma (i.e. fwmh = 2.35 * sigma)
  }

  static void correctELossEffectInAbsorber(TrackParam& param, double eLoss, double sigmaELoss2)
  {
    /// Correct parameters for energy loss and add energy loss fluctuation effect to covariances

    // Get parameter covariances in (X, SlopeX, Y, SlopeY, q*PTot) coordinate system
    TMatrixD newParamCov(param.getCovariances());
    cov2CovP(param.getParameters(), newParamCov);

    // Compute new parameters corrected for energy loss
    double p = param.p();
    double e = TMath::Sqrt(p * p + SMuMass * SMuMass);
    double eCorr = e + eLoss;
    double pCorr = TMath::Sqrt(eCorr * eCorr - SMuMass * SMuMass);
    double nonBendingSlope = param.getNonBendingSlope();
    double bendingSlope = param.getBendingSlope();
    param.setInverseBendingMomentum(param.getCharge() / pCorr *
                                    TMath::Sqrt(1.0 + nonBendingSlope * nonBendingSlope + bendingSlope * bendingSlope) /
                                    TMath::Sqrt(1.0 + bendingSlope * bendingSlope));

    // Add effects of energy loss fluctuation to covariances
    newParamCov(4, 4) += eCorr * eCorr / pCorr / pCorr * sigmaELoss2;

    // Get new parameter covariances in (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
    covP2Cov(param.getParameters(), newParamCov);

    // Set new parameter covariances
    param.setCovariances(newParamCov);
  }

  static void cov2CovP(const TMatrixD& param, TMatrixD& cov)
  {
    /// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, SlopeX, Y, SlopeY, q*PTot)
    /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system

    // charge * total momentum
    double qPTot = TMath::Sqrt(1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0)) /
                   TMath::Sqrt(1. + param(3, 0) * param(3, 0)) / param(4, 0);

    // Jacobian of the opposite transformation
    TMatrixD jacob(5, 5);
    jacob.UnitMatrix();
    jacob(4, 1) = qPTot * param(1, 0) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 3) = -qPTot * param(1, 0) * param(1, 0) * param(3, 0) /
                  (1. + param(3, 0) * param(3, 0)) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 4) = -qPTot / param(4, 0);

    // compute covariances in new coordinate system
    TMatrixD tmp(cov, TMatrixD::kMultTranspose, jacob);
    cov.Mult(jacob, tmp);
  }

  static void covP2Cov(const TMatrixD& param, TMatrixD& covP)
  {
    /// change coordinate system: (X, SlopeX, Y, SlopeY, q*PTot) -> (X, SlopeX, Y, SlopeY, q/Pyz)
    /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system

    // charge * total momentum
    double qPTot = TMath::Sqrt(1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0)) /
                   TMath::Sqrt(1. + param(3, 0) * param(3, 0)) / param(4, 0);

    // Jacobian of the transformation
    TMatrixD jacob(5, 5);
    jacob.UnitMatrix();
    jacob(4, 1) = param(4, 0) * param(1, 0) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 3) = -param(4, 0) * param(1, 0) * param(1, 0) * param(3, 0) /
                  (1. + param(3, 0) * param(3, 0)) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 4) = -param(4, 0) / qPTot;

    // compute covariances in new coordinate system
    TMatrixD tmp(covP, TMatrixD::kMultTranspose, jacob);
    covP.Mult(jacob, tmp);
  }
};

// Conditional MFT prediction using this candidate MCH's momentum. Input detector
// errors are independent; endpoint errors are NOT. Shared absorber energy loss
// uses one common unit-normal pull; MCS is introduced once on the forward path.
struct MftMchJointTransport {
  using Fwd = o2::track::TrackParCovFwd;
  using Conversion = o2::aod::fwdtrackutils::GlobalMuonMatchingPlane;
  struct State {
    Fwd seed, mft, mch;
    double sigmaLossUp = 0.;
  };
  struct Result {
    State state;
    TMatrixD joint{10, 10}, residual{5, 5}, seedCov{5, 5}, seedMchCross{5, 5}, jacobian{10, 10};
    bool hasNominalState = false;
    double chiAll = std::numeric_limits<double>::quiet_NaN(), chiFour = std::numeric_limits<double>::quiet_NaN(), chiXY = std::numeric_limits<double>::quiet_NaN();
    int rankAll = 0, rankFour = 0, rankXY = 0, status = 0;
  };

  static Fwd perturb(const Fwd& track, int index, double amount)
  {
    Fwd out(track);
    auto parameters = out.getParameters();
    parameters[index] += amount;
    out.setParameters(parameters);
    return out;
  }
  static double difference(double a, double b, int index)
  {
    return index % 5 == 2 ? TVector2::Phi_mpi_pi(a - b) : a - b;
  }
  static bool evaluate(Conversion& conversion, const Fwd& mft, const Fwd& mch,
                       double zEnd, double step, bool mchCharge, double pull, State& result,
                       double downstreamPull = std::numeric_limits<double>::quiet_NaN())
  {
    if (mft.getZ() < -90. || mft.getTanl() >= 0. || !std::isfinite(mch.getInvQPt()) || mch.getInvQPt() == 0.) {
      return false;
    }
    auto up = conversion.fwdToMch(mch);
    if (!o2::mch::TrackExtrap::extrapToVertexWithoutBranson(up, mft.getZ())) {
      return false;
    }
    // Reproduce the official unconstrained upstream material chord for sigmaE.
    auto exit = conversion.fwdToMch(mch);
    if (!o2::mch::TrackExtrap::extrapToZ(exit, -505.)) {
      return false;
    }
    auto entrance = exit;
    if (!o2::mch::TrackExtrap::extrapToZ(entrance, -90.)) {
      return false;
    }
    double in[3] = {entrance.getNonBendingCoor(), entrance.getBendingCoor(), -90.};
    double out[3] = {exit.getNonBendingCoor(), exit.getBendingCoor(), -505.};
    double length = 0., f0 = 0., f1 = 0., f2 = 0., rho = 0., loss = 0., variance = 0.;
    if (!MftAbsorberExtrap::getAbsorberCorrectionParam(in, out, exit.p(), length, f0, f1, f2, rho, loss, variance) ||
        !std::isfinite(variance) || variance < 0.) {
      return false;
    }
    result.sigmaLossUp = std::sqrt(variance);
    double energy = std::hypot(up.p(), MftAbsorberExtrap::SMuMass) + pull * result.sigmaLossUp;
    if (!std::isfinite(energy) || energy <= MftAbsorberExtrap::SMuMass) {
      return false;
    }
    double p = std::sqrt((energy - MftAbsorberExtrap::SMuMass) * (energy + MftAbsorberExtrap::SMuMass));
    result.seed = mft;
    double charge = std::copysign(1., mchCharge ? mch.getInvQPt() : mft.getInvQPt());
    result.seed.setInvQPt(charge * std::hypot(1., mft.getTanl()) / p);
    // This covariance is a PROCESS-only forward covariance. Measurement errors
    // are propagated below through the joint Jacobian, never added twice here.
    o2::track::SMatrix55Sym zero;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j <= i; ++j)
        zero(i, j) = 0.;
    result.seed.setCovariances(zero);
    auto forward = conversion.fwdToMch(result.seed);
    MftAbsorberExtrap::NoiseOptions noise;
    noise.lossPull = std::isfinite(downstreamPull) ? downstreamPull : pull;
    noise.includeLossVariance = false;
    noise.includeMCS = true;
    if (!MftAbsorberExtrap::extrapToAbsorberEndWithoutBranson(forward, zEnd, nullptr, step, &noise)) {
      return false;
    }
    result.mft = conversion.mchToFwd(forward);
    auto comparison = conversion.fwdToMch(mch);
    if (!o2::mch::TrackExtrap::extrapToZCov(comparison, zEnd)) {
      return false;
    }
    result.mch = conversion.mchToFwd(comparison);
    for (int i = 0; i < 5; ++i) {
      if (!std::isfinite(result.mft.getParameters()[i]) || !std::isfinite(result.mch.getParameters()[i]))
        return false;
    }
    return true;
  }

  // Use the official matching prescription: directly invert the residual
  // covariance, without a PSD gate or eigenvalue clipping. Preserve the
  // pair correlations already included in the downstream residual covariance.
  template <unsigned int Dimension>
  static bool chiDirect(const TMatrixD& covariance, const std::vector<double>& delta,
                        double& value, int& rank)
  {
    ROOT::Math::SMatrix<double, Dimension> inverse;
    ROOT::Math::SVector<double, Dimension> residual;
    for (unsigned int i = 0; i < Dimension; ++i) {
      residual[i] = delta[i];
      if (!std::isfinite(residual[i]))
        return false;
      for (unsigned int j = 0; j < Dimension; ++j) {
        inverse(i, j) = covariance(i, j);
        if (!std::isfinite(inverse(i, j)))
          return false;
      }
    }
    inverse.Invert(); // As in official matchALL: no PSD/inversion-status gate.
    value = ROOT::Math::Similarity(residual, inverse);
    rank = Dimension;
    return std::isfinite(value);
  }

  static TMatrixD inputCovariance(const Fwd& mft, const Fwd& mch)
  {
    TMatrixD input(10, 10);
    input.Zero();
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j) {
        // MFT momentum is replaced, so its original uncertainty/correlations
        // are deliberately absent. Its measured angular/position marginal stays.
        if (i < 4 && j < 4)
          input(i, j) = mft.getCovariances()(i, j);
        input(i + 5, j + 5) = mch.getCovariances()(i, j);
      }
    return input;
  }
  static bool compute(Conversion& conversion, const Fwd& mft, const Fwd& mch,
                      double zEnd, double step, bool mchCharge, double derivativeScale, Result& result)
  {
    result.status = 1;
    if (!evaluate(conversion, mft, mch, zEnd, step, mchCharge, 0., result.state))
      return false;
    result.hasNominalState = true;
    TMatrixD seedJac(5, 10);
    seedJac.Zero();
    result.jacobian.Zero();
    auto input = inputCovariance(mft, mch);
    const double floors[5] = {1.e-3, 1.e-3, 1.e-5, 1.e-4, 1.e-6};
    for (int column = 0; column < 10; ++column) {
      if (column == 4)
        continue; // Never use MFT's own momentum magnitude.
      int index = column % 5;
      const auto& original = column < 5 ? mft : mch;
      double h = derivativeScale * std::max({floors[index], 1.e-4 * std::abs(original.getParameters()[index]),
                                             .01 * std::sqrt(std::max(0., input(column, column)))});
      if (index == 4)
        h = std::min(h, .1 * std::abs(original.getInvQPt()));
      State plus, minus;
      auto p = perturb(original, index, h), m = perturb(original, index, -h);
      bool ok = column < 5 ? evaluate(conversion, p, mch, zEnd, step, mchCharge, 0., plus) && evaluate(conversion, m, mch, zEnd, step, mchCharge, 0., minus)
                           : evaluate(conversion, mft, p, zEnd, step, mchCharge, 0., plus) && evaluate(conversion, mft, m, zEnd, step, mchCharge, 0., minus);
      if (!ok) {
        result.status = 2;
        return false;
      }
      for (int row = 0; row < 5; ++row) {
        result.jacobian(row, column) = difference(plus.mft.getParameters()[row], minus.mft.getParameters()[row], row) / (2 * h);
        result.jacobian(row + 5, column) = difference(plus.mch.getParameters()[row], minus.mch.getParameters()[row], row) / (2 * h);
        seedJac(row, column) = difference(plus.seed.getParameters()[row], minus.seed.getParameters()[row], row) / (2 * h);
      }
    }
    TMatrixD temporary(result.jacobian, TMatrixD::kMult, input);
    result.joint = TMatrixD(temporary, TMatrixD::kMultTranspose, result.jacobian);
    TMatrixD seedTmp(seedJac, TMatrixD::kMult, input);
    result.seedCov = TMatrixD(seedTmp, TMatrixD::kMultTranspose, seedJac);
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j)
        result.seedMchCross(i, j) = seedTmp(i, j + 5);
    // One absorber-loss nuisance, shared with opposite correction directions.
    // Perfect upstream/downstream correlation is the explicit same-particle
    // model; widths can differ because the two estimated paths differ.
    State plus, minus;
    constexpr double pullStep = 1.e-3;
    if (!evaluate(conversion, mft, mch, zEnd, step, mchCharge, pullStep, plus) ||
        !evaluate(conversion, mft, mch, zEnd, step, mchCharge, -pullStep, minus)) {
      result.status = 3;
      return false;
    }
    double lossDerivative[5], seedDerivative[5];
    for (int i = 0; i < 5; ++i) {
      lossDerivative[i] = difference(plus.mft.getParameters()[i], minus.mft.getParameters()[i], i) / (2 * pullStep);
      seedDerivative[i] = difference(plus.seed.getParameters()[i], minus.seed.getParameters()[i], i) / (2 * pullStep);
    }
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j) {
        result.joint(i, j) += result.state.mft.getCovariances()(i, j) + lossDerivative[i] * lossDerivative[j];
        result.seedCov(i, j) += seedDerivative[i] * seedDerivative[j];
        result.residual(i, j) = result.joint(i, j) + result.joint(i + 5, j + 5) - result.joint(i, j + 5) - result.joint(i + 5, j);
      }
    std::vector<double> delta(5);
    for (int i = 0; i < 5; ++i)
      delta[i] = difference(result.state.mft.getParameters()[i], result.state.mch.getParameters()[i], i);
    if (!chiDirect<5>(result.residual, delta, result.chiAll, result.rankAll) ||
        !chiDirect<4>(result.residual, delta, result.chiFour, result.rankFour) ||
        !chiDirect<2>(result.residual, delta, result.chiXY, result.rankXY)) {
      result.status = 4;
      return false;
    }
    o2::track::SMatrix55Sym mftCov, mchCov, seedCov;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j <= i; ++j) {
        mftCov(i, j) = .5 * (result.joint(i, j) + result.joint(j, i));
        mchCov(i, j) = .5 * (result.joint(i + 5, j + 5) + result.joint(j + 5, i + 5));
        seedCov(i, j) = .5 * (result.seedCov(i, j) + result.seedCov(j, i));
      }
    result.state.mft.setCovariances(mftCov);
    result.state.mch.setCovariances(mchCov);
    result.state.seed.setCovariances(seedCov);
    result.status = 0;
    return true;
  }
  static std::vector<double> flatten(const TMatrixD& matrix)
  {
    std::vector<double> out;
    out.reserve(matrix.GetNrows() * matrix.GetNcols());
    for (int i = 0; i < matrix.GetNrows(); ++i)
      for (int j = 0; j < matrix.GetNcols(); ++j)
        out.push_back(matrix(i, j));
    return out;
  }
};

// Mixed-event alignment at the absorber exit. The independent inputs are the
// correlated (donor MFT, donor MCH) endpoint state and the target MCH endpoint.
struct MftMchMixedAlignment505 {
  using Joint = MftMchJointTransport;
  using Fwd = Joint::Fwd;
  using Conversion = Joint::Conversion;
  using Rotation = o2::analysis::mftmch::alignment::Rotation;
  // No momentum rescaling or charge replacement occurs in this rigid transform.
  static bool evaluate(Conversion& conversion, const Fwd& mft, const Fwd& donor, const Fwd& target,
                       Fwd& out, Rotation* savedRotation = nullptr)
  {
    namespace align = o2::analysis::mftmch::alignment;
    Rotation r;
    if (!align::directionRotation({donor.getPx(), donor.getPy(), donor.getPz()},
                                  {target.getPx(), target.getPy(), target.getPz()}, r))
      return false;
    auto x = align::transformPosition(r, {mft.getX(), mft.getY(), mft.getZ()},
                                      {donor.getX(), donor.getY(), donor.getZ()},
                                      {target.getX(), target.getY(), target.getZ()});
    auto p = align::rotate(r, {mft.getPx(), mft.getPy(), mft.getPz()});
    double pt = std::hypot(p[0], p[1]);
    if (!std::isfinite(pt) || pt <= 0. || !std::isfinite(p[2]) || p[2] >= -1.e-10)
      return false;
    out = mft;
    out.setX(x[0]);
    out.setY(x[1]);
    out.setZ(x[2]);
    out.setPhi(std::atan2(p[1], p[0]));
    out.setTanl(p[2] / pt);
    out.setInvQPt(std::copysign(1. / pt, mft.getInvQPt()));
    o2::track::SMatrix55Sym zero;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j <= i; ++j)
        zero(i, j) = 0.;
    out.setCovariances(zero);
    // A 3D rotation moves off the fixed z plane. Re-intersect using the official
    // mapped-field extrapolator, without reapplying absorber energy loss or MCS.
    if (std::abs(out.getZ() - target.getZ()) > 1.e-9) {
      auto state = conversion.fwdToMch(out);
      if (!o2::mch::TrackExtrap::extrapToZ(state, target.getZ()))
        return false;
      out = conversion.mchToFwd(state);
    }
    for (int i = 0; i < 5; ++i)
      if (!std::isfinite(out.getParameters()[i]))
        return false;
    if (savedRotation)
      *savedRotation = r;
    return true;
  }
  // Status: 20 invalid nominal transform, 21 donor covariance unavailable,
  // 22 numerical Jacobian failed, 23 invalid chi2. No PSD gate is introduced.
  static bool compute(Conversion& conversion, const Joint::Result& donor, const Fwd& target,
                      bool donorCovValid, double scale, Joint::Result& result, Rotation& rotation)
  {
    result.status = 20;
    result.state = donor.state;
    result.state.mch = target;
    if (!evaluate(conversion, donor.state.mft, donor.state.mch, target, result.state.mft, &rotation))
      return false;
    result.hasNominalState = true;
    result.status = 21;
    if (!donorCovValid)
      return false;
    TMatrixD input(15, 15);
    input.Zero();
    for (int i = 0; i < 10; ++i)
      for (int j = 0; j < 10; ++j)
        input(i, j) = donor.joint(i, j);
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j)
        input(i + 10, j + 10) = target.getCovariances()(i, j);
    TMatrixD jac(10, 15);
    jac.Zero();
    const double floors[5] = {1.e-3, 1.e-3, 1.e-5, 1.e-4, 1.e-6};
    for (int col = 0; col < 15; ++col) {
      int ix = col % 5;
      const Fwd& state = col < 5 ? donor.state.mft : (col < 10 ? donor.state.mch : target);
      double h = scale * std::max({floors[ix], 1.e-4 * std::abs(state.getParameters()[ix]), .01 * std::sqrt(std::max(0., input(col, col)))});
      if (ix == 4)
        h = std::min(h, .1 * std::abs(state.getInvQPt()));
      if (!std::isfinite(h) || h <= 0.) {
        result.status = 22;
        return false;
      }
      Fwd p = Joint::perturb(state, ix, h), m = Joint::perturb(state, ix, -h), plus, minus;
      bool ok = col < 5    ? evaluate(conversion, p, donor.state.mch, target, plus) && evaluate(conversion, m, donor.state.mch, target, minus)
                : col < 10 ? evaluate(conversion, donor.state.mft, p, target, plus) && evaluate(conversion, donor.state.mft, m, target, minus)
                           : evaluate(conversion, donor.state.mft, donor.state.mch, p, plus) && evaluate(conversion, donor.state.mft, donor.state.mch, m, minus);
      if (!ok) {
        result.status = 22;
        return false;
      }
      for (int row = 0; row < 5; ++row)
        jac(row, col) = Joint::difference(plus.getParameters()[row], minus.getParameters()[row], row) / (2 * h);
      if (col >= 10)
        jac(col - 5, col) = 1.;
    }
    TMatrixD tmp(jac, TMatrixD::kMult, input);
    result.joint = TMatrixD(tmp, TMatrixD::kMultTranspose, jac);
    std::vector<double> delta(5);
    o2::track::SMatrix55Sym mc, fc;
    for (int i = 0; i < 5; ++i) {
      delta[i] = Joint::difference(result.state.mft.getParameters()[i], target.getParameters()[i], i);
      for (int j = 0; j < 5; ++j)
        result.residual(i, j) = result.joint(i, j) + result.joint(i + 5, j + 5) - result.joint(i, j + 5) - result.joint(i + 5, j);
      for (int j = 0; j <= i; ++j) {
        fc(i, j) = .5 * (result.joint(i, j) + result.joint(j, i));
        mc(i, j) = .5 * (result.joint(i + 5, j + 5) + result.joint(j + 5, i + 5));
      }
    }
    result.state.mft.setCovariances(fc);
    result.state.mch.setCovariances(mc);
    if (!Joint::chiDirect<5>(result.residual, delta, result.chiAll, result.rankAll) ||
        !Joint::chiDirect<4>(result.residual, delta, result.chiFour, result.rankFour) ||
        !Joint::chiDirect<2>(result.residual, delta, result.chiXY, result.rankXY)) {
      result.status = 23;
      return false;
    }
    result.status = 0;
    return true;
  }
};

template <typename Row, typename Output>
inline void copyMixedResidual(const Row& row, Output& output, int groups = -1)
{
  output(row.mchTrackId(), row.mftTrackId(), row.mchCollisionId(), row.mftCollisionId(), row.matchingPlaneZ(), row.bzAtMftCenter(), row.mchX(), row.mchY(), row.mchPhi(), row.mchTanl(), row.mchInvQPt(), row.mftX(), row.mftY(), row.mftPhi(), row.mftTanl(), row.mftInvQPt(), row.deltaX(), row.deltaY(), row.deltaPhi(), row.deltaPhiWrapped(), row.deltaTanl(), row.deltaInvQPt(), row.deltaEta(), row.deltaPt(), row.chi2All(), row.chi2XYPhiTanl(), row.chi2XY(), row.isFinite(), row.absorberStepCm(), row.transportStatus(), row.donorMchTrackId(), row.donorCollisionId(), row.mixedGroupIndex(), groups >= 0 ? groups : row.numMixedGroups(), row.mixDeltaBc(), row.mixDeltaPhi(), row.mixDeltaR(), row.mixDeltaAttemptsRel(), row.mixDeltaZ(), row.targetMatchAttempts(), row.donorMatchAttempts(), row.targetMftMultiplicity(), row.donorMftMultiplicity(), row.targetMchMultiplicity(), row.donorMchMultiplicity(), row.targetCentrality(), row.donorCentrality(), row.targetCbMultiplicity(), row.donorCbMultiplicity(), (groups > 0 ? 1. / groups : row.mixingWeight()), row.targetFrameId(), row.donorFrameId(), row.mixingStreamId());
}

} // namespace o2::analysis::mftmch
// -----------------------------------------------------------------------------
// MftMchInference
// -----------------------------------------------------------------------------
#include "Common/Core/fwdtrackUtilities.h"
#include "Tools/ML/model.h"

#include <set>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>

#include <THnSparse.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <deque>
#include <fstream>
#include <numeric>
#include <optional>
namespace o2::analysis::mftmch
{
struct MatchingInferenceConfig : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> cfgTrainingMode{"cfgTrainingMode", true, "Training-sample tables; does NOT train a model"};
  o2::framework::Configurable<bool> cfgEvaluationMode{"cfgEvaluationMode", false, "Run frozen model, matching, single-muon and same/mixed dimuon sparse histograms"};
  o2::framework::Configurable<bool> cfgSaveSelectedMuonTable{"cfgSaveSelectedMuonTable", false, "Optional selected single-muon table; no candidate table in inference-only mode"};
  o2::framework::Configurable<bool> cfgEnableMlMatching{"cfgEnableMlMatching", true, "Enable ML matching in inference mode"};
  o2::framework::Configurable<bool> cfgEnableChi2Matching{"cfgEnableChi2Matching", true, "Enable chi2 matching in inference mode"};
  o2::framework::Configurable<bool> cfgEnableMchPV{"cfgEnableMchPV", true, "Enable MCH-only PV reference"};
  o2::framework::Configurable<std::string> cfgMlModel{"cfgMlModel", "NN", "NN or XGBoost (frozen Same/Mixed model by default)"};
  o2::framework::Configurable<std::string> cfgMlModelDirectory{"cfgMlModelDirectory", "/Users/syano/mft-mch-study/models/same-mixed-all291-v1", "Directory containing NN/XGBoost .onnx and .json contracts"};
  o2::framework::Configurable<std::string> cfgMlModelFile{"cfgMlModelFile", "", "Optional ONNX override; metadata must describe the same model"};
  o2::framework::Configurable<std::string> cfgMlMetadataFile{"cfgMlMetadataFile", "", "Optional feature/preprocessing metadata JSON override"};
  o2::framework::Configurable<std::vector<std::string>> cfgMlInputFeatures{"cfgMlInputFeatures", std::vector<std::string>{}, "Ordered feature names; empty uses model contract; mismatches are fatal"};
  o2::framework::Configurable<double> cfgMlScoreMin{"cfgMlScoreMin", 0., "Require sigmoid(good logit)>=cut AFTER fixed one-to-one assignment; source score is not a truth probability"};
  o2::framework::Configurable<std::vector<double>> cfgMlScorePtEdges{"cfgMlScorePtEdges", std::vector<double>{}, "Optional MCH PV pT bin edges for ML score cuts; empty uses cfgMlScoreMin"};
  o2::framework::Configurable<std::vector<double>> cfgMlScorePtCuts{"cfgMlScorePtCuts", std::vector<double>{}, "ML score minima for pT bins; size=edges.size()-1; outside edges uses cfgMlScoreMin"};
  o2::framework::Configurable<std::vector<std::string>> cfgMlWorkingPointNames{"cfgMlWorkingPointNames", std::vector<std::string>{}, "Optional named ML dimuon working points from one common uncut assignment"};
  o2::framework::Configurable<std::vector<double>> cfgMlWorkingPointPtEdges{"cfgMlWorkingPointPtEdges", std::vector<double>{}, "Common MCH PV pT bin edges for named working points"};
  o2::framework::Configurable<std::vector<double>> cfgMlWorkingPointCuts{"cfgMlWorkingPointCuts", std::vector<double>{}, "Row-major score minima: working-point index, then pT bin"};
  o2::framework::Configurable<int> cfgMlBatchSize{"cfgMlBatchSize", 4096, "Maximum candidate batch for CPU ONNX inference"};
  o2::framework::Configurable<int> cfgOnnxThreads{"cfgOnnxThreads", 1, "ONNX CPU threads"};
  o2::framework::Configurable<std::string> cfgChi2Method{"cfgChi2Method", "77_all", "77_all,77_four,77_xy,505_all,505_four,505_xy"};
  o2::framework::Configurable<double> cfgChi2Max{"cfgChi2Max", 40., "Maximum chosen chi2 AFTER one-to-one assignment"};
  o2::framework::Configurable<bool> cfgRequireMid{"cfgRequireMid", true, "Require MCH-MID standalone track (trackType=3) for every reference/matched method"};
  o2::framework::Configurable<bool> cfgEnableDimuonMixing{"cfgEnableDimuonMixing", true, "Mix already-reconstructed single muons, separately from MFT-MCH training mixing"};
  o2::framework::Configurable<int> cfgDimuonMixingDepth{"cfgDimuonMixingDepth", 5, "Maximum prior donor events per event; 0 disables"};
  o2::framework::Configurable<int> cfgDimuonPoolMaxEvents{"cfgDimuonPoolMaxEvents", 2000, "Maximum retained prior events across DFs"};
  o2::framework::Configurable<int> cfgDimuonPoolMaxMuons{"cfgDimuonPoolMaxMuons", 200000, "Maximum retained single-muon states across all methods"};
  o2::framework::Configurable<bool> cfgDimuonMixAcrossDF{"cfgDimuonMixAcrossDF", true, "Retain owned single-muon states across DFs; pool resets on run/field changes"};
  o2::framework::Configurable<double> cfgDimuonMixDeltaZ{"cfgDimuonMixDeltaZ", 1., "Maximum absolute PV z difference in cm; negative disables"};
  o2::framework::Configurable<double> cfgDimuonMixMftMultiplicityRel{"cfgDimuonMixMftMultiplicityRel", -1., "Relative MFT multiplicity difference / target; negative disables"};
  o2::framework::Configurable<double> cfgDimuonMixPvMultiplicityRel{"cfgDimuonMixPvMultiplicityRel", -1., "PVMults relative difference / target; requires centrality process; negative disables"};
  o2::framework::Configurable<double> cfgDimuonMixDeltaCentrality{"cfgDimuonMixDeltaCentrality", -1., "FT0C centrality difference; requires centrality process; negative disables"};
  o2::framework::Configurable<int64_t> cfgDimuonMixMinDeltaBC{"cfgDimuonMixMinDeltaBC", 0, "Minimum absolute BC separation between distinct mixed events"};
  o2::framework::Configurable<bool> cfgDimuonMixNormalize{"cfgDimuonMixNormalize", true, "Mixed-pair weight=1/number of compatible donor events (empty selected-method events included)"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonMassAxis{"cfgDimuonMassAxis", std::vector<double>{600, 0., 12.}, "Sparse mass: nbins,min,max GeV/c^2"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPtAxis{"cfgDimuonPtAxis", std::vector<double>{200, 0., 20.}, "Sparse dimuon pT: nbins,min,max GeV/c"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonRapidityAxis{"cfgDimuonRapidityAxis", std::vector<double>{80, -5., -1.}, "Sparse rapidity: nbins,min,max"};
  o2::framework::Configurable<std::vector<double>> cfgSinglePtAxis{"cfgSinglePtAxis", std::vector<double>{200, 0., 20.}, "Single pT axis"};
  o2::framework::Configurable<std::vector<double>> cfgSingleEtaAxis{"cfgSingleEtaAxis", std::vector<double>{100, -5., -1.}, "Single eta axis"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaAxis{"cfgSingleDcaAxis", std::vector<double>{400, -.2, .2}, "Single signed DCA x/y/z in cm"};
  o2::framework::Configurable<std::vector<double>> cfgMlScoreAxis{"cfgMlScoreAxis", std::vector<double>{100, 0., 1.}, "ML selection score axis"};
  o2::framework::Configurable<std::vector<double>> cfgChi2ScoreAxis{"cfgChi2ScoreAxis", std::vector<double>{200, 0., 100.}, "Chi2 score axis (overflow retained)"};
};
struct MatchingInference {
  using Fwd = o2::track::TrackParCovFwd;
  struct Muon {
    int64_t mch = -1, mft = -1;
    int charge = 0;
    double score = NAN, raw = NAN, mchPtPV = NAN;
    TLorentzVector p;
    std::array<double, 3> dca{NAN, NAN, NAN};
    double pdca = NAN, mftChi = NAN, mftClusters = NAN, mchClusters = NAN, midChi = NAN;
  };
  struct Event {
    uint64_t frame = 0, bc = 0;
    int64_t collision = -1;
    int run = -1, mftMult = 0, pvMult = -1;
    double z = 0., cent = NAN;
    std::array<std::vector<Muon>, 3> muons;
  };
  std::unique_ptr<o2::ml::OnnxModel> model;
  std::vector<std::string> featureNames;
  std::vector<float> mean, stddev;
  double logitSign = -1., requiredAbsorberStep = 0.;
  bool requiredMchCharge = true;
  std::string contract, contractContent;
  std::array<std::shared_ptr<THnSparse>, 3> singles, quality;
  std::array<std::array<std::shared_ptr<THnSparse>, 6>, 3> dimuons;
  std::vector<std::array<std::shared_ptr<THnSparse>, 6>> workingDimuons;
  std::vector<std::shared_ptr<THnSparse>> workingScores;
  std::shared_ptr<TH1> status, poolQA;
  std::deque<Event> pool;
  size_t poolMuons = 0;
  int run = -1;
  double field = NAN;
  static std::vector<std::string> supportedFeatures() { return {"MFT77_X", "MFT77_Y", "MFT77_phi", "MFT77_tanLambda", "MFT77_qOverPt", "MCH_PV_pt", "MCH_PV_eta", "MCH_PV_p", "MCH77_X", "MCH77_Y", "MCH77_phi", "MCH77_tanLambda", "MCH77_qOverPt", "MCH505_X", "MCH505_Y", "MCH505_phi", "MCH505_tanLambda", "MCH505_qOverPt", "delta77_X", "delta77_Y", "delta77_phiWrapped", "delta77_tanLambda", "delta77_qOverPt", "delta77_eta", "delta77_pt", "delta505_X", "delta505_Y", "delta505_phiWrapped", "delta505_tanLambda", "delta505_qOverPt", "delta505_eta", "delta505_pt", "MFT77_charge", "MCH77_charge", "MFT505_X", "MFT505_Y", "R77", "R505"}; }
  std::vector<int> featureIndices;
  static o2::framework::AxisSpec axis(const std::vector<double>& v, const char* title)
  {
    if (v.size() != 3 || !std::isfinite(v[0]) || v[0] < 1 || v[0] != std::floor(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]) || v[2] <= v[1])
      LOGF(fatal, "Invalid sparse axis %s", title);
    return {static_cast<int>(v[0]), v[1], v[2], title};
  }
  void init(MatchingInferenceConfig& c, o2::framework::HistogramRegistry& h)
  {
    if (!c.cfgEvaluationMode)
      return;
    if (c.cfgMlScoreMin < 0. || c.cfgMlScoreMin > 1. || !std::isfinite(c.cfgMlScoreMin.value) || c.cfgChi2Max < 0. || !std::isfinite(c.cfgChi2Max.value) || c.cfgMlBatchSize < 1 || c.cfgDimuonMixingDepth < 0 || c.cfgDimuonPoolMaxEvents < 1 || c.cfgDimuonPoolMaxMuons < 1)
      LOGF(fatal, "Invalid inference cuts/pool limits");
    for (double v : {c.cfgDimuonMixDeltaZ.value, c.cfgDimuonMixMftMultiplicityRel.value, c.cfgDimuonMixPvMultiplicityRel.value, c.cfgDimuonMixDeltaCentrality.value})
      if (!std::isfinite(v))
        LOGF(fatal, "Nonfinite dimuon pool setting");
    auto const& scoreEdges = c.cfgMlScorePtEdges.value;
    auto const& scoreCuts = c.cfgMlScorePtCuts.value;
    if (!scoreEdges.empty() || !scoreCuts.empty()) {
      if (scoreEdges.size() < 2 || scoreCuts.size() + 1 != scoreEdges.size())
        LOGF(fatal, "ML pT score cuts require N+1 edges for N cuts");
      for (size_t i = 0; i < scoreEdges.size(); ++i)
        if (!std::isfinite(scoreEdges[i]) || (i && scoreEdges[i] <= scoreEdges[i - 1]))
          LOGF(fatal, "ML score pT edges must be finite and strictly increasing");
      for (double v : scoreCuts)
        if (!std::isfinite(v) || v < 0. || v > 1.)
          LOGF(fatal, "ML pT score minima must lie in [0,1]");
    }
    auto const& wpNames = c.cfgMlWorkingPointNames.value;
    auto const& wpEdges = c.cfgMlWorkingPointPtEdges.value;
    auto const& wpCuts = c.cfgMlWorkingPointCuts.value;
    if (!wpNames.empty()) {
      if (!c.cfgEnableMlMatching || c.cfgMlScoreMin.value != 0. || !c.cfgMlScorePtEdges.value.empty())
        LOGF(fatal, "Named working points require enabled ML with no primary score cut");
      if (wpEdges.size() < 2 || wpCuts.size() != wpNames.size() * (wpEdges.size() - 1))
        LOGF(fatal, "Working point cut matrix has wrong dimensions");
      std::set<std::string> uniqueNames;
      for (auto const& name : wpNames) {
        if (name.empty() || name == "ML" || name == "MCHPV" || name == "Chi2" || !uniqueNames.insert(name).second ||
            name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
          LOGF(fatal, "Invalid or duplicate working point name");
      }
      for (size_t i = 0; i < wpEdges.size(); ++i)
        if (!std::isfinite(wpEdges[i]) || (i && wpEdges[i] <= wpEdges[i - 1]))
          LOGF(fatal, "Invalid working point pT edges");
      for (double cut : wpCuts)
        if (!std::isfinite(cut) || cut < 0. || cut > 1.)
          LOGF(fatal, "Invalid working point score cut");
    } else if (!wpEdges.empty() || !wpCuts.empty())
      LOGF(fatal, "Working point cuts require names");
    if (c.cfgDimuonMixMinDeltaBC < 0 || c.cfgOnnxThreads < 1)
      LOGF(fatal, "Invalid BC gap/ONNX thread count");
    const std::set<std::string> methods{"77_all", "77_four", "77_xy", "505_all", "505_four", "505_xy"};
    if (!methods.count(c.cfgChi2Method.value))
      LOGF(fatal, "Unknown chi2 method");
    if (c.cfgEnableMlMatching) {
      if (c.cfgMlModel.value != "NN" && c.cfgMlModel.value != "XGBoost")
        LOGF(fatal, "Choose NN or XGBoost");
      const std::string metadata = c.cfgMlMetadataFile.value.empty() ? c.cfgMlModelDirectory.value + "/" + c.cfgMlModel.value + ".json" : c.cfgMlMetadataFile.value;
      boost::property_tree::ptree t;
      boost::property_tree::read_json(metadata, t);
      if (t.get<int>("schema") != 1 || t.get<std::string>("transform") != "signed_log1p_finite_then_standardize" || t.get<std::string>("model") != c.cfgMlModel.value)
        LOGF(fatal, "Unsupported/wrong model contract");
      contract = metadata;
      requiredAbsorberStep = t.get<double>("absorber_step_cm");
      requiredMchCharge = t.get<bool>("mch_charge_seed");
      for (auto& v : t.get_child("features"))
        featureNames.push_back(v.second.get_value<std::string>());
      if (!c.cfgMlInputFeatures.value.empty() && c.cfgMlInputFeatures.value != featureNames)
        LOGF(fatal, "JSON feature order differs from trained model contract");
      const auto available = supportedFeatures();
      for (auto& n : featureNames) {
        auto it = std::find(available.begin(), available.end(), n);
        if (it == available.end())
          LOGF(fatal, "Unsupported feature %s", n.c_str());
        featureIndices.push_back(it - available.begin());
      }
      for (auto& v : t.get_child("mean"))
        mean.push_back(v.second.get_value<float>());
      for (auto& v : t.get_child("std"))
        stddev.push_back(v.second.get_value<float>());
      if (mean.size() != 2 * featureNames.size() || stddev.size() != mean.size())
        LOGF(fatal, "Wrong feature preprocessing dimensions");
      for (size_t i = 0; i < mean.size(); ++i)
        if (!std::isfinite(mean[i]) || !std::isfinite(stddev[i]) || stddev[i] <= 0.)
          LOGF(fatal, "Invalid model normalization");
      logitSign = t.get<double>("good_logit_sign");
      if (logitSign != 1. && logitSign != -1.)
        LOGF(fatal, "Invalid score direction");
      model = std::make_unique<o2::ml::OnnxModel>();
      model->initModel(c.cfgMlModelFile.value.empty() ? c.cfgMlModelDirectory.value + "/" + c.cfgMlModel.value + ".onnx" : c.cfgMlModelFile.value, false, c.cfgOnnxThreads);
      const auto shapes = model->getInputShapes();
      if (shapes.size() != 1 || shapes[0].size() != 2 || shapes[0][1] != static_cast<int64_t>(2 * featureNames.size()))
        LOGF(fatal, "ONNX input shape disagrees with model contract");
      std::ostringstream contractStream;
      boost::property_tree::write_json(contractStream, t);
      contractContent = contractStream.str();
      LOGF(info, "Inference model %s contract=%s; %zu features; sigmoid(good logit), NOT correct probability", c.cfgMlModel.value.c_str(), metadata.c_str(), featureNames.size());
    }
    using namespace o2::framework;
    const auto mass = axis(c.cfgDimuonMassAxis, "m_{#mu#mu} (GeV/c^{2})"), pt = axis(c.cfgDimuonPtAxis, "p_{T,#mu#mu} (GeV/c)"), rap = axis(c.cfgDimuonRapidityAxis, "y_{#mu#mu}");
    const std::array<std::string, 3> labels{"MCHPV", "Chi2", "ML"};
    const std::array<std::string, 6> charges{"ULS", "LSpp", "LSmm", "MixULS", "MixLSpp", "MixLSmm"};
    for (int m = 0; m < 3; ++m) {
      if (!enabled(c, m))
        continue;
      for (int q = 0; q < 6; ++q)
        dimuons[m][q] = h.add<THnSparse>(labels[m] + "/" + charges[q], "PV dimuons", HistType::kTHnSparseD, {mass, pt, rap}, true);
      std::vector<AxisSpec> axes{axis(c.cfgSinglePtAxis, "p_{T} (GeV/c)"), axis(c.cfgSingleEtaAxis, "#eta"), {128, 0., 2 * TMath::Pi(), "#varphi"}, {200, 0., 200., "p (GeV/c)"}, axis(c.cfgSingleDcaAxis, "DCA_{x} (cm)"), axis(c.cfgSingleDcaAxis, "DCA_{y} (cm)"), axis(c.cfgSingleDcaAxis, "DCA_{z} (cm)"), {3, -1.5, 1.5, "q"}};
      if (m)
        axes.push_back(m == 1 ? axis(c.cfgChi2ScoreAxis, "chosen #chi^{2}") : axis(c.cfgMlScoreAxis, "ML good-source score"));
      singles[m] = h.add<THnSparse>(labels[m] + "/SingleMuon", "Selected single muons; MCHPV DCA uses no-Branson MCH trajectory", HistType::kTHnSparseD, axes, true);
      std::vector<AxisSpec> qa{axis(c.cfgSinglePtAxis, "p_{T}"), {200, 0., 10., "DCA_{3D} (cm)"}, {200, 0., 2000., "MCH pDCA"}, {200, 0., 100., "MCH-MID #chi^{2}"}, {20, 0., 20., "MCH clusters"}};
      if (m) {
        qa.push_back({200, 0., 100., "MFT fit #chi^{2}"});
        qa.push_back({11, -.5, 10.5, "MFT clusters"});
        qa.push_back(m == 1 ? axis(c.cfgChi2ScoreAxis, "chosen #chi^{2}") : axis(c.cfgMlScoreAxis, "ML good-source score"));
      }
      quality[m] = h.add<THnSparse>(labels[m] + "/SingleQuality", "Selected track quality", HistType::kTHnSparseD, qa, true);
    }
    workingDimuons.resize(c.cfgMlWorkingPointNames.value.size());
    workingScores.resize(c.cfgMlWorkingPointNames.value.size());
    for (size_t w = 0; w < workingDimuons.size(); ++w) {
      for (int q = 0; q < 6; ++q)
        workingDimuons[w][q] = h.add<THnSparse>(c.cfgMlWorkingPointNames.value[w] + "/" + charges[q], "Named ML working-point PV dimuons", HistType::kTHnSparseD, {mass, pt, rap}, true);
      workingScores[w] = h.add<THnSparse>(c.cfgMlWorkingPointNames.value[w] + "/MchPtScore", "Selected MCH PV pT and NN score", HistType::kTHnSparseD,
                                          {axis(c.cfgSinglePtAxis, "MCH p_{T} at PV (GeV/c)"), axis(c.cfgMlScoreAxis, "ML selection score")}, true);
    }
    status = h.add<TH1>("Status", "0=MCH PV failure,1=matched DCA failure,2=ML input pair,3=ML assigned,4=ML accepted,5=chi2 assigned,6=chi2 accepted,7=event,8=pool reset,9=pool eviction", HistType::kTH1D, {{10, -.5, 9.5}});
    poolQA = h.add<TH1>("DimuonMixDonors", "Compatible donor events;N donors;Target events", HistType::kTH1D, {{201, -.5, 200.5}});
  }
  static bool enabled(MatchingInferenceConfig& c, int m) { return m == 0 ? c.cfgEnableMchPV.value : m == 1 ? c.cfgEnableChi2Matching.value
                                                                                                           : c.cfgEnableMlMatching.value; }
  static double sigmoid(double x) { return x >= 0. ? 1. / (1. + std::exp(-x)) : std::exp(x) / (1. + std::exp(x)); }
  // Rectangular Hungarian maximum-weight matching; dummy columns allow unpaired MCH.
  static std::vector<int> assign(const std::vector<std::pair<int64_t, int64_t>>& ids, const std::vector<double>& w)
  {
    std::map<int64_t, int> cm, fm;
    for (auto p : ids) {
      cm.emplace(p.first, 0);
      fm.emplace(p.second, 0);
    }
    int n = 0, k = 0;
    for (auto& p : cm)
      p.second = n++;
    for (auto& p : fm)
      p.second = k++;
    int cols = k + n;
    std::vector<std::vector<double>> cost(n, std::vector<double>(cols, 1e6));
    std::vector<std::vector<int>> lookup(n, std::vector<int>(k, -1));
    for (int i = 0; i < n; ++i)
      for (int j = k; j < cols; ++j)
        cost[i][j] = 0.;
    for (size_t t = 0; t < ids.size(); ++t) {
      int i = cm.at(ids[t].first), j = fm.at(ids[t].second);
      if (std::isfinite(w[t]) && -w[t] < cost[i][j]) {
        cost[i][j] = -w[t];
        lookup[i][j] = t;
      }
    }
    std::vector<double> u(n + 1), v(cols + 1);
    std::vector<int> p(cols + 1), way(cols + 1);
    for (int i = 1; i <= n; ++i) {
      p[0] = i;
      int j0 = 0;
      std::vector<double> minv(cols + 1, 1e20);
      std::vector<bool> used(cols + 1, false);
      do {
        used[j0] = true;
        int i0 = p[j0], j1 = 0;
        double delta = 1e20;
        for (int j = 1; j <= cols; ++j)
          if (!used[j]) {
            double cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
            if (cur < minv[j]) {
              minv[j] = cur;
              way[j] = j0;
            }
            if (minv[j] < delta) {
              delta = minv[j];
              j1 = j;
            }
          }
        for (int j = 0; j <= cols; ++j)
          if (used[j]) {
            u[p[j]] += delta;
            v[j] -= delta;
          } else
            minv[j] -= delta;
        j0 = j1;
      } while (p[j0]);
      do {
        int j1 = way[j0];
        p[j0] = p[j1];
        j0 = j1;
      } while (j0);
    }
    std::vector<int> out;
    for (int j = 1; j <= k; ++j)
      if (p[j] && lookup[p[j] - 1][j - 1] >= 0)
        out.push_back(lookup[p[j] - 1][j - 1]);
    std::sort(out.begin(), out.end());
    return out;
  }
  template <typename Row>
  static std::array<float, 38> features(Row const& a, Row const& b, Fwd const& pv)
  {
    std::array<float, 38> f{};
    f[0] = std::get<11>(a);
    f[1] = std::get<12>(a);
    f[2] = std::get<13>(a);
    f[3] = std::get<14>(a);
    f[4] = std::get<15>(a);
    f[5] = pv.getPt();
    f[6] = pv.getEta();
    f[7] = pv.getP();
    f[8] = std::get<6>(a);
    f[9] = std::get<7>(a);
    f[10] = std::get<8>(a);
    f[11] = std::get<9>(a);
    f[12] = std::get<10>(a);
    f[13] = std::get<6>(b);
    f[14] = std::get<7>(b);
    f[15] = std::get<8>(b);
    f[16] = std::get<9>(b);
    f[17] = std::get<10>(b);
    f[18] = std::get<16>(a);
    f[19] = std::get<17>(a);
    f[20] = std::get<19>(a);
    f[21] = std::get<20>(a);
    f[22] = std::get<21>(a);
    f[23] = std::get<22>(a);
    f[24] = std::get<23>(a);
    f[25] = std::get<16>(b);
    f[26] = std::get<17>(b);
    f[27] = std::get<19>(b);
    f[28] = std::get<20>(b);
    f[29] = std::get<21>(b);
    f[30] = std::get<22>(b);
    f[31] = std::get<23>(b);
    f[32] = std::get<15>(a) > 0. ? 1. : -1.;
    f[33] = std::get<10>(a) > 0. ? 1. : -1.;
    f[34] = std::get<11>(b);
    f[35] = std::get<12>(b);
    f[36] = std::hypot(f[18], f[19]);
    f[37] = std::hypot(f[25], f[26]);
    return f;
  }
  std::vector<float> predict(std::vector<std::array<float, 38>> const& x, MatchingInferenceConfig& c)
  {
    std::vector<float> out;
    out.reserve(x.size());
    const size_t n = featureNames.size();
    for (size_t lo = 0; lo < x.size(); lo += c.cfgMlBatchSize.value) {
      const size_t hi = std::min(x.size(), lo + c.cfgMlBatchSize.value);
      std::vector<float> input((hi - lo) * 2 * n);
      for (size_t i = lo; i < hi; ++i)
        for (size_t j = 0; j < n; ++j) {
          float v = x[i][featureIndices[j]];
          bool finite = std::isfinite(v);
          float z = finite ? std::copysign(std::log1p(std::abs(v)), v) : 0.f;
          input[(i - lo) * 2 * n + j] = (z - mean[j]) / stddev[j];
          input[(i - lo) * 2 * n + n + j] = ((finite ? 1.f : 0.f) - mean[n + j]) / stddev[n + j];
        }
      auto z = model->evalModel<float>(input);
      if (z.size() != hi - lo)
        LOGF(fatal, "ONNX output must be one scalar logit per pair");
      for (auto v : z) {
        if (!std::isfinite(v))
          LOGF(fatal, "Nonfinite ONNX logit");
        out.push_back(logitSign * v);
      }
    }
    return out;
  }
  static void fillPair(std::shared_ptr<THnSparse> const& h, Muon const& a, Muon const& b, double w)
  {
    auto p = a.p + b.p;
    if (!std::isfinite(p.M()) || !std::isfinite(p.Pt()) || !std::isfinite(p.Rapidity()))
      return;
    double x[3]{p.M(), p.Pt(), p.Rapidity()};
    h->Fill(x, w);
  }
  static int signIndex(Muon const& a, Muon const& b) { return a.charge * b.charge < 0 ? 0 : a.charge > 0 ? 1
                                                                                                         : 2; }
  static bool passesWorkingPoint(MatchingInferenceConfig& c, size_t workingPoint, Muon const& u)
  {
    auto const& edges = c.cfgMlWorkingPointPtEdges.value;
    if (!std::isfinite(u.mchPtPV) || u.mchPtPV < edges.front() || u.mchPtPV >= edges.back())
      return false;
    size_t bin = std::upper_bound(edges.begin(), edges.end(), u.mchPtPV) - edges.begin() - 1;
    return u.score >= c.cfgMlWorkingPointCuts.value[workingPoint * (edges.size() - 1) + bin];
  }
  static bool rel(double a, double b, double limit) { return limit < 0. || (a == 0. ? b == 0. : std::abs(a - b) / std::abs(a) <= limit); }
  template <typename Collisions, typename BCs, typename Mchs, typename Mfts, typename Rows, typename Output, typename Conversion>
  void process(MatchingInferenceConfig& c, uint64_t frame, Collisions const& collisions, BCs const& bcs, Mchs const& tracks, Mfts const& mfts, Rows const& upstream, Rows const& downstream, std::unordered_map<int64_t, bool> const& selected, std::unordered_map<int64_t, double> const& cents, std::unordered_map<int64_t, int> const& pvMults, Conversion& conversion, double bz, Output& table)
  {
    if (!c.cfgEvaluationMode || bcs.size() == 0)
      return;
    int currentRun = bcs.begin().runNumber();
    if (run != currentRun || field != bz) {
      pool.clear();
      poolMuons = 0;
      run = currentRun;
      field = bz;
      status->Fill(8);
    }
    if (!c.cfgDimuonMixAcrossDF) {
      pool.clear();
      poolMuons = 0;
    }
    std::map<int64_t, Fwd> pvCache;
    std::map<int64_t, Event> events;
    std::map<int64_t, int> mftMult;
    for (auto& m : mfts)
      if (m.has_collision())
        ++mftMult[m.collisionId()];
    auto baseMuon = [&](auto const& tr) {Muon u;u.mch=tr.globalIndex();u.charge=tr.signed1Pt()>0?1:-1;u.pdca=tr.pDca();u.midChi=tr.chi2MatchMCHMID();u.mchClusters=tr.nClusters();return u; };
    for (auto const& tr : tracks) {
      auto sel = selected.find(tr.globalIndex());
      if (sel == selected.end() || !sel->second || !tr.has_collision() || (c.cfgRequireMid && int(tr.trackType()) != 3))
        continue;
      auto col = collisions.rawIteratorAt(tr.collisionId() - collisions.offset());
      auto mch = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(tr, tr));
      if (!o2::mch::TrackExtrap::extrapToVertex(mch, col.posX(), col.posY(), col.posZ(), col.covXX(), col.covYY())) {
        status->Fill(0);
        continue;
      }
      auto pv = conversion.mchToFwd(mch);
      if (!std::isfinite(pv.getP()) || pv.getP() <= 0.) {
        status->Fill(0);
        continue;
      }
      pvCache.emplace(tr.globalIndex(), pv);
      auto& e = events[tr.collisionId()];
      e.frame = frame;
      e.collision = tr.collisionId();
      e.run = currentRun;
      e.z = col.posZ();
      e.bc = bcs.rawIteratorAt(col.bcId() - bcs.offset()).globalBC();
      e.mftMult = mftMult[tr.collisionId()];
      auto ct = cents.find(tr.collisionId());
      if (ct != cents.end())
        e.cent = ct->second;
      auto mt = pvMults.find(tr.collisionId());
      if (mt != pvMults.end())
        e.pvMult = mt->second;
      if (c.cfgEnableMchPV) {
        auto u = baseMuon(tr);
        u.p.SetXYZM(pv.getPx(), pv.getPy(), pv.getPz(), mMu);
        auto d = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(tr, tr));
        if (o2::mch::TrackExtrap::extrapToVertexWithoutBranson(d, col.posZ())) {
          auto dd = conversion.mchToFwd(d);
          u.dca = {dd.getX() - col.posX(), dd.getY() - col.posY(), dd.getZ() - col.posZ()};
        }
        e.muons[0].push_back(u);
      }
    }
    std::map<int64_t, std::vector<int>> candidates;
    for (size_t i = 0; i < upstream.size(); ++i)
      if (pvCache.count(std::get<0>(upstream[i])) && std::get<2>(upstream[i]) == std::get<3>(upstream[i]))
        candidates[std::get<2>(upstream[i])].push_back(i);
    auto matched = [&](int row, double score, double raw) -> std::optional<Muon> {
      auto const& r = upstream[row];
      auto tr = tracks.rawIteratorAt(std::get<0>(r) - tracks.offset());
      auto mf = mfts.rawIteratorAt(std::get<1>(r) - mfts.offset());
      auto col = collisions.rawIteratorAt(tr.collisionId() - collisions.offset());
      Muon u = baseMuon(tr);
      u.mft = mf.globalIndex();
      u.score = score;
      u.raw = raw;
      u.mchPtPV = pvCache.at(tr.globalIndex()).getPt();
      u.mftChi = mf.chi2();
      uint64_t packed = mf.mftClusterSizesAndTrackFlags();
      u.mftClusters = 0;
      for (int j = 0; j < 10; ++j)
        u.mftClusters += ((packed >> (6 * j)) & 63) != 0;
      o2::track::SMatrix55Sym zero;
      for (int i = 0; i < 5; ++i)
        for (int j = 0; j <= i; ++j)
          zero(i, j) = 0.;
      Fwd native(mf.z(), SMatrix5(mf.x(), mf.y(), mf.phi(), mf.tgl(), mf.signed1Pt()), zero, mf.chi2());
      auto hybrid = o2::aod::fwdtrackutils::refitGlobalMuonCov(pvCache.at(tr.globalIndex()), native);
      std::array<double, 3> dca{NAN, NAN, NAN};
      if (std::abs(bz) > 1e-12) {
        hybrid.propagateToDCAhelix(bz, {col.posX(), col.posY(), col.posZ()}, dca);
      } else {
        double inv = 1. / hybrid.getTgl(), tx = std::cos(hybrid.getPhi()) * inv, ty = std::sin(hybrid.getPhi()) * inv;
        double dz = -((hybrid.getX() - col.posX()) * tx + (hybrid.getY() - col.posY()) * ty + (hybrid.getZ() - col.posZ())) / (1. + tx * tx + ty * ty);
        hybrid.propagateToZlinear(hybrid.getZ() + dz);
        dca = {hybrid.getX() - col.posX(), hybrid.getY() - col.posY(), hybrid.getZ() - col.posZ()};
      }
      if (!std::isfinite(dca[0]) || !std::isfinite(dca[1]) || !std::isfinite(dca[2]) || !std::isfinite(hybrid.getPhi())) {
        status->Fill(1);
        return std::nullopt;
      }
      u.dca = dca;
      const double p = pvCache.at(tr.globalIndex()).getP(), pt = p / std::sqrt(1. + hybrid.getTgl() * hybrid.getTgl());
      u.p.SetXYZM(pt * std::cos(hybrid.getPhi()), pt * std::sin(hybrid.getPhi()), pt * hybrid.getTgl(), mMu);
      return u;
    };
    for (auto const& [col, rows] : candidates) {
      std::vector<std::pair<int64_t, int64_t>> ids;
      for (int i : rows)
        ids.emplace_back(std::get<0>(upstream[i]), std::get<1>(upstream[i]));
      for (int method = 1; method <= 2; ++method) {
        if (!enabled(c, method))
          continue;
        std::vector<double> w(rows.size(), NAN), score(rows.size(), NAN), raw(rows.size(), NAN);
        if (method == 2) {
          std::vector<std::array<float, 38>> x;
          for (int i : rows)
            x.push_back(features(upstream[i], downstream[i], pvCache.at(std::get<0>(upstream[i]))));
          auto z = predict(x, c);
          for (size_t j = 0; j < rows.size(); ++j) {
            raw[j] = z[j];
            score[j] = w[j] = sigmoid(z[j]);
          }
          status->Fill(2, rows.size());
        } else {
          bool down = c.cfgChi2Method.value.rfind("505", 0) == 0;
          std::string name = c.cfgChi2Method.value;
          for (size_t j = 0; j < rows.size(); ++j) {
            auto const& r = down ? downstream[rows[j]] : upstream[rows[j]];
            double chi = name.find("all") != std::string::npos ? std::get<24>(r) : name.find("four") != std::string::npos ? std::get<25>(r)
                                                                                                                          : std::get<26>(r);
            if (std::isfinite(chi) && chi >= 0.) {
              raw[j] = score[j] = chi;
              w[j] = std::exp(-chi / 2.);
            }
          }
        }
        auto selectedRows = assign(ids, w);
        status->Fill(method == 1 ? 5 : 3, selectedRows.size());
        for (int j : selectedRows) {
          double mlCut = c.cfgMlScoreMin.value;
          const double mchPtPV = pvCache.at(std::get<0>(upstream[rows[j]])).getPt();
          auto const& edges = c.cfgMlScorePtEdges.value;
          if (edges.size() >= 2 && mchPtPV >= edges.front() && mchPtPV < edges.back()) {
            const size_t bin = std::upper_bound(edges.begin(), edges.end(), mchPtPV) - edges.begin() - 1;
            mlCut = c.cfgMlScorePtCuts.value[bin];
          }
          if (method == 1 ? score[j] > c.cfgChi2Max.value : score[j] < mlCut)
            continue;
          auto u = matched(rows[j], score[j], raw[j]);
          if (u) {
            events.at(col).muons[method].push_back(*u);
            status->Fill(method == 1 ? 6 : 4);
          }
        }
      }
    }
    for (auto& [col, e] : events) {
      status->Fill(7);
      for (int m = 0; m < 3; ++m) {
        if (!enabled(c, m))
          continue;
        for (auto const& u : e.muons[m]) {
          double x[9]{u.p.Pt(), u.p.Eta(), TVector2::Phi_0_2pi(u.p.Phi()), u.p.P(), u.dca[0], u.dca[1], u.dca[2], double(u.charge), u.score};
          singles[m]->Fill(x);
          double q[8]{u.p.Pt(), std::hypot(std::hypot(u.dca[0], u.dca[1]), u.dca[2]), u.pdca, u.midChi, u.mchClusters, u.mftChi, u.mftClusters, u.score};
          quality[m]->Fill(q);
          if (c.cfgSaveSelectedMuonTable)
            table(frame, col, u.mch, u.mft, static_cast<int8_t>(m), u.score, u.raw, static_cast<int8_t>(u.charge), u.p.Px(), u.p.Py(), u.p.Pz(), u.dca[0], u.dca[1], u.dca[2]);
        }
        for (size_t i = 0; i < e.muons[m].size(); ++i)
          for (size_t j = i + 1; j < e.muons[m].size(); ++j)
            fillPair(dimuons[m][signIndex(e.muons[m][i], e.muons[m][j])], e.muons[m][i], e.muons[m][j], 1.);
      }
      for (size_t w = 0; w < workingDimuons.size(); ++w) {
        auto const& muons = e.muons[2];
        for (size_t i = 0; i < muons.size(); ++i) {
          if (!passesWorkingPoint(c, w, muons[i]))
            continue;
          double point[2]{muons[i].mchPtPV, muons[i].score};
          workingScores[w]->Fill(point);
          for (size_t j = i + 1; j < muons.size(); ++j)
            if (passesWorkingPoint(c, w, muons[j]))
              fillPair(workingDimuons[w][signIndex(muons[i], muons[j])], muons[i], muons[j], 1.);
        }
      }
      if (!c.cfgEnableDimuonMixing || c.cfgDimuonMixingDepth == 0)
        continue;
      std::vector<Event const*> donors;
      for (auto it = pool.rbegin(); it != pool.rend() && donors.size() < static_cast<size_t>(c.cfgDimuonMixingDepth.value); ++it) {
        auto const& d = *it;
        if (e.frame == d.frame && e.collision == d.collision)
          continue;
        auto deltaBC = e.bc > d.bc ? e.bc - d.bc : d.bc - e.bc;
        if (deltaBC < static_cast<uint64_t>(c.cfgDimuonMixMinDeltaBC.value))
          continue;
        if (c.cfgDimuonMixDeltaZ >= 0. && std::abs(e.z - d.z) > c.cfgDimuonMixDeltaZ)
          continue;
        if (!rel(e.mftMult, d.mftMult, c.cfgDimuonMixMftMultiplicityRel))
          continue;
        if (c.cfgDimuonMixPvMultiplicityRel >= 0. && (e.pvMult < 0 || d.pvMult < 0 || !rel(e.pvMult, d.pvMult, c.cfgDimuonMixPvMultiplicityRel)))
          continue;
        if (c.cfgDimuonMixDeltaCentrality >= 0. && (!std::isfinite(e.cent) || !std::isfinite(d.cent) || std::abs(e.cent - d.cent) > c.cfgDimuonMixDeltaCentrality))
          continue;
        donors.push_back(&d);
      }
      poolQA->Fill(donors.size());
      double weight = c.cfgDimuonMixNormalize && donors.size() ? 1. / donors.size() : 1.;
      for (auto* d : donors)
        for (int m = 0; m < 3; ++m)
          if (enabled(c, m))
            for (auto const& u : e.muons[m])
              for (auto const& v : d->muons[m])
                fillPair(dimuons[m][3 + signIndex(u, v)], u, v, weight);
      for (auto* donor : donors)
        for (size_t w = 0; w < workingDimuons.size(); ++w)
          for (auto const& u : e.muons[2])
            if (passesWorkingPoint(c, w, u))
              for (auto const& v : donor->muons[2])
                if (passesWorkingPoint(c, w, v))
                  fillPair(workingDimuons[w][3 + signIndex(u, v)], u, v, weight);
      size_t n = 0;
      for (auto const& v : e.muons)
        n += v.size();
      if (n > static_cast<size_t>(c.cfgDimuonPoolMaxMuons.value))
        continue;
      poolMuons += n;
      pool.push_back(e);
      while (pool.size() > static_cast<size_t>(c.cfgDimuonPoolMaxEvents.value) || poolMuons > static_cast<size_t>(c.cfgDimuonPoolMaxMuons.value)) {
        for (auto const& v : pool.front().muons)
          poolMuons -= v.size();
        pool.pop_front();
        status->Fill(9);
      }
    }
  }
};
} // namespace o2::analysis::mftmch

namespace o2::analysis::mftmch
{
struct match_mft_mch_data_mc {
  MatchingInferenceConfig inferenceConfig;
  MatchingInference inference;
  HistogramRegistry inferenceHistos{"MatchingInference", {}, OutputObjHandlingPolicy::AnalysisObject};
  Produces<aod::SelectedMatchingMuons> selectedMuonTable;
  OutputObj<TNamed> inferenceMetadata{TNamed("InferenceMetadata", "Inference disabled")};

  Produces<o2::aod::RawMixedMatchPlaneResiduals> tableMixedResiduals;
  Configurable<bool> cfgEnableEventMixing{"cfgEnableEventMixing", false, "Save separate mixed-event residuals at both planes"};
  Produces<o2::aod::MftMchMixedAlignments> tableMixedAlignments;
  Configurable<bool> cfgMixingAlign505{"cfgMixingAlign505", false, "Align donor MFT to target MCH position/direction at -505 cm; -77.5 cm remains unaligned"};
  Configurable<int> cfgMixingDepth{"cfgMixingDepth", -1, "Maximum distinct donor collisions per MCH; -1=no limit"};
  Configurable<int64_t> cfgMinDeltaBc{"cfgMinDeltaBc", 3564, "Minimum BC separation of mixed collisions"};
  Configurable<float> cfgMaxDeltaPhi{"cfgMaxDeltaPhi", static_cast<float>(TMath::Pi() / 10.), "Maximum difference of MCH position azimuths in radians; negative disables"};
  Configurable<float> cfgMaxDeltaR{"cfgMaxDeltaR", 10.f, "Maximum donor-target MCH radius difference in cm; negative disables"};
  Configurable<float> cfgMaxDeltaAttempts{"cfgMaxDeltaAttempts", .1f, "Maximum donor-target relative MFT candidate-count difference; negative disables"};
  Configurable<float> cfgMaxDeltaZ{"cfgMaxDeltaZ", 1.f, "Maximum absolute vertex-z difference in cm; negative disables"};
  Configurable<bool> cfgMixingSymmetricCuts{"cfgMixingSymmetricCuts", true, "Use absolute radius/relative-count differences; false preserves official signed cuts"};
  Configurable<float> cfgMixingMaxDeltaCentrality{"cfgMixingMaxDeltaCentrality", -1.f, "Absolute FT0C centrality difference in percentile points; -1 disables"};
  Configurable<double> cfgMixingMaxDeltaMchEtaPV{"cfgMixingMaxDeltaMchEtaPV", -1., "Maximum absolute donor-target MCH momentum-eta difference at PV; negative disables"};
  // Group multiplicity controls to stay within the DPL aggregate-field limit.
  struct : ConfigurableGroup {
    Configurable<double> cfgMixingMaxDeltaMftMultiplicityRel{"cfgMixingMaxDeltaMftMultiplicityRel", .1, "Absolute MFT count difference / target count; -1 disables"};
    Configurable<int> cfgMixingMaxDeltaMchMultiplicity{"cfgMixingMaxDeltaMchMultiplicity", 0, "Absolute MCH count difference; 0 exact, -1 disables"};
    Configurable<double> cfgMixingMaxDeltaCbMultiplicityRel{"cfgMixingMaxDeltaCbMultiplicityRel", -1., "Absolute PVMults multNTracksPV difference / target count; negative disables, zero target requires zero donor"};
    Configurable<bool> cfgMixingCbMultiplicityStrict{"cfgMixingCbMultiplicityStrict", false, "Require PV multiplicity relative difference strictly below the limit"};
    Configurable<bool> cfgMixingVertexZStrict{"cfgMixingVertexZStrict", false, "Require vertex-z difference strictly below the limit"};
    Configurable<bool> cfgMixingMftMultiplicityStrict{"cfgMixingMftMultiplicityStrict", false, "Require MFT multiplicity relative difference strictly below the limit"};
    Configurable<bool> cfgMixingMchRadiusStrict{"cfgMixingMchRadiusStrict", false, "Require MCH position-radius difference strictly below its limit"};
    Configurable<bool> cfgMixingMchEtaStrict{"cfgMixingMchEtaStrict", false, "Require MCH PV eta difference strictly below its limit"};
    Configurable<double> cfgMixingMaxDeltaMchPtRel{"cfgMixingMaxDeltaMchPtRel", -1., "Absolute donor-target MCH pT difference at PV / target pT; negative disables"};
    Configurable<float> cfgMixingMaxDeltaMchMomentumRel{"cfgMixingMaxDeltaMchMomentumRel", -1.f, "Absolute donor-target MCH total momentum difference at PV / target momentum; negative disables"};
    Configurable<float> cfgMixingMaxDeltaMchPhiPV{"cfgMixingMaxDeltaMchPhiPV", -1.f, "Maximum wrapped donor-target MCH momentum-phi difference at PV (rad); negative disables"};
  } mixingMultiplicityCuts;
  struct OwnedMixMft {
    int64_t id = -1, collision = -1;
    MftMchJointTransport::Fwd state;
    int64_t globalIndex() const { return id; }
    int64_t collisionId() const { return collision; }
  };
  struct MixingPoolEntry {
    mixing::Event event;
    MftMchJointTransport::Fwd mch;
    std::vector<OwnedMixMft> mfts;
  };
  struct OwnedMixEvent {
    uint64_t frame = 0;
    int64_t collision = -1;
    size_t stateCount = 0;
    std::vector<MixingPoolEntry> anchors;
  };
  struct : ConfigurableGroup {
    Configurable<bool> cfgMixingAcrossDF{"cfgMixingAcrossDF", false, "Also use owned donor events from previous DFs in this task instance"};
    Configurable<int> cfgMixingPoolMaxEvents{"cfgMixingPoolMaxEvents", 2000, "Maximum retained donor collisions across DFs (FIFO)"};
    Configurable<int> cfgMixingPoolMaxMftStates{"cfgMixingPoolMaxMftStates", 200000, "Maximum retained MFT state copies; oversized events are skipped intact"};
    Configurable<uint64_t> cfgMixingStreamId{"cfgMixingStreamId", 0, "Source stream ID for merged outputs; 0 generates a random nonzero 64-bit ID per task instance"};
    struct State {
      uint64_t stream = 0;
      mixing::History<OwnedMixEvent> history;
      uint64_t frame = 0;
    } state;
  } mixingPersistence;
  Produces<o2::aod::MftMchMixingPoolStates> tableMixingPoolStates;
  std::unordered_map<int64_t, int> mixingCbMultiplicities;
  std::unordered_map<int64_t, double> mixingCentralities;

  Produces<o2::aod::RawMatchPlaneResiduals> tableMatchPlaneResiduals;
  Produces<o2::aod::DownstreamMatchPlaneResiduals> tableDownstreamResiduals;
  Produces<o2::aod::MatchPlaneEventSels> tableMatchPlaneEventSels;
  Produces<o2::aod::MftMchPairCovariances> tablePairCovariances;
  Produces<o2::aod::MftMchMixingDiagnostics> tableMixingDiagnostics;
  Configurable<bool> cfgSaveMixingDiagnostics{"cfgSaveMixingDiagnostics", false, "Save pair covariance matrices for MCHs with mixing candidates"};
  Configurable<std::string> cfgEventSelection{"cfgEventSelection", "all", "Matching-plane output: all, int8 (Run 3 sel8), upc (MCH-BC FV0-veto candidates)"};
  Configurable<bool> cfgFilterBeforeMatching{"cfgFilterBeforeMatching", true, "For int8/upc: reject MCH tracks before propagation; false keeps the previous late-filter path for numerical comparisons"};
  Configurable<int> cfgUpcMinMchMid{"cfgUpcMinMchMid", 1, "Minimum collision-associated MCH-MID tracks in the same corrected BC"};
  Configurable<int> cfgUpcMaxMchMid{"cfgUpcMaxMchMid", -1, "Maximum MCH-MID tracks in the same corrected BC; -1 disables upper cut"};
  Configurable<int> cfgUpcMinMchTracks{"cfgUpcMinMchTracks", 2, "Minimum standalone MCH tracks in the MCH collision"};
  Configurable<int> cfgUpcMinMftTracks{"cfgUpcMinMftTracks", 2, "Minimum MFT tracks in the MCH collision, including tracks without covariance"};
  Configurable<int> cfgUpcMaxMchTracks{"cfgUpcMaxMchTracks", -1, "Maximum standalone MCH tracks in the MCH collision; -1 disables the upper cut"};
  Configurable<bool> cfgUpcRequireEqualTrackCounts{"cfgUpcRequireEqualTrackCounts", false, "Require equal MFT and standalone MCH counts in the MCH collision; MFT count includes tracks without covariance"};
  Configurable<int> cfgUpcFV0BcWindow{"cfgUpcFV0BcWindow", 0, "FIT veto within +/- this many BCs of corrected MCH BC (legacy option name)"};
  Configurable<float> cfgUpcMaxFV0Amplitude{"cfgUpcMaxFV0Amplitude", 100.f, "UDHelpers FV0A amplitude limit; negative disables FV0A"};
  Configurable<float> cfgUpcFV0ValidTime{"cfgUpcFV0ValidTime", 15.f, "UDHelpers maximum absolute FIT time in ns (legacy option name); out-of-time records fail"};
  Configurable<float> cfgUpcMaxFT0AAmplitude{"cfgUpcMaxFT0AAmplitude", -1.f, "UDHelpers FT0A amplitude limit; negative disables"};
  Configurable<float> cfgUpcMaxFT0CAmplitude{"cfgUpcMaxFT0CAmplitude", -1.f, "UDHelpers FT0C amplitude limit; negative disables"};
  Configurable<float> cfgUpcMaxFDDAAmplitude{"cfgUpcMaxFDDAAmplitude", -1.f, "UDHelpers FDDA charge limit; negative disables"};
  Configurable<float> cfgUpcMaxFDDCAmplitude{"cfgUpcMaxFDDCAmplitude", -1.f, "UDHelpers FDDC charge limit; negative disables"};
  Configurable<bool> cfgUpcRequireUnlikeSign{"cfgUpcRequireUnlikeSign", true, "Require at least one positive and one negative MCH"};
  struct : ConfigurableGroup {
    Configurable<bool> cfgRequireMchMftAcceptance{"cfgRequireMchMftAcceptance", true, "Require PV-propagated MCH to cross every nominal MFT disk acceptance envelope"};
    Configurable<bool> cfgRequireMchPVEta{"cfgRequireMchPVEta", false, "Require propagated MCH eta at its collision PV within the configured strict range"};
    Configurable<double> cfgMchPVEtaMin{"cfgMchPVEtaMin", -3.6, "Strict lower limit on propagated MCH eta at its collision PV"};
    Configurable<double> cfgMchPVEtaMax{"cfgMchPVEtaMax", -2.5, "Strict upper limit on propagated MCH eta at its collision PV"};
    Configurable<double> cfgMchMftAcceptanceEtaMin{"cfgMchMftAcceptanceEtaMin", -3.6, "Nominal MFT envelope eta minimum for a vertex at (0,0,0)"};
    Configurable<double> cfgMchMftAcceptanceEtaMax{"cfgMchMftAcceptanceEtaMax", -2.5, "Nominal MFT envelope eta maximum for a vertex at (0,0,0)"};
    Configurable<std::vector<double>> cfgMchMftAcceptanceDiskZ{"cfgMchMftAcceptanceDiskZ", {-46., -49.3, -53.1, -68.4, -76.8}, "Nominal disk z positions (cm), from MFT Geometry.xml"};
    Configurable<bool> cfgPropMuon{"cfgPropMuon", true, "Propagate MCH to its collision PV before evaluating quality-cut eta and pt"};
    Configurable<std::string> cfgMuonCuts{"cfgMuonCuts", "muonQualityCuts", "PWGDQ cut: muonQualityCuts, muonQualityCuts5SigmaPDCA_Run3, muonQualityCuts10SigmaPDCA, or none"};
  } mchAcceptance;
  HistogramRegistry mchAcceptanceQA{"MchMftAcceptance", {}, OutputObjHandlingPolicy::AnalysisObject};

  // This is the nominal eta-defined geometric envelope, not a chip hit-efficiency
  // prediction. Propagate from the ACTUAL collision PV, including its x/y/z;
  // a fixed eta(PV) cut is a separate configurable selection, not this envelope.
  bool crossesNominalMftAcceptance(const o2::mch::TrackParam& mchAtPV)
  {
    for (const double z : mchAcceptance.cfgMchMftAcceptanceDiskZ.value) {
      if (!(mchAtPV.getZ() > z)) {
        return false;
      }
      auto atDisk = mchAtPV;
      // Official field integration; no absorber lies between the PV and MFT.
      // Only central parameters are needed for this geometric selection.
      if (!o2::mch::TrackExtrap::extrapToZ(atDisk, z)) {
        return false;
      }
      const double r = std::hypot(atDisk.getNonBendingCoor(), atDisk.getBendingCoor());
      const double rMin = std::abs(z / std::sinh(mchAcceptance.cfgMchMftAcceptanceEtaMin.value));
      const double rMax = std::abs(z / std::sinh(mchAcceptance.cfgMchMftAcceptanceEtaMax.value));
      if (!std::isfinite(r) || !(r > rMin && r < rMax)) {
        return false;
      }
    }
    return true;
  }
  std::unique_ptr<AnalysisCut> mchQualityCut;
  int8_t matchingEventSelectionMode = 0;
  Configurable<double> cfgMftAbsorberStepCm{"cfgMftAbsorberStepCm", 0., "MFT absorber material-update maximum delta-z (cm): 0=whole absorber, positive=stepwise; ignored at z>=-90"};
  Configurable<bool> cfgMftUseMchMomentum{"cfgMftUseMchMomentum", true, "At z<=-505, replace MFT momentum magnitude using candidate MCH; propagate full joint errors"};
  Configurable<bool> cfgMftUseMchCharge{"cfgMftUseMchCharge", true, "Use reliable MCH charge for momentum-updated MFT; original MFT q/pt is saved"};
  Configurable<bool> cfgRequireMftCovariance{"cfgRequireMftCovariance", true, "Require stored MFT covariance; false admits all accepted MFT tracks and requires cfgComputeMatchingChi2=false"};
  Configurable<bool> cfgComputeMatchingChi2{"cfgComputeMatchingChi2", true, "Compute matching chi2 and joint errors; false saves nominal parameters/residuals and NaN chi2"};
  Configurable<double> cfgPairJacobianScale{"cfgPairJacobianScale", 1., "Finite-difference scale for joint detector-error propagation; positive"};
  o2::aod::fwdtrackutils::GlobalMuonMatchingPlane globalMatchingPlane;
  // Keep buckets across DFs, exactly as global-muon-matcher::mMchTrackInfos.
  std::unordered_map<int64_t, int> globalMchIterationOrder;

  struct : ConfigurableGroup {
    Configurable<int> cfgMuonTaggingNCrossedMftPlanesLow{"cfgMuonTaggingNCrossedMftPlanesLow", 5, ""};
    Configurable<float> cfgMuonTaggingTrackChi2MchUp{"cfgMuonTaggingTrackChi2MchUp", 5.f, ""};
    Configurable<float> cfgMuonTaggingPMchLow{"cfgMuonTaggingPMchLow", 0.0f, ""};
    Configurable<float> cfgMuonTaggingPtMchLow{"cfgMuonTaggingPtMchLow", 0.7f, ""};
    Configurable<float> cfgMuonTaggingEtaMchLow{"cfgMuonTaggingEtaMchLow", -3.6f, ""};
    Configurable<float> cfgMuonTaggingEtaMchUp{"cfgMuonTaggingEtaMchUp", -2.5f, ""};
    Configurable<float> cfgMuonTaggingRabsLow{"cfgMuonTaggingRabsLow", 17.6f, ""};
    Configurable<float> cfgMuonTaggingRabsUp{"cfgMuonTaggingRabsUp", 89.5f, ""};
    Configurable<float> cfgMuonTaggingPdcaUp{"cfgMuonTaggingPdcaUp", 4.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftFrontLow{"cfgMuonTaggingRadiusAtMftFrontLow", 3.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftFrontUp{"cfgMuonTaggingRadiusAtMftFrontUp", 9.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftBackLow{"cfgMuonTaggingRadiusAtMftBackLow", 5.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftBackUp{"cfgMuonTaggingRadiusAtMftBackUp", 12.f, ""};
  } configMuonTagging;

  ////  Variables for matching method
  Configurable<int> fMatchingMethod{"cfgMatchingMethod", 0, ""};

  ////   Variables for selecting muon tracks
  struct : ConfigurableGroup {
    Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
    Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
    Configurable<float> fRabsLow1{"cfgRabsLow1", 17.6f, ""};
    Configurable<float> fRabsUp1{"cfgRabsUp1", 26.5f, ""};
    Configurable<float> fRabsLow2{"cfgRabsLow2", 26.5f, ""};
    Configurable<float> fRabsUp2{"cfgRabsUp2", 89.5f, ""};
    Configurable<float> fPdcaUp1{"cfgPdcaUp1", 594.f, ""};
    Configurable<float> fPdcaUp2{"cfgPdcaUp2", 324.f, ""};
    Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
    Configurable<float> fMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};
  } legacySelection;

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
  Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
  Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  ///    Variables to add preselection for the matching table
  Configurable<float> fPreselectMatchingX{"cfgPreselectMatchingX", 15.f, ""};
  Configurable<float> fPreselectMatchingY{"cfgPreselectMatchingY", 15.f, ""};

  ///    Variables to event mixing criteria
  Configurable<float> fSaveMixedMatchingParamsRate{"cfgSaveMixedMatchingParamsRate", 0.002f, ""};
  Configurable<int> fEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
  Configurable<float> fEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};

  ////   Variables for selecting tag muon
  Configurable<float> fTagMassWindowMin{"cfgTagMassWindowMin", 2.8f, ""};
  Configurable<float> fTagMassWindowMax{"cfgTagMassWindowMax", 3.3f, ""};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  ////   Variables for Tag matching criteria
  Configurable<float> fSigmaXTagMuonCut{"cfgSigmaXTagMuonCut", 1.f, ""};
  Configurable<float> fMeanXTagMuonCut{"cfgMeanXTagMuonCut", 0.f, ""};
  Configurable<float> fSigmaYTagMuonCut{"cfgSigmaYTagMuonCut", 1.f, ""};
  Configurable<float> fMeanYTagMuonCut{"cfgMeanYTagMuonCut", 1.f, ""};

  Configurable<float> fSigmaEtaTagMuonCut{"cfgSigmaEtaTagMuonCut", 0.2f, ""};
  Configurable<float> fMeanEtaTagMuonCut{"cfgMeanEtaTagMuonCut", 0.f, ""};
  Configurable<float> fSigmaPhiTagMuonCut{"cfgSigmaPhiTagMuonCut", 0.2f, ""};
  Configurable<float> fMeanPhiTagMuonCut{"cfgMeanPhiTagMuonCut", 0.f, ""};

  template <typename MUON, typename Collision>
  class FindTagAndProbe
  {
   private:
    o2::dataformats::GlobalFwdTrack muontrack_at_pv[2];

    TLorentzVector mDimuon;
    MUON muontrack1;
    MUON muontrack2;
    Collision collision;
    int tagIdx, probeIdx;

    int16_t mQ;

    inline void fillCovarianceArray(MUON const& muontrack, float cov[15]) const
    {
      cov[0] = muontrack.cXX();
      cov[1] = muontrack.cXY();
      cov[2] = muontrack.cYY();
      cov[3] = muontrack.cPhiX();
      cov[4] = muontrack.cPhiY();
      cov[5] = muontrack.cPhiPhi();
      cov[6] = muontrack.cTglX();
      cov[7] = muontrack.cTglY();
      cov[8] = muontrack.cTglPhi();
      cov[9] = muontrack.cTglTgl();
      cov[10] = muontrack.c1PtX();
      cov[11] = muontrack.c1PtY();
      cov[12] = muontrack.c1PtPhi();
      cov[13] = muontrack.c1PtTgl();
      cov[14] = muontrack.c1Pt21Pt2();
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoPV(MUON const& muontrack) const
    {
      const double mz = muontrack.z();
      const double mchi2 = muontrack.chi2();
      const float mx = muontrack.x();
      const float my = muontrack.y();
      const float mphi = muontrack.phi();
      const float mtgl = muontrack.tgl();
      const float m1pt = muontrack.signed1Pt();

      float cov[15];
      fillCovarianceArray(muontrack, cov);
      SMatrix5 tpars(mx, my, mphi, mtgl, m1pt);
      SMatrix55 tcovs(cov, cov + 15);

      o2::track::TrackParCovFwd parcovmuontrack{mz, tpars, tcovs, mchi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      o2::globaltracking::MatchGlobalFwd mMatching;
      auto mchtrack = mMatching.FwdtoMCH(gtrack);

      o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
      o2::dataformats::GlobalFwdTrack extrap_muontrack;
      extrap_muontrack.setParameters(fwdtrack.getParameters());
      extrap_muontrack.setZ(fwdtrack.getZ());
      extrap_muontrack.setCovariances(fwdtrack.getCovariances());

      return extrap_muontrack;
    }

    inline void setTagAndProbe()
    {
      if (muontrack1.pt() > muontrack2.pt()) {
        tagIdx = 0;
        probeIdx = 1;
      } else {
        tagIdx = 1;
        probeIdx = 0;
      }
    }

   public:
    inline FindTagAndProbe(const MUON& muon1, const MUON& muon2, const Collision& coll)
      : muontrack_at_pv(), mDimuon(), muontrack1(muon1), muontrack2(muon2), collision(coll), tagIdx(-1), probeIdx(-1), mQ(0)
    {
      mQ = muontrack1.sign() + muontrack2.sign();
      setTagAndProbe();
    }

    void calcMuonPairAtPV()
    {
      muontrack_at_pv[0] = propagateMUONtoPV(muontrack1);
      muontrack_at_pv[1] = propagateMUONtoPV(muontrack2);
      TLorentzVector vMuon1, vMuon2;
      vMuon1.SetPtEtaPhiM(muontrack_at_pv[0].getPt(), muontrack_at_pv[0].getEta(), muontrack_at_pv[0].getPhi(), mMu);
      vMuon2.SetPtEtaPhiM(muontrack_at_pv[1].getPt(), muontrack_at_pv[1].getEta(), muontrack_at_pv[1].getPhi(), mMu);
      mDimuon = vMuon1 + vMuon2;
    }
    inline int getTagMuonIndex() const { return tagIdx; }
    inline int getProbeMuonIndex() const { return probeIdx; }
    inline float getMass() const { return mDimuon.M(); }
    inline float getPt() const { return mDimuon.Pt(); }
    inline float getRap() const { return mDimuon.Rapidity(); }
    inline int16_t getCharge() const { return mQ; }
    inline const o2::dataformats::GlobalFwdTrack& getMuonAtPV(int idx) const { return muontrack_at_pv[idx]; }
  }; // end of class FindTagAndProbe

  template <typename MUON, typename MFT, typename Collision>
  class MatchingParamsML
  {
   private:
    MUON muontrack;
    MFT mfttrack;
    Collision collision;

    float mDX, mDY, mDPt, mDPhi, mDEta;
    float mGlobalMuonPtAtDCA, mGlobalMuonEtaAtDCA, mGlobalMuonPhiAtDCA, mGlobalMuonDCAx, mGlobalMuonDCAy, mGlobalMuonQ;
    int mMatchingType;

    o2::field::MagneticField* fieldB;
    o2::globaltracking::MatchGlobalFwd mMatching;

    inline o2::track::TrackParCovFwd propagateMFTtoMatchingPlane()
    {
      double covArr[15]{0.0};
      SMatrix55 tmftcovs(covArr, covArr + 15);

      SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      o2::track::TrackParCovFwd extrap_mfttrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};

      double propVec[3] = {0.};
      float zPlane = 0.f;
      if (mMatchingType == MCH_FIRST_CLUSTER) {
        propVec[0] = muontrack.x() - mfttrack.x();
        propVec[1] = muontrack.y() - mfttrack.y();
        propVec[2] = muontrack.z() - mfttrack.z();
        zPlane = muontrack.z();
      } else if (mMatchingType == END_OF_ABSORBER || mMatchingType == BEGINING_OF_ABSORBER) {
        auto extrap_muontrack = propagateMUONtoMatchingPlane();
        propVec[0] = extrap_muontrack.getX() - mfttrack.x();
        propVec[1] = extrap_muontrack.getY() - mfttrack.y();
        propVec[2] = extrap_muontrack.getZ() - mfttrack.z();
        zPlane = (mMatchingType == END_OF_ABSORBER) ? -505.f : -90.f;
      } else {
        zPlane = mfttrack.z();
      }

      double centerZ[3] = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
      float Bz = fieldB->getBz(centerZ);
      extrap_mfttrack.propagateToZ(zPlane, Bz); // z in cm
      return extrap_mfttrack;
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoMatchingPlane()
    {
      float cov[15] = {
        muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
        muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
        muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
        muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
        muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};

      SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
      SMatrix55 tcovs(cov, cov + 15);
      double chi2 = muontrack.chi2();

      o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      auto mchtrack = mMatching.FwdtoMCH(gtrack);

      if (mMatchingType == MFT_LAST_CLUSTR) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, mfttrack.z());
      } else if (mMatchingType == END_OF_ABSORBER) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, -505.);
      } else if (mMatchingType == BEGINING_OF_ABSORBER) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, -90.);
      }

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);

      o2::dataformats::GlobalFwdTrack extrap_muontrack;
      extrap_muontrack.setParameters(fwdtrack.getParameters());
      extrap_muontrack.setZ(fwdtrack.getZ());
      extrap_muontrack.setCovariances(fwdtrack.getCovariances());
      return extrap_muontrack;
    }

    inline o2::track::TrackParCovFwd propagateMFTtoDCA()
    {
      double covArr[15]{0.0};
      SMatrix55 tmftcovs(covArr, covArr + 15);

      SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      o2::track::TrackParCovFwd extrap_mfttrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};

      double propVec[3] = {};
      propVec[0] = collision.posX() - mfttrack.x();
      propVec[1] = collision.posY() - mfttrack.y();
      propVec[2] = collision.posZ() - mfttrack.z();

      double centerZ[3] = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
      float Bz = fieldB->getBz(centerZ);
      extrap_mfttrack.propagateToZ(collision.posZ(), Bz); // z in cm
      return extrap_mfttrack;
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoPV()
    {
      float cov[15] = {
        muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
        muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
        muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
        muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
        muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};

      SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
      SMatrix55 tcovs(cov, cov + 15);
      double chi2 = muontrack.chi2();

      o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      auto mchtrack = mMatching.FwdtoMCH(gtrack);
      o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
      o2::dataformats::GlobalFwdTrack extrap_muontrack;
      extrap_muontrack.setParameters(fwdtrack.getParameters());
      extrap_muontrack.setZ(fwdtrack.getZ());
      extrap_muontrack.setCovariances(fwdtrack.getCovariances());

      return extrap_muontrack;
    }

   public:
    enum MATCHING_TYPE { MCH_FIRST_CLUSTER,
                         MFT_LAST_CLUSTR,
                         END_OF_ABSORBER,
                         BEGINING_OF_ABSORBER };

    MatchingParamsML(MUON const& muon, MFT const& mft, Collision const& coll, int MType, o2::field::MagneticField* field) : muontrack(muon), mfttrack(mft), collision(coll), mDX(0.f), mDY(0.f), mDPt(0.f), mDPhi(0.f), mDEta(0.f), mGlobalMuonPtAtDCA(0.f), mGlobalMuonEtaAtDCA(0.f), mGlobalMuonPhiAtDCA(0.f), mGlobalMuonDCAx(0.f), mGlobalMuonDCAy(0.f), mGlobalMuonQ(0.f), mMatchingType(MType), fieldB(field) {}
    void calcMatchingParams()
    {
      auto mfttrack_on_matchingP = propagateMFTtoMatchingPlane();
      auto muontrack_on_matchingP = propagateMUONtoMatchingPlane();

      float dphiRaw = mfttrack_on_matchingP.getPhi() - muontrack_on_matchingP.getPhi();
      float dphi = TVector2::Phi_mpi_pi(dphiRaw);
      float deta = mfttrack_on_matchingP.getEta() - muontrack_on_matchingP.getEta();

      mDX = mfttrack_on_matchingP.getX() - muontrack_on_matchingP.getX();
      mDY = mfttrack_on_matchingP.getY() - muontrack_on_matchingP.getY();
      mDPt = mfttrack_on_matchingP.getPt() - muontrack_on_matchingP.getPt();
      mDPhi = dphi;
      mDEta = deta;
    }

    void calcGlobalMuonParams()
    {
      auto mfttrack_at_dca = propagateMFTtoDCA();
      auto muontrack_at_pv = propagateMUONtoPV();

      float momentum = muontrack_at_pv.getP();
      float theta = mfttrack_at_dca.getTheta();
      float phiTrack = mfttrack_at_dca.getPhi();
      float px = momentum * std::sin(theta) * std::cos(phiTrack);
      float py = momentum * std::sin(theta) * std::sin(phiTrack);

      mGlobalMuonQ = muontrack.sign() + mfttrack.sign();
      mGlobalMuonPtAtDCA = std::sqrt(px * px + py * py);
      mGlobalMuonEtaAtDCA = mfttrack_at_dca.getEta();
      mGlobalMuonPhiAtDCA = mfttrack_at_dca.getPhi();
      mGlobalMuonDCAx = mfttrack_at_dca.getX() - collision.posX();
      mGlobalMuonDCAy = mfttrack_at_dca.getY() - collision.posY();
    }

    inline float getDx() const { return mDX; }
    inline float getDy() const { return mDY; }
    inline float getDphi() const { return mDPhi; }
    inline float getDeta() const { return mDEta; }
    inline float getDpt() const { return mDPt; }
    inline float getGMPtAtDCA() const { return mGlobalMuonPtAtDCA; }
    inline float getGMEtaAtDCA() const { return mGlobalMuonEtaAtDCA; }
    inline float getGMPhiAtDCA() const { return mGlobalMuonPhiAtDCA; }
    inline float getGMDcaX() const { return mGlobalMuonDCAx; }
    inline float getGMDcaY() const { return mGlobalMuonDCAy; }
    inline float getGMDcaXY() const { return std::sqrt(mGlobalMuonDCAx * mGlobalMuonDCAx + mGlobalMuonDCAy * mGlobalMuonDCAy); }
    inline int16_t getGMQ() const { return static_cast<int16_t>(mGlobalMuonQ); }

  }; // end of class MatchingParamsML

  template <typename MUON, typename Collisions>
  o2::dataformats::GlobalFwdTrack propagateMUONtoPV(MUON const& muontrack, Collisions const& collisions)
  {
    auto collision = collisions.rawIteratorAt(muontrack.collisionId());
    o2::globaltracking::MatchGlobalFwd mMatching;
    o2::dataformats::GlobalFwdTrack extrap_muontrack;

    SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
    std::vector<float> v1{muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
                          muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
                          muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
                          muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
                          muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    double chi2 = muontrack.chi2();
    o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

    o2::dataformats::GlobalFwdTrack gtrack;
    gtrack.setParameters(tpars);
    gtrack.setZ(parcovmuontrack.getZ());
    gtrack.setCovariances(tcovs);
    auto mchtrack = mMatching.FwdtoMCH(gtrack);

    o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

    auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
    extrap_muontrack.setParameters(fwdtrack.getParameters());
    extrap_muontrack.setZ(fwdtrack.getZ());
    extrap_muontrack.setCovariances(fwdtrack.getCovariances());

    return extrap_muontrack;
  }

  inline bool isGoodTagDimuon(float M)
  {
    return !(M < fTagMassWindowMin || M > fTagMassWindowMax);
  }

  inline bool isGoodTagMatching(float mDX, float mDY, float mDEta, float mDPhi)
  {
    float dxNorm = (mDX - fMeanXTagMuonCut) / (fSigmaXTagMuonCut * 3);
    float dyNorm = (mDY - fMeanYTagMuonCut) / (fSigmaYTagMuonCut * 3);
    float detaNorm = (mDEta - fMeanEtaTagMuonCut) / (fSigmaEtaTagMuonCut * 3);
    float dphiNorm = (mDPhi - fMeanPhiTagMuonCut) / (fSigmaPhiTagMuonCut * 3);

    float rTagXY = dxNorm * dxNorm + dyNorm * dyNorm;
    float rTagEtaPhi = detaNorm * detaNorm + dphiNorm * dphiNorm;

    return (rTagXY < 1.f && rTagEtaPhi > 0.f);
  }

  template <typename MUON, typename MFT>
  bool isCorrectMatching(MUON const& muontrack, MFT const& mfttrack)
  {

    int idmuon = muontrack.mcParticleId();
    int idmft = mfttrack.mcParticleId();

    if (idmuon == -1 || idmft == -1)
      return false;
    if (idmuon != idmft)
      return false;
    else
      return true;
  };

  template <typename MUON>
  bool isGoodMuonQuality(const MUON& muontrack)
  {
    if (!muontrack.has_collision())
      return false;
    if (muontrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
      return false;
    if (muontrack.chi2() > legacySelection.fTrackChi2MchUp)
      return false;
    if (legacySelection.fRabsLow1 > muontrack.rAtAbsorberEnd() || muontrack.rAtAbsorberEnd() > legacySelection.fRabsUp2)
      return false;
    if (muontrack.rAtAbsorberEnd() < legacySelection.fRabsUp1 && legacySelection.fPdcaUp1 < muontrack.pDca())
      return false;
    if (muontrack.rAtAbsorberEnd() > legacySelection.fRabsLow2 && legacySelection.fPdcaUp2 < muontrack.pDca())
      return false;
    return true;
  }

  template <typename MUON>
  bool isGoodMuonKine(const MUON& muontrack)
  {
    if (legacySelection.fEtaMchLow > muontrack.getEta() || muontrack.getEta() > legacySelection.fEtaMchUp)
      return false;
    return true;
  }

  template <typename MFT>
  bool isGoodMFTQuality(const MFT& mfttrack)
  {
    if (!mfttrack.has_collision())
      return false;
    if (mfttrack.chi2() > fTrackChi2MftUp)
      return false;
    if (mfttrack.nClusters() < fTrackNClustMftLow)
      return false;
    return true;
  }

  template <typename MFT>
  bool isGoodMFTKine(const MFT& mfttrack)
  {
    if (fEtaMftLow > mfttrack.getEta() || mfttrack.getEta() > fEtaMftUp)
      return false;
    return true;
  }

  inline bool isPassMatchingPreselection(float Dx, float Dy)
  {
    return !(std::abs(Dx) > fPreselectMatchingX || std::abs(Dy) > fPreselectMatchingY);
  }

  template <typename MUONs, typename Collisions>
  void setMUONs(MUONs const& muontracks, Collisions const& collisions)
  {
    for (const auto& muontrack : muontracks) {
      if (!isGoodMuonQuality(muontrack))
        continue;
      o2::dataformats::GlobalFwdTrack muontrack_at_pv = propagateMUONtoPV(muontrack, collisions);
      if (!isGoodMuonKine(muontrack_at_pv))
        continue;

      auto collision = collisions.rawIteratorAt(muontrack.collisionId());

      bool& has = map_has_muontracks_collisions[muontrack.collisionId()];
      has = true;

      vector<int64_t>& arr_muontracks = map_muontracks[collision.globalIndex()];
      arr_muontracks.push_back(muontrack.globalIndex());
    }
  }

  template <typename MFT, typename Collisions>
  o2::track::TrackParCovFwd PropagateMFTtoDCA(MFT const& mfttrack, Collisions const& collisions, o2::field::MagneticField* field)
  {
    auto collision = collisions.rawIteratorAt(mfttrack.collisionId());
    std::vector<double> mftv1;
    SMatrix55 mftcovs{mftv1.begin(), mftv1.end()};
    SMatrix5 mftpars = {mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt()};
    o2::track::TrackParCovFwd mftpartrack = {mfttrack.z(), mftpars, mftcovs, mfttrack.chi2()};
    double propVec[3] = {fabs(mfttrack.x() - collision.posX()),
                         fabs(mfttrack.y() - collision.posY()),
                         fabs(mfttrack.z() - collision.posZ())};
    double centerZ[3] = {mfttrack.x() - propVec[0] / 2.,
                         mfttrack.y() - propVec[1] / 2.,
                         mfttrack.z() - propVec[2] / 2.};
    float Bz = field->getBz(centerZ);
    mftpartrack.propagateToZ(collision.posZ(), Bz);
    return mftpartrack;
  }

  template <typename MFTs, typename Collisions>
  void setMFTs(MFTs const& mfttracks, Collisions const& collisions, o2::field::MagneticField* field)
  {
    for (const auto& mfttrack : mfttracks) {
      if (!isGoodMFTQuality(mfttrack))
        continue;

      o2::track::TrackParCovFwd mfttrack_at_dca = PropagateMFTtoDCA(mfttrack, collisions, field);
      if (!isGoodMFTKine(mfttrack_at_dca))
        continue;

      auto collision = collisions.rawIteratorAt(mfttrack.collisionId());

      map_vtxz[mfttrack.collisionId()] = collision.posZ();
      map_nmfttrack[mfttrack.collisionId()] += 1;

      bool& has = map_has_mfttracks_collisions[mfttrack.collisionId()];
      has = true;

      vector<int64_t>& arr_mfttracks = map_mfttracks[collision.globalIndex()];
      arr_mfttracks.push_back(mfttrack.globalIndex());
    }
  }

  Produces<o2::aod::MatchParams> tableMatchingParams;
  Produces<o2::aod::TagMatchParams> tableTagMatchingParams;
  Produces<o2::aod::ProbeMatchParams> tableProbeMatchingParams;
  Produces<o2::aod::MixMatchParams> tableMixMatchingParams;
  Produces<o2::aod::MuonPair> tableMuonPair;

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;

  o2::field::MagneticField* fieldB;
  o2::ccdb::CcdbApi ccdbApi;

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    mRunNumber = bc.runNumber();
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(geoPath);
    }
    if (doprocessMatchingPlane || doprocessMatchingPlaneCentrality) {
      // Match both geometry loads and detector-element navigation performed by
      // global-muon-matcher, even with its realignment disabled.
      const std::array<int, 10> elementsPerChamber{4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
      for (int geometryCopy = 0; geometryCopy < 2; ++geometryCopy) {
        auto* geometry = ccdbManager->getForTimeStamp<TGeoManager>(geoPath, bc.timestamp());
        if (!geometry) {
          LOGF(fatal, "Aligned geometry is unavailable at the BC timestamp");
        }
        ccdbManager->clearCache(geoPath);
        auto transformation = o2::mch::geo::transformationFromTGeoManager(*geometry);
        for (int chamber = 1; chamber <= 10; ++chamber) {
          for (int element = 0; element < elementsPerChamber[chamber - 1]; ++element) {
            (void)transformation(100 * chamber + element);
          }
        }
      }
    }
    o2::mch::TrackExtrap::setField();
    fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  }

  void run(o2::framework::ProcessingContext& context)
  {
    mixingPersistence.state.frame = context.services().get<o2::framework::TimingInfo>().timeslice;
  }
  void init(o2::framework::InitContext&)
  {
    if (mchAcceptance.cfgRequireMchPVEta &&
        (!std::isfinite(mchAcceptance.cfgMchPVEtaMin.value) || !std::isfinite(mchAcceptance.cfgMchPVEtaMax.value) ||
         !(mchAcceptance.cfgMchPVEtaMin < mchAcceptance.cfgMchPVEtaMax))) {
      LOGF(fatal, "Invalid MCH PV eta selection range");
    }
    if (mchAcceptance.cfgRequireMchMftAcceptance) {
      if (!std::isfinite(mchAcceptance.cfgMchMftAcceptanceEtaMin.value) || !std::isfinite(mchAcceptance.cfgMchMftAcceptanceEtaMax.value) ||
          !(mchAcceptance.cfgMchMftAcceptanceEtaMin < mchAcceptance.cfgMchMftAcceptanceEtaMax && mchAcceptance.cfgMchMftAcceptanceEtaMax < 0.) || mchAcceptance.cfgMchMftAcceptanceDiskZ.value.empty()) {
        LOGF(fatal, "Invalid nominal MFT acceptance envelope");
      }
      for (const double z : mchAcceptance.cfgMchMftAcceptanceDiskZ.value) {
        if (!std::isfinite(z) || z >= 0.) {
          LOGF(fatal, "MFT acceptance disk z must be finite and negative");
        }
      }
    }
    mchAcceptanceQA.add("status", "MCH acceptance: 0=tested, 1=PV failure, 2=envelope failure, 3=pass;status;MCH tracks", HistType::kTH1F, {{4, -0.5, 3.5}});
    mchAcceptanceQA.add("pvEtaStatus", "MCH PV eta: 0=tested, 1=PV failure, 2=eta failure, 3=pass;status;MCH tracks", HistType::kTH1F, {{4, -0.5, 3.5}});
    LOGF(info, "MCH acceptance: disk envelope=%d, PV eta cut=%d, strict eta range=(%g,%g)",
         int(mchAcceptance.cfgRequireMchMftAcceptance.value), int(mchAcceptance.cfgRequireMchPVEta.value),
         mchAcceptance.cfgMchPVEtaMin.value, mchAcceptance.cfgMchPVEtaMax.value);
    mixingPersistence.state.stream = mixingPersistence.cfgMixingStreamId.value;
    if (mixingPersistence.state.stream == 0) {
      std::random_device random;
      do {
        mixingPersistence.state.stream = (uint64_t(random()) << 32) ^ uint64_t(random());
      } while (mixingPersistence.state.stream == 0);
    }
    LOGF(info, "Mixing source stream ID: %llu", static_cast<unsigned long long>(mixingPersistence.state.stream));
    if (mixingPersistence.cfgMixingPoolMaxEvents <= 0 || mixingPersistence.cfgMixingPoolMaxMftStates <= 0) {
      LOGF(fatal, "Cross-DF pool limits must be positive");
    }
    if (!cfgRequireMftCovariance && cfgComputeMatchingChi2 && !inferenceConfig.cfgEvaluationMode) {
      LOGF(fatal, "All-MFT training without inference requires cfgComputeMatchingChi2=false");
    }
    if (inferenceConfig.cfgEvaluationMode) {
      if (doprocessLegacy || (!doprocessMatchingPlane && !doprocessMatchingPlaneCentrality))
        LOGF(fatal, "Inference requires matching-plane process, no legacy process");
      if (inferenceConfig.cfgEnableChi2Matching && !cfgComputeMatchingChi2)
        LOGF(fatal, "Chi2 inference requires cfgComputeMatchingChi2=true (only actual covariance rows are used)");
      if (!inferenceConfig.cfgTrainingMode && cfgEnableEventMixing)
        LOGF(fatal, "Inference-only: disable MFT-MCH training event mixing; use cfgEnableDimuonMixing");
      if ((inferenceConfig.cfgDimuonMixPvMultiplicityRel >= 0. || inferenceConfig.cfgDimuonMixDeltaCentrality >= 0.) && !doprocessMatchingPlaneCentrality)
        LOGF(fatal, "Requested dimuon centrality/PV-mult pool requires centrality process");
    }
    inference.init(inferenceConfig, inferenceHistos);
    if (inferenceConfig.cfgEvaluationMode) {
      if (inferenceConfig.cfgEnableMlMatching && (!cfgMftUseMchMomentum || cfgMftUseMchCharge.value != inference.requiredMchCharge || cfgMftAbsorberStepCm.value != inference.requiredAbsorberStep))
        LOGF(fatal, "Propagation settings disagree with trained ML feature contract");
      std::ostringstream meta;
      meta << "ML=" << inferenceConfig.cfgMlModel.value << ";contract=" << inference.contract
           << ";MLcut=" << inferenceConfig.cfgMlScoreMin.value << ";chi2=" << inferenceConfig.cfgChi2Method.value
           << ";chi2max=" << inferenceConfig.cfgChi2Max.value << ";MIDrequired=" << inferenceConfig.cfgRequireMid.value
           << ";kinematics=MCH total P at PV + MFT position/direction, MCH charge, then official 3D helix DCA"
           << ";mixDepth=" << inferenceConfig.cfgDimuonMixingDepth.value << ";mixDeltaZ=" << inferenceConfig.cfgDimuonMixDeltaZ.value
           << ";mixNormalize=" << inferenceConfig.cfgDimuonMixNormalize.value << ";source score is NOT correct-match probability\n"
           << inference.contractContent;
      inferenceMetadata->SetTitle(meta.str().c_str());
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel.value) || mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel < -1.) {
      LOGF(fatal, "Invalid mixing Central barrel multiplicity limit");
    }
    if (cfgEnableEventMixing && mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel >= 0. && !doprocessMatchingPlaneCentrality) {
      LOGF(fatal, "Central barrel mixing cut requires processMatchingPlaneCentrality=true and the official PVMults table");
    }
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    ccdbManager->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
    if (cfgEnableEventMixing && !doprocessMatchingPlane && !doprocessMatchingPlaneCentrality) {
      LOGF(fatal, "Event mixing requires processMatchingPlane");
    }
    if (doprocessMatchingPlane && doprocessMatchingPlaneCentrality) {
      LOGF(fatal, "Enable only one matching-plane process");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMftMultiplicityRel.value) || mixingMultiplicityCuts.cfgMixingMaxDeltaMftMultiplicityRel < -1.f ||
        mixingMultiplicityCuts.cfgMixingMaxDeltaMchMultiplicity < -1 || !std::isfinite(cfgMixingMaxDeltaCentrality.value) || cfgMixingMaxDeltaCentrality < -1.f) {
      LOGF(fatal, "Invalid mixing multiplicity/centrality limits");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMchPtRel.value)) {
      LOGF(fatal, "Invalid mixing MCH pT limit");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMchMomentumRel.value)) {
      LOGF(fatal, "Invalid mixing MCH momentum limit");
    }
    if (cfgMixingAlign505 && (!cfgEnableEventMixing || !cfgMftUseMchMomentum)) {
      LOGF(fatal, "Alignment requires event mixing and donor-MCH momentum transport");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV.value) || !std::isfinite(cfgMixingMaxDeltaMchEtaPV.value)) {
      LOGF(fatal, "Invalid mixing MCH direction limits");
    }
    if (cfgEnableEventMixing && cfgMixingMaxDeltaCentrality >= 0.f && !doprocessMatchingPlaneCentrality) {
      LOGF(fatal, "Centrality mixing cut requires processMatchingPlaneCentrality=true and processMatchingPlane=false");
    }
    if (cfgMixingDepth < -1 || cfgMinDeltaBc < 0 || !std::isfinite(cfgMaxDeltaPhi.value) ||
        !std::isfinite(cfgMaxDeltaR.value) || !std::isfinite(cfgMaxDeltaAttempts.value) ||
        !std::isfinite(cfgMaxDeltaZ.value)) {
      LOGF(fatal, "Invalid event-mixing configuration");
    }
    if (doprocessMatchingPlane || doprocessMatchingPlaneCentrality) {
      if (doprocessLegacy) {
        LOGF(fatal, "Enable only one of processLegacy and processMatchingPlane");
      }
      const auto& cutName = mchAcceptance.cfgMuonCuts.value;
      if (cutName != "none" && cutName != "muonQualityCuts" && cutName != "muonQualityCuts5SigmaPDCA_Run3" && cutName != "muonQualityCuts10SigmaPDCA") {
        LOGF(fatal, "Unsupported mchAcceptance.cfgMuonCuts: %s", cutName.c_str());
      }
      if (cutName != "none") {
        mchQualityCut.reset(o2::aod::dqcuts::GetAnalysisCut(cutName.c_str()));
        if (!mchQualityCut) {
          LOGF(fatal, "Cannot load PWGDQ muon cut");
        }
      }
      if (cfgEventSelection.value == "all") {
        matchingEventSelectionMode = 0;
      } else if (cfgEventSelection.value == "int8" || cfgEventSelection.value == "sel8") {
        matchingEventSelectionMode = 1;
      } else if (cfgEventSelection.value == "upc") {
        matchingEventSelectionMode = 2;
      } else {
        LOGF(fatal, "cfgEventSelection must be all, int8/sel8 or upc");
      }
      if (cfgUpcMinMchTracks < 0 || cfgUpcMinMftTracks < 0 || (cfgUpcMaxMchTracks != -1 && cfgUpcMaxMchTracks < cfgUpcMinMchTracks) || cfgUpcMinMchMid < 0 || (cfgUpcMaxMchMid != -1 && cfgUpcMaxMchMid < cfgUpcMinMchMid) ||
          cfgUpcFV0BcWindow < 0 || !std::isfinite(cfgUpcMaxFT0AAmplitude.value) || !std::isfinite(cfgUpcMaxFT0CAmplitude.value) || !std::isfinite(cfgUpcMaxFDDAAmplitude.value) || !std::isfinite(cfgUpcMaxFDDCAmplitude.value) || !std::isfinite(cfgUpcMaxFV0Amplitude.value) ||
          !std::isfinite(cfgUpcFV0ValidTime.value) || cfgUpcFV0ValidTime < 0.f) {
        LOGF(fatal, "Invalid UPC multiplicity, FV0 threshold, BC window or valid-time configuration");
      }
      LOGF(info, "Matching-plane event selection: %s (mode %d); only MCH output is selected, MFT candidates retain the original time window", cfgEventSelection.value.c_str(), matchingEventSelectionMode);
      if (!std::isfinite(cfgMftAbsorberStepCm.value) || cfgMftAbsorberStepCm < 0.) {
        LOGF(fatal, "cfgMftAbsorberStepCm must be finite and nonnegative");
      }
      LOGF(info, "Always save both planes: -77.5 cm official transport; -505 cm material step=%g cm (0=whole absorber); magnetic integration remains official/adaptive", cfgMftAbsorberStepCm.value);
      if (!std::isfinite(cfgPairJacobianScale.value) || cfgPairJacobianScale <= 0.) {
        LOGF(fatal, "cfgPairJacobianScale must be finite and positive");
      }
      o2::conf::ConfigurableParam::setValue("FieldOriginBias.z", "0");
    }
  }

  void processLegacy(MyCollisions const& collisions,
                     MyBCs const& bcs,
                     MyMUONs const& muontracks,
                     MyMFTs const& mfttracks)
  {
    LOG(info) << "Process()  ";
    map_muontracks.clear();
    map_mfttracks.clear();
    map_collisions.clear();
    map_has_muontracks_collisions.clear();
    map_has_mfttracks_collisions.clear();

    initCCDB(bcs.begin());
    setMUONs(muontracks, collisions);
    setMFTs(mfttracks, collisions, fieldB);

    for (auto map_has_muontracks_collision : map_has_muontracks_collisions) {
      auto idmuontrack_collisions = map_has_muontracks_collision.first;
      for (auto map_has_mfttracks_collision : map_has_mfttracks_collisions) {
        auto idmfttrack_collisions = map_has_mfttracks_collision.first;
        if (idmuontrack_collisions != idmfttrack_collisions)
          continue;
        map_collisions[idmfttrack_collisions] = true;
      }
    }

    for (auto const& map_collision : map_collisions) {
      auto const& collision = collisions.rawIteratorAt(map_collision.first);

      for (auto const& imuontrack1 : map_muontracks[map_collision.first]) {
        auto const& muontrack1 = muontracks.rawIteratorAt(imuontrack1);

        for (auto const& imfttrack1 : map_mfttracks[map_collision.first]) {
          auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);

          MatchingParamsML<MyMUON, MyMFT, MyCollision> matching(muontrack1, mfttrack1, collision, fMatchingMethod, fieldB);
          matching.calcMatchingParams();

          if (!isPassMatchingPreselection(matching.getDx(), matching.getDy()))
            continue;

          matching.calcGlobalMuonParams();

          bool isTrue = false;

          tableMatchingParams(matching.getGMPtAtDCA(),
                              matching.getGMEtaAtDCA(),
                              static_cast<int16_t>(matching.getGMQ()),
                              matching.getDpt(),
                              matching.getDx(),
                              matching.getDy(),
                              matching.getDeta(),
                              matching.getDphi(),
                              isTrue);
        }

        for (auto const& map_mfttrack : map_mfttracks) {
          if (map_mfttrack.first == map_collision.first)
            continue;
          if (fabs(map_vtxz[map_mfttrack.first] - map_vtxz[map_collision.first]) > fEventMaxDeltaVtxZ)
            continue;
          if (fabs(map_nmfttrack[map_mfttrack.first] - map_nmfttrack[map_collision.first]) > fEventMaxDeltaNMFT)
            continue;

          for (auto const& imfttrack1 : map_mfttrack.second) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matching(muontrack1, mfttrack1, collision, fMatchingMethod, fieldB);
            matching.calcMatchingParams();
            if (!isPassMatchingPreselection(matching.getDx(), matching.getDy()))
              continue;
            matching.calcGlobalMuonParams();

            bool isTrue = false;

            tableMixMatchingParams(matching.getGMPtAtDCA(),
                                   matching.getGMEtaAtDCA(),
                                   static_cast<int16_t>(matching.getGMQ()),
                                   matching.getDpt(),
                                   matching.getDx(),
                                   matching.getDy(),
                                   matching.getDeta(),
                                   matching.getDphi(),
                                   isTrue);
          }
        }

        for (auto const& imuontrack2 : map_muontracks[map_collision.first]) {

          if (imuontrack1 >= imuontrack2)
            continue;

          auto const& muontrack2 = muontracks.rawIteratorAt(imuontrack2);

          FindTagAndProbe<MyMUON, MyCollision> tagdimuon(muontrack1, muontrack2, collision);
          tagdimuon.calcMuonPairAtPV();
          tableMuonPair(tagdimuon.getCharge(), tagdimuon.getMass(), tagdimuon.getPt(), tagdimuon.getRap());

          if (!isGoodTagDimuon(tagdimuon.getMass()))
            continue;

          auto tagmuontrack = muontrack1;
          auto probemuontrack = muontrack2;

          if (tagdimuon.getTagMuonIndex() == 1) {
            tagmuontrack = muontrack2;
            probemuontrack = muontrack1;
          }

          int nTagMFTCand = 0;
          int nProbeMFTCand = 0;

          int IndexTagMFTCand = -1;
          float tagGMPtAtDCA = 0;
          // float tagGMEtaAtDCA = 0;

          float minimumR = 9999.;
          int minimumIndexProbeMFTCand = -1;

          unordered_map<int, vector<float>> map_tagMatchingParams;
          unordered_map<int, vector<float>> map_probeMatchingParams;

          for (auto const& imfttrack1 : map_mfttracks[map_collision.first]) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matchingTag(tagmuontrack, mfttrack1, collision, fMatchingMethod, fieldB);
            matchingTag.calcMatchingParams();
            matchingTag.calcGlobalMuonParams();
            if (isGoodTagMatching(matchingTag.getDx(), matchingTag.getDy(), matchingTag.getDeta(), matchingTag.getDphi()) &&
                isPassMatchingPreselection(matchingTag.getDx(), matchingTag.getDy())) {
              bool isTrue = false;
              tableTagMatchingParams(matchingTag.getGMPtAtDCA(),
                                     matchingTag.getGMEtaAtDCA(),
                                     matchingTag.getGMQ(),
                                     matchingTag.getDpt(),
                                     matchingTag.getDx(),
                                     matchingTag.getDy(),
                                     matchingTag.getDeta(),
                                     matchingTag.getDphi(),
                                     isTrue);
              IndexTagMFTCand = mfttrack1.globalIndex();
              tagGMPtAtDCA = matchingTag.getGMPtAtDCA();
              // tagGMEtaAtDCA = matchingTag.getGMEtaAtDCA();
              ++nTagMFTCand;
            }
          } // end of loop imfttrack1

          if (nTagMFTCand != 1)
            continue;
          for (auto const& imfttrack1 : map_mfttracks[map_collision.first]) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            if (mfttrack1.globalIndex() == IndexTagMFTCand)
              continue;
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matchingProbe(probemuontrack, mfttrack1, collision, fMatchingMethod, fieldB);
            matchingProbe.calcMatchingParams();
            if (isPassMatchingPreselection(matchingProbe.getDx(), matchingProbe.getDy())) {
              float R = sqrt(matchingProbe.getDx() * matchingProbe.getDx() + matchingProbe.getDy() * matchingProbe.getDy());
              bool isTrue = false;
              matchingProbe.calcGlobalMuonParams();
              vector<float>& probeMatchingParams = map_probeMatchingParams[nProbeMFTCand];
              probeMatchingParams.push_back(tagGMPtAtDCA);
              probeMatchingParams.push_back(matchingProbe.getGMPtAtDCA());
              probeMatchingParams.push_back(matchingProbe.getGMEtaAtDCA());
              probeMatchingParams.push_back(matchingProbe.getGMQ());
              probeMatchingParams.push_back(matchingProbe.getDpt());
              probeMatchingParams.push_back(matchingProbe.getDx());
              probeMatchingParams.push_back(matchingProbe.getDy());
              probeMatchingParams.push_back(matchingProbe.getDeta());
              probeMatchingParams.push_back(matchingProbe.getDphi());
              probeMatchingParams.push_back(isTrue);
              if (R < minimumR) {
                minimumIndexProbeMFTCand = nProbeMFTCand;
                minimumR = R;
              }
              ++nProbeMFTCand;
            }
          } // end of loop imfttrack1

          if (nProbeMFTCand < 1)
            continue;

          if (minimumIndexProbeMFTCand > -1) {
            vector<float>& probeMatchingParams = map_probeMatchingParams[minimumIndexProbeMFTCand];
            tableProbeMatchingParams(probeMatchingParams[0],
                                     probeMatchingParams[1],
                                     probeMatchingParams[2],
                                     static_cast<int16_t>(probeMatchingParams[3]),
                                     probeMatchingParams[4],
                                     probeMatchingParams[5],
                                     probeMatchingParams[6],
                                     probeMatchingParams[7],
                                     probeMatchingParams[8],
                                     static_cast<bool>(probeMatchingParams[9]));
          }

        } // end of loop imuontrack2
      } // end of loop imuontrack1
    } // end of loop map_collision

  } // end of processMC

  void processMatchingPlaneImpl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                aod::BCsWithTimestamps const& bcs,
                                MyMUONs const& muontracks,
                                aod::MFTTracks const& mfttracks,
                                aod::MFTTracksCov const& mftCovs,
                                aod::FV0As const& fv0s,
                                aod::FT0s const& ft0s,
                                aod::FDDs const& fdds, aod::Origins const& origins)
  {
    using o2::aod::fwdtrackutils::getTrackParCovFwd;
    if (!inferenceConfig.cfgTrainingMode && !inferenceConfig.cfgEvaluationMode)
      return;
    auto& history = mixingPersistence.state.history;
    const uint64_t frame = mixingPersistence.state.frame;
    bool poolReset = false;
    // Publish even on early return, so MC caches see every source frame.
    struct PublishOnExit {
      std::function<void()> f;
      ~PublishOnExit() { f(); }
    } publishPool{[&]() {
      std::set<uint64_t> frames;
      for (const auto& event : history.events)
        frames.insert(event.frame);
      tableMixingPoolStates(frame, std::vector<uint64_t>(frames.begin(), frames.end()),
                            history.events.size(), history.states, history.evicted, history.oversized, poolReset, mixingPersistence.state.stream);
    }};
    if (!mixingPersistence.cfgMixingAcrossDF)
      history.clear();
    if (bcs.size() > 0)
      poolReset = history.setRun(bcs.begin().runNumber());
    if (bcs.size() == 0) {
      return;
    }
    const bool filterBeforeMatching = cfgFilterBeforeMatching && (matchingEventSelectionMode != 0 || mchQualityCut != nullptr || mchAcceptance.cfgRequireMchMftAcceptance || mchAcceptance.cfgRequireMchPVEta);

    globalMchIterationOrder.clear();
    for (const auto& mch : muontracks) {
      if (static_cast<int>(mch.trackType()) > 2) {
        globalMchIterationOrder.try_emplace(mch.globalIndex(), 0);
      }
    }

    // Event selection is anchored on the MCH side. UPC multiplicity counts use
    // the original collision ID and all associated MFT tracks before covariance
    // selection. Candidate matching retains the original timing definition.
    // Corrected BC follows PWGUD upcCandProducerMuon (trackTime / LHC bunch spacing + 1).
    // This mode deliberately retains the original collision-association scope;
    // ambiguous/unassigned tracks need a separate BC/time association extension.
    std::unordered_map<int64_t, int> numMchPerCollision;
    std::unordered_map<int64_t, int> numMftPerCollision;
    std::unordered_map<int64_t, int> numPositiveMchPerCollision;
    std::unordered_map<int64_t, int> numNegativeMchPerCollision;
    for (const auto& mch : muontracks) {
      // Global tracks (types 0--2) duplicate standalone MCH tracks: do not count.
      if (static_cast<int>(mch.trackType()) > 2 && mch.has_collision()) {
        ++numMchPerCollision[mch.collisionId()];
        if (std::isfinite(mch.signed1Pt())) {
          if (mch.signed1Pt() > 0.f) {
            ++numPositiveMchPerCollision[mch.collisionId()];
          } else if (mch.signed1Pt() < 0.f) {
            ++numNegativeMchPerCollision[mch.collisionId()];
          }
        }
      }
    }
    for (const auto& mft : mfttracks) {
      if (mft.has_collision()) {
        ++numMftPerCollision[mft.collisionId()];
      }
    }
    std::unordered_map<int64_t, int64_t> mchGlobalBCs;
    std::map<int64_t, int> numMchMidPerBC;
    for (const auto& mch : muontracks) {
      if (static_cast<int>(mch.trackType()) <= 2 || !mch.has_collision() || !std::isfinite(mch.trackTime())) {
        continue;
      }
      auto collision = collisions.rawIteratorAt(mch.collisionId());
      auto bc = bcs.rawIteratorAt(collision.bcId());
      int64_t correctedBC = static_cast<int64_t>(bc.globalBC()) + TMath::FloorNint(mch.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + 1.);
      mchGlobalBCs.emplace(mch.globalIndex(), correctedBC);
      if (mch.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        ++numMchMidPerBC[correctedBC];
      }
    }
    // Apply the official helper to every FIT record; duplicate BCs must all pass.
    // Keep amplitudes separately for the existing FV0 diagnostic columns.
    std::map<int64_t, float> fv0Amplitudes;
    std::map<int64_t, bool> cleanFITPerBC;
    const std::vector<float> fitLimits{cfgUpcMaxFV0Amplitude, cfgUpcMaxFT0AAmplitude,
                                       cfgUpcMaxFT0CAmplitude, cfgUpcMaxFDDAAmplitude, cfgUpcMaxFDDCAmplitude};
    auto recordFIT = [&](int64_t bcId, MatchingFITRecord record) {
      int64_t gbc = bcs.rawIteratorAt(bcId).globalBC();
      bool clean = udhelpers::cleanFIT(record, cfgUpcFV0ValidTime, fitLimits);
      auto [it, inserted] = cleanFITPerBC.emplace(gbc, clean);
      if (!inserted) {
        it->second = it->second && clean;
      }
    };
    if (matchingEventSelectionMode == 2) {
      for (const auto& fv0 : fv0s) {
        recordFIT(fv0.bcId(), MatchingFITRecord{fv0s, ft0s, fdds, fv0.index(), -1, -1});
        auto gbc = static_cast<int64_t>(bcs.rawIteratorAt(fv0.bcId()).globalBC());
        float sum = udhelpers::FV0AmplitudeA(fv0);
        auto [it, inserted] = fv0Amplitudes.emplace(gbc, sum);
        if (!inserted) {
          it->second = std::max(it->second, sum);
        }
      }
      for (const auto& ft0 : ft0s) {
        recordFIT(ft0.bcId(), MatchingFITRecord{fv0s, ft0s, fdds, -1, ft0.index(), -1});
      }
      for (const auto& fdd : fdds) {
        recordFIT(fdd.bcId(), MatchingFITRecord{fv0s, ft0s, fdds, -1, -1, fdd.index()});
      }
    }
    std::unordered_map<int64_t, std::array<double, 4>> mixingMchDirectionsPV;
    const bool mixingNeedsDirectionPV = cfgEnableEventMixing && cfgMixingDepth != 0 &&
                                        (mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV >= 0.f || cfgMixingMaxDeltaMchEtaPV >= 0.f || mixingMultiplicityCuts.cfgMixingMaxDeltaMchMomentumRel >= 0.f || mixingMultiplicityCuts.cfgMixingMaxDeltaMchPtRel >= 0.);
    std::unordered_map<int64_t, bool> selectedMch;
    int numSelectedMch = 0;
    bool qualityFieldReady = false;
    for (const auto& mch : muontracks) {
      if (static_cast<int>(mch.trackType()) <= 2) {
        continue;
      }
      bool sel8 = mch.has_collision() && collisions.rawIteratorAt(mch.collisionId()).sel8();
      auto bcIt = mchGlobalBCs.find(mch.globalIndex());
      bool hasBC = bcIt != mchGlobalBCs.end();
      int64_t globalBC = hasBC ? bcIt->second : 0;
      int numMchMid = hasBC ? numMchMidPerBC[globalBC] : 0;
      bool passFIT = true;
      bool hasFV0 = false;
      float maxFV0Amplitude = 0.f;
      if (hasBC) {
        for (auto it = cleanFITPerBC.lower_bound(globalBC - cfgUpcFV0BcWindow);
             it != cleanFITPerBC.end() && it->first <= globalBC + cfgUpcFV0BcWindow; ++it) {
          passFIT = passFIT && it->second;
        }
        auto first = fv0Amplitudes.lower_bound(globalBC - cfgUpcFV0BcWindow);
        for (auto it = first; it != fv0Amplitudes.end() && it->first <= globalBC + cfgUpcFV0BcWindow; ++it) {
          hasFV0 = true;
          maxFV0Amplitude = std::max(maxFV0Amplitude, it->second);
        }
      }
      int numMch = mch.has_collision() ? numMchPerCollision[mch.collisionId()] : 0;
      int numMft = mch.has_collision() ? numMftPerCollision[mch.collisionId()] : 0;
      bool passTrackMultiplicity = mch.has_collision() && numMch >= cfgUpcMinMchTracks && numMft >= cfgUpcMinMftTracks &&
                                   (cfgUpcMaxMchTracks < 0 || numMch <= cfgUpcMaxMchTracks) &&
                                   (!cfgUpcRequireEqualTrackCounts || numMft == numMch);
      bool passMchUnlikeSign = mch.has_collision() && numPositiveMchPerCollision[mch.collisionId()] > 0 &&
                               numNegativeMchPerCollision[mch.collisionId()] > 0;
      // Inclusive multiplicities: any number of MCH opposite-sign pairs is allowed.
      // Individual J/psi mass tags are attached downstream per participating MCH.
      bool upc = passTrackMultiplicity && (!cfgUpcRequireUnlikeSign || passMchUnlikeSign) && hasBC && numMchMid >= cfgUpcMinMchMid &&
                 (cfgUpcMaxMchMid < 0 || numMchMid <= cfgUpcMaxMchMid) && passFIT;
      bool selected = matchingEventSelectionMode == 0 || (matchingEventSelectionMode == 1 ? sel8 : upc);
      bool passMuonQuality = !mchQualityCut;
      bool propagatedMuon = false;
      float cutEta = std::numeric_limits<float>::quiet_NaN();
      float cutPt = std::numeric_limits<float>::quiet_NaN();
      // Event rejection precedes PV propagation and geometric acceptance.
      if (selected && (mchQualityCut || mchAcceptance.cfgRequireMchMftAcceptance || mchAcceptance.cfgRequireMchPVEta)) {
        cutEta = mch.eta();
        cutPt = mch.pt();
        bool valid = mch.has_collision();
        bool passMftAcceptance = !mchAcceptance.cfgRequireMchMftAcceptance;
        bool passPVEta = !mchAcceptance.cfgRequireMchPVEta;
        if (mchAcceptance.cfgRequireMchPVEta) {
          mchAcceptanceQA.fill(HIST("pvEtaStatus"), 0.);
        }
        if (mchAcceptance.cfgRequireMchMftAcceptance) {
          mchAcceptanceQA.fill(HIST("status"), 0.);
        }
        if (valid && (mchAcceptance.cfgPropMuon || mchAcceptance.cfgRequireMchMftAcceptance || mchAcceptance.cfgRequireMchPVEta)) {
          if (!qualityFieldReady) {
            initCCDB(bcs.begin());
            o2::mch::TrackExtrap::useExtrapV2();
            qualityFieldReady = true;
          }
          const auto collision = collisions.rawIteratorAt(mch.collisionId());
          auto mchAtPV = globalMatchingPlane.fwdToMch(getTrackParCovFwd(mch, mch));
          valid = std::isfinite(collision.posX()) && std::isfinite(collision.posY()) && std::isfinite(collision.posZ()) &&
                  o2::mch::TrackExtrap::extrapToVertex(mchAtPV, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
          if (valid) {
            const auto atPV = globalMatchingPlane.mchToFwd(mchAtPV);
            // Preserve mchAcceptance.cfgPropMuon semantics for quality cuts independently of
            // the mandatory PV propagation used by geometric acceptance.
            if (mchAcceptance.cfgPropMuon) {
              cutEta = atPV.getEta();
              cutPt = atPV.getPt();
            }
            if (mchAcceptance.cfgRequireMchPVEta) {
              const float etaAtPV = atPV.getEta();
              passPVEta = std::isfinite(etaAtPV) &&
                          etaAtPV > static_cast<float>(mchAcceptance.cfgMchPVEtaMin.value) &&
                          etaAtPV < static_cast<float>(mchAcceptance.cfgMchPVEtaMax.value);
            }
            propagatedMuon = true;
            if (mixingNeedsDirectionPV) {
              mixingMchDirectionsPV[mch.globalIndex()] = {atPV.getPhi(), atPV.getEta(), atPV.getP(), atPV.getPt()};
            }
            if (mchAcceptance.cfgRequireMchMftAcceptance) {
              passMftAcceptance = crossesNominalMftAcceptance(mchAtPV);
            }
          }
        }
        if (mchAcceptance.cfgRequireMchMftAcceptance) {
          mchAcceptanceQA.fill(HIST("status"), !valid ? 1. : passMftAcceptance ? 3.
                                                                               : 2.);
        }
        if (mchAcceptance.cfgRequireMchPVEta) {
          mchAcceptanceQA.fill(HIST("pvEtaStatus"), !valid ? 1. : passPVEta ? 3.
                                                                            : 2.);
        }
        if (mchQualityCut) {
          std::array<float, VarManager::kNVars> values{};
          values[VarManager::kEta] = cutEta;
          values[VarManager::kPt] = cutPt;
          values[VarManager::kMuonRAtAbsorberEnd] = mch.rAtAbsorberEnd();
          values[VarManager::kMuonPDca] = mch.pDca();
          values[VarManager::kMuonChi2] = mch.chi2();
          values[VarManager::kMuonChi2MatchMCHMID] = mch.chi2MatchMCHMID();
          const bool validQuality = valid && std::isfinite(cutEta) && std::isfinite(cutPt) &&
                                    std::isfinite(mch.rAtAbsorberEnd()) && std::isfinite(mch.pDca()) &&
                                    std::isfinite(mch.chi2()) && std::isfinite(mch.chi2MatchMCHMID());
          passMuonQuality = validQuality && mchQualityCut->IsSelected(values.data());
        }
        selected = selected && valid && passMuonQuality && passMftAcceptance && passPVEta;
      }
      selectedMch.emplace(mch.globalIndex(), selected);
      numSelectedMch += selected;
      if (inferenceConfig.cfgTrainingMode)
        tableMatchPlaneEventSels(mch.globalIndex(), mch.collisionId(), matchingEventSelectionMode, selected,
                                 hasBC, static_cast<uint64_t>(globalBC), numMchMid, sel8, hasFV0, maxFV0Amplitude,
                                 numMch, numMft, passTrackMultiplicity, passMchUnlikeSign, passFIT, passMuonQuality, propagatedMuon, cutEta, cutPt, cfgMftAbsorberStepCm.value);
    }
    LOGF(info, "Matching-plane event selection mode=%d: %d/%zu MCH tracks selected", matchingEventSelectionMode, numSelectedMch, selectedMch.size());

    // Rejected DFs still publish selection QA, but perform no geometry/field
    // setup, track conversion or propagation in the early-filter path.
    if (filterBeforeMatching && numSelectedMch == 0) {
      LOGF(info, "Matching-plane early selection: no accepted MCH tracks; skip all pair propagation");
      return;
    }
    initCCDB(bcs.begin());
    o2::mch::TrackExtrap::useExtrapV2();
    const std::array<double, 3> centerMft{0., 0., -61.4};
    globalMatchingPlane.mBzAtMftCenter = fieldB->getBz(centerMft.data());
    poolReset = history.setField(fieldB->getFactorSolenoid(), fieldB->getFactorDipole()) || poolReset;

    // Same last-row-wins lookup as global-muon-matcher, including duplicate cov rows.
    std::unordered_map<int64_t, int64_t> covIndexes;
    for (const auto& cov : mftCovs) {
      covIndexes[cov.matchMFTTrackId()] = cov.globalIndex();
    }
    std::unordered_map<int64_t, o2::track::TrackParCovFwd> mftTrackPars;
    for (const auto& mft : mfttracks) {
      // Apply the MFT acceptance to the original track, before either transport.
      // Both matching planes and the data/MC workflows use this candidate set.
      const auto eta = mft.eta();
      if (!std::isfinite(eta) || !(eta > fEtaMftLow && eta < fEtaMftUp)) {
        continue;
      }
      auto covIt = covIndexes.find(mft.globalIndex());
      if (covIt == covIndexes.end() && cfgRequireMftCovariance) {
        continue;
      }
      if (covIt != covIndexes.end()) {
        auto cov = mftCovs.rawIteratorAt(covIt->second);
        mftTrackPars.emplace(mft.globalIndex(), getTrackParCovFwd(mft, cov));
      } else {
        // Only the fitted central parameters are needed. This zero matrix is
        // container storage, never a measured uncertainty or an input to chi2.
        SMatrix5 parameters(mft.x(), mft.y(), mft.phi(), mft.tgl(), mft.signed1Pt());
        o2::track::SMatrix55Sym zero;
        for (int i = 0; i < 5; ++i)
          for (int j = 0; j <= i; ++j)
            zero(i, j) = 0.;
        mftTrackPars.emplace(mft.globalIndex(), o2::track::TrackParCovFwd(mft.z(), parameters, zero, mft.chi2()));
      }
    }
    // Tagging results are unused here. Only the late-filter/reference path
    // replays these calls to preserve the original navigation history.
    if (!filterBeforeMatching) {
      std::vector<int64_t> taggedMuons;
      globalMatchingPlane.getTaggedMuons(collisions, muontracks, taggedMuons, configMuonTagging);
    }
    using PendingRow = std::tuple<int64_t, int64_t, int64_t, int64_t, float, float, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, bool>;
    std::vector<PendingRow> pendingRows;
    std::vector<PendingRow> downstreamRows;
    std::vector<int> downstreamStatuses;
    size_t numFailedMftAbsorberTransport = 0;
    size_t numFailedJointTransport = 0;
    std::array<size_t, static_cast<size_t>(MftAbsorberExtrap::Status::Count)> mftTransportStatus{};
    std::map<int64_t, std::vector<std::pair<float, int64_t>>> candidatesForOutputPropagation;
    std::vector<PendingRow> mixedPendingRows, mixedDownstreamRows;
    std::vector<int> mixedDownstreamStatuses;
    struct MixingDiagnostic {
      int64_t mch, mft;
      float z;
      bool mixed;
      std::vector<double> mchCov, mftCov, residualCov;
      uint64_t donorFrame;
    };
    std::vector<MixingDiagnostic> mixingDiagnostics;
    auto calculatePair = [&](float planeZ, const auto& mch, const auto& mft, bool mixed,
                             const MftMchJointTransport::Fwd* donorMchPar = nullptr, int64_t donorMchId = -1,
                             const MftMchJointTransport::Fwd* ownedMft = nullptr, uint64_t donorSourceFrame = 0) {
      const bool chiPlaneRequested = inferenceConfig.cfgTrainingMode || !inferenceConfig.cfgEvaluationMode ||
                                     ((inferenceConfig.cfgChi2Method.value.rfind("505", 0) == 0) == (planeZ == -505.f));
      const bool calculateChi2 = cfgComputeMatchingChi2 && chiPlaneRequested && (mixed ? cfgRequireMftCovariance.value : covIndexes.count(mft.globalIndex()) > 0);
      const auto mchTrackPar = getTrackParCovFwd(mch, mch);
      const auto it = mftTrackPars.find(mft.globalIndex());
      if (!ownedMft && it == mftTrackPars.end()) {
        LOGF(fatal, "Missing MFT state in pair calculation");
      }
      const auto& mftInput = ownedMft ? *ownedMft : it->second;
      auto& pairPendingRows = mixed ? mixedPendingRows : pendingRows;
      auto& pairDownstreamRows = mixed ? mixedDownstreamRows : downstreamRows;
      auto& pairDownstreamStatuses = mixed ? mixedDownstreamStatuses : downstreamStatuses;
      // Preserve the original per-pair call order. The geometry navigator
      // and field extrapolator have state; do not cache propagated MCH tracks.
      int transportStatus = 0;
      const auto appendFailedDownstream = [&](int status) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        pairDownstreamRows.emplace_back(mch.globalIndex(), mft.globalIndex(), mch.collisionId(), mft.collisionId(),
                                        -505.f, static_cast<float>(globalMatchingPlane.mBzAtMftCenter),
                                        nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
                                        nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, false);
        pairDownstreamStatuses.push_back(status);
      };
      auto mftProp = mftInput;
      auto mchProp = mchTrackPar;
      std::unique_ptr<MftMchJointTransport::Result> joint;
      if (planeZ <= -505.f && cfgMftUseMchMomentum) {
        joint = std::make_unique<MftMchJointTransport::Result>();
        const bool alignMixed = mixed && cfgMixingAlign505;
        if (alignMixed && !donorMchPar) {
          LOGF(fatal, "Missing donor MCH for alignment");
        }
        bool valid = false;
        if (calculateChi2) {
          valid = MftMchJointTransport::compute(globalMatchingPlane, mftInput, alignMixed ? *donorMchPar : mchTrackPar,
                                                planeZ, cfgMftAbsorberStepCm, cfgMftUseMchCharge, cfgPairJacobianScale, *joint);
        } else {
          // Same nominal MCH-momentum-seeded trajectory, without the numeric
          // joint Jacobian, uncertainty construction or chi2 inversion.
          valid = MftMchJointTransport::evaluate(globalMatchingPlane, mftInput, alignMixed ? *donorMchPar : mchTrackPar,
                                                 planeZ, cfgMftAbsorberStepCm, cfgMftUseMchCharge, 0., joint->state);
          joint->hasNominalState = valid;
          joint->status = valid ? 0 : 1;
        }
        if (alignMixed) {
          auto native = std::move(joint);
          joint = std::make_unique<MftMchJointTransport::Result>();
          MftMchMixedAlignment505::Rotation rotation{1., 0., 0., 0., 1., 0., 0., 0., 1.};
          const auto target = globalMatchingPlane.propagateToZMch(mchTrackPar, planeZ);
          if (native->hasNominalState) {
            if (calculateChi2) {
              valid = MftMchMixedAlignment505::compute(globalMatchingPlane, *native, target, valid,
                                                       cfgPairJacobianScale.value, *joint, rotation);
            } else {
              // Reuse the identical central-value transform without constructing
              // Jacobians or interpreting an absent MFT covariance as measured.
              joint->state = native->state;
              joint->state.mch = target;
              valid = MftMchMixedAlignment505::evaluate(globalMatchingPlane, native->state.mft,
                                                        native->state.mch, target, joint->state.mft, &rotation);
              joint->hasNominalState = valid;
              joint->status = valid ? 0 : 20;
            }
          } else {
            valid = false;
            joint->status = 24;
          }
          const auto vectorOf = [](const auto& track) { std::vector<double> v;for(int i=0;i<5;++i)v.push_back(track.getParameters()[i]);return v; };
          std::vector<double> anchorDelta(5), residualBefore(5);
          for (int i = 0; i < 5; ++i) {
            anchorDelta[i] = native->hasNominalState ? MftMchJointTransport::difference(native->state.mch.getParameters()[i], target.getParameters()[i], i) : std::numeric_limits<double>::quiet_NaN();
            residualBefore[i] = native->hasNominalState ? MftMchJointTransport::difference(native->state.mft.getParameters()[i], target.getParameters()[i], i) : std::numeric_limits<double>::quiet_NaN();
          }
          tableMixedAlignments(mch.globalIndex(), mft.globalIndex(), donorMchId, planeZ,
                               joint->hasNominalState, joint->status, native->status,
                               native->hasNominalState ? vectorOf(native->state.mch) : std::vector<double>{}, vectorOf(target),
                               native->hasNominalState ? vectorOf(native->state.mft) : std::vector<double>{}, anchorDelta, residualBefore,
                               std::vector<double>(rotation.begin(), rotation.end()),
                               native->hasNominalState ? native->state.seed.getP() : std::numeric_limits<double>::quiet_NaN(), frame, donorSourceFrame, mixingPersistence.state.stream);
        }
        if (!valid) {
          joint->chiAll = joint->chiFour = joint->chiXY = std::numeric_limits<double>::quiet_NaN();
          joint->rankAll = joint->rankFour = joint->rankXY = 0;
        }
        if (calculateChi2 && inferenceConfig.cfgTrainingMode && !mixed && selectedMch.at(mch.globalIndex())) {
          tablePairCovariances(mch.globalIndex(), mft.globalIndex(), valid, joint->status,
                               joint->hasNominalState ? joint->state.seed.getP() : std::numeric_limits<double>::quiet_NaN(),
                               mftInput.getInvQPt(), cfgMftUseMchCharge.value,
                               joint->rankAll, joint->rankFour, joint->rankXY,
                               valid ? MftMchJointTransport::flatten(joint->joint) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->residual) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->seedCov) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->seedMchCross) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->jacobian) : std::vector<double>{},
                               joint->state.sigmaLossUp, cfgPairJacobianScale.value);
        }
        transportStatus = joint->status;
        if (!valid) {
          ++numFailedJointTransport;
          // Covariance failure must not discard a successful nominal trajectory.
          // All chi2 values/ranks are invalid together, including partial chi2 success.
          joint->chiAll = joint->chiFour = joint->chiXY = std::numeric_limits<double>::quiet_NaN();
          joint->rankAll = joint->rankFour = joint->rankXY = 0;
          if (!joint->hasNominalState) {
            appendFailedDownstream(transportStatus);
            return;
          }
        }
        mftProp = joint->state.mft;
        mchProp = joint->state.mch;
      } else {

        if (planeZ <= -505.f) {
          // Absorber-exit study: use mapped B and downstream material transport.
          auto mftAsMch = globalMatchingPlane.fwdToMch(mftInput);
          MftAbsorberExtrap::TransportInfo transportInfo;
          if (!MftAbsorberExtrap::extrapToAbsorberEndWithoutBranson(mftAsMch, planeZ, &transportInfo, cfgMftAbsorberStepCm)) {
            ++numFailedMftAbsorberTransport;
            ++mftTransportStatus[static_cast<size_t>(transportInfo.status)];
            appendFailedDownstream(100 + static_cast<int>(transportInfo.status));
            return;
          }
          mftProp = globalMatchingPlane.mchToFwd(mftAsMch);
        } else {
          // In particular z=-77.5: preserve the official fixed-Bz propagation.
          mftProp = globalMatchingPlane.propagateToZMft(mftInput, planeZ);
        }
        mchProp = globalMatchingPlane.propagateToZMch(mchTrackPar, planeZ);
      }
      SMatrix5 delta = mftProp.getParameters() - mchProp.getParameters();
      // Save both raw and wrapped phi residuals; all angular chi2 use the wrapped difference.
      const double wrappedPhi = TVector2::Phi_mpi_pi(delta[2]);
      if (calculateChi2 && cfgSaveMixingDiagnostics && selectedMch.at(mch.globalIndex())) {
        MixingDiagnostic diag{mch.globalIndex(), mft.globalIndex(), planeZ, mixed, {}, {}, {}, mixed ? donorSourceFrame : frame};
        for (int r = 0; r < 5; ++r) {
          for (int c = 0; c < 5; ++c) {
            diag.mchCov.push_back(mchProp.getCovariances()(r, c));
            diag.mftCov.push_back(mftProp.getCovariances()(r, c));
            diag.residualCov.push_back(joint ? joint->residual(r, c) : mchProp.getCovariances()(r, c) + (r == c ? mftProp.getCovariances()(r, c) : 0.));
          }
        }
        mixingDiagnostics.push_back(std::move(diag));
      }
      const double nanChi2 = std::numeric_limits<double>::quiet_NaN();
      const double chi2All = !calculateChi2 ? nanChi2 : joint ? joint->chiAll
                                                              : std::get<0>(globalMatchingPlane.mMatchingFunctionMap.at("matchALL")(mchProp, mftProp));
      const double chi2Angles = !calculateChi2 ? nanChi2 : joint ? joint->chiFour
                                                                 : std::get<0>(globalMatchingPlane.mMatchingFunctionMap.at("matchXYPhiTanl")(mchProp, mftProp));
      const double chi2XY = !calculateChi2 ? nanChi2 : joint ? joint->chiXY
                                                             : std::get<0>(globalMatchingPlane.mMatchingFunctionMap.at("matchXY")(mchProp, mftProp));
      if (calculateChi2 && !mixed && !filterBeforeMatching && planeZ == -77.5f) {
        candidatesForOutputPropagation[mch.globalIndex()].emplace_back(static_cast<float>(chi2All), mft.globalIndex());
      }
      bool finite = !calculateChi2 || (std::isfinite(chi2All) && std::isfinite(chi2Angles) && std::isfinite(chi2XY));
      for (int i = 0; i < 5; ++i) {
        finite = finite && std::isfinite(delta[i]);
      }
      auto& outputRows = planeZ == -77.5f ? pairPendingRows : pairDownstreamRows;
      if (planeZ == -505.f) {
        pairDownstreamStatuses.push_back(transportStatus);
      }
      outputRows.emplace_back(mch.globalIndex(), mft.globalIndex(), mch.collisionId(), mft.collisionId(),
                              planeZ, static_cast<float>(globalMatchingPlane.mBzAtMftCenter),
                              mchProp.getX(), mchProp.getY(), mchProp.getPhi(), mchProp.getTanl(), mchProp.getInvQPt(),
                              mftProp.getX(), mftProp.getY(), mftProp.getPhi(), mftProp.getTanl(), mftProp.getInvQPt(),
                              delta[0], delta[1], delta[2], wrappedPhi, delta[3], delta[4],
                              mftProp.getEta() - mchProp.getEta(), mftProp.getPt() - mchProp.getPt(),
                              chi2All, chi2Angles, chi2XY, finite);
    };
    // Two passes preserve the upstream navigation/call order. Both start from
    // the original AO2D states, never from the other plane's propagated state.
    for (const float planeZ : {-77.5f, -505.f}) {
      for (const auto& entry : globalMchIterationOrder) {
        if (filterBeforeMatching && !selectedMch.at(entry.first)) {
          continue;
        }
        // Retain the relative unordered traversal order of the accepted MCH tracks.
        const auto mch = muontracks.rawIteratorAt(entry.first);
        for (const auto& mft : mfttracks) {
          // Match only tracks assigned to the same reconstructed primary vertex.
          // Timing compatibility remains an additional requirement, not an
          // alternative collision association. Apply this before either transport.
          if (!mch.has_collision() || !mft.has_collision() || mch.collisionId() != mft.collisionId()) {
            continue;
          }
          auto it = mftTrackPars.find(mft.globalIndex());
          if (it == mftTrackPars.end() || !globalMatchingPlane.isMftMchTimeCompatible(collisions, bcs, mch, mft)) {
            continue;
          }
          calculatePair(planeZ, mch, mft, false);
        }
      }
    } // plane passes
    {
      LOGF(info, "MFT/MCH joint transport: failed candidate pairs=%zu, enabled=%d", numFailedJointTransport, static_cast<int>(cfgMftUseMchMomentum.value));
      LOGF(info, "MFT absorber transport: stepCm=%g, failedPairs=%zu [invalid=%zu field=%zu geometry=%zu rangeOut=%zu nonfinite=%zu stepLimit=%zu]", cfgMftAbsorberStepCm.value, numFailedMftAbsorberTransport, mftTransportStatus[1], mftTransportStatus[2], mftTransportStatus[3], mftTransportStatus[4], mftTransportStatus[5], mftTransportStatus[6]);
    }
    // The reference path replays propagation to the MFT point for compatibility.
    // The early-filter path leaves this map empty: those states are never saved.
    for (auto& [mchId, candidates] : candidatesForOutputPropagation) {
      std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
      const auto mch = muontracks.rawIteratorAt(mchId);
      const auto mchPar = getTrackParCovFwd(mch, mch);
      for (const auto& candidate : candidates) {
        const auto& mftPar = mftTrackPars.at(candidate.second);
        auto mchAtMft = globalMatchingPlane.fwdToMch(mchPar);
        o2::mch::TrackExtrap::extrapToVertex(mchAtMft, mftPar.getX(), mftPar.getY(), mftPar.getZ(), mftPar.getSigma2X(), mftPar.getSigma2Y());
      }
    }
    if (cfgEnableEventMixing && cfgMixingDepth != 0) {
      using PoolEntry = MixingPoolEntry;
      struct MixedPair {
        int64_t mch, mft, donorMch, donorCollision;
        int group, groups, attemptsTarget, attemptsDonor;
        mixing::Deltas delta;
        int mftTarget, mftDonor, mchTarget, mchDonor;
        double centTarget, centDonor;
        int cbTarget, cbDonor;
        uint64_t donorFrame;
        const OwnedMixMft* ownedMft;
        const MftMchJointTransport::Fwd* donorPar;
      };
      std::vector<PoolEntry> pool;
      // Current-DF donors are followed by explicitly owned historical events.
      // Use only selected donor events and MFTs associated with that donor collision.
      for (const auto& mch : muontracks) {
        if (static_cast<int>(mch.trackType()) <= 2 || !mch.has_collision() || !selectedMch.at(mch.globalIndex())) {
          continue;
        }
        const auto col = collisions.rawIteratorAt(mch.collisionId());
        const auto bc = bcs.rawIteratorAt(col.bcId());
        PoolEntry entry;
        entry.event = {mch.collisionId(), mch.globalIndex(), static_cast<int64_t>(bc.globalBC()), bc.runNumber(), 0,
                       std::atan2(mch.y(), mch.x()), std::hypot(mch.x(), mch.y()), col.posZ()};
        entry.event.frame = frame;
        entry.mch = getTrackParCovFwd(mch, mch);
        if (mixingNeedsDirectionPV) {
          auto direction = mixingMchDirectionsPV.find(mch.globalIndex());
          if (direction == mixingMchDirectionsPV.end()) {
            // Also use PV directions when quality cuts are disabled or mchAcceptance.cfgPropMuon=false.
            // This propagation is exclusively for the mixing pool; same-event selection is unchanged.
            auto mchAtPV = globalMatchingPlane.fwdToMch(getTrackParCovFwd(mch, mch));
            if (!o2::mch::TrackExtrap::extrapToVertex(mchAtPV, col.posX(), col.posY(), col.posZ(), col.covXX(), col.covYY())) {
              continue;
            }
            const auto atPV = globalMatchingPlane.mchToFwd(mchAtPV);
            direction = mixingMchDirectionsPV.emplace(mch.globalIndex(), std::array<double, 4>{atPV.getPhi(), atPV.getEta(), atPV.getP(), atPV.getPt()}).first;
          }
          entry.event.mchPhiPV = direction->second[0];
          entry.event.mchEtaPV = direction->second[1];
          entry.event.mchMomentumPV = direction->second[2];
          entry.event.mchPtPV = direction->second[3];
          if ((mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV >= 0.f && !std::isfinite(entry.event.mchPhiPV)) ||
              (cfgMixingMaxDeltaMchEtaPV >= 0.f && !std::isfinite(entry.event.mchEtaPV))) {
            continue;
          }
        }
        if (auto count = mixingCbMultiplicities.find(mch.collisionId()); count != mixingCbMultiplicities.end()) {
          entry.event.cbMultiplicity = count->second;
        }
        entry.event.mftMultiplicity = numMftPerCollision[mch.collisionId()];
        entry.event.mchMultiplicity = numMchPerCollision[mch.collisionId()];
        if (auto cent = mixingCentralities.find(mch.collisionId()); cent != mixingCentralities.end()) {
          entry.event.centrality = cent->second;
        }
        for (const auto& mft : mfttracks) {
          if (!mft.has_collision() || !globalMatchingPlane.isMftMchTimeCompatible(collisions, bcs, mch, mft)) {
            continue;
          }
          // Official cfgMatchAllTracks=false definition: all time-compatible MFTs,
          // including those without saved covariance, enter the attempts estimator.
          ++entry.event.attempts;
          if (mft.collisionId() == mch.collisionId() && mftTrackPars.count(mft.globalIndex())) {
            entry.mfts.push_back({mft.globalIndex(), mft.collisionId(), mftTrackPars.at(mft.globalIndex())});
          }
        }
        pool.push_back(std::move(entry));
      }
      mixing::Config config{cfgMixingDepth.value, cfgMinDeltaBc.value, cfgMaxDeltaPhi.value, cfgMaxDeltaR.value,
                            cfgMaxDeltaAttempts.value, cfgMaxDeltaZ.value, cfgMixingSymmetricCuts.value,
                            mixingMultiplicityCuts.cfgMixingMaxDeltaMftMultiplicityRel.value, mixingMultiplicityCuts.cfgMixingMaxDeltaMchMultiplicity.value, cfgMixingMaxDeltaCentrality.value,
                            mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV.value, mixingMultiplicityCuts.cfgMixingMaxDeltaMchMomentumRel.value, cfgMixingMaxDeltaMchEtaPV.value,
                            mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel.value, mixingMultiplicityCuts.cfgMixingMftMultiplicityStrict.value,
                            mixingMultiplicityCuts.cfgMixingMchRadiusStrict.value, mixingMultiplicityCuts.cfgMixingMchEtaStrict.value, mixingMultiplicityCuts.cfgMixingMaxDeltaMchPtRel.value,
                            mixingMultiplicityCuts.cfgMixingCbMultiplicityStrict.value, mixingMultiplicityCuts.cfgMixingVertexZStrict.value};
      std::vector<const PoolEntry*> donors;
      for (const auto& entry : pool)
        donors.push_back(&entry);
      if (mixingPersistence.cfgMixingAcrossDF)
        for (const auto& event : history.events)
          for (const auto& entry : event.anchors)
            donors.push_back(&entry);
      std::vector<MixedPair> mixedPairs;
      for (const auto& target : pool) {
        std::set<std::pair<uint64_t, int64_t>> usedCollisions;
        const size_t begin = mixedPairs.size();
        int group = 0;
        for (const auto* donorPtr : donors) {
          const auto& donor = *donorPtr;
          if (config.depth >= 0 && group >= config.depth) {
            break;
          }
          if (donor.mfts.empty() || usedCollisions.count({donor.event.frame, donor.event.collision})) {
            continue;
          }
          mixing::Deltas delta;
          if (!mixing::compatible(target.event, donor.event, config, delta)) {
            continue;
          }
          usedCollisions.insert({donor.event.frame, donor.event.collision});
          for (const auto& mft : donor.mfts) {
            mixedPairs.push_back({target.event.mch, mft.id, donor.event.mch, donor.event.collision,
                                  group, 0, target.event.attempts, donor.event.attempts, delta,
                                  target.event.mftMultiplicity, donor.event.mftMultiplicity, target.event.mchMultiplicity, donor.event.mchMultiplicity,
                                  target.event.centrality, donor.event.centrality, target.event.cbMultiplicity, donor.event.cbMultiplicity, donor.event.frame, &mft, &donor.mch});
          }
          ++group;
        }
        for (size_t i = begin; i < mixedPairs.size(); ++i) {
          mixedPairs[i].groups = group;
        }
      }
      // Separate passes preserve same-event output exactly, including navigator history.
      // Cross-event timing is deliberately NOT reapplied: donor membership was checked above.
      for (const float planeZ : {-77.5f, -505.f}) {
        for (const auto& pair : mixedPairs) {
          const auto mch = muontracks.rawIteratorAt(pair.mch);
          calculatePair(planeZ, mch, *pair.ownedMft, true,
                        cfgMixingAlign505 && planeZ == -505.f ? pair.donorPar : nullptr,
                        pair.donorMch, &pair.ownedMft->state, pair.donorFrame);
        }
      }
      if (mixedPendingRows.size() != mixedPairs.size() || mixedDownstreamRows.size() != mixedPairs.size() || mixedDownstreamStatuses.size() != mixedPairs.size()) {
        LOGF(fatal, "Mixed dual-plane row alignment failure");
      }
      for (size_t i = 0; i < mixedPairs.size(); ++i) {
        const auto& pair = mixedPairs[i];
        for (int plane = 0; plane < 2; ++plane) {
          const auto& row = plane == 0 ? mixedPendingRows[i] : mixedDownstreamRows[i];
          if (std::get<0>(row) != pair.mch || std::get<1>(row) != pair.mft || (frame == pair.donorFrame && std::get<2>(row) == std::get<3>(row))) {
            LOGF(fatal, "Mixed pair identity failure");
          }
          std::apply([&](const auto&... values) { tableMixedResiduals(values..., cfgMftAbsorberStepCm.value,
                                                                      plane == 0 ? 0 : mixedDownstreamStatuses[i], pair.donorMch, pair.donorCollision, pair.group, pair.groups,
                                                                      pair.delta.bc, pair.delta.phi, pair.delta.radius, pair.delta.attempts, pair.delta.z,
                                                                      pair.attemptsTarget, pair.attemptsDonor, pair.mftTarget, pair.mftDonor, pair.mchTarget, pair.mchDonor, pair.centTarget, pair.centDonor, pair.cbTarget, pair.cbDonor, 1. / pair.groups, frame, pair.donorFrame, mixingPersistence.state.stream); }, row);
        }
      }
      if (cfgSaveMixingDiagnostics) {
        std::set<int64_t> targets;
        for (const auto& pair : mixedPairs) {
          targets.insert(pair.mch);
        }
        for (const auto& diag : mixingDiagnostics) {
          if (targets.count(diag.mch)) {
            tableMixingDiagnostics(diag.mch, diag.mft, diag.z, diag.mixed, diag.mchCov, diag.mftCov, diag.residualCov, frame, diag.donorFrame, mixingPersistence.state.stream);
          }
        }
      }
      const size_t crossPairs = std::count_if(mixedPairs.begin(), mixedPairs.end(), [&](const auto& p) { return p.donorFrame != frame; });
      // Copy/retain only after all pointers used above have gone out of use.
      // One FIFO item per donor collision; all its MCH anchors remain available.
      if (mixingPersistence.cfgMixingAcrossDF) {
        std::map<int64_t, OwnedMixEvent> currentEvents;
        for (auto& entry : pool) {
          if (entry.mfts.empty())
            continue;
          auto& event = currentEvents[entry.event.collision];
          event.frame = frame;
          event.collision = entry.event.collision;
          event.stateCount += entry.mfts.size();
          event.anchors.push_back(std::move(entry));
        }
        for (auto& [collision, event] : currentEvents)
          history.push(std::move(event), mixingPersistence.cfgMixingPoolMaxEvents.value, mixingPersistence.cfgMixingPoolMaxMftStates.value);
      }
      LOGF(info, "Cross-DF mixing: frame=%llu crossPairs=%zu retainedEvents=%zu retainedMftStates=%zu evicted=%zu oversized=%zu", static_cast<unsigned long long>(frame), crossPairs, history.events.size(), history.states, history.evicted, history.oversized);
      LOGF(info, "Event mixing: selected MCH=%zu; mixed pairs=%zu; saved rows=%zu; depth=%d", pool.size(), mixedPairs.size(), 2 * mixedPairs.size(), config.depth);
    }
    if (pendingRows.size() != downstreamRows.size() || downstreamRows.size() != downstreamStatuses.size()) {
      LOGF(fatal, "Dual-plane row alignment failure");
    }
    if (inferenceConfig.cfgEvaluationMode) {
      inference.process(inferenceConfig, origins.size() ? origins.begin().dataframeID() : frame, collisions, bcs, muontracks, mfttracks, pendingRows, downstreamRows, selectedMch,
                        mixingCentralities, mixingCbMultiplicities, globalMatchingPlane, globalMatchingPlane.mBzAtMftCenter, selectedMuonTable);
    }
    if (!inferenceConfig.cfgTrainingMode)
      return;
    for (size_t i = 0; i < downstreamRows.size(); ++i) {
      const auto& row = downstreamRows[i];
      if (std::get<0>(row) != std::get<0>(pendingRows[i]) || std::get<1>(row) != std::get<1>(pendingRows[i])) {
        LOGF(fatal, "Dual-plane pair identity/order mismatch");
      }
      if (selectedMch.at(std::get<0>(row))) {
        std::apply([&](const auto&... values) { tableDownstreamResiduals(values..., cfgMftAbsorberStepCm.value, downstreamStatuses[i]); }, row);
      }
    }
    // Publish only after all propagation, as in global-muon-matcher.
    for (const auto& row : pendingRows) {
      // The early path already rejected unselected MCH tracks before propagation.
      // Keep this check for the late-filter/reference path as well.
      if (selectedMch.at(std::get<0>(row))) {
        std::apply([this](const auto&... values) { tableMatchPlaneResiduals(values...); }, row);
      }
    }
  }
  void processMatchingPlane(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                            aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                            aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    processMatchingPlaneImpl(collisions, bcs, muontracks, mfttracks, mftCovs, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneCentrality(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                      aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                                      aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds,
                                      aod::CentFT0Cs const& centralities, aod::PVMults const& pvMultiplicities, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    if (pvMultiplicities.size() != collisions.size()) {
      LOGF(fatal, "PVMults/collision row count mismatch");
    }
    if (centralities.size() != collisions.size()) {
      LOGF(fatal, "FT0C centrality/collision row count mismatch");
    }
    for (const auto& col : collisions) {
      mixingCentralities.emplace(col.globalIndex(), centralities.rawIteratorAt(col.globalIndex() - collisions.offset()).centFT0C());
      mixingCbMultiplicities.emplace(col.globalIndex(), pvMultiplicities.rawIteratorAt(col.globalIndex() - collisions.offset()).multNTracksPV());
    }
    processMatchingPlaneImpl(collisions, bcs, muontracks, mfttracks, mftCovs, fv0s, ft0s, fdds, origins);
  }
  PROCESS_SWITCH(match_mft_mch_data_mc, processMatchingPlaneCentrality, "Matching with FT0C centrality for event mixing", false);
  PROCESS_SWITCH(match_mft_mch_data_mc, processLegacy, "Legacy residual tables and matching-plane conventions", true);
  PROCESS_SWITCH(match_mft_mch_data_mc, processMatchingPlane, "Global-muon-matcher propagation, timing and chi2; no realignment", false);
};

// Separate device: PV propagation must not change matching-plane navigation caches.
struct MchDimuonMassAtPV {
  Produces<aod::MchDimuonsPV> tableMchDimuonsPV;
  Configurable<bool> cfgSaveMchDimuonMass{"cfgSaveMchDimuonMass", true, "Save selected same-collision MCH dimuons using momenta propagated to the primary vertex"};
  Configurable<int> cfgMassBins{"cfgMassBins", 600, "Dimuon mass histogram bins"};
  Configurable<float> cfgMassMax{"cfgMassMax", 12.f, "Dimuon mass display upper edge in GeV/c^2; higher masses kept in overflow and table"};
  Configurable<std::string> ccdbURL{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"MchDimuonMassPV", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::aod::fwdtrackutils::GlobalMuonMatchingPlane conversion;
  int runNumber = -1;

  void init(InitContext&)
  {
    if (!cfgSaveMchDimuonMass)
      return;
    ccdb->setURL(ccdbURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (cfgMassBins < 1 || !std::isfinite(cfgMassMax.value) || cfgMassMax <= 0.f) {
      LOGF(fatal, "Invalid MCH dimuon mass histogram range or bin count");
    }
    const AxisSpec massAxis{cfgMassBins, 0., cfgMassMax, "m_{#mu#mu} at PV (GeV/c^{2})"};
    histos.add("hMassAll", "Selected MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassUnlikeSign", "Unlike-sign MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassPlusPlus", "Positive like-sign MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassMinusMinus", "Negative like-sign MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassVsPt", "Selected MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); p_{T,#mu#mu} (GeV/c)", HistType::kTH2D, {massAxis, {200, 0., 20.}});
    histos.add("hPVPropagation", "MCH PV propagation; status; Tracks", HistType::kTH1D, {{3, 0., 3.}});
    histos.get<TH1>(HIST("hPVPropagation"))->GetXaxis()->SetBinLabel(1, "Attempted");
    histos.get<TH1>(HIST("hPVPropagation"))->GetXaxis()->SetBinLabel(2, "Succeeded");
    histos.get<TH1>(HIST("hPVPropagation"))->GetXaxis()->SetBinLabel(3, "Failed/nonfinite");
    o2::conf::ConfigurableParam::setValue("FieldOriginBias.z", "0");
  }

  void process(aod::Collisions const& collisions,
               aod::BCsWithTimestamps const& bcs,
               MyMUONs const& tracks,
               aod::MatchPlaneEventSels const& selections)
  {
    if (!cfgSaveMchDimuonMass || bcs.size() == 0) {
      return;
    }
    std::map<int64_t, std::vector<int64_t>> selectedByCollision;
    for (const auto& sel : selections) {
      if (sel.isSelected() && sel.mchCollisionId() >= 0) {
        selectedByCollision[sel.mchCollisionId()].push_back(sel.mchTrackId());
      }
    }
    bool hasPair = std::any_of(selectedByCollision.begin(), selectedByCollision.end(), [](const auto& entry) { return entry.second.size() >= 2; });
    if (!hasPair) {
      return;
    }
    const auto& bc = bcs.begin();
    if (runNumber != bc.runNumber()) {
      auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (!grp) {
        LOGF(fatal, "No magnetic field for MCH PV propagation at BC timestamp");
      }
      o2::base::Propagator::initFieldFromGRP(grp);
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->getForTimeStamp<TGeoManager>("GLO/Config/GeometryAligned", bc.timestamp());
      }
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        LOGF(fatal, "No geometry for MCH PV propagation");
      }
      o2::mch::TrackExtrap::setField();
      o2::mch::TrackExtrap::useExtrapV2();
      runNumber = bc.runNumber();
    }
    const double muonMass = TDatabasePDG::Instance()->GetParticle(13)->Mass();
    struct MuonAtPV {
      int64_t id;
      int charge;
      TLorentzVector momentum;
    };
    for (auto& [collisionId, ids] : selectedByCollision) {
      if (ids.size() < 2) {
        continue;
      }
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      const auto collision = collisions.rawIteratorAt(collisionId);
      std::vector<MuonAtPV> muons;
      for (auto id : ids) {
        const auto track = tracks.rawIteratorAt(id);
        if (static_cast<int>(track.trackType()) <= 2 || !track.has_collision() || track.collisionId() != collisionId) {
          LOGF(fatal, "MCH selection contains an invalid track/collision association");
        }
        histos.fill(HIST("hPVPropagation"), 0.5);
        // Start from the ORIGINAL AO2D MCH state, never from matching-plane states.
        auto mch = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(track, track));
        bool success = o2::mch::TrackExtrap::extrapToVertex(mch, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
        if (!success) {
          histos.fill(HIST("hPVPropagation"), 2.5);
          continue;
        }
        auto pv = conversion.mchToFwd(mch);
        double px = pv.getPx(), py = pv.getPy(), pz = pv.getPz();
        if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(pv.getInvQPt()) || pv.getInvQPt() == 0.) {
          histos.fill(HIST("hPVPropagation"), 2.5);
          continue;
        }
        histos.fill(HIST("hPVPropagation"), 1.5);
        TLorentzVector momentum;
        momentum.SetXYZM(px, py, pz, muonMass);
        muons.push_back({id, pv.getInvQPt() > 0. ? 1 : -1, momentum});
      }
      for (size_t i = 0; i < muons.size(); ++i) {
        for (size_t j = i + 1; j < muons.size(); ++j) {
          const auto& a = muons[i];
          const auto& b = muons[j];
          TLorentzVector pair = a.momentum + b.momentum;
          const double mass = pair.M(), pt = pair.Pt(), rapidity = pair.Rapidity();
          if (!std::isfinite(mass) || !std::isfinite(pt) || !std::isfinite(rapidity)) {
            continue;
          }
          histos.fill(HIST("hMassAll"), mass);
          histos.fill(HIST("hMassVsPt"), mass, pt);
          if (a.charge * b.charge < 0) {
            histos.fill(HIST("hMassUnlikeSign"), mass);
          } else if (a.charge > 0) {
            histos.fill(HIST("hMassPlusPlus"), mass);
          } else {
            histos.fill(HIST("hMassMinusMinus"), mass);
          }
          tableMchDimuonsPV(collisionId, a.id, b.id, static_cast<int8_t>(a.charge + b.charge),
                            collision.posX(), collision.posY(), collision.posZ(),
                            a.momentum.Px(), a.momentum.Py(), a.momentum.Pz(),
                            b.momentum.Px(), b.momentum.Py(), b.momentum.Pz(), mass, pt, rapidity);
        }
      }
    }
  }
};

// ML inference/one-to-one assignment is performed externally. A per-AO2D CSV
// selects exact candidates; no MC labels are consumed by this device.
struct MlMatchedDimuonMassAtPV {
  Produces<aod::MlMatchedMuonsPV> muonTable;
  Produces<aod::MlMatchedDimuonsPV> pairTable;
  Configurable<bool> cfgSaveMlDimuons{"cfgSaveMlDimuons", false, "Build PV dimuons from externally ML-selected one-to-one MFT-MCH pairs"};
  Configurable<std::string> cfgMlSelectionFile{"cfgMlSelectionFile", "", "Per-input CSV: frameId,collisionId,mchTrackId,mftTrackId,goodScore (no MC truth)"};
  Configurable<double> cfgMinMlGoodScore{"cfgMinMlGoodScore", 0., "Minimum supplied good-match score; source score is not a calibrated correct-match probability"};
  Configurable<int> cfgMassBins{"cfgMassBins", 600, "Dimuon mass bins"};
  Configurable<float> cfgMassMax{"cfgMassMax", 12.f, "Mass histogram upper edge; all finite masses retained in table"};
  Configurable<std::string> ccdbURL{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"MlMatchedDimuonPV", {}, OutputObjHandlingPolicy::AnalysisObject};
  struct Choice {
    int64_t collision, mch, mft;
    double score;
  };
  std::map<uint64_t, std::vector<Choice>> choices;
  std::set<uint64_t> processedFrames;
  size_t expected = 0, consumed = 0;
  int runNumber = -1;
  double bz = 0.;
  o2::aod::fwdtrackutils::GlobalMuonMatchingPlane conversion;

  void init(InitContext&)
  {
    if (!cfgSaveMlDimuons) {
      return;
    }
    if (cfgMlSelectionFile.value.empty() || cfgMassBins < 1 || !std::isfinite(cfgMassMax.value) || cfgMassMax <= 0.f ||
        !std::isfinite(cfgMinMlGoodScore.value) || cfgMinMlGoodScore < 0. || cfgMinMlGoodScore > 1.) {
      LOGF(fatal, "Invalid ML dimuon configuration: a per-input selection CSV is required");
    }
    std::ifstream input(cfgMlSelectionFile.value);
    if (!input) {
      LOGF(fatal, "Cannot read ML selection file %s", cfgMlSelectionFile.value.c_str());
    }
    std::string line;
    size_t lineNumber = 0;
    std::map<std::pair<uint64_t, int64_t>, std::set<int64_t>> usedMch, usedMft;
    while (std::getline(input, line)) {
      ++lineNumber;
      if (line.empty() || line[0] == '#') {
        continue;
      }
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream row(line);
      uint64_t frame;
      Choice c;
      std::string extra;
      if (!(row >> frame >> c.collision >> c.mch >> c.mft >> c.score) || (row >> extra) ||
          c.collision < 0 || c.mch < 0 || c.mft < 0 || !std::isfinite(c.score) || c.score < 0. || c.score > 1.) {
        LOGF(fatal, "Malformed ML selection CSV at line %zu", lineNumber);
      }
      if (c.score < cfgMinMlGoodScore) {
        continue;
      }
      const auto event = std::make_pair(frame, c.collision);
      if (!usedMch[event].insert(c.mch).second || !usedMft[event].insert(c.mft).second) {
        LOGF(fatal, "ML dimuons require one-to-one selection within each reconstructed event");
      }
      choices[frame].push_back(c);
      ++expected;
    }
    ccdb->setURL(ccdbURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const AxisSpec massAxis{cfgMassBins, 0., cfgMassMax, "m_{#mu#mu} at PV (GeV/c^{2})"};
    histos.add("hMassAll", "ML-selected matched dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassUnlikeSign", "ML-selected unlike-sign dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassPlusPlus", "ML-selected positive like-sign dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassMinusMinus", "ML-selected negative like-sign dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", HistType::kTH1D, {massAxis});
    histos.add("hMassVsPt", "ML-selected dimuons; m_{#mu#mu} (GeV/c^{2}); p_{T,#mu#mu} (GeV/c)", HistType::kTH2D, {massAxis, {200, 0., 20.}});
    histos.add("hMuonStatus", "ML-selected PV reconstruction; status; Tracks", HistType::kTH1D, {{4, 0., 4.}});
  }

  void process(aod::Origins const& origins, aod::Collisions const& collisions,
               aod::BCsWithTimestamps const& bcs, MyMUONs const& mchTracks,
               aod::MFTTracks const& mftTracks, aod::MFTTracksCov const& mftCovs,
               aod::MatchPlaneResiduals const& candidates)
  {
    if (!cfgSaveMlDimuons) {
      return;
    }
    if (origins.size() == 0) {
      LOGF(fatal, "ML dimuon selection needs the original AO2D frame identity");
    }
    const uint64_t frame = origins.begin().dataframeID();
    auto found = choices.find(frame);
    if (found == choices.end()) {
      return;
    }
    if (!processedFrames.insert(frame).second) {
      LOGF(fatal, "Repeated AO2D frame in ML dimuon input");
    }
    std::set<std::tuple<int64_t, int64_t, int64_t>> available;
    for (const auto& row : candidates) {
      if (row.mchCollisionId() == row.mftCollisionId()) {
        available.emplace(row.mchCollisionId(), row.mchTrackId(), row.mftTrackId());
      }
    }
    std::unordered_map<int64_t, int64_t> covIndex;
    for (const auto& cov : mftCovs) {
      covIndex[cov.matchMFTTrackId()] = cov.globalIndex();
    }
    struct Muon {
      Choice choice;
      int charge;
      TLorentzVector momentum;
    };
    std::map<int64_t, std::vector<Muon>> muons;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double mass = TDatabasePDG::Instance()->GetParticle(13)->Mass();
    for (const auto& c : found->second) {
      if (!available.count({c.collision, c.mch, c.mft})) {
        LOGF(fatal, "ML-selected pair is absent from same-collision matching candidates: frame=%llu MCH=%lld MFT=%lld", static_cast<unsigned long long>(frame), static_cast<long long>(c.mch), static_cast<long long>(c.mft));
      }
      const auto collision = collisions.rawIteratorAt(c.collision - collisions.offset());
      const auto mch = mchTracks.rawIteratorAt(c.mch - mchTracks.offset());
      const auto mft = mftTracks.rawIteratorAt(c.mft - mftTracks.offset());
      if (mch.collisionId() != c.collision || mft.collisionId() != c.collision || mch.trackType() <= 2) {
        LOGF(fatal, "Invalid selected MFT-MCH collision or track type");
      }
      const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
        if (!grp) {
          LOGF(fatal, "Missing magnetic field for ML dimuon reconstruction");
        }
        o2::base::Propagator::initFieldFromGRP(grp);
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
          ccdb->getForTimeStamp<TGeoManager>("GLO/Config/GeometryAligned", bc.timestamp());
        }
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
          LOGF(fatal, "Missing geometry for ML dimuons");
        }
        o2::mch::TrackExtrap::setField();
        o2::mch::TrackExtrap::useExtrapV2();
        const double center[3]{0., 0., -61.4};
        auto* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
        bz = field->getBz(center);
        runNumber = bc.runNumber();
      }
      const bool hasCov = covIndex.count(c.mft);
      auto publishFailure = [&](int status) {
        muonTable(frame, c.collision, c.mch, c.mft, c.score, hasCov, status, static_cast<int8_t>(0), nan, nan, nan);
        histos.fill(HIST("hMuonStatus"), status + 0.5);
      };
      ++consumed;
      auto mchAtPV = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(mch, mch));
      if (!o2::mch::TrackExtrap::extrapToVertex(mchAtPV, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY())) {
        publishFailure(1);
        continue;
      }
      const auto momentumAtPV = conversion.mchToFwd(mchAtPV);
      if (!std::isfinite(momentumAtPV.getInvQPt()) || momentumAtPV.getInvQPt() == 0. || !std::isfinite(momentumAtPV.getTgl())) {
        publishFailure(2);
        continue;
      }
      o2::track::SMatrix55Sym zero;
      for (int i = 0; i < 5; ++i)
        for (int j = 0; j <= i; ++j)
          zero(i, j) = 0.;
      o2::track::SMatrix5 parameters(mft.x(), mft.y(), mft.phi(), mft.tgl(), mft.signed1Pt());
      o2::track::TrackParCovFwd mftState(mft.z(), parameters, zero, mft.chi2());
      if (hasCov) {
        auto cov = mftCovs.rawIteratorAt(covIndex.at(c.mft) - mftCovs.offset());
        mftState = o2::aod::fwdtrackutils::getTrackParCovFwd(mft, cov);
      }
      // PWGEM construction: MCH momentum at PV, MFT position/direction at IU,
      // MCH charge. No absorber-end hybrid state is used. Without MFT covariance,
      // zero is only container storage: only central kinematics are published.
      auto refit = o2::aod::fwdtrackutils::refitGlobalMuonCov(momentumAtPV, mftState);
      const auto pv = o2::aod::fwdtrackutils::propagateTrackParCovFwd(refit, static_cast<uint8_t>(0), collision,
                                                                      o2::aod::fwdtrackutils::propagationPoint::kToVertex, -77.5f, bz);
      const double px = pv.getPx(), py = pv.getPy(), pz = pv.getPz();
      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(pv.getInvQPt()) || pv.getInvQPt() == 0.) {
        publishFailure(3);
        continue;
      }
      const int charge = pv.getInvQPt() > 0. ? 1 : -1;
      TLorentzVector momentum;
      momentum.SetXYZM(px, py, pz, mass);
      muonTable(frame, c.collision, c.mch, c.mft, c.score, hasCov, 0, static_cast<int8_t>(charge), px, py, pz);
      histos.fill(HIST("hMuonStatus"), 0.5);
      muons[c.collision].push_back({c, charge, momentum});
    }
    for (const auto& [collisionId, tracks] : muons) {
      for (size_t i = 0; i < tracks.size(); ++i)
        for (size_t j = i + 1; j < tracks.size(); ++j) {
          const auto& a = tracks[i];
          const auto& b = tracks[j];
          const auto pair = a.momentum + b.momentum;
          if (!std::isfinite(pair.M()) || !std::isfinite(pair.Pt()) || !std::isfinite(pair.Rapidity())) {
            continue;
          }
          pairTable(frame, collisionId, a.choice.mch, a.choice.mft, b.choice.mch, b.choice.mft,
                    a.choice.score, b.choice.score, static_cast<int8_t>(a.charge + b.charge),
                    a.momentum.Px(), a.momentum.Py(), a.momentum.Pz(), b.momentum.Px(), b.momentum.Py(), b.momentum.Pz(),
                    pair.M(), pair.Pt(), pair.Rapidity());
          histos.fill(HIST("hMassAll"), pair.M());
          histos.fill(HIST("hMassVsPt"), pair.M(), pair.Pt());
          if (a.charge * b.charge < 0)
            histos.fill(HIST("hMassUnlikeSign"), pair.M());
          else if (a.charge > 0)
            histos.fill(HIST("hMassPlusPlus"), pair.M());
          else
            histos.fill(HIST("hMassMinusMinus"), pair.M());
        }
    }
  }
  void endOfStream(EndOfStreamContext&)
  {
    if (cfgSaveMlDimuons && consumed != expected) {
      LOGF(fatal, "ML selection/input mismatch: consumed %zu of %zu selected muons", consumed, expected);
    }
  }
};

// Join PV mass tags to residual rows downstream, without feeding PV propagation
// back into the matching device. Each candidate inherits the tag of its MCH pair.
struct MchJpsiMatchFlags {
  Produces<aod::MatchPlaneResiduals> residuals;
  Produces<aod::JpsiMatchPlaneResiduals> jpsiResiduals;
  Configurable<double> cfgJpsiMassMin{"cfgJpsiMassMin", 2.8, "Inclusive lower bound on OS MCH-pair PV invariant mass in GeV/c^2"};
  Configurable<double> cfgJpsiMassMax{"cfgJpsiMassMax", 3.3, "Inclusive upper bound on OS MCH-pair PV invariant mass in GeV/c^2"};
  void init(InitContext&)
  {
    if (!std::isfinite(cfgJpsiMassMin.value) || !std::isfinite(cfgJpsiMassMax.value) || cfgJpsiMassMin < 0. || cfgJpsiMassMax <= cfgJpsiMassMin) {
      LOGF(fatal, "Invalid J/psi PV mass window");
    }
  }
  void process(aod::RawMatchPlaneResiduals const& raw,
               aod::MatchPlaneEventSels const& selections,
               aod::MchDimuonsPV const& pairs)
  {
    std::map<int64_t, std::set<int64_t>> eligible;
    for (const auto& sel : selections) {
      if (sel.isSelected() && sel.mchCollisionId() >= 0) {
        eligible[sel.mchCollisionId()].insert(sel.mchTrackId());
      }
    }
    struct TrackTag {
      double mass = std::numeric_limits<double>::quiet_NaN();
      int64_t partner = -1;
      bool representativeInWindow = false;
      std::set<int64_t> jpsiPartners;
    };
    // The key contains the MCH ID: a signal elsewhere in the collision must not
    // tag unrelated MCH tracks. Multiple qualifying partners do not duplicate rows.
    std::map<std::pair<int64_t, int64_t>, TrackTag> tags;
    const double targetMass = TDatabasePDG::Instance()->GetParticle(443)->Mass();
    auto addPartner = [&](int64_t collision, int64_t mch, int64_t partner, double mass) {
      auto& tag = tags[{collision, mch}];
      bool inWindow = mass >= cfgJpsiMassMin && mass <= cfgJpsiMassMax;
      if (inWindow) {
        tag.jpsiPartners.insert(partner);
      }
      // One scalar mass per MFT-MCH row: prioritize an in-window partner, then
      // the mass closest to PDG J/psi. Use partner ID to break exact ties.
      if (tag.partner < 0 || (inWindow && !tag.representativeInWindow) ||
          (inWindow == tag.representativeInWindow &&
           (std::abs(mass - targetMass) < std::abs(tag.mass - targetMass) ||
            (std::abs(mass - targetMass) == std::abs(tag.mass - targetMass) && partner < tag.partner)))) {
        tag.mass = mass;
        tag.partner = partner;
        tag.representativeInWindow = inWindow;
      }
    };
    for (const auto& pair : pairs) {
      auto event = eligible.find(pair.collisionId());
      if (event == eligible.end() || pair.pairCharge() != 0 || !std::isfinite(pair.mass()) ||
          !event->second.count(pair.mchTrackId1()) || !event->second.count(pair.mchTrackId2())) {
        continue;
      }
      addPartner(pair.collisionId(), pair.mchTrackId1(), pair.mchTrackId2(), pair.mass());
      addPartner(pair.collisionId(), pair.mchTrackId2(), pair.mchTrackId1(), pair.mass());
    }
    for (const auto& row : raw) {
      auto found = tags.find({row.mchCollisionId(), row.mchTrackId()});
      bool valid = found != tags.end();
      double mass = valid ? found->second.mass : std::numeric_limits<double>::quiet_NaN();
      bool isJpsi = valid && !found->second.jpsiPartners.empty();
      int64_t partner = valid ? found->second.partner : -1;
      int numPartners = valid ? static_cast<int>(found->second.jpsiPartners.size()) : 0;
      residuals(row.mchTrackId(), row.mftTrackId(), row.mchCollisionId(), row.mftCollisionId(), row.matchingPlaneZ(), row.bzAtMftCenter(), row.mchX(), row.mchY(), row.mchPhi(), row.mchTanl(), row.mchInvQPt(), row.mftX(), row.mftY(), row.mftPhi(), row.mftTanl(), row.mftInvQPt(), row.deltaX(), row.deltaY(), row.deltaPhi(), row.deltaPhiWrapped(), row.deltaTanl(), row.deltaInvQPt(), row.deltaEta(), row.deltaPt(), row.chi2All(), row.chi2XYPhiTanl(), row.chi2XY(), row.isFinite(), mass, valid, isJpsi, partner, numPartners);
      if (isJpsi) {
        jpsiResiduals(row.mchTrackId(), row.mftTrackId(), row.mchCollisionId(), row.mftCollisionId(), row.matchingPlaneZ(), row.bzAtMftCenter(), row.mchX(), row.mchY(), row.mchPhi(), row.mchTanl(), row.mchInvQPt(), row.mftX(), row.mftY(), row.mftPhi(), row.mftTanl(), row.mftInvQPt(), row.deltaX(), row.deltaY(), row.deltaPhi(), row.deltaPhiWrapped(), row.deltaTanl(), row.deltaInvQPt(), row.deltaEta(), row.deltaPt(), row.chi2All(), row.chi2XYPhiTanl(), row.chi2XY(), row.isFinite(), mass, valid, isJpsi, partner, numPartners);
      }
    }
  }
};

struct MixedMatchPlanePublisher {
  Produces<aod::MixedMatchPlaneResiduals> residuals;
  void process(aod::RawMixedMatchPlaneResiduals const& rows)
  {
    for (const auto& row : rows) {
      copyMixedResidual(row, residuals);
    }
  }
};

inline WorkflowSpec makeMatchingWorkflow(ConfigContext const& cfgc, bool modernDefaults = false)
{
  WorkflowSpec workflow{adaptAnalysisTask<match_mft_mch_data_mc>(cfgc, TaskName{"match_mft_mch_data_mc"},
                                                                 modernDefaults ? SetDefaultProcesses{{{"processLegacy", false}, {"processMatchingPlane", true}}} : SetDefaultProcesses{}),
                        adaptAnalysisTask<MchDimuonMassAtPV>(cfgc, TaskName{"mch-dimuon-mass-at-p-v"}),
                        adaptAnalysisTask<MlMatchedDimuonMassAtPV>(cfgc, TaskName{"ml-matched-dimuon-mass-at-p-v"}),
                        adaptAnalysisTask<MchJpsiMatchFlags>(cfgc, TaskName{"mch-jpsi-match-flags"})};
  if (!modernDefaults) {
    workflow.push_back(adaptAnalysisTask<MixedMatchPlanePublisher>(cfgc, TaskName{"mixed-match-plane-publisher"}));
  }
  return workflow;
}

} // namespace o2::analysis::mftmch

// -----------------------------------------------------------------------------
// MftMchMatchAncestry
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cstdint>
#include <deque>
#include <map>
#include <set>
#include <tuple>
#include <vector>
namespace o2::analysis::mftmch::truth
{
struct Node {
  int64_t id = -1;
  int pdg = 0, process = -1;
  bool valid = false, generator = false;
  float vx = 0.f, vy = 0.f, vz = 0.f;
  std::vector<int64_t> mothers;
};
// Completeness refers to the genealogy stored in this AO2D, not omitted MC tracks.
enum SearchStatus : int { Complete = 0,
                          MissingMother = 1,
                          DepthLimit = 2,
                          Cycle = 4,
                          NodeLimit = 8,
                          MissingLabel = 16 };
struct Path {
  std::vector<Node> childToAncestor;
  bool decayOnly = false;
};
struct Search {
  std::map<int64_t, Path> ancestors;
  int status = Complete;
};
// Search all mother branches. A generator genealogy is accepted as decay-like,
// while every transport-produced child must explicitly have the decay process.
// At least one PDG change is required: same-PDG bookkeeping copies are not decays.
template <typename Getter>
Search findAncestors(const Node& descendant, Getter get, int decayProcess, int maxDepth, int maxNodes)
{
  Search out;
  struct State {
    std::vector<Node> path;
    bool clean = true, changed = false;
  };
  std::deque<State> queue;
  queue.push_back({{descendant}, true, false});
  std::set<std::tuple<int64_t, bool, bool>> visited;
  visited.emplace(descendant.id, true, false);
  int expanded = 0;
  while (!queue.empty()) {
    State state = std::move(queue.front());
    queue.pop_front();
    const Node child = state.path.back();
    bool hasMother = std::any_of(child.mothers.begin(), child.mothers.end(), [](int64_t id) { return id != -1; });
    if (static_cast<int>(state.path.size()) - 1 >= maxDepth) {
      if (hasMother) {
        out.status |= DepthLimit;
      }
      continue;
    }
    if (++expanded > maxNodes) {
      out.status |= NodeLimit;
      break;
    }
    auto ids = child.mothers;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    for (int64_t id : ids) {
      if (id == -1) {
        continue;
      }
      if (id < 0) {
        out.status |= MissingMother;
        continue;
      }
      if (std::any_of(state.path.begin(), state.path.end(), [&](const Node& n) { return n.id == id; })) {
        out.status |= Cycle;
        continue;
      }
      Node parent = get(id);
      if (!parent.valid) {
        out.status |= MissingMother;
        continue;
      }
      const bool clean = state.clean && (child.generator || child.process == decayProcess);
      const bool changed = state.changed || parent.pdg != child.pdg;
      State next = state;
      next.path.push_back(parent);
      next.clean = clean;
      next.changed = changed;
      const bool decayOnly = clean && changed;
      auto found = out.ancestors.find(id);
      if (found == out.ancestors.end() || (decayOnly && !found->second.decayOnly) ||
          (decayOnly == found->second.decayOnly && next.path.size() < found->second.childToAncestor.size())) {
        out.ancestors[id] = Path{next.path, decayOnly};
      }
      if (visited.emplace(id, clean, changed).second) {
        queue.push_back(std::move(next));
      }
    }
  }
  return out;
}
struct Relation {
  // -1=unknown, 0=no MFT->MCH-muon ancestry, 1=same particle,
  // 2=MFT ancestor -> MCH muon via decay-only path, 3=other ancestor path.
  int origin = -1, status = Complete, depth = -1;
  bool ancestor = false, decay = false;
  std::vector<Node> parentToChild;
};
inline Relation classify(const Node& mch, const Node& mft, const Search* search)
{
  Relation out;
  if (!mch.valid || !mft.valid) {
    out.status = MissingLabel;
    return out;
  }
  if (mch.id == mft.id) {
    out.origin = 1;
    out.depth = 0;
    out.parentToChild = {mch};
    return out;
  }
  if (mch.pdg != 13 && mch.pdg != -13) {
    out.origin = 0;
    return out;
  }
  if (!search) {
    out.status = MissingMother;
    return out;
  }
  out.status = search->status;
  auto found = search->ancestors.find(mft.id);
  if (found == search->ancestors.end()) {
    out.origin = out.status == Complete ? 0 : -1;
    return out;
  }
  out.ancestor = true;
  out.decay = found->second.decayOnly;
  out.origin = out.decay ? 2 : 3;
  out.parentToChild = found->second.childToAncestor;
  std::reverse(out.parentToChild.begin(), out.parentToChild.end());
  out.depth = static_cast<int>(out.parentToChild.size()) - 1;
  return out;
}
} // namespace o2::analysis::mftmch::truth

#endif // COMMON_CORE_MFTMCHMATCHING_H_
