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
/// \file GlobalMuonMatchingPlane.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_GLOBALMUONMATCHINGPLANE_H_
#define COMMON_CORE_GLOBALMUONMATCHINGPLANE_H_
#include "Common/Core/fwdtrackUtilities.h"

#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackParam.h>
#include <MFTTracking/Constants.h>
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
#endif // COMMON_CORE_GLOBALMUONMATCHINGPLANE_H_
