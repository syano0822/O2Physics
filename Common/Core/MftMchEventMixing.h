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

/// \file MftMchEventMixing.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_MFTMCHEVENTMIXING_H_
#define COMMON_CORE_MFTMCHEVENTMIXING_H_
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
#endif
