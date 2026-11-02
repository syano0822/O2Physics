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
/// \file MftMchMixingAlignment.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_MFTMCHMIXINGALIGNMENT_H_
#define COMMON_CORE_MFTMCHMIXINGALIGNMENT_H_
#include <algorithm>
#include <array>
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
#endif
