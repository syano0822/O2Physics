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

/// \file MftMchMixingHistory.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_MFTMCHMIXINGHISTORY_H_
#define COMMON_CORE_MFTMCHMIXINGHISTORY_H_
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
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
#endif
