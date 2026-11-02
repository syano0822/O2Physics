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

/// \file MftMchMatchAncestry.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_MFTMCHMATCHANCESTRY_H_
#define COMMON_CORE_MFTMCHMATCHANCESTRY_H_
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
#endif
