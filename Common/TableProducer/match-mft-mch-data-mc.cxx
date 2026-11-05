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

#include "Common/Core/MftMchMatching.h"
#include "Common/DataModel/MftMchMatchTruth.h"
#include <TMCProcess.h>
#include <unordered_map>
#include <vector>
#include <Framework/runDataProcessing.h>
#include <limits>

using namespace o2;
using namespace o2::framework;

// This device only annotates final residual rows. No matching/propagation here.
struct MftMchMcTruth {
  Produces<aod::MftMchMatchTruths> truths;
  Produces<aod::MixedMatchPlaneResiduals> mixedResiduals;
  Produces<aod::MixedMftMchMatchTruths> mixedTruths;
  Produces<aod::MftMchMatchRelations> relations;
  Configurable<int> cfgMcAncestryMaxDepth{"cfgMcAncestryMaxDepth", 128, "Maximum stored MC mother-graph depth; truncation is explicitly flagged"};
  Configurable<int> cfgMcAncestryMaxNodes{"cfgMcAncestryMaxNodes", 4096, "Maximum mother-graph states expanded per distinct MCH MC particle"};
  void init(InitContext&)
  {
    if (cfgMcAncestryMaxDepth < 1 || cfgMcAncestryMaxNodes < 1) {
      LOGF(fatal, "MC ancestry limits must be positive");
    }
  }
  using MchLabels = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
  using MftLabels = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  struct ParticleInfo {
    int64_t id = -1, collision = -1;
    int pdg = 0;
    bool valid = false, primary = false;
    float px = std::numeric_limits<float>::quiet_NaN();
    float py = std::numeric_limits<float>::quiet_NaN();
    float pz = std::numeric_limits<float>::quiet_NaN();
  };
  static ParticleInfo particle(int64_t id, aod::McParticles const& particles)
  {
    ParticleInfo out;
    out.id = id;
    const auto offset = static_cast<int64_t>(particles.offset());
    if (id < offset || id - offset >= particles.size()) {
      return out;
    }
    const auto p = particles.rawIteratorAt(id - offset);
    out.valid = true;
    out.collision = p.mcCollisionId();
    out.pdg = p.pdgCode();
    out.px = p.px();
    out.py = p.py();
    out.pz = p.pz();
    out.primary = p.isPhysicalPrimary();
    return out;
  }
  void process(aod::MatchPlaneResiduals const& residuals, MchLabels const& mchTracks,
               MftLabels const& mftTracks, aod::McParticles const& particles)
  {
    using namespace o2::analysis::mftmch::truth;
    // Local caches must not cross DF boundaries: MC IDs are DF scoped.
    std::unordered_map<int64_t, Node> nodes;
    std::unordered_map<int64_t, Search> searches;
    auto node = [&](int64_t id) -> Node {
      auto cached = nodes.find(id);
      if (cached != nodes.end()) {
        return cached->second;
      }
      Node n;
      n.id = id;
      const auto offset = static_cast<int64_t>(particles.offset());
      if (id >= offset && id - offset < particles.size()) {
        const auto p = particles.rawIteratorAt(id - offset);
        n.valid = true;
        n.pdg = p.pdgCode();
        n.process = p.getProcess();
        n.generator = p.producedByGenerator();
        n.vx = p.vx();
        n.vy = p.vy();
        n.vz = p.vz();
        for (auto mother : p.mothersIds()) {
          n.mothers.push_back(mother);
        }
      }
      nodes.emplace(id, n);
      return n;
    };
    auto motherPdgs = [&](const Node& n) {
      std::vector<int> pdgs;
      for (auto id : n.mothers) {
        auto p = node(id);
        pdgs.push_back(p.valid ? p.pdg : 0);
      }
      return pdgs;
    };
    int64_t residualRow = residuals.offset();
    for (const auto& row : residuals) {
      const auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      const auto mft = mftTracks.rawIteratorAt(row.mftTrackId() - mftTracks.offset());
      const auto a = particle(mch.mcParticleId(), particles);
      const auto b = particle(mft.mcParticleId(), particles);
      // Identity criterion of the previous MC task. Masks are saved separately
      // so fake/noise label policy can be studied without altering matching.
      const bool labelsValid = a.valid && b.valid;
      const bool correct = labelsValid && a.id == b.id;
      const bool falseMatch = labelsValid && a.id != b.id;
      const int8_t status = labelsValid ? (correct ? 1 : 0) : -1;
      truths(residualRow, row.mchTrackId(), row.mftTrackId(), a.id, b.id,
             a.valid, b.valid, correct, mch.mcMask(), mft.mcMask(),
             a.collision, b.collision, a.pdg, b.pdg,
             a.px, a.py, a.pz, b.px, b.py, b.pz, a.primary, b.primary, falseMatch, status);
      const Node mchNode = node(a.id), mftNode = node(b.id);
      const Search* search = nullptr;
      if (labelsValid && !correct && (a.pdg == 13 || a.pdg == -13)) {
        auto cached = searches.find(a.id);
        if (cached == searches.end()) {
          cached = searches.emplace(a.id, findAncestors(mchNode, node, static_cast<int>(TMCProcess::kPDecay), cfgMcAncestryMaxDepth, cfgMcAncestryMaxNodes)).first;
        }
        search = &cached->second;
      }
      const auto relation = classify(mchNode, mftNode, search);
      std::vector<int64_t> pathIds;
      std::vector<int> pathPdgs, pathProcesses;
      std::vector<uint8_t> pathGenerator;
      std::vector<float> pathX, pathY, pathZ;
      for (const auto& n : relation.parentToChild) {
        pathIds.push_back(n.id);
        pathPdgs.push_back(n.pdg);
        pathProcesses.push_back(n.process);
        pathGenerator.push_back(n.generator);
        pathX.push_back(n.vx);
        pathY.push_back(n.vy);
        pathZ.push_back(n.vz);
      }
      const float nan = std::numeric_limits<float>::quiet_NaN();
      relations(residualRow++, row.mchTrackId(), row.mftTrackId(),
                static_cast<int8_t>(relation.origin), relation.status, relation.status == Complete,
                a.valid && (a.pdg == 13 || a.pdg == -13), relation.ancestor, relation.decay,
                relation.decay && relation.depth == 1, relation.depth,
                mchNode.process, mftNode.process, mchNode.generator, mftNode.generator,
                mchNode.valid ? mchNode.vx : nan, mchNode.valid ? mchNode.vy : nan, mchNode.valid ? mchNode.vz : nan,
                mftNode.valid ? mftNode.vx : nan, mftNode.valid ? mftNode.vy : nan, mftNode.valid ? mftNode.vz : nan,
                mchNode.mothers, motherPdgs(mchNode), mftNode.mothers, motherPdgs(mftNode),
                pathIds, pathPdgs, pathProcesses, pathGenerator, pathX, pathY, pathZ);
    }
  }
  // Only MC-independent matching is done upstream. This publication guard makes
  // the MC mixed sample strictly different generator collisions and particles.
  // Labels are never forced false to hide a true/unknown pair.
  struct CachedMftTruth {
    ParticleInfo info;
    uint16_t mask = 0;
  };
  // MC particle/collision IDs are only comparable within the same source frame.
  std::map<uint64_t, std::unordered_map<int64_t, CachedMftTruth>> mixedLabelFrames;
  void processMixed(aod::RawMixedMatchPlaneResiduals const& rows, MchLabels const& mchTracks,
                    MftLabels const& mftTracks, aod::McParticles const& particles,
                    aod::MftMchMixingPoolStates const& poolStates)
  {
    if (poolStates.size() == 0) {
      if (rows.size() != 0)
        LOGF(fatal, "Missing mixed-event source-frame metadata");
      mixedLabelFrames.clear();
      return;
    }
    if (poolStates.size() != 1)
      LOGF(fatal, "Expected one pool state per DF");
    const auto state = poolStates.begin();
    const auto frame = state.targetFrameId();
    auto& current = mixedLabelFrames[frame];
    current.clear();
    for (const auto& mft : mftTracks)
      current.emplace(mft.globalIndex(), CachedMftTruth{particle(mft.mcParticleId(), particles), mft.mcMask()});
    auto donorTruth = [&](const auto& row) -> const CachedMftTruth& {
      if (row.targetFrameId() != frame || row.mixingStreamId() != state.mixingStreamId())
        LOGF(fatal, "Mixed source frame/stream mismatch");
      const auto source = mixedLabelFrames.find(row.donorFrameId());
      if (source == mixedLabelFrames.end())
        LOGF(fatal, "Historical MFT label frame is unavailable");
      const auto track = source->second.find(row.mftTrackId());
      if (track == source->second.end())
        LOGF(fatal, "Historical MFT label ID is unavailable");
      return track->second;
    };
    std::set<int64_t> eligible;
    std::map<int64_t, std::set<int>> groups;
    size_t rejectedUnknown = 0, rejectedSameMcEvent = 0;
    int64_t sourceIndex = 0;
    for (const auto& row : rows) {
      const auto sourceRow = sourceIndex++;
      const auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      const auto a = particle(mch.mcParticleId(), particles);
      const auto& b = donorTruth(row).info;
      if (!a.valid || !b.valid || a.collision < 0 || b.collision < 0) {
        ++rejectedUnknown;
        continue;
      }
      if (row.targetFrameId() == row.donorFrameId() && (a.id == b.id || a.collision == b.collision)) {
        ++rejectedSameMcEvent;
        continue;
      }
      eligible.insert(sourceRow);
      groups[row.mchTrackId()].insert(row.mixedGroupIndex());
    }
    int64_t index = 0;
    sourceIndex = 0;
    for (const auto& row : rows) {
      if (!eligible.count(sourceIndex++)) {
        continue;
      }
      const auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      const auto a = particle(mch.mcParticleId(), particles);
      const auto& donor = donorTruth(row);
      const auto& b = donor.info;
      o2::analysis::mftmch::copyMixedResidual(row, mixedResiduals, groups.at(row.mchTrackId()).size());
      mixedTruths(index++, row.mchTrackId(), row.mftTrackId(), a.id, b.id,
                  true, true, false, mch.mcMask(), donor.mask, a.collision, b.collision, a.pdg, b.pdg,
                  a.px, a.py, a.pz, b.px, b.py, b.pz, a.primary, b.primary, true, static_cast<int8_t>(0));
    }
    // Evicted upstream frames may still have been used in the current rows.
    // Prune only after publication, against the exact retained-frame set.
    const auto retained = state.retainedFrames();
    std::set<uint64_t> keep(retained.begin(), retained.end());
    for (auto it = mixedLabelFrames.begin(); it != mixedLabelFrames.end();) {
      if (!keep.count(it->first))
        it = mixedLabelFrames.erase(it);
      else
        ++it;
    }
    LOGF(info, "Mixed MC publication: saved=%lld, unknown rejected=%zu, same-MC-event/particle rejected=%zu, retained label frames=%zu",
         static_cast<long long>(index), rejectedUnknown, rejectedSameMcEvent, mixedLabelFrames.size());
  }
  PROCESS_SWITCH(MftMchMcTruth, processMixed, "Publish mixed rows only for valid, distinct MC collisions/particles", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = o2::analysis::mftmch::makeMatchingWorkflow(cfgc, true);
  workflow.push_back(adaptAnalysisTask<MftMchMcTruth>(cfgc, TaskName{"mft-mch-mc-truth"}));
  return workflow;
}
