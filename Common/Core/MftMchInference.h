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
/// \file MftMchInference.h
/// \brief Shared MFT–MCH matching study utilities and output definitions.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_MFTMCHINFERENCE_H_
#define COMMON_CORE_MFTMCHINFERENCE_H_
#include "Common/Core/fwdtrackUtilities.h"
#include "Tools/ML/model.h"

#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>

#include <THnSparse.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <deque>
#include <fstream>
#include <numeric>
#include <optional>
#include <set>
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
  o2::framework::Configurable<std::string> cfgMlModelDirectory{"cfgMlModelDirectory", ".", "Directory containing NN/XGBoost .onnx and .json contracts"};
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
#endif
