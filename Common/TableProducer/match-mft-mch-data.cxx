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
#include <map>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/MftmchMatchingML.h"
#include "Common/Core/trackUtilities.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DetectorsVertexing/VertexTrackMatcher.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "MathUtils/Primitive2D.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/MatchMFTFT0.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "EventFiltering/Zorro.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

float mMu = TDatabasePDG::Instance()->GetParticle(13)->Mass();

TRandom* rnd = new TRandom();

TLorentzVector muon1LV;
TLorentzVector muon2LV;
TLorentzVector dimuonLV;

TVector3 V1;
TVector3 V2;

namespace o2::aod
{

namespace muon_params
{
DECLARE_SOA_COLUMN(TRACKCHI2, trackChi2, float);
DECLARE_SOA_COLUMN(RABS, rabs, float);
DECLARE_SOA_COLUMN(Q, q, int16_t);
DECLARE_SOA_COLUMN(PT, pt, float);
DECLARE_SOA_COLUMN(ETA, eta, float);
DECLARE_SOA_COLUMN(PHI, phi, float);
DECLARE_SOA_COLUMN(HASMFT, has_mft, bool);
} // namespace muon_params

DECLARE_SOA_TABLE(MUONParams, "AOD", "MUON",
                  muon_params::TRACKCHI2,
                  muon_params::RABS,
                  muon_params::Q,
                  muon_params::PT,
                  muon_params::ETA,
                  muon_params::PHI,
                  muon_params::HASMFT);

namespace mft_params
{
DECLARE_SOA_COLUMN(NCLUST, nclust, int);
DECLARE_SOA_COLUMN(ISCA, isCA, bool);
DECLARE_SOA_COLUMN(TRACKCHI2, trackChi2, float);
DECLARE_SOA_COLUMN(Q, q, int16_t);
DECLARE_SOA_COLUMN(PT_AT_DCA, pt_dca, float);
DECLARE_SOA_COLUMN(ETA_AT_DCA, eta_dca, float);
DECLARE_SOA_COLUMN(PHI_AT_DCA, phi_dca, float);
DECLARE_SOA_COLUMN(DCAx, dcax, float);
DECLARE_SOA_COLUMN(DCAy, dcay, float);
} // namespace mft_params

DECLARE_SOA_TABLE(MFTParams, "AOD", "MFT",
                  mft_params::NCLUST,
                  mft_params::ISCA,
                  mft_params::TRACKCHI2,
                  mft_params::Q,
                  mft_params::PT_AT_DCA,
                  mft_params::ETA_AT_DCA,
                  mft_params::PHI_AT_DCA,
                  mft_params::DCAx,
                  mft_params::DCAy);

namespace matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);
DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchQ, mchq, int16_t);
DECLARE_SOA_COLUMN(MftQ, mftq, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, is_correct, bool);

} // namespace matching_params

DECLARE_SOA_TABLE(MatchParams, "AOD", "MATCHING",
                  matching_params::DeltaPt,
                  matching_params::DeltaEta,
                  matching_params::DeltaPhi,
                  matching_params::DeltaX,
                  matching_params::DeltaY,
                  matching_params::MchPt,
                  matching_params::MchEta,
                  matching_params::MchQ,
                  matching_params::MftQ,
                  matching_params::IsCorrectMatch);

namespace mix_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);
DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchQ, mchq, int16_t);
DECLARE_SOA_COLUMN(MftQ, mftq, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, is_correct, bool);
} // namespace mix_matching_params

DECLARE_SOA_TABLE(MixMatchParams, "AOD", "MIXMATCHING",
                  mix_matching_params::DeltaPt,
                  mix_matching_params::DeltaEta,
                  mix_matching_params::DeltaPhi,
                  mix_matching_params::DeltaX,
                  mix_matching_params::DeltaY,
                  mix_matching_params::MchPt,
                  mix_matching_params::MchEta,
                  mix_matching_params::MchQ,
                  mix_matching_params::MftQ,
                  mix_matching_params::IsCorrectMatch);

namespace tag_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);

DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchQ, mchq, int16_t);
DECLARE_SOA_COLUMN(MftQ, mftq, int16_t);
DECLARE_SOA_COLUMN(IsTaged, isTaged, bool);
DECLARE_SOA_COLUMN(IsCorrectMatch, is_correct, bool);
} // namespace tag_matching_params

DECLARE_SOA_TABLE(TagMatchParams, "AOD", "TAGMATCHING",
                  tag_matching_params::DeltaPt,
                  tag_matching_params::DeltaEta,
                  tag_matching_params::DeltaPhi,
                  tag_matching_params::DeltaX,
                  tag_matching_params::DeltaY,
                  tag_matching_params::MchPt,
                  tag_matching_params::MchEta,
                  tag_matching_params::MchQ,
                  tag_matching_params::MftQ,
                  tag_matching_params::IsCorrectMatch);

namespace probe_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(NMFTCandTagMuon, nTagMFT, int);
DECLARE_SOA_COLUMN(NMFTCandProbeMuon, nProbeMFT, int);
DECLARE_SOA_COLUMN(TagMuonP, tagmuonp, float);

DECLARE_SOA_COLUMN(NClustMFTTracks, nClustMFT, int);
DECLARE_SOA_COLUMN(Chi2MFTTracks, chi2MFT, float);

DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);

DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchQ, mchq, int16_t);
DECLARE_SOA_COLUMN(MftQ, mftq, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, is_correct, bool);
} // namespace probe_matching_params

DECLARE_SOA_TABLE(ProbeMatchParams, "AOD", "PROBEMATCHING",
                  probe_matching_params::NMFTCandProbeMuon,
                  probe_matching_params::TagMuonP,
                  probe_matching_params::DeltaPt,
                  probe_matching_params::DeltaEta,
                  probe_matching_params::DeltaPhi,
                  probe_matching_params::DeltaX,
                  probe_matching_params::DeltaY,
                  probe_matching_params::MchPt,
                  probe_matching_params::MchEta,
                  probe_matching_params::MchQ,
                  probe_matching_params::MftQ,
                  probe_matching_params::IsCorrectMatch);

namespace muon_pair
{
DECLARE_SOA_COLUMN(NMFT, nMft, int);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
} // namespace muon_pair

DECLARE_SOA_TABLE(MuonPair, "AOD", "DIMUON", muon_pair::NMFT, muon_pair::M, muon_pair::Pt, muon_pair::Rap);

} // namespace o2::aod

struct match_mft_mch_data {

  Produces<o2::aod::MatchParams> tableMatchingParams;
  Produces<o2::aod::TagMatchParams> tableTagmatchingParams;
  Produces<o2::aod::ProbeMatchParams> tableProbematchingParams;
  Produces<o2::aod::MixMatchParams> tableMixmatchingParams;
  Produces<o2::aod::MuonPair> tableMuonPairs;
  Produces<o2::aod::MUONParams> muonParams;
  Produces<o2::aod::MFTParams> mftParams;

  HistogramRegistry registry{
    "registry",
    {{"hMchP", "MCH track total momentum (at the first station); p [GeV/c]; Counts", {HistType::kTH1F, {{2000, 0, 200}}}},
     {"hMchCorrP", "MCH track total momentum (propagated to PV); p [GeV/c]; Counts", {HistType::kTH1F, {{2000, 0, 200}}}},
     {"hMassCorrMchPair", "Corrected MCH track pair mass (propagated to PV); m [GeV/c^{2}]; Counts", {HistType::kTH1F, {{1000, 0, 10}}}}}};

  Configurable<bool> fdoMC{"cfgMC", true, ""};

  ////   Variables for selecting muon tracks
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

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
  Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
  Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  ///    Variables to add preselection for the matching table
  Configurable<float> fPreselectMatchingX{"cfgPreselectMatchingX", 999.f, ""};
  Configurable<float> fPreselectMatchingY{"cfgPreselectMatchingY", 999.f, ""};

  ///    Variables to event mixing criteria
  Configurable<float> fSaveMixedMatchingParamsRate{"cfgSaveMixedMatchingParamsRate", 0.002f, ""};
  Configurable<int> fEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
  Configurable<float> fEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};

  ////   Variables for selecting tag muon
  Configurable<float> fTagMassWindowMin{"cfgTagMassWindowMin", 2.8f, ""};
  Configurable<float> fTagMassWindowMax{"cfgTagMassWindowMax", 3.3f, ""};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  o2::globaltracking::MatchGlobalFwd mMatching;
  o2::field::MagneticField* fieldB;

  o2::ccdb::CcdbApi ccdbApi;
  int mRunNumber;

  Configurable<float> fSigmaXTagMuonCut{"cfgSigmaXTagMuonCut", 3.f, ""};
  Configurable<float> fMeanXTagMuonCut{"cfgMeanXTagMuonCut", 0.f, ""};
  Configurable<float> fSigmaYTagMuonCut{"cfgSigmaYTagMuonCut", 1.f, ""};
  Configurable<float> fMeanYTagMuonCut{"cfgMeanYTagMuonCut", 1.f, ""};

  Configurable<float> fSigmaEtaTagMuonCut{"cfgSigmaEtaTagMuonCut", 1.f, ""};
  Configurable<float> fMeanEtaTagMuonCut{"cfgMeanEtaTagMuonCut", 1.f, ""};
  Configurable<float> fSigmaPhiTagMuonCut{"cfgSigmaPhiTagMuonCut", 1.f, ""};
  Configurable<float> fMeanPhiTagMuonCut{"cfgMeanPhiTagMuonCut", 1.f, ""};

  std::unordered_map<int, float> map_vtxZ;
  std::unordered_set<int> bcs_mfttrack;
  std::unordered_map<int, int> nmfttracks;
  std::unordered_map<int, std::vector<int64_t>> map_mfttraks;
  std::unordered_map<int, int> nmuontracks;
  std::unordered_map<int, std::vector<int64_t>> map_muontracks;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
  }

  template <typename T>
  void initCCDB(T const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    std::map<string, string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  }

  enum ProagationPoint { ToVtx,
                         ToDCA };

  template <typename T>
  o2::dataformats::GlobalFwdTrack PropagateMUONTrack(T const& muon, int PropType)
  {

    auto collision = muon.collision();

    o2::dataformats::GlobalFwdTrack propmuon;

    double chi2 = muon.chi2();

    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<float> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                          muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                          muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    if (isGoodMUONTrack(muon)) {
      SMatrix55 tcovs(v1.begin(), v1.end());
      o2::track::TrackParCovFwd muontrack{muon.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack track;
      track.setParameters(tpars);
      track.setZ(muontrack.getZ());
      track.setCovariances(tcovs);
      auto mchTrack = mMatching.FwdtoMCH(track);
      if (PropType == ProagationPoint::ToVtx)
        o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
      else if (PropType == ProagationPoint::ToDCA)
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());

      auto proptrack = mMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());
    }

    v1.clear();
    v1.shrink_to_fit();

    return propmuon;
  }

  template <typename T>
  o2::track::TrackParCovFwd PropagateMFT(T const& mfttrack, int PropType)
  {
    std::vector<double> mftv1;
    SMatrix55 mftcovs{mftv1.begin(), mftv1.end()};
    SMatrix5 mftpars = {mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt()};
    o2::track::TrackParCovFwd mftpartrack = {mfttrack.z(), mftpars, mftcovs, mfttrack.chi2()};
    if (PropType == ProagationPoint::ToDCA) {
      double propVec[3] = {fabs(mfttrack.x() - mfttrack.collision().posX()), fabs(mfttrack.y() - mfttrack.collision().posY()), fabs(mfttrack.z() - mfttrack.collision().posZ())};
      double centerZ[3] = {mfttrack.x() - propVec[0] / 2., mfttrack.y() - propVec[1] / 2., mfttrack.z() - propVec[2] / 2.};
      float Bz = fieldB->getBz(centerZ);
      mftpartrack.propagateToZ(mfttrack.collision().posZ(), Bz);
    }
    return mftpartrack;
  }

  template <typename MFT, typename FWD>
  o2::track::TrackParCovFwd PropagateMFTtoMatchingPlane(MFT const& mfttrack, FWD const& muontrack)
  {
    std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
    double propVec[3] = {muontrack.x() - mfttrack.x(), muontrack.y() - mfttrack.y(), muontrack.z() - mfttrack.z()};
    double centerZ[3] = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
    float Bz = fieldB->getBz(centerZ); // gives error if the propagator is not initFielded
    SMatrix55 tmftcovs(v1.begin(), v1.end());
    SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
    o2::track::TrackParCovFwd extrap_mfttrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};
    extrap_mfttrack.propagateToZ(muontrack.z(), Bz); // z in cm
    return extrap_mfttrack;
  }

  template <typename T>
  bool isGoodMUONTrack(T track)
  {
    if (!track.has_collision())
      return false;
    if (track.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
      return false;
    if (track.chi2() > fTrackChi2MchUp)
      return false;
    if (fRabsLow1 > track.rAtAbsorberEnd() || track.rAtAbsorberEnd() > fRabsUp2)
      return false;
    if (track.rAtAbsorberEnd() < fRabsUp1 && fPdcaUp1 < track.pDca())
      return false;
    if (track.rAtAbsorberEnd() > fRabsLow2 && fPdcaUp2 < track.pDca())
      return false;
    return true;
  }

  template <typename T>
  int selectTagMUON(T track1, T track2)
  {
    if (track1.pt() > track2.pt()) {
      return track1.globalIndex();
    } else {
      return track2.globalIndex();
    }
  }

  template <typename T>
  int selectProbeMUON(T track1, T track2)
  {
    if (track1.pt() < track2.pt()) {
      return track1.globalIndex();
    } else {
      return track2.globalIndex();
    }
  }

  bool isGoodKenematicMUONTrack(o2::dataformats::GlobalFwdTrack track)
  {
    if (fEtaMchLow > track.getEta() || track.getEta() > fEtaMchUp)
      return false;
    return true;
  }

  template <typename T>
  bool isGoodMFTTrack(T track)
  {
    if (!track.has_collision())
      return false;
    if (track.chi2() > fTrackChi2MftUp)
      return false;
    if (track.nClusters() < fTrackNClustMftLow)
      return false;
    return true;
  }

  bool isGoodKenematicMFTTrack(o2::track::TrackParCovFwd track)
  {
    if (fEtaMftLow > track.getEta() || track.getEta() > fEtaMftUp) {
      return false;
    }
    return true;
  }

  template <typename MFTs, typename MUONs>
  bool isCorrectMatching(MFTs const& mft, MUONs const& muon)
  {
    int idmuon = muon.mcParticleId();
    int idmft = mft.mcParticleId();

    if (idmuon == -1 || idmft == -1) {
      return false;
    }
    if (idmuon != idmft) {
      return false;
    } else {
      return true;
    }
  }

  template <typename MFTs>
  void setMFT(MFTs const& mfttracks,
              std::unordered_set<int>& bcs_mfttrack, std::unordered_map<int, std::vector<int64_t>>& map_mfttraks,
              std::unordered_map<int, float>& map_vtxZ, std::unordered_map<int, int>& nmfttracks)
  {
    for (const auto& mfttrack : mfttracks) {

      if (!isGoodMFTTrack(mfttrack)) {
        continue;
      }

      bcs_mfttrack.insert(mfttrack.collisionId());

      std::vector<int64_t>& tracks = map_mfttraks[mfttrack.collisionId()];
      tracks.push_back(mfttrack.globalIndex());

      o2::track::TrackParCovFwd mftpartrack = PropagateMFT(mfttrack, ProagationPoint::ToDCA);

      if (!isGoodKenematicMFTTrack(mftpartrack)) {
        continue;
      }

      map_vtxZ[mfttrack.collisionId()] = mfttrack.collision().posZ();

      float dx = mftpartrack.getX() - mfttrack.collision().posX();
      float dy = mftpartrack.getY() - mfttrack.collision().posY();

      mftParams(mfttrack.nClusters(), mfttrack.isCA(), mfttrack.chi2(), mfttrack.sign(), mftpartrack.getPt(), mftpartrack.getEta(), mftpartrack.getPhi(), dx, dy);
      nmfttracks[mfttrack.collisionId()]++;
    }
  }

  template <typename MUONs>
  void setMUON(MUONs const& muontracks, std::unordered_map<int, int>& nmuontracks, std::unordered_map<int, std::vector<int64_t>>& map_muontracks,
               std::unordered_map<int, std::vector<int64_t>> map_mfttraks)
  {

    for (auto muontrack : muontracks) {

      if (!isGoodMUONTrack(muontrack))
        continue;

      o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);

      if (!isGoodKenematicMUONTrack(propmuonAtPV)) {
        continue;
      }

      std::vector<int64_t>& tracks = map_muontracks[muontrack.collisionId()];
      tracks.push_back(muontrack.globalIndex());

      bool hasMFT = false;

      std::vector<int64_t>& mfttracks = map_mfttraks[muontrack.collisionId()];

      if (mfttracks.size() > 0) {
        hasMFT = true;
      }
      muonParams(muontrack.chi2(), muontrack.rAtAbsorberEnd(), muontrack.sign(), propmuonAtPV.getPt(), propmuonAtPV.getEta(), propmuonAtPV.getPhi(), hasMFT);
      nmuontracks[muontrack.collisionId()]++;
    }
  }

  template <typename MFTs, typename MUONs>
  void calcMatchingParams(MFTs const& mfttrack, MUONs const& muontrack, float& deltaX, float& deltaY, float& deltaPt, float& deltaPhi, float& deltaEta)
  {
    o2::track::TrackParCovFwd mfttrack_at_matching = PropagateMFTtoMatchingPlane(mfttrack, muontrack);
    V1.SetPtEtaPhi(mfttrack_at_matching.getPt(), mfttrack_at_matching.getEta(), mfttrack_at_matching.getPhi());
    V2.SetPtEtaPhi(muontrack.pt(), muontrack.eta(), muontrack.phi());
    deltaX = mfttrack_at_matching.getX() - muontrack.x();
    deltaY = mfttrack_at_matching.getY() - muontrack.y();
    deltaPt = mfttrack_at_matching.getPt() - muontrack.pt();
    deltaPhi = V1.DeltaPhi(V2);
    deltaEta = mfttrack_at_matching.getEta() - muontrack.eta();
  }

  template <typename MFTs, typename MUONs>
  void fillMatchingParams(MFTs const& mfttrack, MUONs const& muontrack)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    calcMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta);
    if (!isPassPreSelection(deltaX, deltaY))
      return;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableMatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                        muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), false);
  }

  template <typename MFTs, typename MUONs>
  void fillMatchingParamsMC(MFTs const& mfttrack, MUONs const& muontrack)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    calcMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta);
    if (!isPassPreSelection(deltaX, deltaY))
      return;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableMatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                        muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), isCorrectMatching(mfttrack, muontrack));
  }

  template <typename MFTs, typename MUONs>
  void fillMixMatchingParams(MFTs const& mfttrack, MUONs const& muontrack)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    calcMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta);
    if (!isPassPreSelection(deltaX, deltaY))
      return;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableMixmatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                           muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), false);
  }

  template <typename MFTs, typename MUONs>
  void fillMixMatchingParamsMC(MFTs const& mfttrack, MUONs const& muontrack)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    calcMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta);
    if (!isPassPreSelection(deltaX, deltaY))
      return;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableMixmatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                           muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), isCorrectMatching(mfttrack, muontrack));
  }

  template <typename MFTs, typename MUONs>
  bool checkTagMatchingParams(MFTs const& mfttrack, MUONs const& muontrack, float& deltaX, float& deltaY, float& deltaPt, float& deltaPhi, float& deltaEta)
  {
    calcMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta);
    if (!isPassPreSelection(deltaX, deltaY))
      return false;
    if (!isGoodTagMuon(deltaX, deltaY, deltaPhi, deltaEta))
      return false;
    return true;
  }

  template <typename MFTs, typename MUONs>
  bool checkProbeMatchingParams(MFTs const& mfttrack, MUONs const& muontrack, float& deltaX, float& deltaY, float& deltaPt, float& deltaPhi, float& deltaEta)
  {
    calcMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta);
    if (!isPassPreSelection(deltaX, deltaY))
      return false;
    if (!isGoodProbeMuon(deltaX, deltaY, deltaPhi, deltaEta))
      return false;
    return true;
  }

  template <typename MFTs, typename MUONs>
  bool fillTagMatchingParams(MFTs const& mfttrack, MUONs const& muontrack)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    if (!checkTagMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta))
      return false;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableTagmatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                           muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), false);
    return true;
  }

  template <typename MFTs, typename MUONs>
  bool fillTagMatchingParamsMC(MFTs const& mfttrack, MUONs const& muontrack)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    if (!checkTagMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta))
      return false;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableTagmatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                           muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), isCorrectMatching(mfttrack, muontrack));
    return true;
  }

  template <typename MFTs, typename MUONs>
  bool fillProbeMatchingParams(MFTs const& mfttrack, MUONs const& muontrack, int nMFTCandsProbeMUON, float tagP)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    if (!checkProbeMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta))
      return false;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableProbematchingParams(nMFTCandsProbeMUON, tagP, deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                             muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), false);
    return true;
  }

  template <typename MFTs, typename MUONs>
  bool fillProbeMatchingParamsMC(MFTs const& mfttrack, MUONs const& muontrack, int nMFTCandsProbeMUON, float tagP)
  {
    float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
    if (!checkProbeMatchingParams(mfttrack, muontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta))
      return false;
    o2::dataformats::GlobalFwdTrack muontrackAtPV = PropagateMUONTrack(muontrack, ProagationPoint::ToVtx);
    tableProbematchingParams(nMFTCandsProbeMUON, tagP, deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                             muontrackAtPV.getPt(), muontrackAtPV.getEta(), muontrack.sign(), mfttrack.sign(), isCorrectMatching(mfttrack, muontrack));
    return true;
  }

  template <typename Colls>
  bool isGoodEventMix(Colls const& coll, int bc1, int bc2)
  {
    if (bc1 == bc2) {
      return false;
    }
    if (fabs(nmfttracks[bc1] - nmfttracks[bc2]) > fEventMaxDeltaNMFT) {
      return false;
    }
    if (fabs(coll.posZ() - map_vtxZ[bc2]) > fEventMaxDeltaVtxZ) {
      return false;
    }
    return true;
  }

  bool isPassPreSelection(float deltaX, float deltaY)
  {
    if (fabs(deltaX) > fPreselectMatchingX) {
      return false;
    }
    if (fabs(deltaY) > fPreselectMatchingY) {
      return false;
    }
    return true;
  }

  bool isGoodTagPair()
  {
    if (fTagMassWindowMin > dimuonLV.M() || dimuonLV.M() > fTagMassWindowMax) {
      return false;
    }
    return true;
  }

  bool isGoodTagMuon(float deltaX, float deltaY, float deltaPhi, float deltaEta)
  {
    float rTagXY = pow((deltaX - fMeanXTagMuonCut) / (fSigmaXTagMuonCut * 3), 2) +
                   pow((deltaY - fMeanYTagMuonCut) / (fSigmaYTagMuonCut * 3), 2);
    float rTagEtaPhi = pow((deltaEta - fMeanEtaTagMuonCut) / (fSigmaEtaTagMuonCut * 3), 2) +
                       pow((deltaPhi - fMeanPhiTagMuonCut) / (fSigmaPhiTagMuonCut * 3), 2);
    if (rTagXY < 1 && rTagEtaPhi) {
      return true;
    }
    return false;
  }

  bool isGoodProbeMuon(float deltaX, float deltaY, float deltaPhi, float deltaEta)
  {
    float rTagXY = pow((deltaX - fMeanXTagMuonCut) / (fSigmaXTagMuonCut * 10), 2) +
                   pow((deltaY - fMeanYTagMuonCut) / (fSigmaYTagMuonCut * 10), 2);
    float rTagEtaPhi = pow((deltaEta - fMeanEtaTagMuonCut) / (fSigmaEtaTagMuonCut * 10), 2) +
                       pow((deltaPhi - fMeanPhiTagMuonCut) / (fSigmaPhiTagMuonCut * 10), 2);
    if (rTagXY < 1 && rTagEtaPhi) {
      return true;
    }
    return false;
  }

  template <typename MUONs, typename MFTs>
  void matchingParams(MUONs const& muontracks, MFTs const& mfttracks)
  {
    for (const auto& map_muontrack : map_muontracks) {
      for (auto imuontrack1 : map_muontrack.second) {
        for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
          fillMatchingParams(mfttracks.rawIteratorAt(imfttrack1), muontracks.rawIteratorAt(imuontrack1));
        } // end of loop imfttrack1
        for (auto bc2 : bcs_mfttrack) {
          if (!isGoodEventMix(muontracks.rawIteratorAt(imuontrack1).collision(), map_muontrack.first, bc2))
            continue;
          for (auto imfttrack1 : map_mfttraks[bc2]) {
            fillMixMatchingParams(mfttracks.rawIteratorAt(imfttrack1), muontracks.rawIteratorAt(imuontrack1));
          }
        } // end of loop bc2
      } // end of loop imuontrack1
    } // end of loop map_muontracks
  } // end of matchingParamsMC

  template <typename MUONs, typename MFTs>
  void matchingParamsMC(MUONs const& muontracks, MFTs const& mfttracks)
  {
    for (const auto& map_muontrack : map_muontracks) {
      for (auto imuontrack1 : map_muontrack.second) {
        for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
          fillMatchingParamsMC(mfttracks.rawIteratorAt(imfttrack1), muontracks.rawIteratorAt(imuontrack1));
        } // end of loop imfttrack1
        for (auto bc2 : bcs_mfttrack) {
          if (!isGoodEventMix(muontracks.rawIteratorAt(imuontrack1).collision(), map_muontrack.first, bc2))
            continue;
          for (auto imfttrack1 : map_mfttraks[bc2]) {
            fillMixMatchingParamsMC(mfttracks.rawIteratorAt(imfttrack1), muontracks.rawIteratorAt(imuontrack1));
          }
        } // end of loop bc2
      } // end of loop imuontrack1
    } // end of loop map_muontracks
  } // end of matchingParamsMC

  template <typename MUONs, typename MFTs>
  void matchingParamsTagAndProbeMC(MUONs const& muontracks, MFTs const& mfttracks)
  {
    for (const auto& map_muontrack : map_muontracks) {
      if (nmfttracks[map_muontrack.first] < 1) {
        continue;
      }
      for (auto imuontrack1 : map_muontrack.second) {
        o2::dataformats::GlobalFwdTrack muontrackAtPV1 = PropagateMUONTrack(muontracks.rawIteratorAt(imuontrack1), ProagationPoint::ToVtx);
        for (auto imuontrack2 : map_muontrack.second) {
          if (imuontrack1 == imuontrack2)
            continue;
          o2::dataformats::GlobalFwdTrack muontrackAtPV2 = PropagateMUONTrack(muontracks.rawIteratorAt(imuontrack2), ProagationPoint::ToVtx);

          if (fabs(muontracks.rawIteratorAt(imuontrack1).sign() + muontracks.rawIteratorAt(imuontrack2).sign()) > 0)
            continue;

          muon1LV.SetPtEtaPhiM(muontrackAtPV1.getPt(), muontrackAtPV1.getEta(), muontrackAtPV1.getPhi(), mMu);
          muon2LV.SetPtEtaPhiM(muontrackAtPV2.getPt(), muontrackAtPV2.getEta(), muontrackAtPV2.getPhi(), mMu);
          dimuonLV = muon1LV + muon2LV;

          tableMuonPairs(nmfttracks[map_muontrack.first], dimuonLV.M(), dimuonLV.Pt(), dimuonLV.Rapidity());

          if (!isGoodTagPair())
            continue;

          auto tagmuontrack = muontracks.rawIteratorAt(selectTagMUON(muontracks.rawIteratorAt(imuontrack1), muontracks.rawIteratorAt(imuontrack2)));

          int nMFTCandsTagMUON = 0;

          for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
            if (!fillTagMatchingParamsMC(mfttracks.rawIteratorAt(imfttrack1), tagmuontrack)) {
              continue;
            }
            nMFTCandsTagMUON++;
          } // end of loop imfttrack1

          if (nMFTCandsTagMUON > 1)
            continue;

          auto probemuontrack = muontracks.rawIteratorAt(selectProbeMUON(muontracks.rawIteratorAt(imuontrack1), muontracks.rawIteratorAt(imuontrack2)));

          int nMFTCandsProbeMUON = 0;

          for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
            float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
            if (!checkProbeMatchingParams(mfttracks.rawIteratorAt(imfttrack1), probemuontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta))
              continue;
            nMFTCandsProbeMUON++;
          }

          for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
            if (!fillProbeMatchingParamsMC(mfttracks.rawIteratorAt(imfttrack1), probemuontrack, nMFTCandsProbeMUON, tagmuontrack.p())) {
              continue;
            }
          } // end of loop imfttrack1
        } // end of loop imuontrack2
      } // end of loop imuontrack1
    } // end of loop map_muontrack
  } // end of matchingParamsTagAndProbeMC

  template <typename MUONs, typename MFTs>
  void matchingParamsTagAndProbe(MUONs const& muontracks, MFTs const& mfttracks)
  {
    for (const auto& map_muontrack : map_muontracks) {
      if (nmfttracks[map_muontrack.first] < 1) {
        continue;
      }
      for (auto imuontrack1 : map_muontrack.second) {
        o2::dataformats::GlobalFwdTrack muontrackAtPV1 = PropagateMUONTrack(muontracks.rawIteratorAt(imuontrack1), ProagationPoint::ToVtx);
        for (auto imuontrack2 : map_muontrack.second) {
          if (imuontrack1 == imuontrack2)
            continue;
          o2::dataformats::GlobalFwdTrack muontrackAtPV2 = PropagateMUONTrack(muontracks.rawIteratorAt(imuontrack2), ProagationPoint::ToVtx);

          if (fabs(muontracks.rawIteratorAt(imuontrack1).sign() + muontracks.rawIteratorAt(imuontrack2).sign()) > 0)
            continue;

          muon1LV.SetPtEtaPhiM(muontrackAtPV1.getPt(), muontrackAtPV1.getEta(), muontrackAtPV1.getPhi(), mMu);
          muon2LV.SetPtEtaPhiM(muontrackAtPV2.getPt(), muontrackAtPV2.getEta(), muontrackAtPV2.getPhi(), mMu);
          dimuonLV = muon1LV + muon2LV;

          tableMuonPairs(nmfttracks[map_muontrack.first], dimuonLV.M(), dimuonLV.Pt(), dimuonLV.Rapidity());

          if (!isGoodTagPair())
            continue;

          auto tagmuontrack = muontracks.rawIteratorAt(selectTagMUON(muontracks.rawIteratorAt(imuontrack1), muontracks.rawIteratorAt(imuontrack2)));

          int nMFTCandsTagMUON = 0;

          for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
            if (!fillTagMatchingParams(mfttracks.rawIteratorAt(imfttrack1), tagmuontrack)) {
              continue;
            }
            nMFTCandsTagMUON++;
          } // end of loop imfttrack1

          if (nMFTCandsTagMUON > 1)
            continue;

          auto probemuontrack = muontracks.rawIteratorAt(selectProbeMUON(muontracks.rawIteratorAt(imuontrack1), muontracks.rawIteratorAt(imuontrack2)));

          int nMFTCandsProbeMUON = 0;

          for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
            float deltaX, deltaY, deltaPt, deltaPhi, deltaEta;
            if (!checkProbeMatchingParams(mfttracks.rawIteratorAt(imfttrack1), probemuontrack, deltaX, deltaY, deltaPt, deltaPhi, deltaEta))
              continue;
            nMFTCandsProbeMUON++;
          }

          for (auto imfttrack1 : map_mfttraks[map_muontrack.first]) {
            if (!fillProbeMatchingParams(mfttracks.rawIteratorAt(imfttrack1), probemuontrack, nMFTCandsProbeMUON, tagmuontrack.p())) {
              continue;
            }
          } // end of loop imfttrack1
        } // end of loop imuontrack2
      } // end of loop imuontrack1
    } // end of loop map_muontrack
  } // end of matchingParamsTagAndProbe

  void process(aod::Collisions const& collisions)
  {
    LOG(info) << "DUMMY PROCESS";
  }

  void processMC(aod::Collisions const& collisions,
                 soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti, aod::MatchedToFT0> const& ebcs,
                 soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::McFwdTrackLabels> const& muontracks,
                 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks)
  {
    map_vtxZ.clear();
    bcs_mfttrack.clear();
    nmfttracks.clear();
    map_mfttraks.clear();
    nmuontracks.clear();
    map_muontracks.clear();
    initCCDB(ebcs.begin());
    setMFT(mfttracks, bcs_mfttrack, map_mfttraks, map_vtxZ, nmfttracks);
    setMUON(muontracks, nmuontracks, map_muontracks, map_mfttraks);
    matchingParamsMC(muontracks, mfttracks);
    matchingParamsTagAndProbeMC(muontracks, mfttracks);
  }

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti, aod::MatchedToFT0> const& ebcs,
                   soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA> const& muontracks,
                   aod::MFTTracks const& mfttracks)
  {
    map_vtxZ.clear();
    bcs_mfttrack.clear();
    nmfttracks.clear();
    map_mfttraks.clear();
    nmuontracks.clear();
    map_muontracks.clear();

    initCCDB(ebcs.begin());
    setMFT(mfttracks, bcs_mfttrack, map_mfttraks, map_vtxZ, nmfttracks);
    setMUON(muontracks, nmuontracks, map_muontracks, map_mfttraks);
    matchingParams(muontracks, mfttracks);
    matchingParamsTagAndProbe(muontracks, mfttracks);
  }

  PROCESS_SWITCH(match_mft_mch_data, processMC, "Produce tables of reconstructed information for MC data", false);
  PROCESS_SWITCH(match_mft_mch_data, processData, "Produce tables of reconstructed information for real data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<match_mft_mch_data>(cfgc)};
}
