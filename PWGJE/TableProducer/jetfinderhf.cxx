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

// jet finder task
//
// Author: Nima Zardoshti

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec hfjetMode = {
    "hfjetMode",
    VariantType::String,
    "",
    {"HF jet finder mode."},
  };
  workflowOptions.push_back(hfjetMode);
}

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTable, typename TrackConstituentTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<TrackConstituentTable> trackConstituents;
  OutputObj<TH1F> hJetPt{"h_jet_pt"};
  OutputObj<TH1F> hD0Pt{"h_D0_pt"};

  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder;

  void init(InitContext const&)
  {
    hJetPt.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100.));
    hD0Pt.setObject(new TH1F("h_D0_pt", "jet p_{T,D};p_{T,D} (GeV/#it{c})",
                             60, 0., 60.));
  }

  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  //need enum as configurable
  enum pdgCode { pdgD0 = 421 };

  Filter trackCuts = (aod::track::pt > 0.15f && aod::track::eta > -0.9f && aod::track::eta < 0.9f);
  Filter seltrack = (aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar);

  void processData(aod::Collision const& collision,
               soa::Filtered<aod::Tracks> const& tracks,
               soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> const& candidates)
  {
    std::cout << "Per Event" << std::endl;
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : tracks) {
        auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
        if (candidate.index0().globalIndex() == track.globalIndex() || candidate.index1().globalIndex() == track.globalIndex()) { //is it global index?
          continue;
        }
        inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
        inputParticles.back().set_user_index(track.globalIndex());
      }
      inputParticles.emplace_back(candidate.px(), candidate.py(), candidate.pz(), candidate.e(RecoDecay::getMassPDG(pdgD0)));
      inputParticles.back().set_user_index(1);

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        isHFJet = false;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == 1 && (candidate.isSelD0() == 1 || candidate.isSelD0bar() == 1)) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.eta(), jet.phi(), jet.pt(),
                    jet.area(), jet.E(), jet.m(), jetFinder.jetR);
          for (const auto& constituent : jet.constituents()) {
            trackConstituents(jetsTable.lastIndex(), constituent.user_index());
          }
          hJetPt->Fill(jet.pt());
          std::cout << "Filling" << std::endl;
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processData, "HF jet finding on data", true);

  // should we use the matching to MC from HFCandidateCreator2Prong here?
  void processMCD(aod::Collision const& collision,
               soa::Filtered<aod::Tracks> const& tracks,
               soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> const& candidates) {
    LOG(debug) << "Per Event MCP";
     std::cout << "Per Event" << std::endl;
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : tracks) {
        auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
        if (candidate.index0().globalIndex() == track.globalIndex() || candidate.index1().globalIndex() == track.globalIndex()) { //is it global index?
          continue;
        }
        inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
        inputParticles.back().set_user_index(track.globalIndex());
      }
      inputParticles.emplace_back(candidate.px(), candidate.py(), candidate.pz(), candidate.e(RecoDecay::getMassPDG(pdgD0)));
      inputParticles.back().set_user_index(1);

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        isHFJet = false;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == 1 && (candidate.isSelD0() == 1 || candidate.isSelD0bar() == 1)) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.eta(), jet.phi(), jet.pt(),
                    jet.area(), jet.E(), jet.m(), jetFinder.jetR);
          for (const auto& constituent : jet.constituents()) {
            trackConstituents(jetsTable.lastIndex(), constituent.user_index());
          }
          hJetPt->Fill(jet.pt());
          std::cout << "Filling" << std::endl;
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCD, "HF jet finding on MC detector level", false);

  // can we get the candidates from another task already?
  void processMCP(aod::McCollision const& collision,
               aod::McParticles const& particles) {
    LOG(debug) << "Per Event MCP";
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    // TODO: candidates should come from HF task!
    std::vector<aod::McParticle> candidates;
    for (auto const &part : particles) {
      if (std::abs(part.pdgCode()) == 421) {
        candidates.push_back(part);
      }
    }

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : particles) {
        // TODO: check what mass to use?
        auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
        // TODO: check use of indices, this should be daughters!
        // TODO: check if D0 decays at particle level
        if (std::find(std::begin(candidate.mothersIds()), std::end(candidate.mothersIds()), track.index()) != std::end(candidate.mothersIds()))
          continue;
        inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
        inputParticles.back().set_user_index(track.globalIndex());
      }
      inputParticles.emplace_back(candidate.px(), candidate.py(), candidate.pz(), candidate.e());
      // TODO: check use of this index, should we write to the table?
      inputParticles.back().set_user_index(1);

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        isHFJet = false;
        for (const auto& constituent : jet.constituents()) {
          // we have only selected D0s above, so enough to test user_index ???
          if (constituent.user_index() == 1) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.eta(), jet.phi(), jet.pt(),
                    jet.area(), jet.E(), jet.m(), jetFinder.jetR);
          for (const auto& constituent : jet.constituents()) {
            trackConstituents(jetsTable.lastIndex(), constituent.user_index());
          }
          hJetPt->Fill(jet.pt());
          std::cout << "Filling" << std::endl;
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCP, "HF jet finding on MC particle level", false);
};

using JetFinderHF = JetFinderHFTask<o2::aod::HFJets, o2::aod::HFJetTrackConstituents>;
using MCParticleLevelJetFinderHF = JetFinderHFTask<o2::aod::MCParticleLevelHFJets, o2::aod::MCParticleLevelHFJetTrackConstituents>;
using MCDetectorLevelJetFinderHF = JetFinderHFTask<o2::aod::MCDetectorLevelHFJets, o2::aod::MCDetectorLevelHFJetTrackConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  auto hfjetMode = cfgc.options().get<std::string>("hfjetMode");
  
  if (hfjetMode.find("data") != std::string::npos)
    tasks.emplace_back(adaptAnalysisTask<JetFinderHF>(cfgc,
      SetDefaultProcesses{{{"processData", true}, {"processMCP", false}, {"processMCD", false}}},
      TaskName{"jet-finder-hf-data"}));

  if (hfjetMode.find("mcp") != std::string::npos)
    tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetFinderHF>(cfgc,
      SetDefaultProcesses{{{"processData", false}, {"processMCP", true}, {"processMCD", false}}},
      TaskName{"jet-finder-hf-mcp"}));

  if (hfjetMode.find("mcd") != std::string::npos)
    tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetFinderHF>(cfgc,
      SetDefaultProcesses{{{"processData", false}, {"processMCP", false}, {"processMCD", true}}},
      TaskName{"jet-finder-hf-mcd"}));

  return WorkflowSpec{tasks};
}
