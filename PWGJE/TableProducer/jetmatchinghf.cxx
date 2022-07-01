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

/// \file jetmatchinghf.cxx
/// \brief Match jets containing the same D0s
///
/// \author Jochen Klein <jochen.klein@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetMatchingHF {
  Produces<aod::MatchedMCParticleLevelJets> jetsMCMatching;
  Produces<aod::MatchedJets> jetsRecMatching;

  void init(InitContext const&)
  {
  }

  void process(
    soa::Filtered<o2::aod::Collisions>::iterator const& collision,
    aod::MCParticleLevelHFJets const& jetsMC,
    aod::MCParticleLevelHFJetTrackConstituents const& constituentsMC,
    aod::McParticles const& particlesMC,
    aod::MCDetectorLevelHFJets const& jetsRec,
    aod::MCDetectorLevelHFJetTrackConstituents const& constituentsRec,
    soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks
    )
  {
    return;
    // match rec to MC
    for (const auto &c : constituentsRec) {
      // TODO: check access to track
      const auto track = tracks.rawIteratorAt(c.track().index());
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for track %d", track.index());
        continue;
      }
      auto mcparticle = track.mcParticle();
      // should use HF candidate instead
      if (TMath::Abs(mcparticle.pdgCode()) != 421)
        continue;
      aod::MCParticleLevelJetTrackConstituent mc_constituent;
      for (const auto &mcc : constituentsMC) {
        if (mcc.track().index() == mcparticle.index()) {
          jetsRecMatching(c.jet().globalIndex(), mc_constituent.jet().globalIndex());
          // mc_constituent = mcc;
          break;
        }
      }
    }

    // match MC to rec
    for (const auto &c : constituentsMC) {
      const auto track = c.track();
      if (TMath::Abs(track.pdgCode()) != 421)
        continue;
      aod::JetTrackConstituent rec_constituent;
      for (const auto &rc : constituentsRec) {
        const auto track = tracks.rawIteratorAt(rc.track().index());
        if (!track.has_mcParticle()) continue;
        if (track.mcParticle().index() == c.track().index()) {
          // rec_constituent = rc;
          jetsMCMatching(c.jet().globalIndex(), rec_constituent.jet().globalIndex());
          break;
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetMatchingHF>(cfgc, TaskName{"jet-matching-hf"})};
}
