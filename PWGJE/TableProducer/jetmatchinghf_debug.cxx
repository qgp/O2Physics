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
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetMatchingHF {
  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using Tracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using HfCandidates = soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>;
  using McParticles = soa::Join<aod::McParticles, aod::HfCandProng2MCGen>;
  using DetectorLevelJets = soa::Join<aod::MCDetectorLevelHFJets, aod::MCDetectorLevelHFJetConstituents>;
  using ParticleLevelJets = soa::Join<aod::MCParticleLevelHFJets, aod::MCParticleLevelHFJetConstituents>;

  Produces<aod::MatchedMCParticleLevelJets> jetsMCMatching;
  Produces<aod::MatchedJets> jetsRecMatching;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<ParticleLevelJets> ParticleLevelJetsPerMcCollision = aod::jet::mcCollisionId;

  void init(InitContext const&)
  {
  }

  void process(Collisions::iterator collision, aod::McCollisions const& mcCollisions,
    ParticleLevelJets const& jetsMC, DetectorLevelJets const& jetsRec,
    HfCandidates hfcandidates, Tracks const& tracks, 
    McParticles const& particlesMC)
  {
    LOGF(debug, "analysing collision %d - %d, MC collision %d", 
         collision.index(), collision.globalIndex(), collision.mcCollisionId());

    for (const auto &t : tracks) {
      LOGF(debug, "track index: %d-%d (coll %d-%d) with mcparticle: %d %d", 
           t.index(), t.globalIndex(), t.collisionId(), t.collision().index(),
           t.has_mcParticle() ? t.mcParticle().index() : -1, t.has_mcParticle() ? t.mcParticle().globalIndex() : -1);
    }

    auto mcps = particlesMC.sliceBy(mcParticlesPerMcCollision, collision.mcCollisionId());
    for (const auto &mcp : mcps) {
      LOGF(debug, "mcparticle index: %d (MC coll %d)", mcp.globalIndex(), mcp.mcCollisionId());
    }

    // match rec to MC
    for (const auto &jet : jetsRec) {
      const auto &coll = jet.collision_as<Collisions>();
      LOGF(info, "jet index: %d (coll %d-%d) with %d tracks, %d HF candidates", 
           jet.index(), jet.collisionId(), coll.index(), jet.tracks().size(), jet.hfcandidates().size());

      const auto &tracks = jet.tracks_as<Tracks>();
      for (const auto &track : tracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for track %d", track.index());
          continue;
        }

        const auto &mcparticle = track.mcParticle_as<McParticles>();
        LOGF(info, "track %d of jet %d in coll %d-%d has mcparticle: %d-%d", 
            track.globalIndex(), jet.globalIndex(), track.collisionId(), track.collision().globalIndex(), track.mcParticleId(), mcparticle.globalIndex());
      }

      for (const auto &cand : hfcandidates) {
        LOGF(info, "candidate %d with prongs: %d, %d", cand.globalIndex(), cand.index0Id(), cand.index1Id());
      }
      for (const auto id : jet.hfcandidatesIds()) {
        const auto &cand = hfcandidates.iteratorAt(id);
        LOGF(info, "candidate %d with prongs: %d, %d", cand.globalIndex(), 0, 0);
      }
      const auto &cands = jet.hfcandidates_as<HfCandidates>();
      for (const auto &cand : cands) {
        const auto &daughter0 = cand.index0_as<Tracks>();
        const auto &daughter1 = cand.index1_as<Tracks>();
        LOGF(info, "MC candidate %d with prongs: %d (MC %d), %d (MC %d)", cand.globalIndex(), 
             daughter0.globalIndex(), daughter0.mcParticle_as<McParticles>().globalIndex(),
             daughter1.globalIndex(), daughter1.mcParticle_as<McParticles>().globalIndex());
      }
    }
      // // should use HF candidate instead
      // if (TMath::Abs(mcparticle.pdgCode()) != 421)
      //   continue;
      // for (const auto &mcc : constituentsMC) {
      //   if (mcc.track().index() == mcparticle.index()) {
      //     LOGF(info, "matched mcd: %d - %d", c.jet().globalIndex(), mcc.jet().globalIndex());
      //     jetsRecMatching(c.jet().globalIndex(), mcc.jet().globalIndex());
      //     break;
      //   }
      // }
    // }

    // match MC to rec
    auto jetsPL = jetsMC.sliceBy(ParticleLevelJetsPerMcCollision, collision.mcCollisionId());
    for (const auto &jet : jetsPL) {
      LOGF(info, "MC jet index: %d (coll %d-%d) with %d tracks, %d HF candidates", 
           jet.index(), jet.mcCollisionId(), jet.mcCollision().index(), jet.tracks().size(), jet.hfcandidates().size());
      for (const auto &cand : jet.hfcandidates_as<McParticles>()) {
        const auto &daughters = cand.daughters_as<McParticles>();
        LOGF(info, "MC candidate %d-%d with prongs: %d, %d", 
             cand.index(), cand.globalIndex(), daughters.iteratorAt(0).globalIndex(), daughters.iteratorAt(1).globalIndex());
      }
    }

    // for (const auto &c : constituentsMC) {
    //   const auto part = c.track();
    //   if (TMath::Abs(part.pdgCode()) != 421)
    //     continue;
    //   // aod::JetTrackConstituent rec_constituent;

    //   for (const auto &rc : constituentsRec) {
    //     const auto track = rc.track_as<soa::Join<aod::Tracks, aod::McTrackLabels>>();
    //     if (!track.has_mcParticle()) continue;
    //     auto mcparticle = track.mcParticle();
    //     if (mcparticle.index() == part.index()) {
    //       LOGF(info, "matched mcp: %d - %d", c.jet().globalIndex(), rc.jet().globalIndex());
    //       jetsMCMatching(c.jet().globalIndex(), rc.jet().globalIndex());
    //       break;
    //     }
    //   }
    // }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetMatchingHF>(cfgc, TaskName{"jet-matching-hf"})};
}
