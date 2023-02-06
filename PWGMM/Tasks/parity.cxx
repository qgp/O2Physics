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

// simple task checking parity in min. bias collisions
//
// Author: Jochen Klein, Magnus Mager

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct ParityTask {
  OutputObj<TH1F> hDet{"h_det"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCut{"trackPtCut", 0.1, "minimum track pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.9, "track eta cut"};

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCut);

  void init(InitContext const&)
  {
    hDet.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100.));
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, aod::Tracks const& tracks)
  {
    LOG(debug) << "Process data!";

    for (const auto &track : tracks) {
        LOG(debug) << "track pt" << track.pt();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<ParityTask>(cfgc));
  return WorkflowSpec{tasks};
}
