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

#include "TMatrixD.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct ParityTask {
  OutputObj<TH1F> hEvSel{"h_evsel"};
  OutputObj<TH1F> hPt{"h_pt"};
  OutputObj<TH2F> hNtracksPt{"h_ntracks_pt"};
  OutputObj<TH1F> hDet{"h_det"};
  OutputObj<TH1F> hNtracks{"h_ntracks"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCut{"trackPtCut", 0.1, "minimum track pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.9, "track eta cut"};

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCut);

  TrackSelection globalTracks;

  void init(InitContext const&)
  {
    // set up global tracks and adjust as necessary
    globalTracks = getGlobalTrackSelection();
    globalTracks.SetEtaRange(-.9, .9);

    hEvSel.setObject(new TH1F("h_evsel", "h_evsel", 10, 0., 10.));
    hPt.setObject(new TH1F("h_pt", "p_{T};p_{T} (GeV/#it{c})",
			    100, 0., 100.));
    hNtracksPt.setObject(new TH2F("h_ntracks_pt", ";n_{trk};p_{T} (GeV/#it{c})", 200, -.5, 199.5, 100, 0., 10.));
    hDet.setObject(new TH1F("h_det", "det", 3, -1.5, 1.5));
    hNtracks.setObject(new TH1F("h_ntracks", ";N_{trk}", 100, 0., 100.));
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>> const& tracks)
  {
    LOG(debug) << "Process data!";
    hEvSel->Fill("all", 1.);

    if (!collision.sel8()) return;
    hEvSel->Fill("sel8", 1.);

    if (tracks.size() > 20) return;
    hEvSel->Fill("tracks", 1.);

    std::vector<double> trackpt;

    for (const auto &track : tracks) {
      LOG(debug) << "track pt" << track.pt();
      if(!globalTracks.IsSelected(track)) continue;
      trackpt.emplace_back(track.pt());
    	hPt->Fill(track.pt());
    }
    hNtracks->Fill(trackpt.size());
    std::sort(trackpt.begin(), trackpt.end(), std::greater<double>());
    int n_track = 1;
    for (const auto & pt : trackpt) {
      hNtracksPt->Fill(n_track, pt);
      n_track++;
    }

    TMatrixD m(3,3);
    for (auto track0 = tracks.begin(); track0 != tracks.end(); ++track0) {
      m[0][0]=(*track0).px();
      m[1][0]=(*track0).py();
      m[2][0]=(*track0).pz();
      for (auto track1 = track0 + 1; track1 != tracks.end(); ++track1) {
        m[0][1]=(*track1).px();
        m[1][1]=(*track1).py();
        m[2][1]=(*track1).pz();
        for (auto track2 = track1 + 1; track2 != tracks.end(); ++track2) {
          m[0][2]=(*track2).px();
          m[1][2]=(*track2).py();
          m[2][2]=(*track2).pz();
          const auto d = m.Determinant();
          hDet->Fill((d > 0) - (d < 0));
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<ParityTask>(cfgc));
  return WorkflowSpec{tasks};
}
