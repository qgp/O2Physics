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

/// \file JetDetivedDataUtilities.h
/// \brief Jet derived data related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_
#define PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_

#include <string>

#include "PWGJE/Core/JetFinder.h"

namespace JetDerivedDataUtilities
{

enum JCollisionSel {
  sel8 = 0,
  sel7 = 1,
  kTVXinEMC = 2
};

template <typename T>
bool selectCollision(T const& collision, int eventSelection = -1)
{
  if (eventSelection == -1) {
    return true;
  }

  return (collision.eventSel() >> eventSelection) & 1;
}

int initialiseEventSelection(std::string eventSelection)
{
  if (eventSelection == "sel8") {
    return JCollisionSel::sel8;
  } else if (eventSelection == "sel7") {
    return JCollisionSel::sel7;
  } else if (eventSelection == "kTVXinEMC") {
    return JCollisionSel::kTVXinEMC;
  }
  return -1;
}

template <typename T>
uint8_t setEventSelectionBit(T const& collision)
{

  uint8_t bit = 0;

  if (collision.sel8()) {
    SETBIT(bit, JCollisionSel::sel8);
  }
  if (collision.sel7()) {
    SETBIT(bit, JCollisionSel::sel7);
  }
  if (collision.alias_bit(kTVXinEMC)) {
    SETBIT(bit, JCollisionSel::kTVXinEMC);
  }

  return bit;
}

enum JTrigSelCh {
  charged = 0
};

enum JTrigSelFull {
  fullHigh = 0,
  fullLow = 1,
  neutralHigh = 2,
  neutralLow = 3,
  gammaHighEMCAL = 4,
  gammaLowEMCAL = 5,
  gammaHighDCAL = 6,
  gammaLowDCAL = 7
};

template <typename T>
bool selectChargedTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.chargedTriggerSel() >> triggerSelection) & 1;
}

int initialiseChargedTriggerSelection(std::string triggerSelection)
{
  if (triggerSelection == "charged") {
    return JTrigSelCh::charged;
  }
  return -1;
}

template <typename T>
uint8_t setChargedTriggerSelectionBit(T const& collision)
{

  uint8_t bit = 0;
  if (collision.hasJetChHighPt()) {
    SETBIT(bit, JTrigSelCh::charged);
  }
  return bit;
}

template <typename T>
bool selectFullTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.fullTriggerSel() >> triggerSelection) & 1;
}

int initialiseFullTriggerSelection(std::string triggerSelection)
{
  if (triggerSelection == "fullHigh") {
    return JTrigSelFull::fullHigh;
  } else if (triggerSelection == "fullLow") {
    return JTrigSelFull::fullLow;
  } else if (triggerSelection == "neutralHigh") {
    return JTrigSelFull::neutralHigh;
  } else if (triggerSelection == "neutralLow") {
    return JTrigSelFull::neutralLow;
  } else if (triggerSelection == "gammaHighEMCAL") {
    return JTrigSelFull::gammaHighEMCAL;
  } else if (triggerSelection == "gammaLowEMCAL") {
    return JTrigSelFull::gammaLowEMCAL;
  } else if (triggerSelection == "gammaHighDCAL") {
    return JTrigSelFull::gammaHighDCAL;
  } else if (triggerSelection == "gammaLowDCAL") {
    return JTrigSelFull::gammaLowDCAL;
  }
  return -1;
}

template <typename T>
uint8_t setFullTriggerSelectionBit(T const& collision)
{
  uint8_t bit = 0;
  if (collision.hasJetFullHighPt()) {
    SETBIT(bit, JTrigSelFull::fullHigh);
  }
  if (collision.hasJetFullLowPt()) {
    SETBIT(bit, JTrigSelFull::fullLow);
  }
  if (collision.hasJetNeutralHighPt()) {
    SETBIT(bit, JTrigSelFull::neutralHigh);
  }
  if (collision.hasJetNeutralLowPt()) {
    SETBIT(bit, JTrigSelFull::neutralLow);
  }
  if (collision.hasGammaHighPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaHighEMCAL);
  }
  if (collision.hasGammaLowPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaLowEMCAL);
  }
  if (collision.hasGammaHighPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaHighDCAL);
  }
  if (collision.hasGammaLowPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaLowDCAL);
  }
  return bit;
}

enum JTrackSel {
  globalTrack = 0,
  qualityTrack = 1,
  hybridTrack = 2
};

template <typename T>
bool selectTrack(T const& track, int trackSelection)
{
  if (trackSelection == -1) {
    return true;
  }
  return (track.trackSel() >> trackSelection) & 1;
}

int initialiseTrackSelection(std::string trackSelection)
{
  if (trackSelection == "globalTracks") {
    return JTrackSel::globalTrack;
  } else if (trackSelection == "QualityTracks") {
    return JTrackSel::qualityTrack;
  } else if (trackSelection == "hybridTracksJE") {
    return JTrackSel::hybridTrack;
  }
  return -1;
}

template <typename T>
uint8_t setTrackSelectionBit(T const& track)
{

  uint8_t bit = 0;

  if (track.isGlobalTrackWoPtEta()) {
    SETBIT(bit, JTrackSel::globalTrack);
  }
  if (track.isQualityTrack()) {
    SETBIT(bit, JTrackSel::qualityTrack);
  }
  if (track.trackCutFlagFb5()) {
    SETBIT(bit, JTrackSel::hybridTrack);
  }

  return bit;
}

template <typename T>
float trackEnergy(T const& track, float mass = JetFinder::mPion)
{
  return std::sqrt((track.p() * track.p()) + (mass * mass));
}

} // namespace JetDerivedDataUtilities

// namespace JetDerivedDataUtilities

#endif // PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_