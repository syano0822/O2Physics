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

/// \file match-mft-mch-data.cxx
/// \brief Produce MFT–MCH matching samples and evaluate frozen matching models for data.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#include "Common/Core/MftMchMatchingWorkflow.h"

#include <Framework/runDataProcessing.h>

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::analysis::mftmch::makeMatchingWorkflow(cfgc);
}
