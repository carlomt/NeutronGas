//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: EventAction.cc 66892 2013-01-17 10:57:59Z gunter $
//

#include "EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleTable.hh"
#include "G4UnitsTable.hh" 
#include "globals.hh"

//#include "g4root.hh"

#ifndef __WITHOUT_ROOT__
#include "TreeManager.hh"
#endif

EventAction::EventAction()
  : G4UserEventAction(),
    fTotalPrimaryTrackLength(0.)
{
}

EventAction::~EventAction()
{;}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
    fTotalPrimaryTrackLength = 0.;
#ifndef __WITHOUT_ROOT__
    TreeManager::Instance()->EventNumber = evt->GetEventID();
#endif
    // auto AnalysisManager = G4AnalysisManager::Instance();

    // Default values (to be reset via /analysis/h1/set command)               
    // G4int nbins = 100;
    // G4double vmin = 0.;
    // G4double vmax = 100.;
    
    // G4int ih = AnalysisManager->CreateH2("h", "Edep (MeV/mm^2) transverse plane",
    // 					  nbins, vmin, vmax,
    // 					  nbins, vmin, vmax );
    // G4cout << "2d histo id: "<< ih <<G4endl;
    // AnalysisManager->SetH2Activation(ih, true);
}

void EventAction::EndOfEventAction(const G4Event*)
{
  // auto eventID = event->GetEventID();
  // G4cout << "---> End of event: " << eventID << G4endl;
  // G4cout << "EventAction::EndOfEventAction total track length: \t"
  // 	 << fTotalPrimaryTrackLength/CLHEP::mm << " mm" << G4endl;
}

