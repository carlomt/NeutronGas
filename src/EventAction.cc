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

EventAction::EventAction()
  : G4UserEventAction(),
    fTotalPrimaryTrackLength(0.)
{
}

EventAction::~EventAction()
{;}

void EventAction::BeginOfEventAction(const G4Event*)
{
    fTotalPrimaryTrackLength = 0.;
}

void EventAction::EndOfEventAction(const G4Event*)
{
  // auto eventID = event->GetEventID();
  // G4cout << "---> End of event: " << eventID << G4endl;
  // G4cout << "EventAction::EndOfEventAction total track length: \t"
  // 	 << fTotalPrimaryTrackLength/CLHEP::mm << " mm" << G4endl;
}

