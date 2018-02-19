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
/// \file electromagnetic/TestEm7/TestEm7.cc
/// \brief Main program of the electromagnetic/TestEm7 example
//
// $Id: TestEm7.cc 82280 2014-06-13 14:45:31Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "EventAction.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
 #include "G4UIExecutive.hh"
#endif

#ifndef __WITHOUT_ROOT__
#include "TreeManager.hh"
#endif

#include "StepMax.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
int main(int argc,char** argv) {

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  G4int nThreads = G4Threading::G4GetNumberOfCores();
  if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
  runManager->SetNumberOfThreads(nThreads);  
#else
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  G4RunManager* runManager = new G4RunManager;
#endif

#ifndef __WITHOUT_ROOT__  
  TreeManager *treeManager=TreeManager::Instance();
#endif  
  
  //set mandatory initialization classes
  //
  DetectorConstruction*   det  = new DetectorConstruction();
  PhysicsList*            phys = new PhysicsList();
  
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(phys);

  runManager->SetUserInitialization(new ActionInitialization(det));    
  
  // //set user action classes
  // //
  // PrimaryGeneratorAction* kin   = new PrimaryGeneratorAction(det);
  // RunAction*              run   = new RunAction(det,phys,kin);
  // EventAction*            ev    = new EventAction();
  // TrackingAction*         track = new TrackingAction(det,run);
  // SteppingAction*         step  = new SteppingAction(det,run,ev);
  
  // runManager->SetUserAction(kin);
  // runManager->SetUserAction(ev); 
  // runManager->SetUserAction(run); 
  // runManager->SetUserAction(track);  
  // runManager->SetUserAction(step);

  // G4double stepMax = phys->GetStepMaxProcess()->GetMaxStep();
  // G4cout << " max step: " << stepMax/CLHEP::mm << G4endl;

  // Get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc!=1)
    {
      // batch mode
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);    
    }
  else           //define visualization and UI terminal for interactive mode
    {
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
      
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      UI->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
#endif
      
#ifdef G4VIS_USE
      delete visManager;
#endif
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
