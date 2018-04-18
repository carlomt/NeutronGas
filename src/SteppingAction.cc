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
/// \file electromagnetic/TestEm7/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 101250 2016-11-10 08:54:02Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "Randomize.hh"
#include "EventAction.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "g4root.hh"
#include "G4Ions.hh"

#ifndef __WITHOUT_ROOT__
#include "TreeManager.hh"
#include "TreeManager2.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* RuAct, EventAction* EvAct)
  :G4UserSteppingAction(),fDetector(det), fRunAction(RuAct), fEventAction(EvAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  
  // G4double edep = step->GetTotalEnergyDeposit();
  // if (edep <= 0.) return;

  G4StepPoint* prePoint  = step->GetPreStepPoint();
  G4StepPoint* postPoint = step->GetPostStepPoint();



  
  // G4int copyNb = prePoint->GetTouchableHandle()->GetCopyNumber();
  // if (copyNb > 0) { fRunAction->FillTallyEdep(copyNb-1, edep); }

  // G4double niel = step->GetNonIonizingEnergyDeposit();
  // fRunAction->FillEdep(edep, niel);
  
  // if (step->GetTrack()->GetTrackID() == 1) {
  //   fRunAction->AddPrimaryStep();
    /*
    G4cout << step->GetTrack()->GetMaterial()->GetName()
           << "  E1= " << step->GetPreStepPoint()->GetKineticEnergy()
           << "  E2= " << step->GetPostStepPoint()->GetKineticEnergy()
           << " Edep= " << edep 
           << " Q= " << step->GetTrack()->GetDynamicParticle()->GetCharge()
           << " Qp= " << step->GetPostStepPoint()->GetCharge()
           << G4endl;
    */
    // G4cout << "SteppingAction::UserSteppingAction step length     \t"
    // 	   // << G4BestUnit(step->GetStepLength(),"Length") << G4endl;
    // 	   << step->GetStepLength()/CLHEP::mm << " mm " << G4endl;

    // fRunAction->AddThisTotalRange(step->GetStepLength());
    // fEventAction->AddPrimaryTrackLength(step->GetStepLength());
  // } 

  //Bragg curve
  //        
  G4double xmax = fDetector->GetAbsorSizeX();
  G4double ymax = fDetector->GetAbsorSizeYZ();
  G4double zmax = fDetector->GetAbsorSizeYZ();
   
  G4double x1 = prePoint->GetPosition().x() ;// + xmax*0.5;
  G4double x2 = postPoint->GetPosition().x() ;// + xmax*0.5;
  G4double y1 = prePoint->GetPosition().y() ;//+ ymax*0.5;
  G4double y2 = postPoint->GetPosition().y();//+ ymax*0.5;
  G4double z1 = prePoint->GetPosition().z() ;//+ zmax*0.5;
  G4double z2 = postPoint->GetPosition().z();//+ zmax*0.5;

#ifndef __WITHOUT_ROOT__
  const G4Track* thisSteppingTrack = step->GetTrack();
	TreeManager2* treeManager2 = TreeManager2::Instance();
	treeManager2->trackID = thisSteppingTrack->GetTrackID();
	treeManager2->parentID = thisSteppingTrack->GetParentID();
	treeManager2->PDGencoding = thisSteppingTrack->GetDefinition()->GetPDGEncoding();
	treeManager2->ParticleName = thisSteppingTrack->GetDefinition()->GetParticleName();
	treeManager2->ParticleType = thisSteppingTrack->GetDefinition()->GetParticleType();
	treeManager2->A = thisSteppingTrack->GetDefinition()->GetAtomicMass();
	treeManager2->Z = thisSteppingTrack->GetDefinition()->GetAtomicNumber();
	treeManager2->x2 = x2;	
	treeManager2->y2 = y2;	
	treeManager2->z2 = z2;
	treeManager2->x1 = x1;	
	treeManager2->y1 = y1;	
	treeManager2->z1 = z1;			
	treeManager2->px = thisSteppingTrack->GetMomentum().x()/CLHEP::MeV;	
	treeManager2->py = thisSteppingTrack->GetMomentum().y()/CLHEP::MeV;	
	treeManager2->pz = thisSteppingTrack->GetMomentum().z()/CLHEP::MeV;
	treeManager2->Ek = thisSteppingTrack->GetKineticEnergy()/CLHEP::MeV;
	treeManager2->EDeposited = step->GetTotalEnergyDeposit()/CLHEP::MeV;
	treeManager2->Fill();
	//	treeManager2->Clear();
#endif
  
  // if(x1 >= 0.0 && x2 <= xmax)
    {  
      G4double x  = x1 + G4UniformRand()*(x2-x1);
      // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      // analysisManager->FillH1(1, x, edep);  
      // analysisManager->FillH1(2, x, edep);

      // if (step->GetTrack()->GetTrackID() == 1)
      // 	{
      // if((y1 >= 0.0 && y2 <= ymax) && (z1 >= 0.0 && z2 <= zmax))
      	{
	  G4double y  = y1 + G4UniformRand()*(y2-y1);
	  G4double z  = z1 + G4UniformRand()*(z2-z1);
	  // analysisManager->FillH2(0, y1/CLHEP::mm, z1/CLHEP::mm, edep/CLHEP::keV);
	  // analysisManager->FillH2(1, x1/CLHEP::mm, y1/CLHEP::mm, edep/CLHEP::keV);
      // 	  G4double dx  = (x2-x1);
      // 	  G4double dy  = (y2-y1);
      // 	  G4double dz  = (z2-z1);
      // 	  G4double ss = std::sqrt(dx*dx + dy*dy + dz*dz);
      // 	  G4cout << "SteppingAction::UserSteppingAction step manually calc \t"
      // 		 << ss/CLHEP::mm << " mm " << G4endl;
      // 	  G4cout <<" prePoint:  ("<<x1<<", "<<y1<<", "<<z1<<")" <<G4endl;
      // 	  G4cout <<" postPoint: ("<<x2<<", "<<y2<<", "<<z2<<")" <<G4endl;
      // // 	  fRunAction->AddThisTotalRange(ss);
      // 	}
	}
    }

    //secondaries
    //
    const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
#ifndef __WITHOUT_ROOT__
      TreeManager::Instance()->nSec = secondary->size();
#endif
    for (size_t lp=0; lp<(*secondary).size(); lp++)
      {
	const G4Track* thisTrack = (*secondary)[lp];
	G4ParticleDefinition* particle = thisTrack->GetDefinition();
	G4int PDGenc =  particle->GetPDGEncoding();
	G4int trackID =  thisTrack->GetTrackID();
	G4int parentID = thisTrack->GetParentID();
	G4String name   = particle->GetParticleName();
	G4String type   = particle->GetParticleType();
	G4double energy = (*secondary)[lp]->GetKineticEnergy()/CLHEP::MeV;
	G4int A = particle->GetAtomicMass();
	G4int Z =  particle->GetAtomicNumber();
	G4LorentzVector fourMomentum((*secondary)[lp]->GetMomentum(), (*secondary)[lp]->GetTotalEnergy());
	G4double px = fourMomentum.x()/CLHEP::MeV;
	G4double py = fourMomentum.y()/CLHEP::MeV;
	G4double pz = fourMomentum.z()/CLHEP::MeV;
	if (type == "nucleus")
        {
	    G4double ExcitationE = ((G4Ions*)particle)->GetExcitationEnergy()/CLHEP::eV;
        }
	// G4AnalysisManager* man = G4AnalysisManager::Instance();
	// man->FillNtupleIColumn(0, trackID);
	// man->FillNtupleIColumn(1, parentID);
	// man->FillNtupleIColumn(2, -99);
	// man->FillNtupleDColumn(3, energy);
	// man->AddNtupleRow();
#ifndef __WITHOUT_ROOT__
	TreeManager* treeManager = TreeManager::Instance();
	treeManager->trackID = trackID;
	treeManager->parentID = parentID;
	treeManager->PDGencoding = PDGenc;
	treeManager->ParticleName = name;
	treeManager->ParticleType = type;
	treeManager->A = A;
	treeManager->Z = Z;
	treeManager->x = x2;	
	treeManager->y = y2;	
	treeManager->z = z2;		
	treeManager->px = px;	
	treeManager->py = py;	
	treeManager->pz = pz;	
	treeManager->Ek = energy;
	treeManager->Fill();
	treeManager->Clear();
#endif
       }//end loop on secondaries
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


