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
/// \file electromagnetic/TestEm7/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 101390 2016-11-16 14:13:20Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh" 

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
    fMagField(nullptr), 
    fLAbsor(nullptr),
    fLWorld(nullptr)
{
  // default parameter values
  fAbsorSizeX = fAbsorSizeYZ = 50*cm;
  fWorldSizeX = fWorldSizeYZ = 1.2*fAbsorSizeX;
  
  fTallyNumber = 0;
  for (G4int j=0; j<kMaxTally; j++) {
    fTallySize[j] = fTallyPosition[j] = G4ThreeVector(0.,0.,0.);
    fTallyMass[j]     = 0.;
    fLTally[j]        = nullptr; 
  }
    
  DefineMaterials();

  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //add new units for pressure and temperature
  //
  G4double torr = 1.33322e-3*bar;
  new G4UnitDefinition( "torr", "torr", "Pressure", torr);
  // //
  // // define Elements
  // //
  // G4Element* H = new G4Element("Hydrogen", "H", 1, 1.008*g/mole);
  // G4Element* N = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  // G4Element* O = new G4Element("Oxygen"  , "O", 8, 16.00*g/mole);
  // // G4Element* F = new G4Element("Fluorine","F", 9, 18.998*g/mole);
  // // G4Element* S = new G4Element("Sulfur", "S", 16, 32.065*g/mole);
  // //
  // // define Materials.
  // //
  G4double density, temperature, pressure;
  // G4int    ncomponents, natoms;
  // G4double fractionmass;
 
  // G4Material* H2O = 
  //   new G4Material("Water", density= 1.0*g/cm3, ncomponents=2);
  // H2O->AddElement(H, natoms=2);
  // H2O->AddElement(O, natoms=1);
  // H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  // // In this line both G4_WATER and Water_1.05 will be constructed
  // G4NistManager::Instance()->
  //   BuildMaterialWithNewDensity("Water_1.05","G4_WATER",1.05*g/cm3);

  // G4Material* Air = 
  //   new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  // Air->AddElement(N, fractionmass=0.7);
  // Air->AddElement(O, fractionmass=0.3);

  // density     = 1.e-5*g/cm3;
  // pressure    = 2.e-2*bar;
  // temperature = STP_Temperature;  // From PhysicalConstants.h .
  // G4Material* vac = new G4Material( "TechVacuum", density, 1, kStateGas, temperature, pressure );
  // vac->AddMaterial( Air, 1. );

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* vacuum = 
    new G4Material("Galactic", 1, 1.008*g/mole,density,
                   kStateGas,temperature,pressure);

  // temperature = 273*kelvin;
  // pressure    = 75*torr;
  // G4Material* SF6 = DefineSF6( pressure, temperature );
    
  //default materials
  fAbsorMaterial = vacuum;
  fWorldMaterial = vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::DefineSF6(const G4double pressure, const G4double temperature)
{
  G4Element* F = new G4Element("Fluorine","F", 9, 18.998*g/mole);
  G4Element* S = new G4Element("Sulfur", "S", 16, 32.065*g/mole);

  G4double STP_density = 6.43e-3 *g/cm3;
  G4double density     = STP_density * pressure/STP_Pressure * STP_Temperature/temperature;

  G4Material* SF6 = new G4Material( "SF6", density, 2, kStateGas, temperature, pressure );  
  SF6->AddElement(S, 1);
  SF6->AddElement(F, 6);

  fAbsorMaterial = SF6;
  
  return SF6;
}

G4Material* DetectorConstruction::DefineMixture(const G4double pressure, const G4double temperature)
{
  G4Element* H = new G4Element("Hydrogen", "H", 1,  1.008*g/mole);
  G4Element* C = new G4Element("Carbon",   "C", 6, 12.011*g/mole);
  G4Element* F = new G4Element("Fluorine", "F", 9, 18.998*g/mole);

  // densities at STP
  G4double rho_CF4 = 3.72 *kg/m3;
  G4double rho_CHF3 = 2.946 *kg/m3;
  G4double rho_C4H10 = 2.51 *kg/m3;

  G4double reference_pressure = 1.*bar;
  G4double reference_temperature = 288.15*kelvin;
  
  G4Material* CF4 = new G4Material( "Tetrafluoromethane", rho_CF4, 2, kStateGas, reference_temperature, reference_pressure );
  CF4->AddElement(C,  1);
  CF4->AddElement(F,  4);  

  G4Material* CHF3 = new G4Material( "Fluoroform",rho_CHF3, 3, kStateGas, reference_temperature, reference_pressure );
  CHF3->AddElement(C,  1);
  CHF3->AddElement(H,  1);
  CHF3->AddElement(F,  3);  

  G4Material* C4H10 = new G4Material( "Isobutane",rho_C4H10, 2, kStateGas, reference_temperature, reference_pressure);
  C4H10->AddElement(C,  4);
  C4H10->AddElement(H, 10);

  G4double fraction_CF4 = 70.*perCent;
  G4double fraction_CHF3 = 28.*perCent;
  G4double fraction_C4H10 = 2.*perCent;
   
  // // scale densities with partial pressures and temperatures
  // G4double density_CF4     = fraction_CF4   * pressure/STP_Pressure * STP_Temperature/temperature;
  // G4double density_CHF3    = fraction_CHF3  * pressure/STP_Pressure * STP_Temperature/temperature;
  // G4double density_C4H10   = fraction_C4H10 * pressure/STP_Pressure * STP_Temperature/temperature;  

  // G4double density     = density_CF4  + density_CHF3 +density_C4H10;
  G4double density     = ( rho_CF4*fraction_CF4 + rho_CHF3*fraction_CHF3 + rho_C4H10*fraction_C4H10 )
    *pressure/reference_pressure * reference_temperature/temperature;
  
  G4Material* Mixture = new G4Material( "CF4-CHF3-C4H10", density, 3, kStateGas, temperature, pressure);
  Mixture->AddMaterial(CF4,   fraction_CF4);
  Mixture->AddMaterial(CHF3,  fraction_CHF3);
  Mixture->AddMaterial(C4H10, fraction_C4H10);

  fAbsorMaterial = Mixture;

  return Mixture;
}

G4Material* DetectorConstruction::DefineHeCF4(const G4double pressure, const G4double temperature, const G4double fraction_He, const G4double fraction_CF4)
{
  G4Element* He = new G4Element("Helium",  "He", 2, 4.0026*g/mole);
  G4Element* C  = new G4Element("Carbon",   "C", 6, 12.011*g/mole);
  G4Element* F  = new G4Element("Fluorine", "F", 9, 18.998*g/mole);

  // densities at STP
  G4double rho_CF4 = 3.72 *kg/m3;
  G4double rho_He  = 0.1786 *kg/m3;

  G4double reference_pressure = 1.*bar;
  G4double reference_temperature = 288.15*kelvin;

  G4double density     = ( rho_CF4*fraction_CF4 + rho_He*fraction_He)
    *pressure/reference_pressure * reference_temperature/temperature;

  G4Material* Mixture = new G4Material( "HeCF4", density, 3, kStateGas, temperature, pressure);
  Mixture->AddElement(He, fraction_He);
  Mixture->AddElement(C,  fraction_CF4*1./5);
  Mixture->AddElement(F,  fraction_CF4*4./5);

  fAbsorMaterial = Mixture;
    
  return Mixture;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // World
  //
  G4Box*
  sWorld = new G4Box("World",                                      //name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);   //dimensions

  fLWorld = new G4LogicalVolume(sWorld,                        //shape
                                fWorldMaterial,                //material
                                "World");                      //name

  G4VPhysicalVolume*                                   
  pWorld = new G4PVPlacement(0,                           //no rotation
                             G4ThreeVector(0.,0.,0.),     //at (0,0,0)
                             fLWorld,                     //logical volume
                             "World",                     //name
                             0,                           //mother  volume
                             false,                       //no boolean operation
                             0);                          //copy number

  // CECILIA

  G4NistManager* man = G4NistManager::Instance();

  
  // outer PVC box
  G4double Box1SizeX = 30.3 *cm;
  G4double Box1SizeY = 30.3 *cm;
  G4double Box1SizeZ = 8 *cm;
  G4double Box1ThicknessX = 2*cm;
  G4double Box1ThicknessY = 2*cm;
  G4double Box1ThicknessZ = 2*cm;
  G4Material* Box1Material  = man->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
  G4Box*  Box1 = new G4Box("Box1", Box1SizeX/2,Box1SizeY/2,Box1SizeZ/2);
  G4LogicalVolume* LBox1 = new G4LogicalVolume(Box1,Box1Material,"Box1");                     
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),LBox1, "Box1",fLWorld,false,0); 

  // gas between PVC box and GEM stuff
  G4double Box2SizeX = Box1SizeX-2*Box1ThicknessX;
  G4double Box2SizeY = Box1SizeY-2*Box1ThicknessY;
  G4double Box2SizeZ = Box1SizeZ-2*Box1ThicknessZ;
  G4Box*  Box2 = new G4Box("Box2", Box2SizeX/2,Box2SizeY/2,Box2SizeZ/2);
  G4LogicalVolume* LBox2 = new G4LogicalVolume(Box2,fAbsorMaterial,"Box2");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),LBox2, "Box2",LBox1,false,0);
  
  // GEM FR4 frame+PCB
  G4double Box3SizeX = 18*cm;
  G4double Box3SizeY = 18*cm;
  G4double Box3SizeZ = 2.5*cm;
  G4double Box3Z = -Box2SizeZ/2.+Box3SizeZ/2.;

  G4int ncomponents, natoms;
  G4String symbol;
  G4double a, z;
  G4double abundance, fractionmass;
  //Epoxy 
  G4double density = 1.2*g/cm3;
  G4Material* Epoxy = new G4Material("Epoxy" , density, ncomponents=2);
  G4Element* H = new G4Element("Hydrogen", "H", 1,  1.008*g/mole);
  G4Element* C = new G4Element("Carbon",   "C", 6, 12.011*g/mole);
  Epoxy->AddElement(H, natoms=2);
  Epoxy->AddElement(C, natoms=2);
  //FR4 (Glass + Epoxy)
  density = 1.86*g/cm3;
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);
  G4Material* FR4 = new G4Material("FR4"  , density, ncomponents=2);
  FR4->AddMaterial(SiO2, fractionmass=0.528);
  FR4->AddMaterial(Epoxy, fractionmass=0.472);
  G4Material* Box3Material  = FR4;
  G4Box*  Box3 = new G4Box("Box3", Box3SizeX/2,Box3SizeY/2,Box3SizeZ/2);
  G4LogicalVolume* LBox3 = new G4LogicalVolume(Box3,Box3Material,"Box3");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,Box3Z),LBox3, "Box3",LBox2,false,0);

  // GEM (for now single GEM)
  G4double Box4SizeX = 10 *cm;
  G4double Box4SizeY = 10 *cm;
  G4double Box4SizeZ = 0.5 *cm;
  G4double Box4Z = Box3Z+Box3SizeZ/2.+Box4SizeZ/2.;
  G4Material* Box4Material  = man->FindOrBuildMaterial("G4_KAPTON");
  G4Box*  Box4 = new G4Box("Box4", Box4SizeX/2,Box4SizeY/2,Box4SizeZ/2);
  G4LogicalVolume* LBox4 = new G4LogicalVolume(Box4,Box4Material,"Box4");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,Box4Z),LBox4, "Box5",LBox2,false,0);  


  /*G4Box*
  sAbsor = new G4Box("Absorber",                                 //name
                   fAbsorSizeX/2,fAbsorSizeYZ/2,fAbsorSizeYZ/2); //dimensions
                                                                 
  fLAbsor = new G4LogicalVolume(sAbsor,                   //shape
                                fAbsorMaterial,           //material
                                "Absorber");              //name

                              
  new G4PVPlacement(0,                           //no rotation
                    G4ThreeVector(0.,0.,0.),     //at (0,0,0)
                    fLAbsor,                     //logical volume
                    "Absorber",                  //name
                    fLWorld,                     //mother  volume
                    false,                       //no boolean operation
                    0);                          //copy number
  */
  //
  // Tallies (optional)
  //
  if (fTallyNumber > 0) {
    for (G4int j=0; j<fTallyNumber; ++j) {
            
       G4Box* sTally = new G4Box("Tally",
                   fTallySize[j].x()/2,fTallySize[j].y()/2,fTallySize[j].z()/2);
                      
       fLTally[j] = new G4LogicalVolume(sTally,fAbsorMaterial,"Tally");
           
       new G4PVPlacement(0,                        //no rotation
                         fTallyPosition[j],        //position
                         fLTally[j],               //logical volume
                         "Tally",                  //name
                         fLAbsor,                  //mother  volume
                         false,                    //no boolean operation
                         j+1);                     //copy number
       
      fTallyMass[j] = fTallySize[j].x()*fTallySize[j].y()*fTallySize[j].z()
               *(fAbsorMaterial->GetDensity());
    }               
  } 

  PrintParameters();
    
  //
  //always return the World volume
  //  
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() const
{
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(fAbsorSizeX,"Length")
         << " of " << fAbsorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  
  if (fTallyNumber > 0) {
    G4cout << "---> There are " << fTallyNumber << " tallies : " << G4endl;    
    for (G4int j=0; j<fTallyNumber; ++j) {
      G4cout << "fTally " << j << ": "
             << fAbsorMaterial->GetName()
             << ",  mass = " << G4BestUnit(fTallyMass[j],"Mass")
             << " size = "   << G4BestUnit(fTallySize[j],"Length")           
             << " position = " << G4BestUnit(fTallyPosition[j],"Length")
             << G4endl;
    }                 
    G4cout << "\n---------------------------------------------------------\n";
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double value)
{
  fAbsorSizeX = value; 
  fWorldSizeX = 1.2*fAbsorSizeX;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeYZ(G4double value)
{
  fAbsorSizeYZ = value; 
  fWorldSizeYZ = 1.2*fAbsorSizeYZ;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial && pttoMaterial != fAbsorMaterial) {
    // change target material everywhere
    fAbsorMaterial = pttoMaterial;
    for (G4int j=0; j<fTallyNumber; ++j) {
      if(fLTally[j]) { 
        fLTally[j]->SetMaterial(pttoMaterial); 
        fTallyMass[j] = fTallySize[j].x()*fTallySize[j].y()*fTallySize[j].z()
          *(pttoMaterial->GetDensity());
      }
    } 
    if(fLAbsor) {
      fLAbsor->SetMaterial(fAbsorMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial && pttoMaterial != fWorldMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLWorld) {
      fLWorld->SetMaterial(fAbsorMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if (fMagField) delete fMagField;        //delete the existing magn field  

  if (fieldValue!=0.)                        // create a new one if non nul
    {
      fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
      fieldMgr->SetDetectorField(fMagField);
      fieldMgr->CreateChordFinder(fMagField);
    }
   else
    {
      fMagField = nullptr;
      fieldMgr->SetDetectorField(fMagField);
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallyNumber(G4int value)
{
  if(value >= 0 && value <kMaxTally) {
    fTallyNumber = value;
  } else {
    G4cout << "### DetectorConstruction::SetTallyNumber WARNING: wrong tally "
           << "number " << value << " is ignored" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallySize(G4int j, const G4ThreeVector& value)
{
  if(j >= 0 && j < kMaxTally) {
    fTallySize[j] = value;
  } else {
    G4cout << "### DetectorConstruction::SetTallyNumber WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
  } 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallyPosition(G4int j, const G4ThreeVector& value)
{
  if(j >= 0 && j < kMaxTally) {
    fTallyPosition[j] = value; 
  } else {
    G4cout << "### DetectorConstruction::SetTallyPosition WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
  } 
}  

G4double DetectorConstruction::GetTallyMass(G4int j) const
{
  if(j >= 0 && j < kMaxTally) {
    return fTallyMass[j];
  } else {
    G4cout << "### DetectorConstruction::GetTallyMass WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
    return 0.0;
  } 
}

const G4LogicalVolume* DetectorConstruction::GetLogicalTally(G4int j) const 
{
  if(j >= 0 && j < kMaxTally) {
    return fLTally[j];
  } else {
    G4cout << "### DetectorConstruction::GetLOgicalTally WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
    return nullptr;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
