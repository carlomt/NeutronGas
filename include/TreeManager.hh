#ifndef TreeManager_h
#define TreeManager_h

#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "TreeManagerMessenger.hh"

#include "G4ThreadLocalSingleton.hh"

class TreeManager
{
  // private:
  // static const Int_t  MaxNum = 2000000;
  // static const Int_t  MaxNumPro = 1000;
  friend class G4ThreadLocalSingleton<TreeManager>;

public:
  static  TreeManager* Instance();
  virtual ~TreeManager();

  inline TTree* GetTree() const { return fTree; };
  inline std::string GetOutputFileName() const { return fName;};
  Int_t Fill();
  Bool_t Write();
  void Clear();
  void SetOutputFilename(const std::string filename) {fName = filename; fIsActive=true;};
  void SetPhysicsList(const std::string physicsList) {fPhysicsList = physicsList;};
  void OpenFile();
  void Print();
  inline bool IsActive(){return fIsActive;};
  // inline bool IsActive(){return messenger->IsActive();;};  
  // Declaration of data types for the DataTree
  // TString G4version;
  Int_t	Run;
  Int_t EventNumber;
  // TString PrimaryName;
  // Double_t PrimaryE;
  // //  Double_t xsec;
  Int_t nSec;
  Int_t PDGencoding;
  TString ParticleName;
  TString ParticleType;
  Double_t Ek;
  Int_t A;
  Int_t Z;
  // Double_t ExcitationE;
  // Bool_t IsShortLived;  
  // Bool_t IsStable;
  // Double_t LifeTime;
  // Int_t NDecayChannels;
  // Double_t theta;
  // Double_t phi;
  TString PhysicsList;
  // Double_t b;

  // Double_t hot_px;
  // Double_t hot_py;
  // Double_t hot_pz;
  // Double_t hot_Eecc;
  // Double_t hot_EeccCorr;
  // Int_t hot_A;
  // Int_t hot_Z;

  Double_t px;
  Double_t py;
  Double_t pz;
  Double_t x;
  Double_t y;
  Double_t z;
  // Double_t cm_px;
  // Double_t cm_py;
  // Double_t cm_pz;
  // Double_t beta;

  Int_t trackID;
  Int_t parentID;
  
private:
  TreeManager();
  // static TreeManager* fInstance;
    
  TTree*                fTree;
  TFile*                fFile;
  TreeManagerMessenger* messenger;
  // static G4ThreadLocal std::string fName;  
  // static G4ThreadLocal bool        fIsActive;
  static std::string fName;  
  static bool        fIsActive;  
  static std::string fPhysicsList;
};

#endif
