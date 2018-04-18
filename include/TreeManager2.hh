#ifndef TreeManager2_h
#define TreeManager2_h

#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "TreeManagerMessenger2.hh"

#include "G4ThreadLocalSingleton.hh"

class TreeManager2
{
  friend class G4ThreadLocalSingleton<TreeManager2>;
  
public:
  static  TreeManager2* Instance();
  virtual ~TreeManager2();
  
  inline TTree* GetTree() const { return fTree; };
  inline std::string GetOutputFileName() const { return fName;};
  Int_t Fill();
  Bool_t Write();
  void Clear();
  void SetOutputFilename(const std::string filename) {fName = filename; fIsActive=true;};
  void OpenFile();
  void Print();
  inline bool IsActive(){return fIsActive;};

  // Declaration of data types for the DataTree
  Int_t	Run;
  Int_t EventNumber;
  Int_t PDGencoding;
  TString ParticleName;
  TString ParticleType;
  Double_t Ek;
  Int_t A;
  Int_t Z;
  Double_t px;
  Double_t py;
  Double_t pz;
  Double_t x1;
  Double_t y1;
  Double_t z1;
  Double_t x2;
  Double_t y2;
  Double_t z2;  
  Int_t trackID;
  Int_t parentID;
  Double_t EDeposited;
  
private:
  TreeManager2();
    
  TTree*                fTree;
  TFile*                fFile;
  TreeManagerMessenger2* messenger;
  static std::string fName;  
  static bool        fIsActive;  
};

#endif
