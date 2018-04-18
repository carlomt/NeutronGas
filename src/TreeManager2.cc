#include "G4AutoLock.hh"

#include "TreeManager2.hh"
#include "TreeManagerMessenger2.hh"
#include <vector>

#include "G4Threading.hh"

std::string    TreeManager2::fName = "tm-out";
bool TreeManager2::fIsActive = false;

namespace { 
  //Mutex to acquire access to singleton instance check/creation
  G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
  //Mutex to acquire accss to histograms creation/access
  //It is also used to control all operations related to histos 
  //File writing and check analysis
  G4Mutex dataManipulationMutex = G4MUTEX_INITIALIZER;
}

TreeManager2* TreeManager2::Instance()
{
  static G4ThreadLocalSingleton<TreeManager2> theInstance;
  return theInstance.Instance();
}


TreeManager2::TreeManager2() :
  fFile(0)//,
{
  messenger = new TreeManagerMessenger2(this);

  //Booking of histograms has to be protected.
  //In addition there are issues with ROOT that is 
  //heavily non thread-safe. In particular I/O related operations
  //are not thread safe. To avoid problems let's mutex everything
  //here
  G4AutoLock l(&dataManipulationMutex);


  if(G4Threading::G4GetThreadId()>=0)
    {
  OpenFile();
  
  fTree = new TTree("tree2","name2");

  fTree->Branch("Run",&Run,"Run/I");
  fTree->Branch("EventNumber", &EventNumber, "EventNumber/I");
  fTree->Branch("PDGencoding", &PDGencoding, "PDGencoding/I");  
  fTree->Branch("ParticleName",&ParticleName);
  fTree->Branch("ParticleType",&ParticleType);
  fTree->Branch("Ek",&Ek,"Ek/D");
  fTree->Branch("A",&A,"A/I");
  fTree->Branch("Z",&Z,"Z/I");
  fTree->Branch("x1",&x1,"x1/D");
  fTree->Branch("y1",&y1,"y1/D");
  fTree->Branch("z1",&z1,"z1/D");
  fTree->Branch("x2",&x2,"x2/D");
  fTree->Branch("y2",&y2,"y2/D");
  fTree->Branch("z2",&z2,"z2/D");
  fTree->Branch("px",&px,"px/D");
  fTree->Branch("py",&py,"py/D");
  fTree->Branch("pz",&pz,"pz/D");  
  fTree->Branch("trackID",&trackID,"trackID/I");
  fTree->Branch("parentID",&parentID,"parentID/I");
  
  fTree->Branch("EDeposited", &EDeposited,"EDeposited/D");          
  // fTree->Branch("DepositionX",&DepositionX);        
  // fTree->Branch("DepositionY",&DepositionY);        
  // fTree->Branch("DepositionZ",&DepositionZ);        
  
  this->Clear();
    }
}

TreeManager2::~TreeManager2()
{
}

Int_t TreeManager2::Fill()
{
    if(G4Threading::G4GetThreadId()>=0)
      fTree->Fill();
}

void TreeManager2::Print()
{
  if(fTree)
    {
      fTree->Print();
    }
  else
    {
      std::cout<<"ERROR TreeManager2::Print"<<std::endl;
    }
}

void TreeManager2::OpenFile()
{
  std::string TLfName =  fName + "_"+std::to_string(G4Threading::G4GetThreadId()) +".root"; ;  
  G4cout<<"... open Root tree file : "<<TLfName;
  fFile = new TFile(TLfName.c_str(),"RECREATE");
  G4cout<<" - done"<<G4endl;
}

Bool_t TreeManager2::Write()
{
  G4AutoLock l(&dataManipulationMutex);
  if(G4Threading::G4GetThreadId()<0)
    return false;
  if (!fFile->IsOpen())
    {
      G4Exception("TreeManager2::Write()","Hadr02",FatalException,
		  "Trying to write on a ROOT file which is not open");
      return false;
    }
  G4cout<<"... trying to write Root file : "<<fName;
  fFile->cd();
      
  fTree->Write();
  fFile->Close();
  G4cout<<" - done"<<G4endl;
  return true;
}

void TreeManager2::Clear()
{
  // nSec = -1;
  // trackID = -1;
  // parentID = -1;
  // PDGencoding = -99;
  // Ek = -1.;
  // ParticleName = "";
  // ParticleType = "";
  // A = -1;
  // Z = -1;
  // px = -99.;
  // py = -99.;
  // pz = -99.;
  // EDeposited=-99;  
  // EDeposited.clear();     
  // DepositionX.clear();    
  // DepositionY.clear();    
  // DepositionZ.clear();    
  
}
