#include "G4AutoLock.hh"

#include "TreeManager.hh"
#include "TreeManagerMessenger.hh"
#include <vector>

#include "G4Threading.hh"

// TreeManager* TreeManager::fInstance = NULL;
// G4ThreadLocal std::string    TreeManager::fName = "tm-out";
// G4ThreadLocal bool TreeManager::fIsActive = false;
std::string    TreeManager::fName = "tm-out";
bool TreeManager::fIsActive = false;
std::string    TreeManager::fPhysicsList = "NONE";

namespace { 
  //Mutex to acquire access to singleton instance check/creation
  G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
  //Mutex to acquire accss to histograms creation/access
  //It is also used to control all operations related to histos 
  //File writing and check analysis
  G4Mutex dataManipulationMutex = G4MUTEX_INITIALIZER;
}

TreeManager* TreeManager::Instance()
{
  static G4ThreadLocalSingleton<TreeManager> theInstance;
  return theInstance.Instance();
  
    // G4AutoLock l(&instanceMutex);
    // if (fInstance == NULL)
    //   fInstance = new TreeManager();
    // return fInstance;
}


TreeManager::TreeManager() :
  fFile(0)//,
  // fName("tm-out")
  // fIsActive(false)
{
  messenger = new TreeManagerMessenger(this);
  
  // if(fInstance) 
  //   {
  //     return;
  //   }

  // fName="out";

  //Booking of histograms has to be protected.
  //In addition there are issues with ROOT that is 
  //heavily non thread-safe. In particular I/O related operations
  //are not thread safe. To avoid problems let's mutex everything
  //here
  G4AutoLock l(&dataManipulationMutex);


  if(G4Threading::G4GetThreadId()>=0)
    {
  OpenFile();
  
  // fInstance = this;
  fTree = new TTree("tree","name");

  // fTree->Branch("G4version",&G4version);
  fTree->Branch("Run",&Run,"Run/I");
  fTree->Branch("EventNumber", &EventNumber, "EventNumber/I");
  // fTree->Branch("PrimaryName",&PrimaryName);  
  // fTree->Branch("PrimaryE",&PrimaryE,"PrimaryE/D");
  // //  fTree->Branch("xsec",&xsec,"xsec/D");  
  fTree->Branch("nSec", &nSec, "nSec/I");
  fTree->Branch("PDGencoding", &PDGencoding, "PDGencoding/I");  
  fTree->Branch("ParticleName",&ParticleName);
  fTree->Branch("ParticleType",&ParticleType);
  fTree->Branch("Ek",&Ek,"Ek/D");
  fTree->Branch("A",&A,"A/I");
  fTree->Branch("Z",&Z,"Z/I");
  // fTree->Branch("ExcitationE",&ExcitationE,"ExcitationE/D");
  // fTree->Branch("IsShortLived",&IsShortLived,"IsShortLived/O");  
  // fTree->Branch("IsStable",&IsStable,"IsStable/O");
  // fTree->Branch("LifeTime",&LifeTime,"LifeTime/D");
  // fTree->Branch("NDecayChannels",&NDecayChannels,"NDecayChannels/I");
  // fTree->Branch("theta",&theta,"theta/D");
  // fTree->Branch("phi",&phi,"phi/D");
  fTree->Branch("PhysicsList",&PhysicsList);
  // fTree->Branch("b",&b,"b/D");
  // fTree->Branch("hot_px",&hot_px,"hot_px/D");
  // fTree->Branch("hot_py",&hot_py,"hot_py/D");
  // fTree->Branch("hot_pz",&hot_pz,"hot_pz/D");
  // fTree->Branch("hot_Eecc",&hot_Eecc,"hot_Eecc/D");
  // fTree->Branch("hot_EeccCorr",&hot_EeccCorr,"hot_EeccCorr/D");
  // fTree->Branch("hot_A",&hot_A,"hot_A/D");
  // fTree->Branch("hot_Z",&hot_Z,"hot_Z/D");
  fTree->Branch("x",&x,"x/D");
  fTree->Branch("y",&y,"y/D");
  fTree->Branch("z",&z,"z/D");
  fTree->Branch("px",&px,"px/D");
  fTree->Branch("py",&py,"py/D");
  fTree->Branch("pz",&pz,"pz/D");  
  // fTree->Branch("cm_px",&cm_px,"cm_px/D");
  // fTree->Branch("cm_py",&cm_py,"cm_py/D");
  // fTree->Branch("cm_pz",&cm_pz,"cm_pz/D");
  // fTree->Branch("beta",&beta,"beta/D");
  fTree->Branch("trackID",&trackID,"trackID/I");
  fTree->Branch("parentID",&parentID,"parentID/I");  
  this->Clear();
    }
}

TreeManager::~TreeManager()
{
  // G4AutoLock l(&dataManipulationMutex);
  // if (fFile->IsOpen())
  //   fFile->Close();
  // //No need to mutex, this is a real singleton.
  // if(messenger)
  //   delete messenger;
  // if(fTree)
  //   delete fTree;
  // if(fFile)
  //   delete fFile;

  // G4cout<<"TreeManager::~TreeManager G4Threading::G4GetThreadId():"<<G4Threading::G4GetThreadId()<<G4endl;  
}

Int_t TreeManager::Fill()
{
  // G4cout<<"\n\nTreeManager::Fill G4Threading::G4GetThreadId():"<<G4Threading::G4GetThreadId()<<G4endl;  
    // G4AutoLock l(&dataManipulationMutex);
    PhysicsList = fPhysicsList;
    if(G4Threading::G4GetThreadId()>=0)
      fTree->Fill();
}

void TreeManager::Print()
{
  if(fTree)
    {
      fTree->Print();
    }
  else
    {
      std::cout<<"ERROR TreeManager::Print"<<std::endl;
    }
}

void TreeManager::OpenFile()
{
  // G4AutoLock l(&dataManipulationMutex);
  
  // std::string filename = fName + "_"+std::to_string(G4Threading::G4GetThreadId()) +".root"; ;
  // fName = messenger->GetOutFilename()+ "_"+std::to_string(G4Threading::G4GetThreadId()) +".root"; ;
  std::string TLfName =  fName + "_"+std::to_string(G4Threading::G4GetThreadId()) +".root"; ;  
  G4cout<<"... open Root tree file : "<<TLfName;
  fFile = new TFile(TLfName.c_str(),"RECREATE");
  G4cout<<" - done"<<G4endl;
}

Bool_t TreeManager::Write()
{
  // TString filename = fName;
  // filename += ".root";
  // TFile* file = new TFile(filename,"RECREATE");
  G4AutoLock l(&dataManipulationMutex);
  if(G4Threading::G4GetThreadId()<0)
    return false;
  // if(fName=="NONE") //file not created at all: e.g. for a vis-only execution
  // if(this->IsActive())
  //   {
  //     G4cout<<"TreeManager::Write() Warning: not writing"<<G4endl;
  //     return false;
  //   }

  if (!fFile->IsOpen())
    {
      G4Exception("TreeManager::Write()","Hadr02",FatalException,
		  "Trying to write on a ROOT file which is not open");
      return false;
    }
  // std::string filename = fName + "_"+std::to_string(G4Threading::G4GetThreadId()) +".root"; ;  
  G4cout<<"... trying to write Root file : "<<fName;
  fFile->cd();
      
  fTree->Write();
  // fFile->Write();
  // G4cout<<" - done"<<G4endl;
  
  // G4cout<<"... write Root file : "<<filename;
  fFile->Close();
  G4cout<<" - done"<<G4endl;
  
  //  delete file;
  return true;
}

void TreeManager::Clear()
{
  nSec = -1;
  // Run=0;
  // ExcitationE=-1;
  trackID = -1;
  parentID = -1;
  PDGencoding = -99;
  Ek = -1.;
  ParticleName = "";
  ParticleType = "";
  A = -1;
  Z = -1;
  px = -99.;
  py = -99.;
  pz = -99.;
  
}

// void TreeManager::SetOutputFilename(const std::string fname)
// {
//   // G4AutoLock l(&dataManipulationMutex);
//   // if(G4Threading::G4GetThreadId()>0)
//     // {
//       //
//   G4cout<<"TreeManager::SetOutputFilename G4Threading::G4GetThreadId():"<<G4Threading::G4GetThreadId()<<G4endl;
	
//       // fName = fname;

//       // fName += "_"+std::to_string(G4Threading::G4GetThreadId());
//       // fName +=".root";
//       // G4cout<<"... open Root tree file : "<<fName;
//       // fFile = new TFile(fName,"RECREATE");
//       // G4cout<<" - done"<<G4endl;
//     // }
//   // else
//   //   {
//   //     G4cout<<"TreeManager::SetOutputFilename G4Threading::G4GetThreadId():"<<G4Threading::G4GetThreadId()<<G4endl;
//   //   }

//     fIsActive = true;
//     fName = fname;
// }


