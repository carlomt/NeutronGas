#include "TreeManagerMessenger2.hh"
#include "TreeManager2.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

TreeManagerMessenger2::TreeManagerMessenger2(TreeManager2* aTreeManager)
  : treeManager(aTreeManager)
{
  treeManager2Dir = new G4UIdirectory("/treeManager2/");
  selectFileName = new G4UIcmdWithAString("/treeManager2/setFileName",this);
  selectFileName->SetParameterName("setIonInelasticProcess",false);
  selectFileName->SetToBeBroadcasted(false);
}

TreeManagerMessenger2::~TreeManagerMessenger2()
{
  delete selectFileName;
  delete treeManager2Dir;
}

void TreeManagerMessenger2::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == selectFileName )
    {
      G4cout<<"TreeManagerMessenger2: setting output filename "<<newValue<<G4endl;
      treeManager->SetOutputFilename(newValue);
      // TreeManager::Instance()->PhysicsList = newValue;
    } 
  else
    {
      G4cerr << "***** Command is not found !!! " << newValue << G4endl;
    }

}
