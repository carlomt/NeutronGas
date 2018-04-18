#ifndef TreeManagerMessenger2_h
#define TreeManagerMessenger2_h

#include "G4UImessenger.hh"

class TreeManager2;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

class TreeManagerMessenger2: public G4UImessenger {

public:
  TreeManagerMessenger2(TreeManager2*);
  virtual ~TreeManagerMessenger2();
  void SetNewValue(G4UIcommand*, G4String);

private:
  TreeManager2* treeManager;
  G4UIdirectory* treeManager2Dir;
  G4UIcmdWithAString* selectFileName;
  
};

#endif
