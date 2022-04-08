
#ifndef REST_APPLICATION_H
#define REST_APPLICATION_H

#include "SimulationManager.h"

class G4VisManager;
class G4UIExecutive;
class G4RunManager;

class Application {
   public:
    void Run(int argc, char** argv);
};

#endif  // REST_APPLICATION_H
