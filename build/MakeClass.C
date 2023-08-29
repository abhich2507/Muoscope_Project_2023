#include <TFile.h>
#include <TTree.h>

void MakeClass()
{ TFile f("RPC.root");
TTree* Default; 
f.GetObject("Default", Default);
Default->MakeClass("MyClass");
 
}
