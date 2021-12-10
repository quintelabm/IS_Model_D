#include "IS_Model.h"

using namespace std;

int main(){
  //IS_Model m(1,1);
//  IS_Model* m = new IS_Model();
  IS_Model* model = new IS_Model();
  //model->setSaveFiles(1);
  //model->setSimulationCase(1);
  model->solve();//solve
  return 0;
}
