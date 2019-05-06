#include "JToyFlowInputs.h"

void drawvn(){	
	gSystem->Load("JToyFlowInputs_cxx.so");
	JToyFlowInputs *jflowinputs = new JToyFlowInputs();
	jflowinputs->Load();
	jflowinputs->pgr_v[1]->Draw();
}
