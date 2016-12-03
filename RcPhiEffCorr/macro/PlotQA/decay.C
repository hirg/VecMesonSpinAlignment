#include "TPythia6.h"

void decay()
{
TPythia6* py;
py = TPythia6::Instance();
py->Pylist(12);
}
