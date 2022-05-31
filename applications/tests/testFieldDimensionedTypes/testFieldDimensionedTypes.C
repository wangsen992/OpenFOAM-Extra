#include "fvCFD.H"
#include "dimensionedTypes.H"
#include "Field.H"


int main(int argc, char *argv[])
{
    List<dimensionedScalar> dimScalar
    (
      2, dimensionedScalar(dimLength, 1)
    );
    Info << "Test Over." << endl;
    return 0;
}
