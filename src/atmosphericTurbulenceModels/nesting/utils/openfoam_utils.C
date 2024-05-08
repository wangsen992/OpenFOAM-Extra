#include "nesting_utils.H"
#include "polyMesh.H"
// compute cells close to a patch by directly calculating the minimun distance
// of each cell to the patchFaces
// This is a slow method, scales with O(N x M) with N faces and M cells
Foam::labelList getPatchCloseCells(const Foam::polyMesh& mesh, Foam::string patchName, double dis_lim)
{
  // certain distance
  Foam::polyPatch patch(mesh.boundaryMesh()[patchName]);
  Foam::pointField patchFaceCenters(patch.faceCentres());
  Foam::vectorField patchFaceNormals(patch.faceNormals());

  // Find cells cut by rays from face centers
  // by limite the eucleadian distance to set boundarys

  Foam::pointField meshCellCenters(mesh.cellCentres());
  Foam::DynamicList<Foam::label, 10> layerCellsInd;
  for(int i = 0; i != meshCellCenters.size(); i++)
  {
    Foam::tmp<Foam::vectorField> tDisField = (meshCellCenters[i] - patchFaceCenters);
    Foam::scalar dis = min(mag(tDisField));
    if (dis < dis_lim)
    {
      layerCellsInd.append(i);
    }
  }

  Foam::labelList layerCells(layerCellsInd);
  return layerCells;
}

