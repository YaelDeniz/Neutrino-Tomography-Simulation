#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoShape.h"
#include "TGeoMatrix.h"

void createSphericalShell() {
    // Create a TGeoManager instance
    TGeoManager *geom = new TGeoManager("world", "Example with spherical shell and segment");

    // Define materials
    TGeoMaterial *matAir = new TGeoMaterial("Air", 0, 0, 0);
    TGeoMaterial *matSteel = new TGeoMaterial("Steel", 55.85, 26, 7.87);
    TGeoMaterial *matAluminum = new TGeoMaterial("Aluminum", 26.98, 13, 2.7);

    // Define media
    TGeoMedium *medAir = new TGeoMedium("Air", 1, matAir);
    TGeoMedium *medSteel = new TGeoMedium("Steel", 2, matSteel);
    TGeoMedium *medAluminum = new TGeoMedium("Aluminum", 3, matAluminum);

    // Create the top world volume
    TGeoVolume *world = geom->MakeBox("world", medAir, 200, 200, 200);
    geom->SetTopVolume(world);

    // Create the spherical shell
    TGeoVolume *shell = geom->MakeSphere("shell", medSteel, 90, 100);
    shell->SetTransparency(50);  // Set transparency to 50%
    world->AddNode(shell, 1);

    // Create a spherical segment (theta from 0 to 45 degrees)
    TGeoVolume *segment = geom->MakeSphere("segment", medAluminum, 90, 100, 0, 45);
    segment->SetTransparency(70);  // Set transparency to 70%
    // Create a rotation matrix for the segment
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(90);
    TGeoCombiTrans *trans = new TGeoCombiTrans(0, 0, 0, rot);

    // Add the segment to the shell
    shell->AddNodeOverlap(segment, 1, trans);

    // Close the geometry
    geom->CloseGeometry();

    // Draw the geometry
    world->Draw("ogl");
}

int main() {
    createSphericalShell();
    return 0;
}

