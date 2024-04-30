// Gmsh project created on Sat Mar 23 18:15:58 2024
// REFERENCE: Gmsh GEO tutorial 13 https://gmsh.info/doc/texinfo/gmsh.html#t13

Merge "brain-h3.0.stl";

// Create a topology
CreateTopology; // AGGIUNTO

// Classify ("color") the surfaces
DefineConstant[
    // Angle between two triangles above which an edge is considered as sharp
    angle = {90, Min 0, Max 180, Step 1,
      Name "Parameters/Angle for surface detection"},
    // For complex geometries, patches can be too complex, too elongated or too
    // large to be parametrized; setting the following option will force the
    // creation of patches that are amenable to reparametrization:
    forceParametrizablePatches = {1, Choices{0,1},
      Name "Parameters/Create surfaces guaranteed to be parametrizable"},
    // For open surfaces include the boundary edges in the classification process:
    includeBoundary = 1,
    // Force curves to be split on given angle:
    curveAngle = 180
  ];
  ClassifySurfaces{angle * Pi/180, includeBoundary, forceParametrizablePatches,
                   curveAngle * Pi / 180};

// Create a geometry
CreateGeometry;

// Create a volume
Surface Loop(1) = Surface{:};
Volume(1) = {1};

// Create a physical group for the volume
Physical Volume("Brain") = {1};

// Specify element sizes
funny = DefineNumber[0, Choices{0,1},
    Name "Parameters/Apply funny mesh size field?" ];
  
  Field[1] = MathEval;
  If(funny)
    Field[1].F = "2*Sin((x+y)/5) + 3";
  Else
    Field[1].F = "4";
  EndIf
  Background Field = 1;

Save "brain-h3.0.msh";

