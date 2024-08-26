// Define parameters for the ellipse
a = 1.0;    // Major axis length
b = 0.5;    // Minor axis length
refinement = 0.1; // Refinement level, smaller values create finer meshes

// Define points
Point(1) = {a, 0, 0, refinement};
Point(2) = {0, b, 0, refinement};
Point(3) = {-a, 0, 0, refinement};
Point(4) = {0, -b, 0, refinement};

// Define ellipse arcs
Ellipse(1) = {1, 0, 2};
Ellipse(2) = {2, 0, 3};
Ellipse(3) = {3, 0, 4};
Ellipse(4) = {4, 0, 1};

// Create a closed curve (Loop) and define the surface
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Mesh the surface
Mesh 2;

// Save the mesh to a file (optional, you can also do this from the command line)
Save "ellipse.msh";
