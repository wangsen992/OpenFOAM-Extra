dimensionedVector f
(
  "f", 
  dimensionSet(0, 0, -1, 0, 0), 
  vector(0, 0, 1.00 * pow(10, -4))
);

dimensionedVector Ug
(
  "Ug",
  dimensionSet(0, 1, -1, 0, 0),
  vector(10, 0, 0)
);

volVectorField rhofU_Ug("rhofU_Ug", rho * ( f ^ (U - Ug)));
