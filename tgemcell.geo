// All length quantities are in cm.
rI = 0.0025; // inner hole radius
rO = 0.0035; // outer hole radius
rE = 0.0040; // etch mount(= cupper hole radius > rO(hole radius on dielectric layer)) 
tC = 0.0005; // copper thickness
tD = 0.0050; // dielectric thickness
p = 0.0140; // the "pitch", or distance between GEM holes
dz12 = 0.2;
dz23 = 0.2;
dzU = 0.2;
dzP = 0.2;
// Fourth input expression for Point() to control mesh sizes at the points.
mshSizeDrift = 0.005;
// mesh sizes at points on circles.
mshSizeEtchHole = 0.0008;
mshSizeDielInnerHole = 0.0010;
mshSizeDielOuterHole = 0.0015;
// mesh size at points on circle centers
mshSizeHoleCenter = 0.0025;

// z-coordinate of planes
zP = 0;
z3 = zP + dzP;
z2 = zP + dzP + dz12;
z1 = zP + dzP + dz12 + dz23;
zU = zP + dzP + dz12 + dz23 + dzU;
bdrHalfX = p/2;
bdrHalfY = p*Sqrt(3)/2;

// ---------------------------------------------------------------------------------------------------------------
// GEM1
// ---------------------------------------------------------------------------------------------------------------
/*---------------------------------- Point Definitions -------------------------------------*/
// points on upper copper plate(top and bottom plane)
// hole1 (bottom left quarter hole)
p_g1_uc_h1_t1 = newp; Point(p_g1_uc_h1_t1) = {-bdrHalfX, -bdrHalfY,       z1 + tD/2 + tC, mshSizeHoleCenter};
p_g1_uc_h1_t2 = newp; Point(p_g1_uc_h1_t2) = {-bdrHalfX, -bdrHalfY + rE,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h1_t3 = newp; Point(p_g1_uc_h1_t3) = {-bdrHalfX + rE, -bdrHalfY,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h1_b1 = newp; Point(p_g1_uc_h1_b1) = {-bdrHalfX, -bdrHalfY,       z1 + tD/2, mshSizeHoleCenter};
p_g1_uc_h1_b2 = newp; Point(p_g1_uc_h1_b2) = {-bdrHalfX, -bdrHalfY + rE,  z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h1_b3 = newp; Point(p_g1_uc_h1_b3) = {-bdrHalfX + rE, -bdrHalfY,  z1 + tD/2, mshSizeEtchHole};
// hole2 (bottom right quarter hole)
p_g1_uc_h2_t1 = newp; Point(p_g1_uc_h2_t1) = {bdrHalfX, -bdrHalfY,        z1 + tD/2 + tC, mshSizeHoleCenter};
p_g1_uc_h2_t2 = newp; Point(p_g1_uc_h2_t2) = {bdrHalfX - rE, -bdrHalfY,   z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h2_t3 = newp; Point(p_g1_uc_h2_t3) = {bdrHalfX, -bdrHalfY + rE,   z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h2_b1 = newp; Point(p_g1_uc_h2_b1) = {bdrHalfX, -bdrHalfY,        z1 + tD/2, mshSizeHoleCenter};
p_g1_uc_h2_b2 = newp; Point(p_g1_uc_h2_b2) = {bdrHalfX - rE, -bdrHalfY,   z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h2_b3 = newp; Point(p_g1_uc_h2_b3) = {bdrHalfX, -bdrHalfY + rE,   z1 + tD/2, mshSizeEtchHole};
// hole3 (top right quarter hole)
p_g1_uc_h3_t1 = newp; Point(p_g1_uc_h3_t1) = {bdrHalfX, bdrHalfY,       z1 + tD/2 + tC, mshSizeHoleCenter};
p_g1_uc_h3_t2 = newp; Point(p_g1_uc_h3_t2) = {bdrHalfX, bdrHalfY - rE,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h3_t3 = newp; Point(p_g1_uc_h3_t3) = {bdrHalfX - rE, bdrHalfY,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h3_b1 = newp; Point(p_g1_uc_h3_b1) = {bdrHalfX, bdrHalfY,       z1 + tD/2, mshSizeHoleCenter};
p_g1_uc_h3_b2 = newp; Point(p_g1_uc_h3_b2) = {bdrHalfX, bdrHalfY - rE,  z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h3_b3 = newp; Point(p_g1_uc_h3_b3) = {bdrHalfX - rE, bdrHalfY,  z1 + tD/2, mshSizeEtchHole};
// hole4 (top left quarter hole)
p_g1_uc_h4_t1 = newp; Point(p_g1_uc_h4_t1) = {-bdrHalfX, bdrHalfY,       z1 + tD/2 + tC, mshSizeHoleCenter};
p_g1_uc_h4_t2 = newp; Point(p_g1_uc_h4_t2) = {-bdrHalfX + rE, bdrHalfY,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h4_t3 = newp; Point(p_g1_uc_h4_t3) = {-bdrHalfX, bdrHalfY - rE,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h4_b1 = newp; Point(p_g1_uc_h4_b1) = {-bdrHalfX, bdrHalfY,       z1 + tD/2, mshSizeHoleCenter};
p_g1_uc_h4_b2 = newp; Point(p_g1_uc_h4_b2) = {-bdrHalfX + rE, bdrHalfY,  z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h4_b3 = newp; Point(p_g1_uc_h4_b3) = {-bdrHalfX, bdrHalfY - rE,  z1 + tD/2, mshSizeEtchHole};
// hole5 (center circular hole)
p_g1_uc_h5_t1 = newp; Point(p_g1_uc_h5_t1) = {0, 0,   z1 + tD/2 + tC, mshSizeHoleCenter};
p_g1_uc_h5_t2 = newp; Point(p_g1_uc_h5_t2) = {rE, 0,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h5_t3 = newp; Point(p_g1_uc_h5_t3) = {0, rE,  z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h5_t4 = newp; Point(p_g1_uc_h5_t4) = {-rE, 0, z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h5_t5 = newp; Point(p_g1_uc_h5_t5) = {0, -rE, z1 + tD/2 + tC, mshSizeEtchHole};
p_g1_uc_h5_b1 = newp; Point(p_g1_uc_h5_b1) = {0, 0,   z1 + tD/2, mshSizeHoleCenter};
p_g1_uc_h5_b2 = newp; Point(p_g1_uc_h5_b2) = {rE, 0,  z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h5_b3 = newp; Point(p_g1_uc_h5_b3) = {0, rE,  z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h5_b4 = newp; Point(p_g1_uc_h5_b4) = {-rE, 0, z1 + tD/2, mshSizeEtchHole};
p_g1_uc_h5_b5 = newp; Point(p_g1_uc_h5_b5) = {0, -rE, z1 + tD/2, mshSizeEtchHole};

// points on lower copper plate (top and bottom plane)
p_g1_lc_h1_t1 = newp; Point(p_g1_lc_h1_t1) = {-bdrHalfX, -bdrHalfY,       z1 - tD/2, mshSizeHoleCenter};
p_g1_lc_h1_t2 = newp; Point(p_g1_lc_h1_t2) = {-bdrHalfX, -bdrHalfY + rE,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h1_t3 = newp; Point(p_g1_lc_h1_t3) = {-bdrHalfX + rE, -bdrHalfY,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h1_b1 = newp; Point(p_g1_lc_h1_b1) = {-bdrHalfX, -bdrHalfY,       z1 - (tD/2 + tC), mshSizeHoleCenter};
p_g1_lc_h1_b2 = newp; Point(p_g1_lc_h1_b2) = {-bdrHalfX, -bdrHalfY + rE,  z1 - (tD/2 + tC), mshSizeEtchHole};
p_g1_lc_h1_b3 = newp; Point(p_g1_lc_h1_b3) = {-bdrHalfX + rE, -bdrHalfY,  z1 - (tD/2 + tC), mshSizeEtchHole};
// hole2 (bottom right quarter hole)
p_g1_lc_h2_t1 = newp; Point(p_g1_lc_h2_t1) = {bdrHalfX, -bdrHalfY,        z1 - tD/2, mshSizeHoleCenter};
p_g1_lc_h2_t2 = newp; Point(p_g1_lc_h2_t2) = {bdrHalfX - rE, -bdrHalfY,   z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h2_t3 = newp; Point(p_g1_lc_h2_t3) = {bdrHalfX, -bdrHalfY + rE,   z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h2_b1 = newp; Point(p_g1_lc_h2_b1) = {bdrHalfX, -bdrHalfY,        z1 - (tD/2 + tC), mshSizeHoleCenter};
p_g1_lc_h2_b2 = newp; Point(p_g1_lc_h2_b2) = {bdrHalfX - rE, -bdrHalfY,   z1 - (tD/2 + tC), mshSizeEtchHole};
p_g1_lc_h2_b3 = newp; Point(p_g1_lc_h2_b3) = {bdrHalfX, -bdrHalfY + rE,   z1 - (tD/2 + tC), mshSizeEtchHole};
// hole3 (top right quarter hole)
p_g1_lc_h3_t1 = newp; Point(p_g1_lc_h3_t1) = {bdrHalfX, bdrHalfY,       z1 - tD/2, mshSizeHoleCenter};
p_g1_lc_h3_t2 = newp; Point(p_g1_lc_h3_t2) = {bdrHalfX, bdrHalfY - rE,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h3_t3 = newp; Point(p_g1_lc_h3_t3) = {bdrHalfX - rE, bdrHalfY,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h3_b1 = newp; Point(p_g1_lc_h3_b1) = {bdrHalfX, bdrHalfY,       z1 - (tD/2 + tC), mshSizeHoleCenter};
p_g1_lc_h3_b2 = newp; Point(p_g1_lc_h3_b2) = {bdrHalfX, bdrHalfY - rE,  z1 - (tD/2 + tC), mshSizeEtchHole};
p_g1_lc_h3_b3 = newp; Point(p_g1_lc_h3_b3) = {bdrHalfX - rE, bdrHalfY,  z1 - (tD/2 + tC), mshSizeEtchHole};
// hole4 (top left quarter hole)
p_g1_lc_h4_t1 = newp; Point(p_g1_lc_h4_t1) = {-bdrHalfX, bdrHalfY,       z1 - tD/2, mshSizeHoleCenter};
p_g1_lc_h4_t2 = newp; Point(p_g1_lc_h4_t2) = {-bdrHalfX + rE, bdrHalfY,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h4_t3 = newp; Point(p_g1_lc_h4_t3) = {-bdrHalfX, bdrHalfY - rE,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h4_b1 = newp; Point(p_g1_lc_h4_b1) = {-bdrHalfX, bdrHalfY,       z1 - (tD/2 + tC), mshSizeHoleCenter};
p_g1_lc_h4_b2 = newp; Point(p_g1_lc_h4_b2) = {-bdrHalfX + rE, bdrHalfY,  z1 - (tD/2 + tC), mshSizeEtchHole};
p_g1_lc_h4_b3 = newp; Point(p_g1_lc_h4_b3) = {-bdrHalfX, bdrHalfY - rE,  z1 - (tD/2 + tC), mshSizeEtchHole};
// hole5 (center circular quarter hole)
p_g1_lc_h5_t1 = newp; Point(p_g1_lc_h5_t1) = {0, 0,   z1 - tD/2, mshSizeHoleCenter};
p_g1_lc_h5_t2 = newp; Point(p_g1_lc_h5_t2) = {rE, 0,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h5_t3 = newp; Point(p_g1_lc_h5_t3) = {0, rE,  z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h5_t4 = newp; Point(p_g1_lc_h5_t4) = {-rE, 0, z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h5_t5 = newp; Point(p_g1_lc_h5_t5) = {0, -rE, z1 - tD/2, mshSizeEtchHole};
p_g1_lc_h5_b1 = newp; Point(p_g1_lc_h5_b1) = {0, 0,   z1 - tD/2 - tC, mshSizeHoleCenter};
p_g1_lc_h5_b2 = newp; Point(p_g1_lc_h5_b2) = {rE, 0,  z1 - tD/2 - tC, mshSizeEtchHole};
p_g1_lc_h5_b3 = newp; Point(p_g1_lc_h5_b3) = {0, rE,  z1 - tD/2 - tC, mshSizeEtchHole};
p_g1_lc_h5_b4 = newp; Point(p_g1_lc_h5_b4) = {-rE, 0, z1 - tD/2 - tC, mshSizeEtchHole};
p_g1_lc_h5_b5 = newp; Point(p_g1_lc_h5_b5) = {0, -rE, z1 - tD/2 - tC, mshSizeEtchHole};

// points on upper dielectric layer(top, middel and bottom plane)
// hole1 (bottom left quarter hole)
p_g1_diel_h1_t1 = p_g1_uc_h1_b1;
p_g1_diel_h1_t2 = newp; Point(p_g1_diel_h1_t2) = {-bdrHalfX, -bdrHalfY + rO,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h1_t3 = newp; Point(p_g1_diel_h1_t3) = {-bdrHalfX + rO, -bdrHalfY,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h1_m1 = newp; Point(p_g1_diel_h1_m1) = {-bdrHalfX, -bdrHalfY,       z1, mshSizeHoleCenter};
p_g1_diel_h1_m2 = newp; Point(p_g1_diel_h1_m2) = {-bdrHalfX, -bdrHalfY + rI,  z1, mshSizeDielInnerHole};
p_g1_diel_h1_m3 = newp; Point(p_g1_diel_h1_m3) = {-bdrHalfX + rI, -bdrHalfY,  z1, mshSizeDielInnerHole};
p_g1_diel_h1_b1 = p_g1_lc_h1_t1;
p_g1_diel_h1_b2 = newp; Point(p_g1_diel_h1_b2) = {-bdrHalfX, -bdrHalfY + rO,  z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h1_b3 = newp; Point(p_g1_diel_h1_b3) = {-bdrHalfX + rO, -bdrHalfY,  z1 - tD/2, mshSizeDielOuterHole};
// hole2 (bottom right quarter hole)
p_g1_diel_h2_t1 = p_g1_uc_h2_b1;
p_g1_diel_h2_t2 = newp; Point(p_g1_diel_h2_t2) = {bdrHalfX - rO, -bdrHalfY,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h2_t3 = newp; Point(p_g1_diel_h2_t3) = {bdrHalfX, -bdrHalfY + rO,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h2_m1 = newp; Point(p_g1_diel_h2_m1) = {bdrHalfX, -bdrHalfY,       z1, mshSizeHoleCenter};
p_g1_diel_h2_m2 = newp; Point(p_g1_diel_h2_m2) = {bdrHalfX - rI, -bdrHalfY,  z1, mshSizeDielInnerHole};
p_g1_diel_h2_m3 = newp; Point(p_g1_diel_h2_m3) = {bdrHalfX, -bdrHalfY + rI,  z1, mshSizeDielInnerHole};
p_g1_diel_h2_b1 = p_g1_lc_h2_t1;
p_g1_diel_h2_b2 = newp; Point(p_g1_diel_h2_b2) = {bdrHalfX - rO, -bdrHalfY,  z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h2_b3 = newp; Point(p_g1_diel_h2_b3) = {bdrHalfX, -bdrHalfY + rO,  z1 - tD/2, mshSizeDielOuterHole};
// hole3 (top right quarter hole)
p_g1_diel_h3_t1 = p_g1_uc_h3_b1;
p_g1_diel_h3_t2 = newp; Point(p_g1_diel_h3_t2) = {bdrHalfX, bdrHalfY - rO,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h3_t3 = newp; Point(p_g1_diel_h3_t3) = {bdrHalfX - rO, bdrHalfY,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h3_m1 = newp; Point(p_g1_diel_h3_m1) = {bdrHalfX, bdrHalfY,       z1, mshSizeHoleCenter};
p_g1_diel_h3_m2 = newp; Point(p_g1_diel_h3_m2) = {bdrHalfX, bdrHalfY - rI,  z1, mshSizeDielOuterHole};
p_g1_diel_h3_m3 = newp; Point(p_g1_diel_h3_m3) = {bdrHalfX - rI, bdrHalfY,  z1, mshSizeDielOuterHole};
p_g1_diel_h3_b1 = p_g1_lc_h3_t1;
p_g1_diel_h3_b2 = newp; Point(p_g1_diel_h3_b2) = {bdrHalfX, bdrHalfY - rO,  z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h3_b3 = newp; Point(p_g1_diel_h3_b3) = {bdrHalfX - rO, bdrHalfY,  z1 - tD/2, mshSizeDielOuterHole};
// hole4 (top left quarter hole)
p_g1_diel_h4_t1 = p_g1_uc_h4_b1;
p_g1_diel_h4_t2 = newp; Point(p_g1_diel_h4_t2) = {-bdrHalfX + rO, bdrHalfY,    z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h4_t3 = newp; Point(p_g1_diel_h4_t3) = {-bdrHalfX, bdrHalfY - rO,    z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h4_m1 = newp; Point(p_g1_diel_h4_m1) = {-bdrHalfX, bdrHalfY,         z1, mshSizeHoleCenter};
p_g1_diel_h4_m2 = newp; Point(p_g1_diel_h4_m2) = {-bdrHalfX + rI, bdrHalfY,    z1, mshSizeDielInnerHole};
p_g1_diel_h4_m3 = newp; Point(p_g1_diel_h4_m3) = {-bdrHalfX, bdrHalfY - rI ,   z1, mshSizeDielInnerHole};
p_g1_diel_h4_b1 = p_g1_lc_h4_t1;
p_g1_diel_h4_b2 = newp; Point(p_g1_diel_h4_b2) = {-bdrHalfX + rO, bdrHalfY,    z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h4_b3 = newp; Point(p_g1_diel_h4_b3) = {-bdrHalfX, bdrHalfY - rO,    z1 - tD/2, mshSizeDielOuterHole};
// hole5 (center circular quarter hole)
p_g1_diel_h5_t1 = p_g1_uc_h5_b1;
p_g1_diel_h5_t2 = newp; Point(p_g1_diel_h5_t2) = {rO, 0,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_t3 = newp; Point(p_g1_diel_h5_t3) = {0, rO,  z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_t4 = newp; Point(p_g1_diel_h5_t4) = {-rO, 0, z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_t5 = newp; Point(p_g1_diel_h5_t5) = {0, -rO, z1 + tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_m1 = newp; Point(p_g1_diel_h5_m1) = {0, 0,   z1, mshSizeHoleCenter};
p_g1_diel_h5_m2 = newp; Point(p_g1_diel_h5_m2) = {rI, 0,  z1, mshSizeDielInnerHole};
p_g1_diel_h5_m3 = newp; Point(p_g1_diel_h5_m3) = {0, rI,  z1, mshSizeDielInnerHole};
p_g1_diel_h5_m4 = newp; Point(p_g1_diel_h5_m4) = {-rI, 0, z1, mshSizeDielInnerHole};
p_g1_diel_h5_m5 = newp; Point(p_g1_diel_h5_m5) = {0, -rI, z1, mshSizeDielInnerHole};
p_g1_diel_h5_b1 = p_g1_lc_h5_t1;
p_g1_diel_h5_b2 = newp; Point(p_g1_diel_h5_b2) = {rO, 0,  z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_b3 = newp; Point(p_g1_diel_h5_b3) = {0, rO,  z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_b4 = newp; Point(p_g1_diel_h5_b4) = {-rO, 0, z1 - tD/2, mshSizeDielOuterHole};
p_g1_diel_h5_b5 = newp; Point(p_g1_diel_h5_b5) = {0, -rO, z1 - tD/2, mshSizeDielOuterHole};
/*---------------------------------- Line Definitions -------------------------------------*/
// lines on upper plane
// on the top plane
l_g1_uc_h1_t1 = newl; Circle(l_g1_uc_h1_t1) = {p_g1_uc_h1_t2, p_g1_uc_h1_t1, p_g1_uc_h1_t3}; // hole 1
l_g1_uc_h12_t1 = newl; Line(l_g1_uc_h12_t1) = {p_g1_uc_h1_t3, p_g1_uc_h2_t2}; // between hole 1 and hole 2
l_g1_uc_h2_t1 = newl; Circle(l_g1_uc_h2_t1) = {p_g1_uc_h2_t2, p_g1_uc_h2_t1, p_g1_uc_h2_t3}; // hole 2
l_g1_uc_h23_t1 = newl; Line(l_g1_uc_h23_t1) = {p_g1_uc_h2_t3, p_g1_uc_h3_t2}; // between hole 2 and hole 3
l_g1_uc_h3_t1 = newl; Circle(l_g1_uc_h3_t1) = {p_g1_uc_h3_t2, p_g1_uc_h3_t1, p_g1_uc_h3_t3}; // hole 3
l_g1_uc_h34_t1 = newl; Line(l_g1_uc_h34_t1) = {p_g1_uc_h3_t3, p_g1_uc_h4_t2}; // between hole 3 and hole 4
l_g1_uc_h4_t1 = newl; Circle(l_g1_uc_h4_t1) = {p_g1_uc_h4_t2, p_g1_uc_h4_t1, p_g1_uc_h4_t3}; // hole 4
l_g1_uc_h41_t1 = newl; Line(l_g1_uc_h41_t1) = {p_g1_uc_h4_t3, p_g1_uc_h1_t2}; // between hole 4 and hole 1
l_g1_uc_h5_t1 = newl; Circle(l_g1_uc_h5_t1) = {p_g1_uc_h5_t2, p_g1_uc_h5_t1, p_g1_uc_h5_t3}; // hole 5
l_g1_uc_h5_t2 = newl; Circle(l_g1_uc_h5_t2) = {p_g1_uc_h5_t3, p_g1_uc_h5_t1, p_g1_uc_h5_t4};
l_g1_uc_h5_t3 = newl; Circle(l_g1_uc_h5_t3) = {p_g1_uc_h5_t4, p_g1_uc_h5_t1, p_g1_uc_h5_t5};
l_g1_uc_h5_t4 = newl; Circle(l_g1_uc_h5_t4) = {p_g1_uc_h5_t5, p_g1_uc_h5_t1, p_g1_uc_h5_t2};
// on the bottom plane
l_g1_uc_h1_b1 = newl; Circle(l_g1_uc_h1_b1) = {p_g1_uc_h1_b2, p_g1_uc_h1_b1, p_g1_uc_h1_b3}; // hole 1
l_g1_uc_h12_b1 = newl; Line(l_g1_uc_h12_b1) = {p_g1_uc_h1_b3, p_g1_uc_h2_b2}; // between hole 1 and hole 2
l_g1_uc_h2_b1 = newl; Circle(l_g1_uc_h2_b1) = {p_g1_uc_h2_b2, p_g1_uc_h2_b1, p_g1_uc_h2_b3}; // hole 2
l_g1_uc_h23_b1 = newl; Line(l_g1_uc_h23_b1) = {p_g1_uc_h2_b3, p_g1_uc_h3_b2}; // between hole 2 and hole 3
l_g1_uc_h3_b1 = newl; Circle(l_g1_uc_h3_b1) = {p_g1_uc_h3_b2, p_g1_uc_h3_b1, p_g1_uc_h3_b3}; // hole 3
l_g1_uc_h34_b1 = newl; Line(l_g1_uc_h34_b1) = {p_g1_uc_h3_b3, p_g1_uc_h4_b2}; // between hole 3 and hole 4
l_g1_uc_h4_b1 = newl; Circle(l_g1_uc_h4_b1) = {p_g1_uc_h4_b2, p_g1_uc_h4_b1, p_g1_uc_h4_b3}; // hole 4
l_g1_uc_h41_b1 = newl; Line(l_g1_uc_h41_b1) = {p_g1_uc_h4_b3, p_g1_uc_h1_b2}; // between hole 4 and hole 1
l_g1_uc_h5_b1 = newl; Circle(l_g1_uc_h5_b1) = {p_g1_uc_h5_b2, p_g1_uc_h5_b1, p_g1_uc_h5_b3}; // hole 5
l_g1_uc_h5_b2 = newl; Circle(l_g1_uc_h5_b2) = {p_g1_uc_h5_b3, p_g1_uc_h5_b1, p_g1_uc_h5_b4};
l_g1_uc_h5_b3 = newl; Circle(l_g1_uc_h5_b3) = {p_g1_uc_h5_b4, p_g1_uc_h5_b1, p_g1_uc_h5_b5};
l_g1_uc_h5_b4 = newl; Circle(l_g1_uc_h5_b4) = {p_g1_uc_h5_b5, p_g1_uc_h5_b1, p_g1_uc_h5_b2};
// connecting top and bottom
l_g1_uc_h1_tb1 = newl; Line(l_g1_uc_h1_tb1) = {p_g1_uc_h1_t1, p_g1_uc_h1_b1}; // hole 1
l_g1_uc_h1_tb2 = newl; Line(l_g1_uc_h1_tb2) = {p_g1_uc_h1_t2, p_g1_uc_h1_b2}; 
l_g1_uc_h1_tb3 = newl; Line(l_g1_uc_h1_tb3) = {p_g1_uc_h1_t3, p_g1_uc_h1_b3};
l_g1_uc_h2_tb1 = newl; Line(l_g1_uc_h2_tb1) = {p_g1_uc_h2_t1, p_g1_uc_h2_b1}; // hole 2
l_g1_uc_h2_tb2 = newl; Line(l_g1_uc_h2_tb2) = {p_g1_uc_h2_t2, p_g1_uc_h2_b2};
l_g1_uc_h2_tb3 = newl; Line(l_g1_uc_h2_tb3) = {p_g1_uc_h2_t3, p_g1_uc_h2_b3};
l_g1_uc_h3_tb1 = newl; Line(l_g1_uc_h3_tb1) = {p_g1_uc_h3_t1, p_g1_uc_h3_b1}; // hole 3
l_g1_uc_h3_tb2 = newl; Line(l_g1_uc_h3_tb2) = {p_g1_uc_h3_t2, p_g1_uc_h3_b2};
l_g1_uc_h3_tb3 = newl; Line(l_g1_uc_h3_tb3) = {p_g1_uc_h3_t3, p_g1_uc_h3_b3};
l_g1_uc_h4_tb1 = newl; Line(l_g1_uc_h4_tb1) = {p_g1_uc_h4_t1, p_g1_uc_h4_b1}; // hole 4
l_g1_uc_h4_tb2 = newl; Line(l_g1_uc_h4_tb2) = {p_g1_uc_h4_t2, p_g1_uc_h4_b2};
l_g1_uc_h4_tb3 = newl; Line(l_g1_uc_h4_tb3) = {p_g1_uc_h4_t3, p_g1_uc_h4_b3};
l_g1_uc_h5_tb1 = newl; Line(l_g1_uc_h5_tb1) = {p_g1_uc_h5_t1, p_g1_uc_h5_b1}; // hole 5
l_g1_uc_h5_tb2 = newl; Line(l_g1_uc_h5_tb2) = {p_g1_uc_h5_t2, p_g1_uc_h5_b2};
l_g1_uc_h5_tb3 = newl; Line(l_g1_uc_h5_tb3) = {p_g1_uc_h5_t3, p_g1_uc_h5_b3};
l_g1_uc_h5_tb4 = newl; Line(l_g1_uc_h5_tb4) = {p_g1_uc_h5_t4, p_g1_uc_h5_b4};
l_g1_uc_h5_tb5 = newl; Line(l_g1_uc_h5_tb5) = {p_g1_uc_h5_t5, p_g1_uc_h5_b5};

// lines on lower plane
// on the top plane
l_g1_lc_h1_t1 = newl; Circle(l_g1_lc_h1_t1) = {p_g1_lc_h1_t2, p_g1_lc_h1_t1, p_g1_lc_h1_t3}; // hole 1
l_g1_lc_h12_t1 = newl; Line(l_g1_lc_h12_t1) = {p_g1_lc_h1_t3, p_g1_lc_h2_t2}; // between hole 1 and hole 2
l_g1_lc_h2_t1 = newl; Circle(l_g1_lc_h2_t1) = {p_g1_lc_h2_t2, p_g1_lc_h2_t1, p_g1_lc_h2_t3}; // hole 2
l_g1_lc_h23_t1 = newl; Line(l_g1_lc_h23_t1) = {p_g1_lc_h2_t3, p_g1_lc_h3_t2}; // between hole 2 and hole 3
l_g1_lc_h3_t1 = newl; Circle(l_g1_lc_h3_t1) = {p_g1_lc_h3_t2, p_g1_lc_h3_t1, p_g1_lc_h3_t3}; // hole 3
l_g1_lc_h34_t1 = newl; Line(l_g1_lc_h34_t1) = {p_g1_lc_h3_t3, p_g1_lc_h4_t2}; // between hole 3 and hole 4
l_g1_lc_h4_t1 = newl; Circle(l_g1_lc_h4_t1) = {p_g1_lc_h4_t2, p_g1_lc_h4_t1, p_g1_lc_h4_t3}; // hole 4
l_g1_lc_h41_t1 = newl; Line(l_g1_lc_h41_t1) = {p_g1_lc_h4_t3, p_g1_lc_h1_t2}; // between hole 4 and hole 1
l_g1_lc_h5_t1 = newl; Circle(l_g1_lc_h5_t1) = {p_g1_lc_h5_t2, p_g1_lc_h5_t1, p_g1_lc_h5_t3}; // hole 5
l_g1_lc_h5_t2 = newl; Circle(l_g1_lc_h5_t2) = {p_g1_lc_h5_t3, p_g1_lc_h5_t1, p_g1_lc_h5_t4};
l_g1_lc_h5_t3 = newl; Circle(l_g1_lc_h5_t3) = {p_g1_lc_h5_t4, p_g1_lc_h5_t1, p_g1_lc_h5_t5};
l_g1_lc_h5_t4 = newl; Circle(l_g1_lc_h5_t4) = {p_g1_lc_h5_t5, p_g1_lc_h5_t1, p_g1_lc_h5_t2};
// on the bottom plane
l_g1_lc_h1_b1 = newl; Circle(l_g1_lc_h1_b1) = {p_g1_lc_h1_b2, p_g1_lc_h1_b1, p_g1_lc_h1_b3}; // hole 1
l_g1_lc_h12_b1 = newl; Line(l_g1_lc_h12_b1) = {p_g1_lc_h1_b3, p_g1_lc_h2_b2}; // between hole 1 and hole 2
l_g1_lc_h2_b1 = newl; Circle(l_g1_lc_h2_b1) = {p_g1_lc_h2_b2, p_g1_lc_h2_b1, p_g1_lc_h2_b3}; // hole 2
l_g1_lc_h23_b1 = newl; Line(l_g1_lc_h23_b1) = {p_g1_lc_h2_b3, p_g1_lc_h3_b2}; // between hole 2 and hole 3
l_g1_lc_h3_b1 = newl; Circle(l_g1_lc_h3_b1) = {p_g1_lc_h3_b2, p_g1_lc_h3_b1, p_g1_lc_h3_b3}; // hole 3
l_g1_lc_h34_b1 = newl; Line(l_g1_lc_h34_b1) = {p_g1_lc_h3_b3, p_g1_lc_h4_b2}; // between hole 3 and hole 4
l_g1_lc_h4_b1 = newl; Circle(l_g1_lc_h4_b1) = {p_g1_lc_h4_b2, p_g1_lc_h4_b1, p_g1_lc_h4_b3}; // hole 4
l_g1_lc_h41_b1 = newl; Line(l_g1_lc_h41_b1) = {p_g1_lc_h4_b3, p_g1_lc_h1_b2}; // between hole 4 and hole 1
l_g1_lc_h5_b1 = newl; Circle(l_g1_lc_h5_b1) = {p_g1_lc_h5_b2, p_g1_lc_h5_b1, p_g1_lc_h5_b3}; // hole 5
l_g1_lc_h5_b2 = newl; Circle(l_g1_lc_h5_b2) = {p_g1_lc_h5_b3, p_g1_lc_h5_b1, p_g1_lc_h5_b4};
l_g1_lc_h5_b3 = newl; Circle(l_g1_lc_h5_b3) = {p_g1_lc_h5_b4, p_g1_lc_h5_b1, p_g1_lc_h5_b5};
l_g1_lc_h5_b4 = newl; Circle(l_g1_lc_h5_b4) = {p_g1_lc_h5_b5, p_g1_lc_h5_b1, p_g1_lc_h5_b2};
// connecting top and bottom
l_g1_lc_h1_tb1 = newl; Line(l_g1_lc_h1_tb1) = {p_g1_lc_h1_t1, p_g1_lc_h1_b1}; // hole 1
l_g1_lc_h1_tb2 = newl; Line(l_g1_lc_h1_tb2) = {p_g1_lc_h1_t2, p_g1_lc_h1_b2}; 
l_g1_lc_h1_tb3 = newl; Line(l_g1_lc_h1_tb3) = {p_g1_lc_h1_t3, p_g1_lc_h1_b3};
l_g1_lc_h2_tb1 = newl; Line(l_g1_lc_h2_tb1) = {p_g1_lc_h2_t1, p_g1_lc_h2_b1}; // hole 2
l_g1_lc_h2_tb2 = newl; Line(l_g1_lc_h2_tb2) = {p_g1_lc_h2_t2, p_g1_lc_h2_b2};
l_g1_lc_h2_tb3 = newl; Line(l_g1_lc_h2_tb3) = {p_g1_lc_h2_t3, p_g1_lc_h2_b3};
l_g1_lc_h3_tb1 = newl; Line(l_g1_lc_h3_tb1) = {p_g1_lc_h3_t1, p_g1_lc_h3_b1}; // hole 3
l_g1_lc_h3_tb2 = newl; Line(l_g1_lc_h3_tb2) = {p_g1_lc_h3_t2, p_g1_lc_h3_b2};
l_g1_lc_h3_tb3 = newl; Line(l_g1_lc_h3_tb3) = {p_g1_lc_h3_t3, p_g1_lc_h3_b3};
l_g1_lc_h4_tb1 = newl; Line(l_g1_lc_h4_tb1) = {p_g1_lc_h4_t1, p_g1_lc_h4_b1}; // hole 4
l_g1_lc_h4_tb2 = newl; Line(l_g1_lc_h4_tb2) = {p_g1_lc_h4_t2, p_g1_lc_h4_b2};
l_g1_lc_h4_tb3 = newl; Line(l_g1_lc_h4_tb3) = {p_g1_lc_h4_t3, p_g1_lc_h4_b3};
l_g1_lc_h5_tb1 = newl; Line(l_g1_lc_h5_tb1) = {p_g1_lc_h5_t1, p_g1_lc_h5_b1}; // hole 5
l_g1_lc_h5_tb2 = newl; Line(l_g1_lc_h5_tb2) = {p_g1_lc_h5_t2, p_g1_lc_h5_b2};
l_g1_lc_h5_tb3 = newl; Line(l_g1_lc_h5_tb3) = {p_g1_lc_h5_t3, p_g1_lc_h5_b3};
l_g1_lc_h5_tb4 = newl; Line(l_g1_lc_h5_tb4) = {p_g1_lc_h5_t4, p_g1_lc_h5_b4};
l_g1_lc_h5_tb5 = newl; Line(l_g1_lc_h5_tb5) = {p_g1_lc_h5_t5, p_g1_lc_h5_b5};

// lines on the dielectric layer
// on the top plane
l_g1_diel_h1_t1 = newl; Circle(l_g1_diel_h1_t1) = {p_g1_diel_h1_t2, p_g1_diel_h1_t1, p_g1_diel_h1_t3}; // hole 1
l_g1_diel_h1_t2 = newl; Line(l_g1_diel_h1_t2) = {p_g1_diel_h1_t2, p_g1_uc_h1_b2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h1_t3 = newl; Line(l_g1_diel_h1_t3) = {p_g1_diel_h1_t3, p_g1_uc_h1_b3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h2_t1 = newl; Circle(l_g1_diel_h2_t1) = {p_g1_diel_h2_t2, p_g1_diel_h2_t1, p_g1_diel_h2_t3}; // hole 2
l_g1_diel_h2_t2 = newl; Line(l_g1_diel_h2_t2) = {p_g1_diel_h2_t2, p_g1_uc_h2_b2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h2_t3 = newl; Line(l_g1_diel_h2_t3) = {p_g1_diel_h2_t3, p_g1_uc_h2_b3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h3_t1 = newl; Circle(l_g1_diel_h3_t1) = {p_g1_diel_h3_t2, p_g1_diel_h3_t1, p_g1_diel_h3_t3}; // hole 3
l_g1_diel_h3_t2 = newl; Line(l_g1_diel_h3_t2) = {p_g1_diel_h3_t2, p_g1_uc_h3_b2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h3_t3 = newl; Line(l_g1_diel_h3_t3) = {p_g1_diel_h3_t3, p_g1_uc_h3_b3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h4_t1 = newl; Circle(l_g1_diel_h4_t1) = {p_g1_diel_h4_t2, p_g1_diel_h4_t1, p_g1_diel_h4_t3}; // hole 4
l_g1_diel_h4_t2 = newl; Line(l_g1_diel_h4_t2) = {p_g1_diel_h4_t2, p_g1_uc_h4_b2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h4_t3 = newl; Line(l_g1_diel_h4_t3) = {p_g1_diel_h4_t3, p_g1_uc_h4_b3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_t1 = newl; Circle(l_g1_diel_h5_t1) = {p_g1_diel_h5_t2, p_g1_diel_h5_t1, p_g1_diel_h5_t3}; // hole 5
l_g1_diel_h5_t2 = newl; Circle(l_g1_diel_h5_t2) = {p_g1_diel_h5_t3, p_g1_diel_h5_t1, p_g1_diel_h5_t4};
l_g1_diel_h5_t3 = newl; Circle(l_g1_diel_h5_t3) = {p_g1_diel_h5_t4, p_g1_diel_h5_t1, p_g1_diel_h5_t5};
l_g1_diel_h5_t4 = newl; Circle(l_g1_diel_h5_t4) = {p_g1_diel_h5_t5, p_g1_diel_h5_t1, p_g1_diel_h5_t2};
l_g1_diel_h5_t5 = newl; Line(l_g1_diel_h5_t5) = {p_g1_diel_h5_t2, p_g1_uc_h5_b2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_t6 = newl; Line(l_g1_diel_h5_t6) = {p_g1_diel_h5_t3, p_g1_uc_h5_b3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_t7 = newl; Line(l_g1_diel_h5_t7) = {p_g1_diel_h5_t4, p_g1_uc_h5_b4}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_t8 = newl; Line(l_g1_diel_h5_t8) = {p_g1_diel_h5_t5, p_g1_uc_h5_b5}; // this line is on the surface exposed to gas by etching
// connecting the top and middel plane
l_g1_diel_h1_tm1 = newl; Line(l_g1_diel_h1_tm1) = {p_g1_diel_h1_t1, p_g1_diel_h1_m1};
l_g1_diel_h1_tm2 = newl; Line(l_g1_diel_h1_tm2) = {p_g1_diel_h1_t2, p_g1_diel_h1_m2};
l_g1_diel_h1_tm3 = newl; Line(l_g1_diel_h1_tm3) = {p_g1_diel_h1_t3, p_g1_diel_h1_m3};
l_g1_diel_h2_tm1 = newl; Line(l_g1_diel_h2_tm1) = {p_g1_diel_h2_t1, p_g1_diel_h2_m1};
l_g1_diel_h2_tm2 = newl; Line(l_g1_diel_h2_tm2) = {p_g1_diel_h2_t2, p_g1_diel_h2_m2};
l_g1_diel_h2_tm3 = newl; Line(l_g1_diel_h2_tm3) = {p_g1_diel_h2_t3, p_g1_diel_h2_m3};
l_g1_diel_h3_tm1 = newl; Line(l_g1_diel_h3_tm1) = {p_g1_diel_h3_t1, p_g1_diel_h3_m1};
l_g1_diel_h3_tm2 = newl; Line(l_g1_diel_h3_tm2) = {p_g1_diel_h3_t2, p_g1_diel_h3_m2};
l_g1_diel_h3_tm3 = newl; Line(l_g1_diel_h3_tm3) = {p_g1_diel_h3_t3, p_g1_diel_h3_m3};
l_g1_diel_h4_tm1 = newl; Line(l_g1_diel_h4_tm1) = {p_g1_diel_h4_t1, p_g1_diel_h4_m1};
l_g1_diel_h4_tm2 = newl; Line(l_g1_diel_h4_tm2) = {p_g1_diel_h4_t2, p_g1_diel_h4_m2};
l_g1_diel_h4_tm3 = newl; Line(l_g1_diel_h4_tm3) = {p_g1_diel_h4_t3, p_g1_diel_h4_m3};
l_g1_diel_h5_tm1 = newl; Line(l_g1_diel_h5_tm1) = {p_g1_diel_h5_t1, p_g1_diel_h5_m1};
l_g1_diel_h5_tm2 = newl; Line(l_g1_diel_h5_tm2) = {p_g1_diel_h5_t2, p_g1_diel_h5_m2};
l_g1_diel_h5_tm3 = newl; Line(l_g1_diel_h5_tm3) = {p_g1_diel_h5_t3, p_g1_diel_h5_m3};
l_g1_diel_h5_tm4 = newl; Line(l_g1_diel_h5_tm4) = {p_g1_diel_h5_t4, p_g1_diel_h5_m4};
l_g1_diel_h5_tm5 = newl; Line(l_g1_diel_h5_tm5) = {p_g1_diel_h5_t5, p_g1_diel_h5_m5};
// on the middle plane
l_g1_diel_h1_m1 = newl; Circle(l_g1_diel_h1_m1) = {p_g1_diel_h1_m2, p_g1_diel_h1_m1, p_g1_diel_h1_m3}; // hole 1
l_g1_diel_h12_m1 = newl; Line(l_g1_diel_h12_m1) = {p_g1_diel_h1_m3, p_g1_diel_h2_m2}; // between hole 1 and hole 2
l_g1_diel_h2_m1 = newl; Circle(l_g1_diel_h2_m1) = {p_g1_diel_h2_m2, p_g1_diel_h2_m1, p_g1_diel_h2_m3}; // hole 2
l_g1_diel_h23_m1 = newl; Line(l_g1_diel_h23_m1) = {p_g1_diel_h2_m3, p_g1_diel_h3_m2}; // between hole 2 and hole 3
l_g1_diel_h3_m1 = newl; Circle(l_g1_diel_h3_m1) = {p_g1_diel_h3_m2, p_g1_diel_h3_m1, p_g1_diel_h3_m3}; // hole 3
l_g1_diel_h34_m1 = newl; Line(l_g1_diel_h34_m1) = {p_g1_diel_h3_m3, p_g1_diel_h4_m2}; // between hole 3 and hole 4
l_g1_diel_h4_m1 = newl; Circle(l_g1_diel_h4_m1) = {p_g1_diel_h4_m2, p_g1_diel_h4_m1, p_g1_diel_h4_m3}; // hole 4
l_g1_diel_h41_m1 = newl; Line(l_g1_diel_h41_m1) = {p_g1_diel_h4_m3, p_g1_diel_h1_m2}; // between hole 4 and hole 1
l_g1_diel_h5_m1 = newl; Circle(l_g1_diel_h5_m1) = {p_g1_diel_h5_m2, p_g1_diel_h5_m1, p_g1_diel_h5_m3}; // hole 5
l_g1_diel_h5_m2 = newl; Circle(l_g1_diel_h5_m2) = {p_g1_diel_h5_m3, p_g1_diel_h5_m1, p_g1_diel_h5_m4};
l_g1_diel_h5_m3 = newl; Circle(l_g1_diel_h5_m3) = {p_g1_diel_h5_m4, p_g1_diel_h5_m1, p_g1_diel_h5_m5};
l_g1_diel_h5_m4 = newl; Circle(l_g1_diel_h5_m4) = {p_g1_diel_h5_m5, p_g1_diel_h5_m1, p_g1_diel_h5_m2};
// connecting the middle and bottom plane
l_g1_diel_h1_mb1 = newl; Line(l_g1_diel_h1_mb1) = {p_g1_diel_h1_m1, p_g1_diel_h1_b1};
l_g1_diel_h1_mb2 = newl; Line(l_g1_diel_h1_mb2) = {p_g1_diel_h1_m2, p_g1_diel_h1_b2};
l_g1_diel_h1_mb3 = newl; Line(l_g1_diel_h1_mb3) = {p_g1_diel_h1_m3, p_g1_diel_h1_b3};
l_g1_diel_h2_mb1 = newl; Line(l_g1_diel_h2_mb1) = {p_g1_diel_h2_m1, p_g1_diel_h2_b1};
l_g1_diel_h2_mb2 = newl; Line(l_g1_diel_h2_mb2) = {p_g1_diel_h2_m2, p_g1_diel_h2_b2};
l_g1_diel_h2_mb3 = newl; Line(l_g1_diel_h2_mb3) = {p_g1_diel_h2_m3, p_g1_diel_h2_b3};
l_g1_diel_h3_mb1 = newl; Line(l_g1_diel_h3_mb1) = {p_g1_diel_h3_m1, p_g1_diel_h3_b1};
l_g1_diel_h3_mb2 = newl; Line(l_g1_diel_h3_mb2) = {p_g1_diel_h3_m2, p_g1_diel_h3_b2};
l_g1_diel_h3_mb3 = newl; Line(l_g1_diel_h3_mb3) = {p_g1_diel_h3_m3, p_g1_diel_h3_b3};
l_g1_diel_h4_mb1 = newl; Line(l_g1_diel_h4_mb1) = {p_g1_diel_h4_m1, p_g1_diel_h4_b1};
l_g1_diel_h4_mb2 = newl; Line(l_g1_diel_h4_mb2) = {p_g1_diel_h4_m2, p_g1_diel_h4_b2};
l_g1_diel_h4_mb3 = newl; Line(l_g1_diel_h4_mb3) = {p_g1_diel_h4_m3, p_g1_diel_h4_b3};
l_g1_diel_h5_mb1 = newl; Line(l_g1_diel_h5_mb1) = {p_g1_diel_h5_m1, p_g1_diel_h5_b1};
l_g1_diel_h5_mb2 = newl; Line(l_g1_diel_h5_mb2) = {p_g1_diel_h5_m2, p_g1_diel_h5_b2};
l_g1_diel_h5_mb3 = newl; Line(l_g1_diel_h5_mb3) = {p_g1_diel_h5_m3, p_g1_diel_h5_b3};
l_g1_diel_h5_mb4 = newl; Line(l_g1_diel_h5_mb4) = {p_g1_diel_h5_m4, p_g1_diel_h5_b4};
l_g1_diel_h5_mb5 = newl; Line(l_g1_diel_h5_mb5) = {p_g1_diel_h5_m5, p_g1_diel_h5_b5};
// on the bottom plane
l_g1_diel_h1_b1 = newl; Circle(l_g1_diel_h1_b1) = {p_g1_diel_h1_b2, p_g1_diel_h1_b1, p_g1_diel_h1_b3}; // hole 1
l_g1_diel_h1_b2 = newl; Line(l_g1_diel_h1_b2) = {p_g1_diel_h1_b2, p_g1_lc_h1_t2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h1_b3 = newl; Line(l_g1_diel_h1_b3) = {p_g1_diel_h1_b3, p_g1_lc_h1_t3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h2_b1 = newl; Circle(l_g1_diel_h2_b1) = {p_g1_diel_h2_b2, p_g1_diel_h2_b1, p_g1_diel_h2_b3}; // hole 2
l_g1_diel_h2_b2 = newl; Line(l_g1_diel_h2_b2) = {p_g1_diel_h2_b2, p_g1_lc_h2_t2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h2_b3 = newl; Line(l_g1_diel_h2_b3) = {p_g1_diel_h2_b3, p_g1_lc_h2_t3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h3_b1 = newl; Circle(l_g1_diel_h3_b1) = {p_g1_diel_h3_b2, p_g1_diel_h3_b1, p_g1_diel_h3_b3}; // hole 3
l_g1_diel_h3_b2 = newl; Line(l_g1_diel_h3_b2) = {p_g1_diel_h3_b2, p_g1_lc_h3_t2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h3_b3 = newl; Line(l_g1_diel_h3_b3) = {p_g1_diel_h3_b3, p_g1_lc_h3_t3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h4_b1 = newl; Circle(l_g1_diel_h4_b1) = {p_g1_diel_h4_b2, p_g1_diel_h4_b1, p_g1_diel_h4_b3}; // hole 4
l_g1_diel_h4_b2 = newl; Line(l_g1_diel_h4_b2) = {p_g1_diel_h4_b2, p_g1_lc_h4_t2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h4_b3 = newl; Line(l_g1_diel_h4_b3) = {p_g1_diel_h4_b3, p_g1_lc_h4_t3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_b1 = newl; Circle(l_g1_diel_h5_b1) = {p_g1_diel_h5_b2, p_g1_diel_h5_b1, p_g1_diel_h5_b3}; // hole 5
l_g1_diel_h5_b2 = newl; Circle(l_g1_diel_h5_b2) = {p_g1_diel_h5_b3, p_g1_diel_h5_b1, p_g1_diel_h5_b4};
l_g1_diel_h5_b3 = newl; Circle(l_g1_diel_h5_b3) = {p_g1_diel_h5_b4, p_g1_diel_h5_b1, p_g1_diel_h5_b5};
l_g1_diel_h5_b4 = newl; Circle(l_g1_diel_h5_b4) = {p_g1_diel_h5_b5, p_g1_diel_h5_b1, p_g1_diel_h5_b2};
l_g1_diel_h5_b5 = newl; Line(l_g1_diel_h5_b5) = {p_g1_diel_h5_b2, p_g1_lc_h5_t2}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_b6 = newl; Line(l_g1_diel_h5_b6) = {p_g1_diel_h5_b3, p_g1_lc_h5_t3}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_b7 = newl; Line(l_g1_diel_h5_b7) = {p_g1_diel_h5_b4, p_g1_lc_h5_t4}; // this line is on the surface exposed to gas by etching
l_g1_diel_h5_b8 = newl; Line(l_g1_diel_h5_b8) = {p_g1_diel_h5_b5, p_g1_lc_h5_t5}; // this line is on the surface exposed to gas by etching

/*---------------------------------- Curve Loop and Surface Definitions -------------------------------------*/
// top plane of upper copper plate
cl_g1_uc_t1 = newll; Curve Loop(cl_g1_uc_t1) = {l_g1_uc_h1_t1, l_g1_uc_h12_t1, l_g1_uc_h2_t1, l_g1_uc_h23_t1,
    l_g1_uc_h3_t1, l_g1_uc_h34_t1, l_g1_uc_h4_t1, l_g1_uc_h41_t1};
cl_g1_uc_t2 = newll; Curve Loop(cl_g1_uc_t2) = {l_g1_uc_h5_t1, l_g1_uc_h5_t2, l_g1_uc_h5_t3, l_g1_uc_h5_t4}; // curve loop for hole 5
// planes between top and bottom plane
cl_g1_uc_h1_tb1 = newll; Curve Loop(cl_g1_uc_h1_tb1) = {-l_g1_uc_h1_tb2, l_g1_uc_h1_t1, l_g1_uc_h1_tb3, -l_g1_uc_h1_b1}; // cylinderical surface on hole 1
cl_g1_uc_h12_tb1 = newll; Curve Loop(cl_g1_uc_h12_tb1) = {-l_g1_uc_h1_tb3, l_g1_uc_h12_t1, l_g1_uc_h2_tb2, -l_g1_uc_h12_b1}; // plane surface bwt hole 1 and hole  2
cl_g1_uc_h2_tb1 = newll; Curve Loop(cl_g1_uc_h2_tb1) = {-l_g1_uc_h2_tb2, l_g1_uc_h2_t1, l_g1_uc_h2_tb3, -l_g1_uc_h2_b1}; // cylinderical surface on hole 2
cl_g1_uc_h23_tb1 = newll; Curve Loop(cl_g1_uc_h23_tb1) = {-l_g1_uc_h2_tb3, l_g1_uc_h23_t1, l_g1_uc_h3_tb2, -l_g1_uc_h23_b1}; // plane surface bwt hole 2 and hole  3
cl_g1_uc_h3_tb1 = newll; Curve Loop(cl_g1_uc_h3_tb1) = {-l_g1_uc_h3_tb2, l_g1_uc_h3_t1, l_g1_uc_h3_tb3, -l_g1_uc_h3_b1}; // cylinderical surface on hole 3
cl_g1_uc_h34_tb1 = newll; Curve Loop(cl_g1_uc_h34_tb1) = {-l_g1_uc_h3_tb3, l_g1_uc_h34_t1, l_g1_uc_h4_tb2, -l_g1_uc_h34_b1}; // plane surface bwt hole 3 and hole  4
cl_g1_uc_h4_tb1 = newll; Curve Loop(cl_g1_uc_h4_tb1) = {-l_g1_uc_h4_tb2, l_g1_uc_h4_t1, l_g1_uc_h4_tb3, -l_g1_uc_h4_b1}; // cylinderical surface on hole 4
cl_g1_uc_h41_tb1 = newll; Curve Loop(cl_g1_uc_h41_tb1) = {-l_g1_uc_h4_tb3, l_g1_uc_h41_t1, l_g1_uc_h1_tb2, -l_g1_uc_h41_b1}; // plane surface bwt hole 4 and hole  1
cl_g1_uc_h5_tb1 = newll; Curve Loop(cl_g1_uc_h5_tb1) = {-l_g1_uc_h5_tb2, l_g1_uc_h5_t1, l_g1_uc_h5_tb3, -l_g1_uc_h5_b1}; // cylinderical surface on hole 5
cl_g1_uc_h5_tb2 = newll; Curve Loop(cl_g1_uc_h5_tb2) = {-l_g1_uc_h5_tb3, l_g1_uc_h5_t2, l_g1_uc_h5_tb4, -l_g1_uc_h5_b2};
cl_g1_uc_h5_tb3 = newll; Curve Loop(cl_g1_uc_h5_tb3) = {-l_g1_uc_h5_tb4, l_g1_uc_h5_t3, l_g1_uc_h5_tb5, -l_g1_uc_h5_b3};
cl_g1_uc_h5_tb4 = newll; Curve Loop(cl_g1_uc_h5_tb4) = {-l_g1_uc_h5_tb5, l_g1_uc_h5_t4, l_g1_uc_h5_tb2, -l_g1_uc_h5_b4};
// bottom plane of upper copper plate
cl_g1_uc_b1 = newll; Curve Loop(cl_g1_uc_b1) = {l_g1_uc_h1_b1, l_g1_uc_h12_b1, l_g1_uc_h2_b1, l_g1_uc_h23_b1,
    l_g1_uc_h3_b1, l_g1_uc_h34_b1, l_g1_uc_h4_b1, l_g1_uc_h41_b1};
cl_g1_uc_b2 = newll; Curve Loop(cl_g1_uc_b2) = {l_g1_uc_h5_b1, l_g1_uc_h5_b2, l_g1_uc_h5_b3, l_g1_uc_h5_b4}; // curve loop for hole 5
// surface definitions of upper plate
s_g1_uc_t1 = news; Plane Surface(s_g1_uc_t1) = {cl_g1_uc_t1, cl_g1_uc_t2};
s_g1_uc_h1_tb1 = news; Surface(s_g1_uc_h1_tb1) = {cl_g1_uc_h1_tb1}; // cylinderical surface on hole 1
s_g1_uc_h12_tb1 = news; Plane Surface(s_g1_uc_h12_tb1) = {cl_g1_uc_h12_tb1}; // plane surface bwt hole 1 and hole  2
s_g1_uc_h2_tb1 = news; Surface(s_g1_uc_h2_tb1) = {cl_g1_uc_h2_tb1}; // cylinderical surface on hole 2
s_g1_uc_h23_tb1 = news; Plane Surface(s_g1_uc_h23_tb1) = {cl_g1_uc_h23_tb1}; // plane surface bwt hole 2 and hole  3
s_g1_uc_h3_tb1 = news; Surface(s_g1_uc_h3_tb1) = {cl_g1_uc_h3_tb1}; // cylinderical surface on hole 3
s_g1_uc_h34_tb1 = news; Plane Surface(s_g1_uc_h34_tb1) = {cl_g1_uc_h34_tb1}; // plane surface bwt hole 3 and hole  4
s_g1_uc_h4_tb1 = news; Surface(s_g1_uc_h4_tb1) = {cl_g1_uc_h4_tb1}; // cylinderical surface on hole 4
s_g1_uc_h41_tb1 = news; Plane Surface(s_g1_uc_h41_tb1) = {cl_g1_uc_h41_tb1}; // plane surface bwt hole 4 and hole  1
s_g1_uc_h5_tb1 = news; Surface(s_g1_uc_h5_tb1) = {cl_g1_uc_h5_tb1}; // cylinderical surface on hole 5
s_g1_uc_h5_tb2 = news; Surface(s_g1_uc_h5_tb2) = {cl_g1_uc_h5_tb2};
s_g1_uc_h5_tb3 = news; Surface(s_g1_uc_h5_tb3) = {cl_g1_uc_h5_tb3};
s_g1_uc_h5_tb4 = news; Surface(s_g1_uc_h5_tb4) = {cl_g1_uc_h5_tb4};
s_g1_uc_b1 = news; Plane Surface(s_g1_uc_b1) = {cl_g1_uc_b1, cl_g1_uc_b2};

// curve loop definitions of upper copper plate
// top plane of upper copper plate
cl_g1_lc_t1 = newll; Curve Loop(cl_g1_lc_t1) = {l_g1_lc_h1_t1, l_g1_lc_h12_t1, l_g1_lc_h2_t1, l_g1_lc_h23_t1,
    l_g1_lc_h3_t1, l_g1_lc_h34_t1, l_g1_lc_h4_t1, l_g1_lc_h41_t1};
cl_g1_lc_t2 = newll; Curve Loop(cl_g1_lc_t2) = {l_g1_lc_h5_t1, l_g1_lc_h5_t2, l_g1_lc_h5_t3, l_g1_lc_h5_t4}; // curve loop for hole 5
// planes between top and bottom plane
cl_g1_lc_h1_tb1 = newll; Curve Loop(cl_g1_lc_h1_tb1) = {-l_g1_lc_h1_tb2, l_g1_lc_h1_t1, l_g1_lc_h1_tb3, -l_g1_lc_h1_b1}; // cylinderical surface on hole 1
cl_g1_lc_h12_tb1 = newll; Curve Loop(cl_g1_lc_h12_tb1) = {-l_g1_lc_h1_tb3, l_g1_lc_h12_t1, l_g1_lc_h2_tb2, -l_g1_lc_h12_b1}; // plane surface bwt hole 1 and hole  2
cl_g1_lc_h2_tb1 = newll; Curve Loop(cl_g1_lc_h2_tb1) = {-l_g1_lc_h2_tb2, l_g1_lc_h2_t1, l_g1_lc_h2_tb3, -l_g1_lc_h2_b1}; // cylinderical surface on hole 2
cl_g1_lc_h23_tb1 = newll; Curve Loop(cl_g1_lc_h23_tb1) = {-l_g1_lc_h2_tb3, l_g1_lc_h23_t1, l_g1_lc_h3_tb2, -l_g1_lc_h23_b1}; // plane surface bwt hole 2 and hole  3
cl_g1_lc_h3_tb1 = newll; Curve Loop(cl_g1_lc_h3_tb1) = {-l_g1_lc_h3_tb2, l_g1_lc_h3_t1, l_g1_lc_h3_tb3, -l_g1_lc_h3_b1}; // cylinderical surface on hole 3
cl_g1_lc_h34_tb1 = newll; Curve Loop(cl_g1_lc_h34_tb1) = {-l_g1_lc_h3_tb3, l_g1_lc_h34_t1, l_g1_lc_h4_tb2, -l_g1_lc_h34_b1}; // plane surface bwt hole 3 and hole  4
cl_g1_lc_h4_tb1 = newll; Curve Loop(cl_g1_lc_h4_tb1) = {-l_g1_lc_h4_tb2, l_g1_lc_h4_t1, l_g1_lc_h4_tb3, -l_g1_lc_h4_b1}; // cylinderical surface on hole 4
cl_g1_lc_h41_tb1 = newll; Curve Loop(cl_g1_lc_h41_tb1) = {-l_g1_lc_h4_tb3, l_g1_lc_h41_t1, l_g1_lc_h1_tb2, -l_g1_lc_h41_b1}; // plane surface bwt hole 4 and hole  1
cl_g1_lc_h5_tb1 = newll; Curve Loop(cl_g1_lc_h5_tb1) = {-l_g1_lc_h5_tb2, l_g1_lc_h5_t1, l_g1_lc_h5_tb3, -l_g1_lc_h5_b1}; // cylinderical surface on hole 5
cl_g1_lc_h5_tb2 = newll; Curve Loop(cl_g1_lc_h5_tb2) = {-l_g1_lc_h5_tb3, l_g1_lc_h5_t2, l_g1_lc_h5_tb4, -l_g1_lc_h5_b2};
cl_g1_lc_h5_tb3 = newll; Curve Loop(cl_g1_lc_h5_tb3) = {-l_g1_lc_h5_tb4, l_g1_lc_h5_t3, l_g1_lc_h5_tb5, -l_g1_lc_h5_b3};
cl_g1_lc_h5_tb4 = newll; Curve Loop(cl_g1_lc_h5_tb4) = {-l_g1_lc_h5_tb5, l_g1_lc_h5_t4, l_g1_lc_h5_tb2, -l_g1_lc_h5_b4};
// bottom plane of upper copper plate
cl_g1_lc_b1 = newll; Curve Loop(cl_g1_lc_b1) = {l_g1_lc_h1_b1, l_g1_lc_h12_b1, l_g1_lc_h2_b1, l_g1_lc_h23_b1,
    l_g1_lc_h3_b1, l_g1_lc_h34_b1, l_g1_lc_h4_b1, l_g1_lc_h41_b1};
cl_g1_lc_b2 = newll; Curve Loop(cl_g1_lc_b2) = {l_g1_lc_h5_b1, l_g1_lc_h5_b2, l_g1_lc_h5_b3, l_g1_lc_h5_b4}; // curve loop for hole 5
// surface definitions of upper plate
s_g1_lc_t1 = news; Plane Surface(s_g1_lc_t1) = {cl_g1_lc_t1, cl_g1_lc_t2};
s_g1_lc_h1_tb1 = news; Surface(s_g1_lc_h1_tb1) = {cl_g1_lc_h1_tb1}; // cylinderical surface on hole 1
s_g1_lc_h12_tb1 = news; Plane Surface(s_g1_lc_h12_tb1) = {cl_g1_lc_h12_tb1}; // plane surface bwt hole 1 and hole  2
s_g1_lc_h2_tb1 = news; Surface(s_g1_lc_h2_tb1) = {cl_g1_lc_h2_tb1}; // cylinderical surface on hole 2
s_g1_lc_h23_tb1 = news; Plane Surface(s_g1_lc_h23_tb1) = {cl_g1_lc_h23_tb1}; // plane surface bwt hole 2 and hole  3
s_g1_lc_h3_tb1 = news; Surface(s_g1_lc_h3_tb1) = {cl_g1_lc_h3_tb1}; // cylinderical surface on hole 3
s_g1_lc_h34_tb1 = news; Plane Surface(s_g1_lc_h34_tb1) = {cl_g1_lc_h34_tb1}; // plane surface bwt hole 3 and hole  4
s_g1_lc_h4_tb1 = news; Surface(s_g1_lc_h4_tb1) = {cl_g1_lc_h4_tb1}; // cylinderical surface on hole 4
s_g1_lc_h41_tb1 = news; Plane Surface(s_g1_lc_h41_tb1) = {cl_g1_lc_h41_tb1}; // plane surface bwt hole 4 and hole  1
s_g1_lc_h5_tb1 = news; Surface(s_g1_lc_h5_tb1) = {cl_g1_lc_h5_tb1}; // cylinderical surface on hole 5
s_g1_lc_h5_tb2 = news; Surface(s_g1_lc_h5_tb2) = {cl_g1_lc_h5_tb2};
s_g1_lc_h5_tb3 = news; Surface(s_g1_lc_h5_tb3) = {cl_g1_lc_h5_tb3};
s_g1_lc_h5_tb4 = news; Surface(s_g1_lc_h5_tb4) = {cl_g1_lc_h5_tb4};
s_g1_lc_b1 = news; Plane Surface(s_g1_lc_b1) = {cl_g1_lc_b1, cl_g1_lc_b2};

// curve loop definitions of dielectric layer
// top plane surfaces
// a surface exposed to gas by etching
cl_g1_diel_h1_t1 = newll; Curve Loop(cl_g1_diel_h1_t1) = {-l_g1_diel_h1_t2, l_g1_diel_h1_t1, l_g1_diel_h1_t3, -l_g1_uc_h1_b1}; 
cl_g1_diel_h2_t1 = newll; Curve Loop(cl_g1_diel_h2_t1) = {-l_g1_diel_h2_t2, l_g1_diel_h2_t1, l_g1_diel_h2_t3, -l_g1_uc_h2_b1};
cl_g1_diel_h3_t1 = newll; Curve Loop(cl_g1_diel_h3_t1) = {-l_g1_diel_h3_t2, l_g1_diel_h3_t1, l_g1_diel_h3_t3, -l_g1_uc_h3_b1};
cl_g1_diel_h4_t1 = newll; Curve Loop(cl_g1_diel_h4_t1) = {-l_g1_diel_h4_t2, l_g1_diel_h4_t1, l_g1_diel_h4_t3, -l_g1_uc_h4_b1};
cl_g1_diel_h5_t1 = newll; Curve Loop(cl_g1_diel_h5_t1) = {-l_g1_diel_h5_t5, l_g1_diel_h5_t1, l_g1_diel_h5_t6, -l_g1_uc_h5_b1};
cl_g1_diel_h5_t2 = newll; Curve Loop(cl_g1_diel_h5_t2) = {-l_g1_diel_h5_t6, l_g1_diel_h5_t2, l_g1_diel_h5_t7, -l_g1_uc_h5_b2};
cl_g1_diel_h5_t3 = newll; Curve Loop(cl_g1_diel_h5_t3) = {-l_g1_diel_h5_t7, l_g1_diel_h5_t3, l_g1_diel_h5_t8, -l_g1_uc_h5_b3};
cl_g1_diel_h5_t4 = newll; Curve Loop(cl_g1_diel_h5_t4) = {-l_g1_diel_h5_t8, l_g1_diel_h5_t4, l_g1_diel_h5_t5, -l_g1_uc_h5_b4};
// planes between top and middle plane
cl_g1_diel_h1_tm1 = newll; Curve Loop(cl_g1_diel_h1_tm1) = {l_g1_diel_h1_tm2, l_g1_diel_h1_m1, -l_g1_diel_h1_tm3, -l_g1_diel_h1_t1};
cl_g1_diel_h1_tm1 = newll; Curve Loop(cl_g1_diel_h1_tm1) = {l_g1_diel_h1_tm2, l_g1_diel_h1_m1, -l_g1_diel_h1_tm3, -l_g1_diel_h1_t1};
cl_g1_diel_h2_tm1 = newll; Curve Loop(cl_g1_diel_h2_tm1) = {l_g1_diel_h2_tm2, l_g1_diel_h2_m1, -l_g1_diel_h2_tm3, -l_g1_diel_h2_t1};
cl_g1_diel_h3_tm1 = newll; Curve Loop(cl_g1_diel_h3_tm1) = {l_g1_diel_h3_tm2, l_g1_diel_h3_m1, -l_g1_diel_h3_tm3, -l_g1_diel_h3_t1};
cl_g1_diel_h4_tm1 = newll; Curve Loop(cl_g1_diel_h4_tm1) = {l_g1_diel_h4_tm2, l_g1_diel_h4_m1, -l_g1_diel_h4_tm3, -l_g1_diel_h4_t1};
cl_g1_diel_h5_tm1 = newll; Curve Loop(cl_g1_diel_h5_tm1) = {l_g1_diel_h5_tm2, l_g1_diel_h5_m1, -l_g1_diel_h5_tm3, -l_g1_diel_h5_t1};
cl_g1_diel_h5_tm2 = newll; Curve Loop(cl_g1_diel_h5_tm2) = {l_g1_diel_h5_tm3, l_g1_diel_h5_m2, -l_g1_diel_h5_tm4, -l_g1_diel_h5_t2};
cl_g1_diel_h5_tm3 = newll; Curve Loop(cl_g1_diel_h5_tm3) = {l_g1_diel_h5_tm4, l_g1_diel_h5_m3, -l_g1_diel_h5_tm5, -l_g1_diel_h5_t3};
cl_g1_diel_h5_tm4 = newll; Curve Loop(cl_g1_diel_h5_tm4) = {l_g1_diel_h5_tm5, l_g1_diel_h5_m4, -l_g1_diel_h5_tm2, -l_g1_diel_h5_t4};
// planes between middle and bottom plane
cl_g1_diel_h1_mb1 = newll; Curve Loop(cl_g1_diel_h1_mb1) = {l_g1_diel_h1_mb2, l_g1_diel_h1_b1, -l_g1_diel_h1_mb3, -l_g1_diel_h1_m1};
cl_g1_diel_h2_mb1 = newll; Curve Loop(cl_g1_diel_h2_mb1) = {l_g1_diel_h2_mb2, l_g1_diel_h2_b1, -l_g1_diel_h2_mb3, -l_g1_diel_h2_m1};
cl_g1_diel_h3_mb1 = newll; Curve Loop(cl_g1_diel_h3_mb1) = {l_g1_diel_h3_mb2, l_g1_diel_h3_b1, -l_g1_diel_h3_mb3, -l_g1_diel_h3_m1};
cl_g1_diel_h4_mb1 = newll; Curve Loop(cl_g1_diel_h4_mb1) = {l_g1_diel_h4_mb2, l_g1_diel_h4_b1, -l_g1_diel_h4_mb3, -l_g1_diel_h4_m1};
cl_g1_diel_h5_mb1 = newll; Curve Loop(cl_g1_diel_h5_mb1) = {l_g1_diel_h5_mb2, l_g1_diel_h5_b1, -l_g1_diel_h5_mb3, -l_g1_diel_h5_m1};
cl_g1_diel_h5_mb2 = newll; Curve Loop(cl_g1_diel_h5_mb2) = {l_g1_diel_h5_mb3, l_g1_diel_h5_b2, -l_g1_diel_h5_mb4, -l_g1_diel_h5_m2};
cl_g1_diel_h5_mb3 = newll; Curve Loop(cl_g1_diel_h5_mb3) = {l_g1_diel_h5_mb4, l_g1_diel_h5_b3, -l_g1_diel_h5_mb5, -l_g1_diel_h5_m3};
cl_g1_diel_h5_mb4 = newll; Curve Loop(cl_g1_diel_h5_mb4) = {l_g1_diel_h5_mb5, l_g1_diel_h5_b4, -l_g1_diel_h5_mb2, -l_g1_diel_h5_m4};
// planes between top and bottom plane
cl_g1_diel_h12_tb1 = newll; Curve Loop(cl_g1_diel_h12_tb1) = {-l_g1_lc_h12_t1, -l_g1_diel_h1_b3, -l_g1_diel_h1_mb3, -l_g1_diel_h1_tm3, l_g1_diel_h1_t3,
    l_g1_uc_h12_b1, -l_g1_diel_h2_t2, l_g1_diel_h2_tm2, l_g1_diel_h2_mb2, l_g1_diel_h2_b2};
cl_g1_diel_h23_tb1 = newll; Curve Loop(cl_g1_diel_h23_tb1) = {-l_g1_lc_h23_t1, -l_g1_diel_h2_b3, -l_g1_diel_h2_mb3, -l_g1_diel_h2_tm3, l_g1_diel_h2_t3,
    l_g1_uc_h23_b1, -l_g1_diel_h3_t2, l_g1_diel_h3_tm2, l_g1_diel_h3_mb2, l_g1_diel_h3_b2};
cl_g1_diel_h34_tb1 = newll; Curve Loop(cl_g1_diel_h34_tb1) = {-l_g1_lc_h34_t1, -l_g1_diel_h3_b3, -l_g1_diel_h3_mb3, -l_g1_diel_h3_tm3, l_g1_diel_h3_t3,
    l_g1_uc_h34_b1, -l_g1_diel_h4_t2, l_g1_diel_h4_tm2, l_g1_diel_h4_mb2, l_g1_diel_h4_b2};
cl_g1_diel_h41_tb1 = newll; Curve Loop(cl_g1_diel_h41_tb1) = {-l_g1_lc_h41_t1, -l_g1_diel_h4_b3, -l_g1_diel_h4_mb3, -l_g1_diel_h4_tm3, l_g1_diel_h4_t3,
    l_g1_uc_h41_b1, -l_g1_diel_h1_t2, l_g1_diel_h1_tm2, l_g1_diel_h1_mb2, l_g1_diel_h1_b2};
// bottom plane surfaces
// a surface exposed to gas by etching
cl_g1_diel_h1_b1 = newll; Curve Loop(cl_g1_diel_h1_b1) = {-l_g1_diel_h1_b2, l_g1_diel_h1_b1, l_g1_diel_h1_b3, -l_g1_lc_h1_t1}; 
cl_g1_diel_h2_b1 = newll; Curve Loop(cl_g1_diel_h2_b1) = {-l_g1_diel_h2_b2, l_g1_diel_h2_b1, l_g1_diel_h2_b3, -l_g1_lc_h2_t1};
cl_g1_diel_h3_b1 = newll; Curve Loop(cl_g1_diel_h3_b1) = {-l_g1_diel_h3_b2, l_g1_diel_h3_b1, l_g1_diel_h3_b3, -l_g1_lc_h3_t1};
cl_g1_diel_h4_b1 = newll; Curve Loop(cl_g1_diel_h4_b1) = {-l_g1_diel_h4_b2, l_g1_diel_h4_b1, l_g1_diel_h4_b3, -l_g1_lc_h4_t1};
cl_g1_diel_h5_b1 = newll; Curve Loop(cl_g1_diel_h5_b1) = {-l_g1_diel_h5_b5, l_g1_diel_h5_b1, l_g1_diel_h5_b6, -l_g1_lc_h5_t1};
cl_g1_diel_h5_b2 = newll; Curve Loop(cl_g1_diel_h5_b2) = {-l_g1_diel_h5_b6, l_g1_diel_h5_b2, l_g1_diel_h5_b7, -l_g1_lc_h5_t2};
cl_g1_diel_h5_b3 = newll; Curve Loop(cl_g1_diel_h5_b3) = {-l_g1_diel_h5_b7, l_g1_diel_h5_b3, l_g1_diel_h5_b8, -l_g1_lc_h5_t3};
cl_g1_diel_h5_b4 = newll; Curve Loop(cl_g1_diel_h5_b4) = {-l_g1_diel_h5_b8, l_g1_diel_h5_b4, l_g1_diel_h5_b5, -l_g1_lc_h5_t4};

// surface definitions of dielectric layer
s_g1_diel_h1_t1 = news; Plane Surface(s_g1_diel_h1_t1) = {cl_g1_diel_h1_t1};
s_g1_diel_h2_t1 = news; Plane Surface(s_g1_diel_h2_t1) = {cl_g1_diel_h2_t1};
s_g1_diel_h3_t1 = news; Plane Surface(s_g1_diel_h3_t1) = {cl_g1_diel_h3_t1};
s_g1_diel_h4_t1 = news; Plane Surface(s_g1_diel_h4_t1) = {cl_g1_diel_h4_t1};
s_g1_diel_h5_t1 = news; Plane Surface(s_g1_diel_h5_t1) = {cl_g1_diel_h5_t1};
s_g1_diel_h5_t2 = news; Plane Surface(s_g1_diel_h5_t2) = {cl_g1_diel_h5_t2};
s_g1_diel_h5_t3 = news; Plane Surface(s_g1_diel_h5_t3) = {cl_g1_diel_h5_t3};
s_g1_diel_h5_t4 = news; Plane Surface(s_g1_diel_h5_t4) = {cl_g1_diel_h5_t4};
s_g1_diel_h1_tm1 = news; Surface(s_g1_diel_h1_tm1) = {cl_g1_diel_h1_tm1};
s_g1_diel_h1_tm1 = news; Surface(s_g1_diel_h1_tm1) = {cl_g1_diel_h1_tm1};
s_g1_diel_h2_tm1 = news; Surface(s_g1_diel_h2_tm1) = {cl_g1_diel_h2_tm1};
s_g1_diel_h3_tm1 = news; Surface(s_g1_diel_h3_tm1) = {cl_g1_diel_h3_tm1};
s_g1_diel_h4_tm1 = news; Surface(s_g1_diel_h4_tm1) = {cl_g1_diel_h4_tm1};
s_g1_diel_h5_tm1 = news; Surface(s_g1_diel_h5_tm1) = {cl_g1_diel_h5_tm1};
s_g1_diel_h5_tm2 = news; Surface(s_g1_diel_h5_tm2) = {cl_g1_diel_h5_tm2};
s_g1_diel_h5_tm3 = news; Surface(s_g1_diel_h5_tm3) = {cl_g1_diel_h5_tm3};
s_g1_diel_h5_tm4 = news; Surface(s_g1_diel_h5_tm4) = {cl_g1_diel_h5_tm4};
s_g1_diel_h1_mb1 = news; Surface(s_g1_diel_h1_mb1) = {cl_g1_diel_h1_mb1};
s_g1_diel_h2_mb1 = news; Surface(s_g1_diel_h2_mb1) = {cl_g1_diel_h2_mb1};
s_g1_diel_h3_mb1 = news; Surface(s_g1_diel_h3_mb1) = {cl_g1_diel_h3_mb1};
s_g1_diel_h4_mb1 = news; Surface(s_g1_diel_h4_mb1) = {cl_g1_diel_h4_mb1};
s_g1_diel_h5_mb1 = news; Surface(s_g1_diel_h5_mb1) = {cl_g1_diel_h5_mb1};
s_g1_diel_h5_mb2 = news; Surface(s_g1_diel_h5_mb2) = {cl_g1_diel_h5_mb2};
s_g1_diel_h5_mb3 = news; Surface(s_g1_diel_h5_mb3) = {cl_g1_diel_h5_mb3};
s_g1_diel_h5_mb4 = news; Surface(s_g1_diel_h5_mb4) = {cl_g1_diel_h5_mb4};
s_g1_diel_h12_tb1 = news; Plane Surface(s_g1_diel_h12_tb1) = {cl_g1_diel_h12_tb1};
s_g1_diel_h23_tb1 = news; Plane Surface(s_g1_diel_h23_tb1) = {cl_g1_diel_h23_tb1};
s_g1_diel_h34_tb1 = news; Plane Surface(s_g1_diel_h34_tb1) = {cl_g1_diel_h34_tb1};
s_g1_diel_h41_tb1 = news; Plane Surface(s_g1_diel_h41_tb1) = {cl_g1_diel_h41_tb1};
s_g1_diel_h1_b1 = news; Plane Surface(s_g1_diel_h1_b1) = {cl_g1_diel_h1_b1};
s_g1_diel_h2_b1 = news; Plane Surface(s_g1_diel_h2_b1) = {cl_g1_diel_h2_b1};
s_g1_diel_h3_b1 = news; Plane Surface(s_g1_diel_h3_b1) = {cl_g1_diel_h3_b1};
s_g1_diel_h4_b1 = news; Plane Surface(s_g1_diel_h4_b1) = {cl_g1_diel_h4_b1};
s_g1_diel_h5_b1 = news; Plane Surface(s_g1_diel_h5_b1) = {cl_g1_diel_h5_b1};
s_g1_diel_h5_b2 = news; Plane Surface(s_g1_diel_h5_b2) = {cl_g1_diel_h5_b2};
s_g1_diel_h5_b3 = news; Plane Surface(s_g1_diel_h5_b3) = {cl_g1_diel_h5_b3};
s_g1_diel_h5_b4 = news; Plane Surface(s_g1_diel_h5_b4) = {cl_g1_diel_h5_b4};
/*---------------------------------- Surface Loop and Volume Definitions -------------------------------------*/
// surface loop definitions of upper plate
sl_g1_uc = newsl; Surface Loop(sl_g1_uc) = {s_g1_uc_t1,
    s_g1_uc_h1_tb1, s_g1_uc_h2_tb1, s_g1_uc_h3_tb1, s_g1_uc_h4_tb1,
    s_g1_uc_h12_tb1, s_g1_uc_h23_tb1, s_g1_uc_h34_tb1, s_g1_uc_h41_tb1,
    s_g1_uc_h5_tb1, s_g1_uc_h5_tb2, s_g1_uc_h5_tb3, s_g1_uc_h5_tb4,
    s_g1_uc_b1};
// surface loop definitions of upper plate
sl_g1_lc = newsl; Surface Loop(sl_g1_lc) = {s_g1_lc_t1,
    s_g1_lc_h1_tb1, s_g1_lc_h2_tb1, s_g1_lc_h3_tb1, s_g1_lc_h4_tb1,
    s_g1_lc_h12_tb1, s_g1_lc_h23_tb1, s_g1_lc_h34_tb1, s_g1_lc_h41_tb1,
    s_g1_lc_h5_tb1, s_g1_lc_h5_tb2, s_g1_lc_h5_tb3, s_g1_lc_h5_tb4,
    s_g1_lc_b1};
// surface loop definitions of dielectric layer
sl_g1_diel = newsl; Surface Loop(sl_g1_diel) = {
    s_g1_uc_b1,
    s_g1_diel_h1_t1, s_g1_diel_h2_t1, s_g1_diel_h3_t1, s_g1_diel_h4_t1,
    s_g1_diel_h5_t1, s_g1_diel_h5_t2, s_g1_diel_h5_t3, s_g1_diel_h5_t4,
    s_g1_diel_h1_tm1, s_g1_diel_h2_tm1, s_g1_diel_h3_tm1, s_g1_diel_h4_tm1,
    s_g1_diel_h5_tm1, s_g1_diel_h5_tm2, s_g1_diel_h5_tm3, s_g1_diel_h5_tm4,
    s_g1_diel_h1_mb1, s_g1_diel_h2_mb1, s_g1_diel_h3_mb1, s_g1_diel_h4_mb1,
    s_g1_diel_h5_mb1, s_g1_diel_h5_mb2, s_g1_diel_h5_mb3, s_g1_diel_h5_mb4,
    s_g1_diel_h12_tb1, s_g1_diel_h23_tb1, s_g1_diel_h34_tb1, s_g1_diel_h41_tb1,
    s_g1_diel_h1_b1, s_g1_diel_h2_b1, s_g1_diel_h3_b1, s_g1_diel_h4_b1,
    s_g1_diel_h5_b1, s_g1_diel_h5_b2, s_g1_diel_h5_b3, s_g1_diel_h5_b4,
    s_g1_lc_t1
};
v_g1_uc = newv; Volume(v_g1_uc) = {sl_g1_uc};
v_g1_lc = newv; Volume(v_g1_lc) = {sl_g1_lc};
v_g1_diel = newv; Volume(v_g1_diel) = {sl_g1_diel};
// ---------------------------------------------------------------------------------------------------------------
// GEM2
// ---------------------------------------------------------------------------------------------------------------
/*---------------------------------- Point Definitions -------------------------------------*/
// points on upper copper plate(top and bottom plane)
// hole1 (bottom left quarter hole)
p_g2_uc_h1_t1 = newp; Point(p_g2_uc_h1_t1) = {-bdrHalfX, -bdrHalfY,       z2 + tD/2 + tC, mshSizeHoleCenter};
p_g2_uc_h1_t2 = newp; Point(p_g2_uc_h1_t2) = {-bdrHalfX, -bdrHalfY + rE,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h1_t3 = newp; Point(p_g2_uc_h1_t3) = {-bdrHalfX + rE, -bdrHalfY,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h1_b1 = newp; Point(p_g2_uc_h1_b1) = {-bdrHalfX, -bdrHalfY,       z2 + tD/2, mshSizeHoleCenter};
p_g2_uc_h1_b2 = newp; Point(p_g2_uc_h1_b2) = {-bdrHalfX, -bdrHalfY + rE,  z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h1_b3 = newp; Point(p_g2_uc_h1_b3) = {-bdrHalfX + rE, -bdrHalfY,  z2 + tD/2, mshSizeEtchHole};
// hole2 (bottom right quarter hole)
p_g2_uc_h2_t1 = newp; Point(p_g2_uc_h2_t1) = {bdrHalfX, -bdrHalfY,        z2 + tD/2 + tC, mshSizeHoleCenter};
p_g2_uc_h2_t2 = newp; Point(p_g2_uc_h2_t2) = {bdrHalfX - rE, -bdrHalfY,   z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h2_t3 = newp; Point(p_g2_uc_h2_t3) = {bdrHalfX, -bdrHalfY + rE,   z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h2_b1 = newp; Point(p_g2_uc_h2_b1) = {bdrHalfX, -bdrHalfY,        z2 + tD/2, mshSizeHoleCenter};
p_g2_uc_h2_b2 = newp; Point(p_g2_uc_h2_b2) = {bdrHalfX - rE, -bdrHalfY,   z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h2_b3 = newp; Point(p_g2_uc_h2_b3) = {bdrHalfX, -bdrHalfY + rE,   z2 + tD/2, mshSizeEtchHole};
// hole3 (top right quarter hole)
p_g2_uc_h3_t1 = newp; Point(p_g2_uc_h3_t1) = {bdrHalfX, bdrHalfY,       z2 + tD/2 + tC, mshSizeHoleCenter};
p_g2_uc_h3_t2 = newp; Point(p_g2_uc_h3_t2) = {bdrHalfX, bdrHalfY - rE,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h3_t3 = newp; Point(p_g2_uc_h3_t3) = {bdrHalfX - rE, bdrHalfY,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h3_b1 = newp; Point(p_g2_uc_h3_b1) = {bdrHalfX, bdrHalfY,       z2 + tD/2, mshSizeHoleCenter};
p_g2_uc_h3_b2 = newp; Point(p_g2_uc_h3_b2) = {bdrHalfX, bdrHalfY - rE,  z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h3_b3 = newp; Point(p_g2_uc_h3_b3) = {bdrHalfX - rE, bdrHalfY,  z2 + tD/2, mshSizeEtchHole};
// hole4 (top left quarter hole)
p_g2_uc_h4_t1 = newp; Point(p_g2_uc_h4_t1) = {-bdrHalfX, bdrHalfY,       z2 + tD/2 + tC, mshSizeHoleCenter};
p_g2_uc_h4_t2 = newp; Point(p_g2_uc_h4_t2) = {-bdrHalfX + rE, bdrHalfY,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h4_t3 = newp; Point(p_g2_uc_h4_t3) = {-bdrHalfX, bdrHalfY - rE,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h4_b1 = newp; Point(p_g2_uc_h4_b1) = {-bdrHalfX, bdrHalfY,       z2 + tD/2, mshSizeHoleCenter};
p_g2_uc_h4_b2 = newp; Point(p_g2_uc_h4_b2) = {-bdrHalfX + rE, bdrHalfY,  z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h4_b3 = newp; Point(p_g2_uc_h4_b3) = {-bdrHalfX, bdrHalfY - rE,  z2 + tD/2, mshSizeEtchHole};
// hole5 (center circular hole)
p_g2_uc_h5_t1 = newp; Point(p_g2_uc_h5_t1) = {0, 0,   z2 + tD/2 + tC, mshSizeHoleCenter};
p_g2_uc_h5_t2 = newp; Point(p_g2_uc_h5_t2) = {rE, 0,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h5_t3 = newp; Point(p_g2_uc_h5_t3) = {0, rE,  z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h5_t4 = newp; Point(p_g2_uc_h5_t4) = {-rE, 0, z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h5_t5 = newp; Point(p_g2_uc_h5_t5) = {0, -rE, z2 + tD/2 + tC, mshSizeEtchHole};
p_g2_uc_h5_b1 = newp; Point(p_g2_uc_h5_b1) = {0, 0,   z2 + tD/2, mshSizeHoleCenter};
p_g2_uc_h5_b2 = newp; Point(p_g2_uc_h5_b2) = {rE, 0,  z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h5_b3 = newp; Point(p_g2_uc_h5_b3) = {0, rE,  z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h5_b4 = newp; Point(p_g2_uc_h5_b4) = {-rE, 0, z2 + tD/2, mshSizeEtchHole};
p_g2_uc_h5_b5 = newp; Point(p_g2_uc_h5_b5) = {0, -rE, z2 + tD/2, mshSizeEtchHole};

// points on lower copper plate (top and bottom plane)
p_g2_lc_h1_t1 = newp; Point(p_g2_lc_h1_t1) = {-bdrHalfX, -bdrHalfY,       z2 - tD/2, mshSizeHoleCenter};
p_g2_lc_h1_t2 = newp; Point(p_g2_lc_h1_t2) = {-bdrHalfX, -bdrHalfY + rE,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h1_t3 = newp; Point(p_g2_lc_h1_t3) = {-bdrHalfX + rE, -bdrHalfY,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h1_b1 = newp; Point(p_g2_lc_h1_b1) = {-bdrHalfX, -bdrHalfY,       z2 - (tD/2 + tC), mshSizeHoleCenter};
p_g2_lc_h1_b2 = newp; Point(p_g2_lc_h1_b2) = {-bdrHalfX, -bdrHalfY + rE,  z2 - (tD/2 + tC), mshSizeEtchHole};
p_g2_lc_h1_b3 = newp; Point(p_g2_lc_h1_b3) = {-bdrHalfX + rE, -bdrHalfY,  z2 - (tD/2 + tC), mshSizeEtchHole};
// hole2 (bottom right quarter hole)
p_g2_lc_h2_t1 = newp; Point(p_g2_lc_h2_t1) = {bdrHalfX, -bdrHalfY,        z2 - tD/2, mshSizeHoleCenter};
p_g2_lc_h2_t2 = newp; Point(p_g2_lc_h2_t2) = {bdrHalfX - rE, -bdrHalfY,   z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h2_t3 = newp; Point(p_g2_lc_h2_t3) = {bdrHalfX, -bdrHalfY + rE,   z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h2_b1 = newp; Point(p_g2_lc_h2_b1) = {bdrHalfX, -bdrHalfY,        z2 - (tD/2 + tC), mshSizeHoleCenter};
p_g2_lc_h2_b2 = newp; Point(p_g2_lc_h2_b2) = {bdrHalfX - rE, -bdrHalfY,   z2 - (tD/2 + tC), mshSizeEtchHole};
p_g2_lc_h2_b3 = newp; Point(p_g2_lc_h2_b3) = {bdrHalfX, -bdrHalfY + rE,   z2 - (tD/2 + tC), mshSizeEtchHole};
// hole3 (top right quarter hole)
p_g2_lc_h3_t1 = newp; Point(p_g2_lc_h3_t1) = {bdrHalfX, bdrHalfY,       z2 - tD/2, mshSizeHoleCenter};
p_g2_lc_h3_t2 = newp; Point(p_g2_lc_h3_t2) = {bdrHalfX, bdrHalfY - rE,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h3_t3 = newp; Point(p_g2_lc_h3_t3) = {bdrHalfX - rE, bdrHalfY,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h3_b1 = newp; Point(p_g2_lc_h3_b1) = {bdrHalfX, bdrHalfY,       z2 - (tD/2 + tC), mshSizeHoleCenter};
p_g2_lc_h3_b2 = newp; Point(p_g2_lc_h3_b2) = {bdrHalfX, bdrHalfY - rE,  z2 - (tD/2 + tC), mshSizeEtchHole};
p_g2_lc_h3_b3 = newp; Point(p_g2_lc_h3_b3) = {bdrHalfX - rE, bdrHalfY,  z2 - (tD/2 + tC), mshSizeEtchHole};
// hole4 (top left quarter hole)
p_g2_lc_h4_t1 = newp; Point(p_g2_lc_h4_t1) = {-bdrHalfX, bdrHalfY,       z2 - tD/2, mshSizeHoleCenter};
p_g2_lc_h4_t2 = newp; Point(p_g2_lc_h4_t2) = {-bdrHalfX + rE, bdrHalfY,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h4_t3 = newp; Point(p_g2_lc_h4_t3) = {-bdrHalfX, bdrHalfY - rE,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h4_b1 = newp; Point(p_g2_lc_h4_b1) = {-bdrHalfX, bdrHalfY,       z2 - (tD/2 + tC), mshSizeHoleCenter};
p_g2_lc_h4_b2 = newp; Point(p_g2_lc_h4_b2) = {-bdrHalfX + rE, bdrHalfY,  z2 - (tD/2 + tC), mshSizeEtchHole};
p_g2_lc_h4_b3 = newp; Point(p_g2_lc_h4_b3) = {-bdrHalfX, bdrHalfY - rE,  z2 - (tD/2 + tC), mshSizeEtchHole};
// hole5 (center circular quarter hole)
p_g2_lc_h5_t1 = newp; Point(p_g2_lc_h5_t1) = {0, 0,   z2 - tD/2, mshSizeHoleCenter};
p_g2_lc_h5_t2 = newp; Point(p_g2_lc_h5_t2) = {rE, 0,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h5_t3 = newp; Point(p_g2_lc_h5_t3) = {0, rE,  z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h5_t4 = newp; Point(p_g2_lc_h5_t4) = {-rE, 0, z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h5_t5 = newp; Point(p_g2_lc_h5_t5) = {0, -rE, z2 - tD/2, mshSizeEtchHole};
p_g2_lc_h5_b1 = newp; Point(p_g2_lc_h5_b1) = {0, 0,   z2 - tD/2 - tC, mshSizeHoleCenter};
p_g2_lc_h5_b2 = newp; Point(p_g2_lc_h5_b2) = {rE, 0,  z2 - tD/2 - tC, mshSizeEtchHole};
p_g2_lc_h5_b3 = newp; Point(p_g2_lc_h5_b3) = {0, rE,  z2 - tD/2 - tC, mshSizeEtchHole};
p_g2_lc_h5_b4 = newp; Point(p_g2_lc_h5_b4) = {-rE, 0, z2 - tD/2 - tC, mshSizeEtchHole};
p_g2_lc_h5_b5 = newp; Point(p_g2_lc_h5_b5) = {0, -rE, z2 - tD/2 - tC, mshSizeEtchHole};

// points on upper dielectric layer(top, middel and bottom plane)
// hole1 (bottom left quarter hole)
p_g2_diel_h1_t1 = p_g2_uc_h1_b1;
p_g2_diel_h1_t2 = newp; Point(p_g2_diel_h1_t2) = {-bdrHalfX, -bdrHalfY + rO,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h1_t3 = newp; Point(p_g2_diel_h1_t3) = {-bdrHalfX + rO, -bdrHalfY,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h1_m1 = newp; Point(p_g2_diel_h1_m1) = {-bdrHalfX, -bdrHalfY,       z2, mshSizeHoleCenter};
p_g2_diel_h1_m2 = newp; Point(p_g2_diel_h1_m2) = {-bdrHalfX, -bdrHalfY + rI,  z2, mshSizeDielInnerHole};
p_g2_diel_h1_m3 = newp; Point(p_g2_diel_h1_m3) = {-bdrHalfX + rI, -bdrHalfY,  z2, mshSizeDielInnerHole};
p_g2_diel_h1_b1 = p_g2_lc_h1_t1;
p_g2_diel_h1_b2 = newp; Point(p_g2_diel_h1_b2) = {-bdrHalfX, -bdrHalfY + rO,  z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h1_b3 = newp; Point(p_g2_diel_h1_b3) = {-bdrHalfX + rO, -bdrHalfY,  z2 - tD/2, mshSizeDielOuterHole};
// hole2 (bottom right quarter hole)
p_g2_diel_h2_t1 = p_g2_uc_h2_b1;
p_g2_diel_h2_t2 = newp; Point(p_g2_diel_h2_t2) = {bdrHalfX - rO, -bdrHalfY,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h2_t3 = newp; Point(p_g2_diel_h2_t3) = {bdrHalfX, -bdrHalfY + rO,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h2_m1 = newp; Point(p_g2_diel_h2_m1) = {bdrHalfX, -bdrHalfY,       z2, mshSizeHoleCenter};
p_g2_diel_h2_m2 = newp; Point(p_g2_diel_h2_m2) = {bdrHalfX - rI, -bdrHalfY,  z2, mshSizeDielInnerHole};
p_g2_diel_h2_m3 = newp; Point(p_g2_diel_h2_m3) = {bdrHalfX, -bdrHalfY + rI,  z2, mshSizeDielInnerHole};
p_g2_diel_h2_b1 = p_g2_lc_h2_t1;
p_g2_diel_h2_b2 = newp; Point(p_g2_diel_h2_b2) = {bdrHalfX - rO, -bdrHalfY,  z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h2_b3 = newp; Point(p_g2_diel_h2_b3) = {bdrHalfX, -bdrHalfY + rO,  z2 - tD/2, mshSizeDielOuterHole};
// hole3 (top right quarter hole)
p_g2_diel_h3_t1 = p_g2_uc_h3_b1;
p_g2_diel_h3_t2 = newp; Point(p_g2_diel_h3_t2) = {bdrHalfX, bdrHalfY - rO,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h3_t3 = newp; Point(p_g2_diel_h3_t3) = {bdrHalfX - rO, bdrHalfY,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h3_m1 = newp; Point(p_g2_diel_h3_m1) = {bdrHalfX, bdrHalfY,       z2, mshSizeHoleCenter};
p_g2_diel_h3_m2 = newp; Point(p_g2_diel_h3_m2) = {bdrHalfX, bdrHalfY - rI,  z2, mshSizeDielOuterHole};
p_g2_diel_h3_m3 = newp; Point(p_g2_diel_h3_m3) = {bdrHalfX - rI, bdrHalfY,  z2, mshSizeDielOuterHole};
p_g2_diel_h3_b1 = p_g2_lc_h3_t1;
p_g2_diel_h3_b2 = newp; Point(p_g2_diel_h3_b2) = {bdrHalfX, bdrHalfY - rO,  z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h3_b3 = newp; Point(p_g2_diel_h3_b3) = {bdrHalfX - rO, bdrHalfY,  z2 - tD/2, mshSizeDielOuterHole};
// hole4 (top left quarter hole)
p_g2_diel_h4_t1 = p_g2_uc_h4_b1;
p_g2_diel_h4_t2 = newp; Point(p_g2_diel_h4_t2) = {-bdrHalfX + rO, bdrHalfY,    z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h4_t3 = newp; Point(p_g2_diel_h4_t3) = {-bdrHalfX, bdrHalfY - rO,    z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h4_m1 = newp; Point(p_g2_diel_h4_m1) = {-bdrHalfX, bdrHalfY,         z2, mshSizeHoleCenter};
p_g2_diel_h4_m2 = newp; Point(p_g2_diel_h4_m2) = {-bdrHalfX + rI, bdrHalfY,    z2, mshSizeDielInnerHole};
p_g2_diel_h4_m3 = newp; Point(p_g2_diel_h4_m3) = {-bdrHalfX, bdrHalfY - rI ,   z2, mshSizeDielInnerHole};
p_g2_diel_h4_b1 = p_g2_lc_h4_t1;
p_g2_diel_h4_b2 = newp; Point(p_g2_diel_h4_b2) = {-bdrHalfX + rO, bdrHalfY,    z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h4_b3 = newp; Point(p_g2_diel_h4_b3) = {-bdrHalfX, bdrHalfY - rO,    z2 - tD/2, mshSizeDielOuterHole};
// hole5 (center circular quarter hole)
p_g2_diel_h5_t1 = p_g2_uc_h5_b1;
p_g2_diel_h5_t2 = newp; Point(p_g2_diel_h5_t2) = {rO, 0,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_t3 = newp; Point(p_g2_diel_h5_t3) = {0, rO,  z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_t4 = newp; Point(p_g2_diel_h5_t4) = {-rO, 0, z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_t5 = newp; Point(p_g2_diel_h5_t5) = {0, -rO, z2 + tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_m1 = newp; Point(p_g2_diel_h5_m1) = {0, 0,   z2, mshSizeHoleCenter};
p_g2_diel_h5_m2 = newp; Point(p_g2_diel_h5_m2) = {rI, 0,  z2, mshSizeDielInnerHole};
p_g2_diel_h5_m3 = newp; Point(p_g2_diel_h5_m3) = {0, rI,  z2, mshSizeDielInnerHole};
p_g2_diel_h5_m4 = newp; Point(p_g2_diel_h5_m4) = {-rI, 0, z2, mshSizeDielInnerHole};
p_g2_diel_h5_m5 = newp; Point(p_g2_diel_h5_m5) = {0, -rI, z2, mshSizeDielInnerHole};
p_g2_diel_h5_b1 = p_g2_lc_h5_t1;
p_g2_diel_h5_b2 = newp; Point(p_g2_diel_h5_b2) = {rO, 0,  z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_b3 = newp; Point(p_g2_diel_h5_b3) = {0, rO,  z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_b4 = newp; Point(p_g2_diel_h5_b4) = {-rO, 0, z2 - tD/2, mshSizeDielOuterHole};
p_g2_diel_h5_b5 = newp; Point(p_g2_diel_h5_b5) = {0, -rO, z2 - tD/2, mshSizeDielOuterHole};
/*---------------------------------- Line Definitions -------------------------------------*/
// lines on upper plane
// on the top plane
l_g2_uc_h1_t1 = newl; Circle(l_g2_uc_h1_t1) = {p_g2_uc_h1_t2, p_g2_uc_h1_t1, p_g2_uc_h1_t3}; // hole 1
l_g2_uc_h12_t1 = newl; Line(l_g2_uc_h12_t1) = {p_g2_uc_h1_t3, p_g2_uc_h2_t2}; // between hole 1 and hole 2
l_g2_uc_h2_t1 = newl; Circle(l_g2_uc_h2_t1) = {p_g2_uc_h2_t2, p_g2_uc_h2_t1, p_g2_uc_h2_t3}; // hole 2
l_g2_uc_h23_t1 = newl; Line(l_g2_uc_h23_t1) = {p_g2_uc_h2_t3, p_g2_uc_h3_t2}; // between hole 2 and hole 3
l_g2_uc_h3_t1 = newl; Circle(l_g2_uc_h3_t1) = {p_g2_uc_h3_t2, p_g2_uc_h3_t1, p_g2_uc_h3_t3}; // hole 3
l_g2_uc_h34_t1 = newl; Line(l_g2_uc_h34_t1) = {p_g2_uc_h3_t3, p_g2_uc_h4_t2}; // between hole 3 and hole 4
l_g2_uc_h4_t1 = newl; Circle(l_g2_uc_h4_t1) = {p_g2_uc_h4_t2, p_g2_uc_h4_t1, p_g2_uc_h4_t3}; // hole 4
l_g2_uc_h41_t1 = newl; Line(l_g2_uc_h41_t1) = {p_g2_uc_h4_t3, p_g2_uc_h1_t2}; // between hole 4 and hole 1
l_g2_uc_h5_t1 = newl; Circle(l_g2_uc_h5_t1) = {p_g2_uc_h5_t2, p_g2_uc_h5_t1, p_g2_uc_h5_t3}; // hole 5
l_g2_uc_h5_t2 = newl; Circle(l_g2_uc_h5_t2) = {p_g2_uc_h5_t3, p_g2_uc_h5_t1, p_g2_uc_h5_t4};
l_g2_uc_h5_t3 = newl; Circle(l_g2_uc_h5_t3) = {p_g2_uc_h5_t4, p_g2_uc_h5_t1, p_g2_uc_h5_t5};
l_g2_uc_h5_t4 = newl; Circle(l_g2_uc_h5_t4) = {p_g2_uc_h5_t5, p_g2_uc_h5_t1, p_g2_uc_h5_t2};
// on the bottom plane
l_g2_uc_h1_b1 = newl; Circle(l_g2_uc_h1_b1) = {p_g2_uc_h1_b2, p_g2_uc_h1_b1, p_g2_uc_h1_b3}; // hole 1
l_g2_uc_h12_b1 = newl; Line(l_g2_uc_h12_b1) = {p_g2_uc_h1_b3, p_g2_uc_h2_b2}; // between hole 1 and hole 2
l_g2_uc_h2_b1 = newl; Circle(l_g2_uc_h2_b1) = {p_g2_uc_h2_b2, p_g2_uc_h2_b1, p_g2_uc_h2_b3}; // hole 2
l_g2_uc_h23_b1 = newl; Line(l_g2_uc_h23_b1) = {p_g2_uc_h2_b3, p_g2_uc_h3_b2}; // between hole 2 and hole 3
l_g2_uc_h3_b1 = newl; Circle(l_g2_uc_h3_b1) = {p_g2_uc_h3_b2, p_g2_uc_h3_b1, p_g2_uc_h3_b3}; // hole 3
l_g2_uc_h34_b1 = newl; Line(l_g2_uc_h34_b1) = {p_g2_uc_h3_b3, p_g2_uc_h4_b2}; // between hole 3 and hole 4
l_g2_uc_h4_b1 = newl; Circle(l_g2_uc_h4_b1) = {p_g2_uc_h4_b2, p_g2_uc_h4_b1, p_g2_uc_h4_b3}; // hole 4
l_g2_uc_h41_b1 = newl; Line(l_g2_uc_h41_b1) = {p_g2_uc_h4_b3, p_g2_uc_h1_b2}; // between hole 4 and hole 1
l_g2_uc_h5_b1 = newl; Circle(l_g2_uc_h5_b1) = {p_g2_uc_h5_b2, p_g2_uc_h5_b1, p_g2_uc_h5_b3}; // hole 5
l_g2_uc_h5_b2 = newl; Circle(l_g2_uc_h5_b2) = {p_g2_uc_h5_b3, p_g2_uc_h5_b1, p_g2_uc_h5_b4};
l_g2_uc_h5_b3 = newl; Circle(l_g2_uc_h5_b3) = {p_g2_uc_h5_b4, p_g2_uc_h5_b1, p_g2_uc_h5_b5};
l_g2_uc_h5_b4 = newl; Circle(l_g2_uc_h5_b4) = {p_g2_uc_h5_b5, p_g2_uc_h5_b1, p_g2_uc_h5_b2};
// connecting top and bottom
l_g2_uc_h1_tb1 = newl; Line(l_g2_uc_h1_tb1) = {p_g2_uc_h1_t1, p_g2_uc_h1_b1}; // hole 1
l_g2_uc_h1_tb2 = newl; Line(l_g2_uc_h1_tb2) = {p_g2_uc_h1_t2, p_g2_uc_h1_b2}; 
l_g2_uc_h1_tb3 = newl; Line(l_g2_uc_h1_tb3) = {p_g2_uc_h1_t3, p_g2_uc_h1_b3};
l_g2_uc_h2_tb1 = newl; Line(l_g2_uc_h2_tb1) = {p_g2_uc_h2_t1, p_g2_uc_h2_b1}; // hole 2
l_g2_uc_h2_tb2 = newl; Line(l_g2_uc_h2_tb2) = {p_g2_uc_h2_t2, p_g2_uc_h2_b2};
l_g2_uc_h2_tb3 = newl; Line(l_g2_uc_h2_tb3) = {p_g2_uc_h2_t3, p_g2_uc_h2_b3};
l_g2_uc_h3_tb1 = newl; Line(l_g2_uc_h3_tb1) = {p_g2_uc_h3_t1, p_g2_uc_h3_b1}; // hole 3
l_g2_uc_h3_tb2 = newl; Line(l_g2_uc_h3_tb2) = {p_g2_uc_h3_t2, p_g2_uc_h3_b2};
l_g2_uc_h3_tb3 = newl; Line(l_g2_uc_h3_tb3) = {p_g2_uc_h3_t3, p_g2_uc_h3_b3};
l_g2_uc_h4_tb1 = newl; Line(l_g2_uc_h4_tb1) = {p_g2_uc_h4_t1, p_g2_uc_h4_b1}; // hole 4
l_g2_uc_h4_tb2 = newl; Line(l_g2_uc_h4_tb2) = {p_g2_uc_h4_t2, p_g2_uc_h4_b2};
l_g2_uc_h4_tb3 = newl; Line(l_g2_uc_h4_tb3) = {p_g2_uc_h4_t3, p_g2_uc_h4_b3};
l_g2_uc_h5_tb1 = newl; Line(l_g2_uc_h5_tb1) = {p_g2_uc_h5_t1, p_g2_uc_h5_b1}; // hole 5
l_g2_uc_h5_tb2 = newl; Line(l_g2_uc_h5_tb2) = {p_g2_uc_h5_t2, p_g2_uc_h5_b2};
l_g2_uc_h5_tb3 = newl; Line(l_g2_uc_h5_tb3) = {p_g2_uc_h5_t3, p_g2_uc_h5_b3};
l_g2_uc_h5_tb4 = newl; Line(l_g2_uc_h5_tb4) = {p_g2_uc_h5_t4, p_g2_uc_h5_b4};
l_g2_uc_h5_tb5 = newl; Line(l_g2_uc_h5_tb5) = {p_g2_uc_h5_t5, p_g2_uc_h5_b5};

// lines on lower plane
// on the top plane
l_g2_lc_h1_t1 = newl; Circle(l_g2_lc_h1_t1) = {p_g2_lc_h1_t2, p_g2_lc_h1_t1, p_g2_lc_h1_t3}; // hole 1
l_g2_lc_h12_t1 = newl; Line(l_g2_lc_h12_t1) = {p_g2_lc_h1_t3, p_g2_lc_h2_t2}; // between hole 1 and hole 2
l_g2_lc_h2_t1 = newl; Circle(l_g2_lc_h2_t1) = {p_g2_lc_h2_t2, p_g2_lc_h2_t1, p_g2_lc_h2_t3}; // hole 2
l_g2_lc_h23_t1 = newl; Line(l_g2_lc_h23_t1) = {p_g2_lc_h2_t3, p_g2_lc_h3_t2}; // between hole 2 and hole 3
l_g2_lc_h3_t1 = newl; Circle(l_g2_lc_h3_t1) = {p_g2_lc_h3_t2, p_g2_lc_h3_t1, p_g2_lc_h3_t3}; // hole 3
l_g2_lc_h34_t1 = newl; Line(l_g2_lc_h34_t1) = {p_g2_lc_h3_t3, p_g2_lc_h4_t2}; // between hole 3 and hole 4
l_g2_lc_h4_t1 = newl; Circle(l_g2_lc_h4_t1) = {p_g2_lc_h4_t2, p_g2_lc_h4_t1, p_g2_lc_h4_t3}; // hole 4
l_g2_lc_h41_t1 = newl; Line(l_g2_lc_h41_t1) = {p_g2_lc_h4_t3, p_g2_lc_h1_t2}; // between hole 4 and hole 1
l_g2_lc_h5_t1 = newl; Circle(l_g2_lc_h5_t1) = {p_g2_lc_h5_t2, p_g2_lc_h5_t1, p_g2_lc_h5_t3}; // hole 5
l_g2_lc_h5_t2 = newl; Circle(l_g2_lc_h5_t2) = {p_g2_lc_h5_t3, p_g2_lc_h5_t1, p_g2_lc_h5_t4};
l_g2_lc_h5_t3 = newl; Circle(l_g2_lc_h5_t3) = {p_g2_lc_h5_t4, p_g2_lc_h5_t1, p_g2_lc_h5_t5};
l_g2_lc_h5_t4 = newl; Circle(l_g2_lc_h5_t4) = {p_g2_lc_h5_t5, p_g2_lc_h5_t1, p_g2_lc_h5_t2};
// on the bottom plane
l_g2_lc_h1_b1 = newl; Circle(l_g2_lc_h1_b1) = {p_g2_lc_h1_b2, p_g2_lc_h1_b1, p_g2_lc_h1_b3}; // hole 1
l_g2_lc_h12_b1 = newl; Line(l_g2_lc_h12_b1) = {p_g2_lc_h1_b3, p_g2_lc_h2_b2}; // between hole 1 and hole 2
l_g2_lc_h2_b1 = newl; Circle(l_g2_lc_h2_b1) = {p_g2_lc_h2_b2, p_g2_lc_h2_b1, p_g2_lc_h2_b3}; // hole 2
l_g2_lc_h23_b1 = newl; Line(l_g2_lc_h23_b1) = {p_g2_lc_h2_b3, p_g2_lc_h3_b2}; // between hole 2 and hole 3
l_g2_lc_h3_b1 = newl; Circle(l_g2_lc_h3_b1) = {p_g2_lc_h3_b2, p_g2_lc_h3_b1, p_g2_lc_h3_b3}; // hole 3
l_g2_lc_h34_b1 = newl; Line(l_g2_lc_h34_b1) = {p_g2_lc_h3_b3, p_g2_lc_h4_b2}; // between hole 3 and hole 4
l_g2_lc_h4_b1 = newl; Circle(l_g2_lc_h4_b1) = {p_g2_lc_h4_b2, p_g2_lc_h4_b1, p_g2_lc_h4_b3}; // hole 4
l_g2_lc_h41_b1 = newl; Line(l_g2_lc_h41_b1) = {p_g2_lc_h4_b3, p_g2_lc_h1_b2}; // between hole 4 and hole 1
l_g2_lc_h5_b1 = newl; Circle(l_g2_lc_h5_b1) = {p_g2_lc_h5_b2, p_g2_lc_h5_b1, p_g2_lc_h5_b3}; // hole 5
l_g2_lc_h5_b2 = newl; Circle(l_g2_lc_h5_b2) = {p_g2_lc_h5_b3, p_g2_lc_h5_b1, p_g2_lc_h5_b4};
l_g2_lc_h5_b3 = newl; Circle(l_g2_lc_h5_b3) = {p_g2_lc_h5_b4, p_g2_lc_h5_b1, p_g2_lc_h5_b5};
l_g2_lc_h5_b4 = newl; Circle(l_g2_lc_h5_b4) = {p_g2_lc_h5_b5, p_g2_lc_h5_b1, p_g2_lc_h5_b2};
// connecting top and bottom
l_g2_lc_h1_tb1 = newl; Line(l_g2_lc_h1_tb1) = {p_g2_lc_h1_t1, p_g2_lc_h1_b1}; // hole 1
l_g2_lc_h1_tb2 = newl; Line(l_g2_lc_h1_tb2) = {p_g2_lc_h1_t2, p_g2_lc_h1_b2}; 
l_g2_lc_h1_tb3 = newl; Line(l_g2_lc_h1_tb3) = {p_g2_lc_h1_t3, p_g2_lc_h1_b3};
l_g2_lc_h2_tb1 = newl; Line(l_g2_lc_h2_tb1) = {p_g2_lc_h2_t1, p_g2_lc_h2_b1}; // hole 2
l_g2_lc_h2_tb2 = newl; Line(l_g2_lc_h2_tb2) = {p_g2_lc_h2_t2, p_g2_lc_h2_b2};
l_g2_lc_h2_tb3 = newl; Line(l_g2_lc_h2_tb3) = {p_g2_lc_h2_t3, p_g2_lc_h2_b3};
l_g2_lc_h3_tb1 = newl; Line(l_g2_lc_h3_tb1) = {p_g2_lc_h3_t1, p_g2_lc_h3_b1}; // hole 3
l_g2_lc_h3_tb2 = newl; Line(l_g2_lc_h3_tb2) = {p_g2_lc_h3_t2, p_g2_lc_h3_b2};
l_g2_lc_h3_tb3 = newl; Line(l_g2_lc_h3_tb3) = {p_g2_lc_h3_t3, p_g2_lc_h3_b3};
l_g2_lc_h4_tb1 = newl; Line(l_g2_lc_h4_tb1) = {p_g2_lc_h4_t1, p_g2_lc_h4_b1}; // hole 4
l_g2_lc_h4_tb2 = newl; Line(l_g2_lc_h4_tb2) = {p_g2_lc_h4_t2, p_g2_lc_h4_b2};
l_g2_lc_h4_tb3 = newl; Line(l_g2_lc_h4_tb3) = {p_g2_lc_h4_t3, p_g2_lc_h4_b3};
l_g2_lc_h5_tb1 = newl; Line(l_g2_lc_h5_tb1) = {p_g2_lc_h5_t1, p_g2_lc_h5_b1}; // hole 5
l_g2_lc_h5_tb2 = newl; Line(l_g2_lc_h5_tb2) = {p_g2_lc_h5_t2, p_g2_lc_h5_b2};
l_g2_lc_h5_tb3 = newl; Line(l_g2_lc_h5_tb3) = {p_g2_lc_h5_t3, p_g2_lc_h5_b3};
l_g2_lc_h5_tb4 = newl; Line(l_g2_lc_h5_tb4) = {p_g2_lc_h5_t4, p_g2_lc_h5_b4};
l_g2_lc_h5_tb5 = newl; Line(l_g2_lc_h5_tb5) = {p_g2_lc_h5_t5, p_g2_lc_h5_b5};

// lines on the dielectric layer
// on the top plane
l_g2_diel_h1_t1 = newl; Circle(l_g2_diel_h1_t1) = {p_g2_diel_h1_t2, p_g2_diel_h1_t1, p_g2_diel_h1_t3}; // hole 1
l_g2_diel_h1_t2 = newl; Line(l_g2_diel_h1_t2) = {p_g2_diel_h1_t2, p_g2_uc_h1_b2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h1_t3 = newl; Line(l_g2_diel_h1_t3) = {p_g2_diel_h1_t3, p_g2_uc_h1_b3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h2_t1 = newl; Circle(l_g2_diel_h2_t1) = {p_g2_diel_h2_t2, p_g2_diel_h2_t1, p_g2_diel_h2_t3}; // hole 2
l_g2_diel_h2_t2 = newl; Line(l_g2_diel_h2_t2) = {p_g2_diel_h2_t2, p_g2_uc_h2_b2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h2_t3 = newl; Line(l_g2_diel_h2_t3) = {p_g2_diel_h2_t3, p_g2_uc_h2_b3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h3_t1 = newl; Circle(l_g2_diel_h3_t1) = {p_g2_diel_h3_t2, p_g2_diel_h3_t1, p_g2_diel_h3_t3}; // hole 3
l_g2_diel_h3_t2 = newl; Line(l_g2_diel_h3_t2) = {p_g2_diel_h3_t2, p_g2_uc_h3_b2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h3_t3 = newl; Line(l_g2_diel_h3_t3) = {p_g2_diel_h3_t3, p_g2_uc_h3_b3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h4_t1 = newl; Circle(l_g2_diel_h4_t1) = {p_g2_diel_h4_t2, p_g2_diel_h4_t1, p_g2_diel_h4_t3}; // hole 4
l_g2_diel_h4_t2 = newl; Line(l_g2_diel_h4_t2) = {p_g2_diel_h4_t2, p_g2_uc_h4_b2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h4_t3 = newl; Line(l_g2_diel_h4_t3) = {p_g2_diel_h4_t3, p_g2_uc_h4_b3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_t1 = newl; Circle(l_g2_diel_h5_t1) = {p_g2_diel_h5_t2, p_g2_diel_h5_t1, p_g2_diel_h5_t3}; // hole 5
l_g2_diel_h5_t2 = newl; Circle(l_g2_diel_h5_t2) = {p_g2_diel_h5_t3, p_g2_diel_h5_t1, p_g2_diel_h5_t4};
l_g2_diel_h5_t3 = newl; Circle(l_g2_diel_h5_t3) = {p_g2_diel_h5_t4, p_g2_diel_h5_t1, p_g2_diel_h5_t5};
l_g2_diel_h5_t4 = newl; Circle(l_g2_diel_h5_t4) = {p_g2_diel_h5_t5, p_g2_diel_h5_t1, p_g2_diel_h5_t2};
l_g2_diel_h5_t5 = newl; Line(l_g2_diel_h5_t5) = {p_g2_diel_h5_t2, p_g2_uc_h5_b2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_t6 = newl; Line(l_g2_diel_h5_t6) = {p_g2_diel_h5_t3, p_g2_uc_h5_b3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_t7 = newl; Line(l_g2_diel_h5_t7) = {p_g2_diel_h5_t4, p_g2_uc_h5_b4}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_t8 = newl; Line(l_g2_diel_h5_t8) = {p_g2_diel_h5_t5, p_g2_uc_h5_b5}; // this line is on the surface exposed to gas by etching
// connecting the top and middel plane
l_g2_diel_h1_tm1 = newl; Line(l_g2_diel_h1_tm1) = {p_g2_diel_h1_t1, p_g2_diel_h1_m1};
l_g2_diel_h1_tm2 = newl; Line(l_g2_diel_h1_tm2) = {p_g2_diel_h1_t2, p_g2_diel_h1_m2};
l_g2_diel_h1_tm3 = newl; Line(l_g2_diel_h1_tm3) = {p_g2_diel_h1_t3, p_g2_diel_h1_m3};
l_g2_diel_h2_tm1 = newl; Line(l_g2_diel_h2_tm1) = {p_g2_diel_h2_t1, p_g2_diel_h2_m1};
l_g2_diel_h2_tm2 = newl; Line(l_g2_diel_h2_tm2) = {p_g2_diel_h2_t2, p_g2_diel_h2_m2};
l_g2_diel_h2_tm3 = newl; Line(l_g2_diel_h2_tm3) = {p_g2_diel_h2_t3, p_g2_diel_h2_m3};
l_g2_diel_h3_tm1 = newl; Line(l_g2_diel_h3_tm1) = {p_g2_diel_h3_t1, p_g2_diel_h3_m1};
l_g2_diel_h3_tm2 = newl; Line(l_g2_diel_h3_tm2) = {p_g2_diel_h3_t2, p_g2_diel_h3_m2};
l_g2_diel_h3_tm3 = newl; Line(l_g2_diel_h3_tm3) = {p_g2_diel_h3_t3, p_g2_diel_h3_m3};
l_g2_diel_h4_tm1 = newl; Line(l_g2_diel_h4_tm1) = {p_g2_diel_h4_t1, p_g2_diel_h4_m1};
l_g2_diel_h4_tm2 = newl; Line(l_g2_diel_h4_tm2) = {p_g2_diel_h4_t2, p_g2_diel_h4_m2};
l_g2_diel_h4_tm3 = newl; Line(l_g2_diel_h4_tm3) = {p_g2_diel_h4_t3, p_g2_diel_h4_m3};
l_g2_diel_h5_tm1 = newl; Line(l_g2_diel_h5_tm1) = {p_g2_diel_h5_t1, p_g2_diel_h5_m1};
l_g2_diel_h5_tm2 = newl; Line(l_g2_diel_h5_tm2) = {p_g2_diel_h5_t2, p_g2_diel_h5_m2};
l_g2_diel_h5_tm3 = newl; Line(l_g2_diel_h5_tm3) = {p_g2_diel_h5_t3, p_g2_diel_h5_m3};
l_g2_diel_h5_tm4 = newl; Line(l_g2_diel_h5_tm4) = {p_g2_diel_h5_t4, p_g2_diel_h5_m4};
l_g2_diel_h5_tm5 = newl; Line(l_g2_diel_h5_tm5) = {p_g2_diel_h5_t5, p_g2_diel_h5_m5};
// on the middle plane
l_g2_diel_h1_m1 = newl; Circle(l_g2_diel_h1_m1) = {p_g2_diel_h1_m2, p_g2_diel_h1_m1, p_g2_diel_h1_m3}; // hole 1
l_g2_diel_h12_m1 = newl; Line(l_g2_diel_h12_m1) = {p_g2_diel_h1_m3, p_g2_diel_h2_m2}; // between hole 1 and hole 2
l_g2_diel_h2_m1 = newl; Circle(l_g2_diel_h2_m1) = {p_g2_diel_h2_m2, p_g2_diel_h2_m1, p_g2_diel_h2_m3}; // hole 2
l_g2_diel_h23_m1 = newl; Line(l_g2_diel_h23_m1) = {p_g2_diel_h2_m3, p_g2_diel_h3_m2}; // between hole 2 and hole 3
l_g2_diel_h3_m1 = newl; Circle(l_g2_diel_h3_m1) = {p_g2_diel_h3_m2, p_g2_diel_h3_m1, p_g2_diel_h3_m3}; // hole 3
l_g2_diel_h34_m1 = newl; Line(l_g2_diel_h34_m1) = {p_g2_diel_h3_m3, p_g2_diel_h4_m2}; // between hole 3 and hole 4
l_g2_diel_h4_m1 = newl; Circle(l_g2_diel_h4_m1) = {p_g2_diel_h4_m2, p_g2_diel_h4_m1, p_g2_diel_h4_m3}; // hole 4
l_g2_diel_h41_m1 = newl; Line(l_g2_diel_h41_m1) = {p_g2_diel_h4_m3, p_g2_diel_h1_m2}; // between hole 4 and hole 1
l_g2_diel_h5_m1 = newl; Circle(l_g2_diel_h5_m1) = {p_g2_diel_h5_m2, p_g2_diel_h5_m1, p_g2_diel_h5_m3}; // hole 5
l_g2_diel_h5_m2 = newl; Circle(l_g2_diel_h5_m2) = {p_g2_diel_h5_m3, p_g2_diel_h5_m1, p_g2_diel_h5_m4};
l_g2_diel_h5_m3 = newl; Circle(l_g2_diel_h5_m3) = {p_g2_diel_h5_m4, p_g2_diel_h5_m1, p_g2_diel_h5_m5};
l_g2_diel_h5_m4 = newl; Circle(l_g2_diel_h5_m4) = {p_g2_diel_h5_m5, p_g2_diel_h5_m1, p_g2_diel_h5_m2};
// connecting the middle and bottom plane
l_g2_diel_h1_mb1 = newl; Line(l_g2_diel_h1_mb1) = {p_g2_diel_h1_m1, p_g2_diel_h1_b1};
l_g2_diel_h1_mb2 = newl; Line(l_g2_diel_h1_mb2) = {p_g2_diel_h1_m2, p_g2_diel_h1_b2};
l_g2_diel_h1_mb3 = newl; Line(l_g2_diel_h1_mb3) = {p_g2_diel_h1_m3, p_g2_diel_h1_b3};
l_g2_diel_h2_mb1 = newl; Line(l_g2_diel_h2_mb1) = {p_g2_diel_h2_m1, p_g2_diel_h2_b1};
l_g2_diel_h2_mb2 = newl; Line(l_g2_diel_h2_mb2) = {p_g2_diel_h2_m2, p_g2_diel_h2_b2};
l_g2_diel_h2_mb3 = newl; Line(l_g2_diel_h2_mb3) = {p_g2_diel_h2_m3, p_g2_diel_h2_b3};
l_g2_diel_h3_mb1 = newl; Line(l_g2_diel_h3_mb1) = {p_g2_diel_h3_m1, p_g2_diel_h3_b1};
l_g2_diel_h3_mb2 = newl; Line(l_g2_diel_h3_mb2) = {p_g2_diel_h3_m2, p_g2_diel_h3_b2};
l_g2_diel_h3_mb3 = newl; Line(l_g2_diel_h3_mb3) = {p_g2_diel_h3_m3, p_g2_diel_h3_b3};
l_g2_diel_h4_mb1 = newl; Line(l_g2_diel_h4_mb1) = {p_g2_diel_h4_m1, p_g2_diel_h4_b1};
l_g2_diel_h4_mb2 = newl; Line(l_g2_diel_h4_mb2) = {p_g2_diel_h4_m2, p_g2_diel_h4_b2};
l_g2_diel_h4_mb3 = newl; Line(l_g2_diel_h4_mb3) = {p_g2_diel_h4_m3, p_g2_diel_h4_b3};
l_g2_diel_h5_mb1 = newl; Line(l_g2_diel_h5_mb1) = {p_g2_diel_h5_m1, p_g2_diel_h5_b1};
l_g2_diel_h5_mb2 = newl; Line(l_g2_diel_h5_mb2) = {p_g2_diel_h5_m2, p_g2_diel_h5_b2};
l_g2_diel_h5_mb3 = newl; Line(l_g2_diel_h5_mb3) = {p_g2_diel_h5_m3, p_g2_diel_h5_b3};
l_g2_diel_h5_mb4 = newl; Line(l_g2_diel_h5_mb4) = {p_g2_diel_h5_m4, p_g2_diel_h5_b4};
l_g2_diel_h5_mb5 = newl; Line(l_g2_diel_h5_mb5) = {p_g2_diel_h5_m5, p_g2_diel_h5_b5};
// on the bottom plane
l_g2_diel_h1_b1 = newl; Circle(l_g2_diel_h1_b1) = {p_g2_diel_h1_b2, p_g2_diel_h1_b1, p_g2_diel_h1_b3}; // hole 1
l_g2_diel_h1_b2 = newl; Line(l_g2_diel_h1_b2) = {p_g2_diel_h1_b2, p_g2_lc_h1_t2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h1_b3 = newl; Line(l_g2_diel_h1_b3) = {p_g2_diel_h1_b3, p_g2_lc_h1_t3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h2_b1 = newl; Circle(l_g2_diel_h2_b1) = {p_g2_diel_h2_b2, p_g2_diel_h2_b1, p_g2_diel_h2_b3}; // hole 2
l_g2_diel_h2_b2 = newl; Line(l_g2_diel_h2_b2) = {p_g2_diel_h2_b2, p_g2_lc_h2_t2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h2_b3 = newl; Line(l_g2_diel_h2_b3) = {p_g2_diel_h2_b3, p_g2_lc_h2_t3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h3_b1 = newl; Circle(l_g2_diel_h3_b1) = {p_g2_diel_h3_b2, p_g2_diel_h3_b1, p_g2_diel_h3_b3}; // hole 3
l_g2_diel_h3_b2 = newl; Line(l_g2_diel_h3_b2) = {p_g2_diel_h3_b2, p_g2_lc_h3_t2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h3_b3 = newl; Line(l_g2_diel_h3_b3) = {p_g2_diel_h3_b3, p_g2_lc_h3_t3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h4_b1 = newl; Circle(l_g2_diel_h4_b1) = {p_g2_diel_h4_b2, p_g2_diel_h4_b1, p_g2_diel_h4_b3}; // hole 4
l_g2_diel_h4_b2 = newl; Line(l_g2_diel_h4_b2) = {p_g2_diel_h4_b2, p_g2_lc_h4_t2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h4_b3 = newl; Line(l_g2_diel_h4_b3) = {p_g2_diel_h4_b3, p_g2_lc_h4_t3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_b1 = newl; Circle(l_g2_diel_h5_b1) = {p_g2_diel_h5_b2, p_g2_diel_h5_b1, p_g2_diel_h5_b3}; // hole 5
l_g2_diel_h5_b2 = newl; Circle(l_g2_diel_h5_b2) = {p_g2_diel_h5_b3, p_g2_diel_h5_b1, p_g2_diel_h5_b4};
l_g2_diel_h5_b3 = newl; Circle(l_g2_diel_h5_b3) = {p_g2_diel_h5_b4, p_g2_diel_h5_b1, p_g2_diel_h5_b5};
l_g2_diel_h5_b4 = newl; Circle(l_g2_diel_h5_b4) = {p_g2_diel_h5_b5, p_g2_diel_h5_b1, p_g2_diel_h5_b2};
l_g2_diel_h5_b5 = newl; Line(l_g2_diel_h5_b5) = {p_g2_diel_h5_b2, p_g2_lc_h5_t2}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_b6 = newl; Line(l_g2_diel_h5_b6) = {p_g2_diel_h5_b3, p_g2_lc_h5_t3}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_b7 = newl; Line(l_g2_diel_h5_b7) = {p_g2_diel_h5_b4, p_g2_lc_h5_t4}; // this line is on the surface exposed to gas by etching
l_g2_diel_h5_b8 = newl; Line(l_g2_diel_h5_b8) = {p_g2_diel_h5_b5, p_g2_lc_h5_t5}; // this line is on the surface exposed to gas by etching

/*---------------------------------- Curve Loop and Surface Definitions -------------------------------------*/
// top plane of upper copper plate
cl_g2_uc_t1 = newll; Curve Loop(cl_g2_uc_t1) = {l_g2_uc_h1_t1, l_g2_uc_h12_t1, l_g2_uc_h2_t1, l_g2_uc_h23_t1,
    l_g2_uc_h3_t1, l_g2_uc_h34_t1, l_g2_uc_h4_t1, l_g2_uc_h41_t1};
cl_g2_uc_t2 = newll; Curve Loop(cl_g2_uc_t2) = {l_g2_uc_h5_t1, l_g2_uc_h5_t2, l_g2_uc_h5_t3, l_g2_uc_h5_t4}; // curve loop for hole 5
// planes between top and bottom plane
cl_g2_uc_h1_tb1 = newll; Curve Loop(cl_g2_uc_h1_tb1) = {-l_g2_uc_h1_tb2, l_g2_uc_h1_t1, l_g2_uc_h1_tb3, -l_g2_uc_h1_b1}; // cylinderical surface on hole 1
cl_g2_uc_h12_tb1 = newll; Curve Loop(cl_g2_uc_h12_tb1) = {-l_g2_uc_h1_tb3, l_g2_uc_h12_t1, l_g2_uc_h2_tb2, -l_g2_uc_h12_b1}; // plane surface bwt hole 1 and hole  2
cl_g2_uc_h2_tb1 = newll; Curve Loop(cl_g2_uc_h2_tb1) = {-l_g2_uc_h2_tb2, l_g2_uc_h2_t1, l_g2_uc_h2_tb3, -l_g2_uc_h2_b1}; // cylinderical surface on hole 2
cl_g2_uc_h23_tb1 = newll; Curve Loop(cl_g2_uc_h23_tb1) = {-l_g2_uc_h2_tb3, l_g2_uc_h23_t1, l_g2_uc_h3_tb2, -l_g2_uc_h23_b1}; // plane surface bwt hole 2 and hole  3
cl_g2_uc_h3_tb1 = newll; Curve Loop(cl_g2_uc_h3_tb1) = {-l_g2_uc_h3_tb2, l_g2_uc_h3_t1, l_g2_uc_h3_tb3, -l_g2_uc_h3_b1}; // cylinderical surface on hole 3
cl_g2_uc_h34_tb1 = newll; Curve Loop(cl_g2_uc_h34_tb1) = {-l_g2_uc_h3_tb3, l_g2_uc_h34_t1, l_g2_uc_h4_tb2, -l_g2_uc_h34_b1}; // plane surface bwt hole 3 and hole  4
cl_g2_uc_h4_tb1 = newll; Curve Loop(cl_g2_uc_h4_tb1) = {-l_g2_uc_h4_tb2, l_g2_uc_h4_t1, l_g2_uc_h4_tb3, -l_g2_uc_h4_b1}; // cylinderical surface on hole 4
cl_g2_uc_h41_tb1 = newll; Curve Loop(cl_g2_uc_h41_tb1) = {-l_g2_uc_h4_tb3, l_g2_uc_h41_t1, l_g2_uc_h1_tb2, -l_g2_uc_h41_b1}; // plane surface bwt hole 4 and hole  1
cl_g2_uc_h5_tb1 = newll; Curve Loop(cl_g2_uc_h5_tb1) = {-l_g2_uc_h5_tb2, l_g2_uc_h5_t1, l_g2_uc_h5_tb3, -l_g2_uc_h5_b1}; // cylinderical surface on hole 5
cl_g2_uc_h5_tb2 = newll; Curve Loop(cl_g2_uc_h5_tb2) = {-l_g2_uc_h5_tb3, l_g2_uc_h5_t2, l_g2_uc_h5_tb4, -l_g2_uc_h5_b2};
cl_g2_uc_h5_tb3 = newll; Curve Loop(cl_g2_uc_h5_tb3) = {-l_g2_uc_h5_tb4, l_g2_uc_h5_t3, l_g2_uc_h5_tb5, -l_g2_uc_h5_b3};
cl_g2_uc_h5_tb4 = newll; Curve Loop(cl_g2_uc_h5_tb4) = {-l_g2_uc_h5_tb5, l_g2_uc_h5_t4, l_g2_uc_h5_tb2, -l_g2_uc_h5_b4};
// bottom plane of upper copper plate
cl_g2_uc_b1 = newll; Curve Loop(cl_g2_uc_b1) = {l_g2_uc_h1_b1, l_g2_uc_h12_b1, l_g2_uc_h2_b1, l_g2_uc_h23_b1,
    l_g2_uc_h3_b1, l_g2_uc_h34_b1, l_g2_uc_h4_b1, l_g2_uc_h41_b1};
cl_g2_uc_b2 = newll; Curve Loop(cl_g2_uc_b2) = {l_g2_uc_h5_b1, l_g2_uc_h5_b2, l_g2_uc_h5_b3, l_g2_uc_h5_b4}; // curve loop for hole 5
// surface definitions of upper plate
s_g2_uc_t1 = news; Plane Surface(s_g2_uc_t1) = {cl_g2_uc_t1, cl_g2_uc_t2};
s_g2_uc_h1_tb1 = news; Surface(s_g2_uc_h1_tb1) = {cl_g2_uc_h1_tb1}; // cylinderical surface on hole 1
s_g2_uc_h12_tb1 = news; Plane Surface(s_g2_uc_h12_tb1) = {cl_g2_uc_h12_tb1}; // plane surface bwt hole 1 and hole  2
s_g2_uc_h2_tb1 = news; Surface(s_g2_uc_h2_tb1) = {cl_g2_uc_h2_tb1}; // cylinderical surface on hole 2
s_g2_uc_h23_tb1 = news; Plane Surface(s_g2_uc_h23_tb1) = {cl_g2_uc_h23_tb1}; // plane surface bwt hole 2 and hole  3
s_g2_uc_h3_tb1 = news; Surface(s_g2_uc_h3_tb1) = {cl_g2_uc_h3_tb1}; // cylinderical surface on hole 3
s_g2_uc_h34_tb1 = news; Plane Surface(s_g2_uc_h34_tb1) = {cl_g2_uc_h34_tb1}; // plane surface bwt hole 3 and hole  4
s_g2_uc_h4_tb1 = news; Surface(s_g2_uc_h4_tb1) = {cl_g2_uc_h4_tb1}; // cylinderical surface on hole 4
s_g2_uc_h41_tb1 = news; Plane Surface(s_g2_uc_h41_tb1) = {cl_g2_uc_h41_tb1}; // plane surface bwt hole 4 and hole  1
s_g2_uc_h5_tb1 = news; Surface(s_g2_uc_h5_tb1) = {cl_g2_uc_h5_tb1}; // cylinderical surface on hole 5
s_g2_uc_h5_tb2 = news; Surface(s_g2_uc_h5_tb2) = {cl_g2_uc_h5_tb2};
s_g2_uc_h5_tb3 = news; Surface(s_g2_uc_h5_tb3) = {cl_g2_uc_h5_tb3};
s_g2_uc_h5_tb4 = news; Surface(s_g2_uc_h5_tb4) = {cl_g2_uc_h5_tb4};
s_g2_uc_b1 = news; Plane Surface(s_g2_uc_b1) = {cl_g2_uc_b1, cl_g2_uc_b2};

// curve loop definitions of upper copper plate
// top plane of upper copper plate
cl_g2_lc_t1 = newll; Curve Loop(cl_g2_lc_t1) = {l_g2_lc_h1_t1, l_g2_lc_h12_t1, l_g2_lc_h2_t1, l_g2_lc_h23_t1,
    l_g2_lc_h3_t1, l_g2_lc_h34_t1, l_g2_lc_h4_t1, l_g2_lc_h41_t1};
cl_g2_lc_t2 = newll; Curve Loop(cl_g2_lc_t2) = {l_g2_lc_h5_t1, l_g2_lc_h5_t2, l_g2_lc_h5_t3, l_g2_lc_h5_t4}; // curve loop for hole 5
// planes between top and bottom plane
cl_g2_lc_h1_tb1 = newll; Curve Loop(cl_g2_lc_h1_tb1) = {-l_g2_lc_h1_tb2, l_g2_lc_h1_t1, l_g2_lc_h1_tb3, -l_g2_lc_h1_b1}; // cylinderical surface on hole 1
cl_g2_lc_h12_tb1 = newll; Curve Loop(cl_g2_lc_h12_tb1) = {-l_g2_lc_h1_tb3, l_g2_lc_h12_t1, l_g2_lc_h2_tb2, -l_g2_lc_h12_b1}; // plane surface bwt hole 1 and hole  2
cl_g2_lc_h2_tb1 = newll; Curve Loop(cl_g2_lc_h2_tb1) = {-l_g2_lc_h2_tb2, l_g2_lc_h2_t1, l_g2_lc_h2_tb3, -l_g2_lc_h2_b1}; // cylinderical surface on hole 2
cl_g2_lc_h23_tb1 = newll; Curve Loop(cl_g2_lc_h23_tb1) = {-l_g2_lc_h2_tb3, l_g2_lc_h23_t1, l_g2_lc_h3_tb2, -l_g2_lc_h23_b1}; // plane surface bwt hole 2 and hole  3
cl_g2_lc_h3_tb1 = newll; Curve Loop(cl_g2_lc_h3_tb1) = {-l_g2_lc_h3_tb2, l_g2_lc_h3_t1, l_g2_lc_h3_tb3, -l_g2_lc_h3_b1}; // cylinderical surface on hole 3
cl_g2_lc_h34_tb1 = newll; Curve Loop(cl_g2_lc_h34_tb1) = {-l_g2_lc_h3_tb3, l_g2_lc_h34_t1, l_g2_lc_h4_tb2, -l_g2_lc_h34_b1}; // plane surface bwt hole 3 and hole  4
cl_g2_lc_h4_tb1 = newll; Curve Loop(cl_g2_lc_h4_tb1) = {-l_g2_lc_h4_tb2, l_g2_lc_h4_t1, l_g2_lc_h4_tb3, -l_g2_lc_h4_b1}; // cylinderical surface on hole 4
cl_g2_lc_h41_tb1 = newll; Curve Loop(cl_g2_lc_h41_tb1) = {-l_g2_lc_h4_tb3, l_g2_lc_h41_t1, l_g2_lc_h1_tb2, -l_g2_lc_h41_b1}; // plane surface bwt hole 4 and hole  1
cl_g2_lc_h5_tb1 = newll; Curve Loop(cl_g2_lc_h5_tb1) = {-l_g2_lc_h5_tb2, l_g2_lc_h5_t1, l_g2_lc_h5_tb3, -l_g2_lc_h5_b1}; // cylinderical surface on hole 5
cl_g2_lc_h5_tb2 = newll; Curve Loop(cl_g2_lc_h5_tb2) = {-l_g2_lc_h5_tb3, l_g2_lc_h5_t2, l_g2_lc_h5_tb4, -l_g2_lc_h5_b2};
cl_g2_lc_h5_tb3 = newll; Curve Loop(cl_g2_lc_h5_tb3) = {-l_g2_lc_h5_tb4, l_g2_lc_h5_t3, l_g2_lc_h5_tb5, -l_g2_lc_h5_b3};
cl_g2_lc_h5_tb4 = newll; Curve Loop(cl_g2_lc_h5_tb4) = {-l_g2_lc_h5_tb5, l_g2_lc_h5_t4, l_g2_lc_h5_tb2, -l_g2_lc_h5_b4};
// bottom plane of upper copper plate
cl_g2_lc_b1 = newll; Curve Loop(cl_g2_lc_b1) = {l_g2_lc_h1_b1, l_g2_lc_h12_b1, l_g2_lc_h2_b1, l_g2_lc_h23_b1,
    l_g2_lc_h3_b1, l_g2_lc_h34_b1, l_g2_lc_h4_b1, l_g2_lc_h41_b1};
cl_g2_lc_b2 = newll; Curve Loop(cl_g2_lc_b2) = {l_g2_lc_h5_b1, l_g2_lc_h5_b2, l_g2_lc_h5_b3, l_g2_lc_h5_b4}; // curve loop for hole 5
// surface definitions of upper plate
s_g2_lc_t1 = news; Plane Surface(s_g2_lc_t1) = {cl_g2_lc_t1, cl_g2_lc_t2};
s_g2_lc_h1_tb1 = news; Surface(s_g2_lc_h1_tb1) = {cl_g2_lc_h1_tb1}; // cylinderical surface on hole 1
s_g2_lc_h12_tb1 = news; Plane Surface(s_g2_lc_h12_tb1) = {cl_g2_lc_h12_tb1}; // plane surface bwt hole 1 and hole  2
s_g2_lc_h2_tb1 = news; Surface(s_g2_lc_h2_tb1) = {cl_g2_lc_h2_tb1}; // cylinderical surface on hole 2
s_g2_lc_h23_tb1 = news; Plane Surface(s_g2_lc_h23_tb1) = {cl_g2_lc_h23_tb1}; // plane surface bwt hole 2 and hole  3
s_g2_lc_h3_tb1 = news; Surface(s_g2_lc_h3_tb1) = {cl_g2_lc_h3_tb1}; // cylinderical surface on hole 3
s_g2_lc_h34_tb1 = news; Plane Surface(s_g2_lc_h34_tb1) = {cl_g2_lc_h34_tb1}; // plane surface bwt hole 3 and hole  4
s_g2_lc_h4_tb1 = news; Surface(s_g2_lc_h4_tb1) = {cl_g2_lc_h4_tb1}; // cylinderical surface on hole 4
s_g2_lc_h41_tb1 = news; Plane Surface(s_g2_lc_h41_tb1) = {cl_g2_lc_h41_tb1}; // plane surface bwt hole 4 and hole  1
s_g2_lc_h5_tb1 = news; Surface(s_g2_lc_h5_tb1) = {cl_g2_lc_h5_tb1}; // cylinderical surface on hole 5
s_g2_lc_h5_tb2 = news; Surface(s_g2_lc_h5_tb2) = {cl_g2_lc_h5_tb2};
s_g2_lc_h5_tb3 = news; Surface(s_g2_lc_h5_tb3) = {cl_g2_lc_h5_tb3};
s_g2_lc_h5_tb4 = news; Surface(s_g2_lc_h5_tb4) = {cl_g2_lc_h5_tb4};
s_g2_lc_b1 = news; Plane Surface(s_g2_lc_b1) = {cl_g2_lc_b1, cl_g2_lc_b2};

// curve loop definitions of dielectric layer
// top plane surfaces
// a surface exposed to gas by etching
cl_g2_diel_h1_t1 = newll; Curve Loop(cl_g2_diel_h1_t1) = {-l_g2_diel_h1_t2, l_g2_diel_h1_t1, l_g2_diel_h1_t3, -l_g2_uc_h1_b1}; 
cl_g2_diel_h2_t1 = newll; Curve Loop(cl_g2_diel_h2_t1) = {-l_g2_diel_h2_t2, l_g2_diel_h2_t1, l_g2_diel_h2_t3, -l_g2_uc_h2_b1};
cl_g2_diel_h3_t1 = newll; Curve Loop(cl_g2_diel_h3_t1) = {-l_g2_diel_h3_t2, l_g2_diel_h3_t1, l_g2_diel_h3_t3, -l_g2_uc_h3_b1};
cl_g2_diel_h4_t1 = newll; Curve Loop(cl_g2_diel_h4_t1) = {-l_g2_diel_h4_t2, l_g2_diel_h4_t1, l_g2_diel_h4_t3, -l_g2_uc_h4_b1};
cl_g2_diel_h5_t1 = newll; Curve Loop(cl_g2_diel_h5_t1) = {-l_g2_diel_h5_t5, l_g2_diel_h5_t1, l_g2_diel_h5_t6, -l_g2_uc_h5_b1};
cl_g2_diel_h5_t2 = newll; Curve Loop(cl_g2_diel_h5_t2) = {-l_g2_diel_h5_t6, l_g2_diel_h5_t2, l_g2_diel_h5_t7, -l_g2_uc_h5_b2};
cl_g2_diel_h5_t3 = newll; Curve Loop(cl_g2_diel_h5_t3) = {-l_g2_diel_h5_t7, l_g2_diel_h5_t3, l_g2_diel_h5_t8, -l_g2_uc_h5_b3};
cl_g2_diel_h5_t4 = newll; Curve Loop(cl_g2_diel_h5_t4) = {-l_g2_diel_h5_t8, l_g2_diel_h5_t4, l_g2_diel_h5_t5, -l_g2_uc_h5_b4};
// planes between top and middle plane
cl_g2_diel_h1_tm1 = newll; Curve Loop(cl_g2_diel_h1_tm1) = {l_g2_diel_h1_tm2, l_g2_diel_h1_m1, -l_g2_diel_h1_tm3, -l_g2_diel_h1_t1};
cl_g2_diel_h1_tm1 = newll; Curve Loop(cl_g2_diel_h1_tm1) = {l_g2_diel_h1_tm2, l_g2_diel_h1_m1, -l_g2_diel_h1_tm3, -l_g2_diel_h1_t1};
cl_g2_diel_h2_tm1 = newll; Curve Loop(cl_g2_diel_h2_tm1) = {l_g2_diel_h2_tm2, l_g2_diel_h2_m1, -l_g2_diel_h2_tm3, -l_g2_diel_h2_t1};
cl_g2_diel_h3_tm1 = newll; Curve Loop(cl_g2_diel_h3_tm1) = {l_g2_diel_h3_tm2, l_g2_diel_h3_m1, -l_g2_diel_h3_tm3, -l_g2_diel_h3_t1};
cl_g2_diel_h4_tm1 = newll; Curve Loop(cl_g2_diel_h4_tm1) = {l_g2_diel_h4_tm2, l_g2_diel_h4_m1, -l_g2_diel_h4_tm3, -l_g2_diel_h4_t1};
cl_g2_diel_h5_tm1 = newll; Curve Loop(cl_g2_diel_h5_tm1) = {l_g2_diel_h5_tm2, l_g2_diel_h5_m1, -l_g2_diel_h5_tm3, -l_g2_diel_h5_t1};
cl_g2_diel_h5_tm2 = newll; Curve Loop(cl_g2_diel_h5_tm2) = {l_g2_diel_h5_tm3, l_g2_diel_h5_m2, -l_g2_diel_h5_tm4, -l_g2_diel_h5_t2};
cl_g2_diel_h5_tm3 = newll; Curve Loop(cl_g2_diel_h5_tm3) = {l_g2_diel_h5_tm4, l_g2_diel_h5_m3, -l_g2_diel_h5_tm5, -l_g2_diel_h5_t3};
cl_g2_diel_h5_tm4 = newll; Curve Loop(cl_g2_diel_h5_tm4) = {l_g2_diel_h5_tm5, l_g2_diel_h5_m4, -l_g2_diel_h5_tm2, -l_g2_diel_h5_t4};
// planes between middle and bottom plane
cl_g2_diel_h1_mb1 = newll; Curve Loop(cl_g2_diel_h1_mb1) = {l_g2_diel_h1_mb2, l_g2_diel_h1_b1, -l_g2_diel_h1_mb3, -l_g2_diel_h1_m1};
cl_g2_diel_h2_mb1 = newll; Curve Loop(cl_g2_diel_h2_mb1) = {l_g2_diel_h2_mb2, l_g2_diel_h2_b1, -l_g2_diel_h2_mb3, -l_g2_diel_h2_m1};
cl_g2_diel_h3_mb1 = newll; Curve Loop(cl_g2_diel_h3_mb1) = {l_g2_diel_h3_mb2, l_g2_diel_h3_b1, -l_g2_diel_h3_mb3, -l_g2_diel_h3_m1};
cl_g2_diel_h4_mb1 = newll; Curve Loop(cl_g2_diel_h4_mb1) = {l_g2_diel_h4_mb2, l_g2_diel_h4_b1, -l_g2_diel_h4_mb3, -l_g2_diel_h4_m1};
cl_g2_diel_h5_mb1 = newll; Curve Loop(cl_g2_diel_h5_mb1) = {l_g2_diel_h5_mb2, l_g2_diel_h5_b1, -l_g2_diel_h5_mb3, -l_g2_diel_h5_m1};
cl_g2_diel_h5_mb2 = newll; Curve Loop(cl_g2_diel_h5_mb2) = {l_g2_diel_h5_mb3, l_g2_diel_h5_b2, -l_g2_diel_h5_mb4, -l_g2_diel_h5_m2};
cl_g2_diel_h5_mb3 = newll; Curve Loop(cl_g2_diel_h5_mb3) = {l_g2_diel_h5_mb4, l_g2_diel_h5_b3, -l_g2_diel_h5_mb5, -l_g2_diel_h5_m3};
cl_g2_diel_h5_mb4 = newll; Curve Loop(cl_g2_diel_h5_mb4) = {l_g2_diel_h5_mb5, l_g2_diel_h5_b4, -l_g2_diel_h5_mb2, -l_g2_diel_h5_m4};
// planes between top and bottom plane
cl_g2_diel_h12_tb1 = newll; Curve Loop(cl_g2_diel_h12_tb1) = {-l_g2_lc_h12_t1, -l_g2_diel_h1_b3, -l_g2_diel_h1_mb3, -l_g2_diel_h1_tm3, l_g2_diel_h1_t3,
    l_g2_uc_h12_b1, -l_g2_diel_h2_t2, l_g2_diel_h2_tm2, l_g2_diel_h2_mb2, l_g2_diel_h2_b2};
cl_g2_diel_h23_tb1 = newll; Curve Loop(cl_g2_diel_h23_tb1) = {-l_g2_lc_h23_t1, -l_g2_diel_h2_b3, -l_g2_diel_h2_mb3, -l_g2_diel_h2_tm3, l_g2_diel_h2_t3,
    l_g2_uc_h23_b1, -l_g2_diel_h3_t2, l_g2_diel_h3_tm2, l_g2_diel_h3_mb2, l_g2_diel_h3_b2};
cl_g2_diel_h34_tb1 = newll; Curve Loop(cl_g2_diel_h34_tb1) = {-l_g2_lc_h34_t1, -l_g2_diel_h3_b3, -l_g2_diel_h3_mb3, -l_g2_diel_h3_tm3, l_g2_diel_h3_t3,
    l_g2_uc_h34_b1, -l_g2_diel_h4_t2, l_g2_diel_h4_tm2, l_g2_diel_h4_mb2, l_g2_diel_h4_b2};
cl_g2_diel_h41_tb1 = newll; Curve Loop(cl_g2_diel_h41_tb1) = {-l_g2_lc_h41_t1, -l_g2_diel_h4_b3, -l_g2_diel_h4_mb3, -l_g2_diel_h4_tm3, l_g2_diel_h4_t3,
    l_g2_uc_h41_b1, -l_g2_diel_h1_t2, l_g2_diel_h1_tm2, l_g2_diel_h1_mb2, l_g2_diel_h1_b2};
// bottom plane surfaces
// a surface exposed to gas by etching
cl_g2_diel_h1_b1 = newll; Curve Loop(cl_g2_diel_h1_b1) = {-l_g2_diel_h1_b2, l_g2_diel_h1_b1, l_g2_diel_h1_b3, -l_g2_lc_h1_t1}; 
cl_g2_diel_h2_b1 = newll; Curve Loop(cl_g2_diel_h2_b1) = {-l_g2_diel_h2_b2, l_g2_diel_h2_b1, l_g2_diel_h2_b3, -l_g2_lc_h2_t1};
cl_g2_diel_h3_b1 = newll; Curve Loop(cl_g2_diel_h3_b1) = {-l_g2_diel_h3_b2, l_g2_diel_h3_b1, l_g2_diel_h3_b3, -l_g2_lc_h3_t1};
cl_g2_diel_h4_b1 = newll; Curve Loop(cl_g2_diel_h4_b1) = {-l_g2_diel_h4_b2, l_g2_diel_h4_b1, l_g2_diel_h4_b3, -l_g2_lc_h4_t1};
cl_g2_diel_h5_b1 = newll; Curve Loop(cl_g2_diel_h5_b1) = {-l_g2_diel_h5_b5, l_g2_diel_h5_b1, l_g2_diel_h5_b6, -l_g2_lc_h5_t1};
cl_g2_diel_h5_b2 = newll; Curve Loop(cl_g2_diel_h5_b2) = {-l_g2_diel_h5_b6, l_g2_diel_h5_b2, l_g2_diel_h5_b7, -l_g2_lc_h5_t2};
cl_g2_diel_h5_b3 = newll; Curve Loop(cl_g2_diel_h5_b3) = {-l_g2_diel_h5_b7, l_g2_diel_h5_b3, l_g2_diel_h5_b8, -l_g2_lc_h5_t3};
cl_g2_diel_h5_b4 = newll; Curve Loop(cl_g2_diel_h5_b4) = {-l_g2_diel_h5_b8, l_g2_diel_h5_b4, l_g2_diel_h5_b5, -l_g2_lc_h5_t4};

// surface definitions of dielectric layer
s_g2_diel_h1_t1 = news; Plane Surface(s_g2_diel_h1_t1) = {cl_g2_diel_h1_t1};
s_g2_diel_h2_t1 = news; Plane Surface(s_g2_diel_h2_t1) = {cl_g2_diel_h2_t1};
s_g2_diel_h3_t1 = news; Plane Surface(s_g2_diel_h3_t1) = {cl_g2_diel_h3_t1};
s_g2_diel_h4_t1 = news; Plane Surface(s_g2_diel_h4_t1) = {cl_g2_diel_h4_t1};
s_g2_diel_h5_t1 = news; Plane Surface(s_g2_diel_h5_t1) = {cl_g2_diel_h5_t1};
s_g2_diel_h5_t2 = news; Plane Surface(s_g2_diel_h5_t2) = {cl_g2_diel_h5_t2};
s_g2_diel_h5_t3 = news; Plane Surface(s_g2_diel_h5_t3) = {cl_g2_diel_h5_t3};
s_g2_diel_h5_t4 = news; Plane Surface(s_g2_diel_h5_t4) = {cl_g2_diel_h5_t4};
s_g2_diel_h1_tm1 = news; Surface(s_g2_diel_h1_tm1) = {cl_g2_diel_h1_tm1};
s_g2_diel_h1_tm1 = news; Surface(s_g2_diel_h1_tm1) = {cl_g2_diel_h1_tm1};
s_g2_diel_h2_tm1 = news; Surface(s_g2_diel_h2_tm1) = {cl_g2_diel_h2_tm1};
s_g2_diel_h3_tm1 = news; Surface(s_g2_diel_h3_tm1) = {cl_g2_diel_h3_tm1};
s_g2_diel_h4_tm1 = news; Surface(s_g2_diel_h4_tm1) = {cl_g2_diel_h4_tm1};
s_g2_diel_h5_tm1 = news; Surface(s_g2_diel_h5_tm1) = {cl_g2_diel_h5_tm1};
s_g2_diel_h5_tm2 = news; Surface(s_g2_diel_h5_tm2) = {cl_g2_diel_h5_tm2};
s_g2_diel_h5_tm3 = news; Surface(s_g2_diel_h5_tm3) = {cl_g2_diel_h5_tm3};
s_g2_diel_h5_tm4 = news; Surface(s_g2_diel_h5_tm4) = {cl_g2_diel_h5_tm4};
s_g2_diel_h1_mb1 = news; Surface(s_g2_diel_h1_mb1) = {cl_g2_diel_h1_mb1};
s_g2_diel_h2_mb1 = news; Surface(s_g2_diel_h2_mb1) = {cl_g2_diel_h2_mb1};
s_g2_diel_h3_mb1 = news; Surface(s_g2_diel_h3_mb1) = {cl_g2_diel_h3_mb1};
s_g2_diel_h4_mb1 = news; Surface(s_g2_diel_h4_mb1) = {cl_g2_diel_h4_mb1};
s_g2_diel_h5_mb1 = news; Surface(s_g2_diel_h5_mb1) = {cl_g2_diel_h5_mb1};
s_g2_diel_h5_mb2 = news; Surface(s_g2_diel_h5_mb2) = {cl_g2_diel_h5_mb2};
s_g2_diel_h5_mb3 = news; Surface(s_g2_diel_h5_mb3) = {cl_g2_diel_h5_mb3};
s_g2_diel_h5_mb4 = news; Surface(s_g2_diel_h5_mb4) = {cl_g2_diel_h5_mb4};
s_g2_diel_h12_tb1 = news; Plane Surface(s_g2_diel_h12_tb1) = {cl_g2_diel_h12_tb1};
s_g2_diel_h23_tb1 = news; Plane Surface(s_g2_diel_h23_tb1) = {cl_g2_diel_h23_tb1};
s_g2_diel_h34_tb1 = news; Plane Surface(s_g2_diel_h34_tb1) = {cl_g2_diel_h34_tb1};
s_g2_diel_h41_tb1 = news; Plane Surface(s_g2_diel_h41_tb1) = {cl_g2_diel_h41_tb1};
s_g2_diel_h1_b1 = news; Plane Surface(s_g2_diel_h1_b1) = {cl_g2_diel_h1_b1};
s_g2_diel_h2_b1 = news; Plane Surface(s_g2_diel_h2_b1) = {cl_g2_diel_h2_b1};
s_g2_diel_h3_b1 = news; Plane Surface(s_g2_diel_h3_b1) = {cl_g2_diel_h3_b1};
s_g2_diel_h4_b1 = news; Plane Surface(s_g2_diel_h4_b1) = {cl_g2_diel_h4_b1};
s_g2_diel_h5_b1 = news; Plane Surface(s_g2_diel_h5_b1) = {cl_g2_diel_h5_b1};
s_g2_diel_h5_b2 = news; Plane Surface(s_g2_diel_h5_b2) = {cl_g2_diel_h5_b2};
s_g2_diel_h5_b3 = news; Plane Surface(s_g2_diel_h5_b3) = {cl_g2_diel_h5_b3};
s_g2_diel_h5_b4 = news; Plane Surface(s_g2_diel_h5_b4) = {cl_g2_diel_h5_b4};
/*---------------------------------- Surface Loop and Volume Definitions -------------------------------------*/
// surface loop definitions of upper plate
sl_g2_uc = newsl; Surface Loop(sl_g2_uc) = {s_g2_uc_t1,
    s_g2_uc_h1_tb1, s_g2_uc_h2_tb1, s_g2_uc_h3_tb1, s_g2_uc_h4_tb1,
    s_g2_uc_h12_tb1, s_g2_uc_h23_tb1, s_g2_uc_h34_tb1, s_g2_uc_h41_tb1,
    s_g2_uc_h5_tb1, s_g2_uc_h5_tb2, s_g2_uc_h5_tb3, s_g2_uc_h5_tb4,
    s_g2_uc_b1};
// surface loop definitions of upper plate
sl_g2_lc = newsl; Surface Loop(sl_g2_lc) = {s_g2_lc_t1,
    s_g2_lc_h1_tb1, s_g2_lc_h2_tb1, s_g2_lc_h3_tb1, s_g2_lc_h4_tb1,
    s_g2_lc_h12_tb1, s_g2_lc_h23_tb1, s_g2_lc_h34_tb1, s_g2_lc_h41_tb1,
    s_g2_lc_h5_tb1, s_g2_lc_h5_tb2, s_g2_lc_h5_tb3, s_g2_lc_h5_tb4,
    s_g2_lc_b1};
// surface loop definitions of dielectric layer
sl_g2_diel = newsl; Surface Loop(sl_g2_diel) = {
    s_g2_uc_b1,
    s_g2_diel_h1_t1, s_g2_diel_h2_t1, s_g2_diel_h3_t1, s_g2_diel_h4_t1,
    s_g2_diel_h5_t1, s_g2_diel_h5_t2, s_g2_diel_h5_t3, s_g2_diel_h5_t4,
    s_g2_diel_h1_tm1, s_g2_diel_h2_tm1, s_g2_diel_h3_tm1, s_g2_diel_h4_tm1,
    s_g2_diel_h5_tm1, s_g2_diel_h5_tm2, s_g2_diel_h5_tm3, s_g2_diel_h5_tm4,
    s_g2_diel_h1_mb1, s_g2_diel_h2_mb1, s_g2_diel_h3_mb1, s_g2_diel_h4_mb1,
    s_g2_diel_h5_mb1, s_g2_diel_h5_mb2, s_g2_diel_h5_mb3, s_g2_diel_h5_mb4,
    s_g2_diel_h12_tb1, s_g2_diel_h23_tb1, s_g2_diel_h34_tb1, s_g2_diel_h41_tb1,
    s_g2_diel_h1_b1, s_g2_diel_h2_b1, s_g2_diel_h3_b1, s_g2_diel_h4_b1,
    s_g2_diel_h5_b1, s_g2_diel_h5_b2, s_g2_diel_h5_b3, s_g2_diel_h5_b4,
    s_g2_lc_t1
};
v_g2_uc = newv; Volume(v_g2_uc) = {sl_g2_uc};
v_g2_lc = newv; Volume(v_g2_lc) = {sl_g2_lc};
v_g2_diel = newv; Volume(v_g2_diel) = {sl_g2_diel};
// ---------------------------------------------------------------------------------------------------------------
// GEM3
// ---------------------------------------------------------------------------------------------------------------
/*---------------------------------- Point Definitions -------------------------------------*/
// points on upper copper plate(top and bottom plane)
// hole1 (bottom left quarter hole)
p_g3_uc_h1_t1 = newp; Point(p_g3_uc_h1_t1) = {-bdrHalfX, -bdrHalfY,       z3 + tD/2 + tC, mshSizeHoleCenter};
p_g3_uc_h1_t2 = newp; Point(p_g3_uc_h1_t2) = {-bdrHalfX, -bdrHalfY + rE,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h1_t3 = newp; Point(p_g3_uc_h1_t3) = {-bdrHalfX + rE, -bdrHalfY,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h1_b1 = newp; Point(p_g3_uc_h1_b1) = {-bdrHalfX, -bdrHalfY,       z3 + tD/2, mshSizeHoleCenter};
p_g3_uc_h1_b2 = newp; Point(p_g3_uc_h1_b2) = {-bdrHalfX, -bdrHalfY + rE,  z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h1_b3 = newp; Point(p_g3_uc_h1_b3) = {-bdrHalfX + rE, -bdrHalfY,  z3 + tD/2, mshSizeEtchHole};
// hole2 (bottom right quarter hole)
p_g3_uc_h2_t1 = newp; Point(p_g3_uc_h2_t1) = {bdrHalfX, -bdrHalfY,        z3 + tD/2 + tC, mshSizeHoleCenter};
p_g3_uc_h2_t2 = newp; Point(p_g3_uc_h2_t2) = {bdrHalfX - rE, -bdrHalfY,   z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h2_t3 = newp; Point(p_g3_uc_h2_t3) = {bdrHalfX, -bdrHalfY + rE,   z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h2_b1 = newp; Point(p_g3_uc_h2_b1) = {bdrHalfX, -bdrHalfY,        z3 + tD/2, mshSizeHoleCenter};
p_g3_uc_h2_b2 = newp; Point(p_g3_uc_h2_b2) = {bdrHalfX - rE, -bdrHalfY,   z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h2_b3 = newp; Point(p_g3_uc_h2_b3) = {bdrHalfX, -bdrHalfY + rE,   z3 + tD/2, mshSizeEtchHole};
// hole3 (top right quarter hole)
p_g3_uc_h3_t1 = newp; Point(p_g3_uc_h3_t1) = {bdrHalfX, bdrHalfY,       z3 + tD/2 + tC, mshSizeHoleCenter};
p_g3_uc_h3_t2 = newp; Point(p_g3_uc_h3_t2) = {bdrHalfX, bdrHalfY - rE,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h3_t3 = newp; Point(p_g3_uc_h3_t3) = {bdrHalfX - rE, bdrHalfY,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h3_b1 = newp; Point(p_g3_uc_h3_b1) = {bdrHalfX, bdrHalfY,       z3 + tD/2, mshSizeHoleCenter};
p_g3_uc_h3_b2 = newp; Point(p_g3_uc_h3_b2) = {bdrHalfX, bdrHalfY - rE,  z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h3_b3 = newp; Point(p_g3_uc_h3_b3) = {bdrHalfX - rE, bdrHalfY,  z3 + tD/2, mshSizeEtchHole};
// hole4 (top left quarter hole)
p_g3_uc_h4_t1 = newp; Point(p_g3_uc_h4_t1) = {-bdrHalfX, bdrHalfY,       z3 + tD/2 + tC, mshSizeHoleCenter};
p_g3_uc_h4_t2 = newp; Point(p_g3_uc_h4_t2) = {-bdrHalfX + rE, bdrHalfY,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h4_t3 = newp; Point(p_g3_uc_h4_t3) = {-bdrHalfX, bdrHalfY - rE,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h4_b1 = newp; Point(p_g3_uc_h4_b1) = {-bdrHalfX, bdrHalfY,       z3 + tD/2, mshSizeHoleCenter};
p_g3_uc_h4_b2 = newp; Point(p_g3_uc_h4_b2) = {-bdrHalfX + rE, bdrHalfY,  z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h4_b3 = newp; Point(p_g3_uc_h4_b3) = {-bdrHalfX, bdrHalfY - rE,  z3 + tD/2, mshSizeEtchHole};
// hole5 (center circular hole)
p_g3_uc_h5_t1 = newp; Point(p_g3_uc_h5_t1) = {0, 0,   z3 + tD/2 + tC, mshSizeHoleCenter};
p_g3_uc_h5_t2 = newp; Point(p_g3_uc_h5_t2) = {rE, 0,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h5_t3 = newp; Point(p_g3_uc_h5_t3) = {0, rE,  z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h5_t4 = newp; Point(p_g3_uc_h5_t4) = {-rE, 0, z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h5_t5 = newp; Point(p_g3_uc_h5_t5) = {0, -rE, z3 + tD/2 + tC, mshSizeEtchHole};
p_g3_uc_h5_b1 = newp; Point(p_g3_uc_h5_b1) = {0, 0,   z3 + tD/2, mshSizeHoleCenter};
p_g3_uc_h5_b2 = newp; Point(p_g3_uc_h5_b2) = {rE, 0,  z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h5_b3 = newp; Point(p_g3_uc_h5_b3) = {0, rE,  z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h5_b4 = newp; Point(p_g3_uc_h5_b4) = {-rE, 0, z3 + tD/2, mshSizeEtchHole};
p_g3_uc_h5_b5 = newp; Point(p_g3_uc_h5_b5) = {0, -rE, z3 + tD/2, mshSizeEtchHole};

// points on lower copper plate (top and bottom plane)
p_g3_lc_h1_t1 = newp; Point(p_g3_lc_h1_t1) = {-bdrHalfX, -bdrHalfY,       z3 - tD/2, mshSizeHoleCenter};
p_g3_lc_h1_t2 = newp; Point(p_g3_lc_h1_t2) = {-bdrHalfX, -bdrHalfY + rE,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h1_t3 = newp; Point(p_g3_lc_h1_t3) = {-bdrHalfX + rE, -bdrHalfY,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h1_b1 = newp; Point(p_g3_lc_h1_b1) = {-bdrHalfX, -bdrHalfY,       z3 - (tD/2 + tC), mshSizeHoleCenter};
p_g3_lc_h1_b2 = newp; Point(p_g3_lc_h1_b2) = {-bdrHalfX, -bdrHalfY + rE,  z3 - (tD/2 + tC), mshSizeEtchHole};
p_g3_lc_h1_b3 = newp; Point(p_g3_lc_h1_b3) = {-bdrHalfX + rE, -bdrHalfY,  z3 - (tD/2 + tC), mshSizeEtchHole};
// hole2 (bottom right quarter hole)
p_g3_lc_h2_t1 = newp; Point(p_g3_lc_h2_t1) = {bdrHalfX, -bdrHalfY,        z3 - tD/2, mshSizeHoleCenter};
p_g3_lc_h2_t2 = newp; Point(p_g3_lc_h2_t2) = {bdrHalfX - rE, -bdrHalfY,   z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h2_t3 = newp; Point(p_g3_lc_h2_t3) = {bdrHalfX, -bdrHalfY + rE,   z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h2_b1 = newp; Point(p_g3_lc_h2_b1) = {bdrHalfX, -bdrHalfY,        z3 - (tD/2 + tC), mshSizeHoleCenter};
p_g3_lc_h2_b2 = newp; Point(p_g3_lc_h2_b2) = {bdrHalfX - rE, -bdrHalfY,   z3 - (tD/2 + tC), mshSizeEtchHole};
p_g3_lc_h2_b3 = newp; Point(p_g3_lc_h2_b3) = {bdrHalfX, -bdrHalfY + rE,   z3 - (tD/2 + tC), mshSizeEtchHole};
// hole3 (top right quarter hole)
p_g3_lc_h3_t1 = newp; Point(p_g3_lc_h3_t1) = {bdrHalfX, bdrHalfY,       z3 - tD/2, mshSizeHoleCenter};
p_g3_lc_h3_t2 = newp; Point(p_g3_lc_h3_t2) = {bdrHalfX, bdrHalfY - rE,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h3_t3 = newp; Point(p_g3_lc_h3_t3) = {bdrHalfX - rE, bdrHalfY,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h3_b1 = newp; Point(p_g3_lc_h3_b1) = {bdrHalfX, bdrHalfY,       z3 - (tD/2 + tC), mshSizeHoleCenter};
p_g3_lc_h3_b2 = newp; Point(p_g3_lc_h3_b2) = {bdrHalfX, bdrHalfY - rE,  z3 - (tD/2 + tC), mshSizeEtchHole};
p_g3_lc_h3_b3 = newp; Point(p_g3_lc_h3_b3) = {bdrHalfX - rE, bdrHalfY,  z3 - (tD/2 + tC), mshSizeEtchHole};
// hole4 (top left quarter hole)
p_g3_lc_h4_t1 = newp; Point(p_g3_lc_h4_t1) = {-bdrHalfX, bdrHalfY,       z3 - tD/2, mshSizeHoleCenter};
p_g3_lc_h4_t2 = newp; Point(p_g3_lc_h4_t2) = {-bdrHalfX + rE, bdrHalfY,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h4_t3 = newp; Point(p_g3_lc_h4_t3) = {-bdrHalfX, bdrHalfY - rE,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h4_b1 = newp; Point(p_g3_lc_h4_b1) = {-bdrHalfX, bdrHalfY,       z3 - (tD/2 + tC), mshSizeHoleCenter};
p_g3_lc_h4_b2 = newp; Point(p_g3_lc_h4_b2) = {-bdrHalfX + rE, bdrHalfY,  z3 - (tD/2 + tC), mshSizeEtchHole};
p_g3_lc_h4_b3 = newp; Point(p_g3_lc_h4_b3) = {-bdrHalfX, bdrHalfY - rE,  z3 - (tD/2 + tC), mshSizeEtchHole};
// hole5 (center circular quarter hole)
p_g3_lc_h5_t1 = newp; Point(p_g3_lc_h5_t1) = {0, 0,   z3 - tD/2, mshSizeHoleCenter};
p_g3_lc_h5_t2 = newp; Point(p_g3_lc_h5_t2) = {rE, 0,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h5_t3 = newp; Point(p_g3_lc_h5_t3) = {0, rE,  z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h5_t4 = newp; Point(p_g3_lc_h5_t4) = {-rE, 0, z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h5_t5 = newp; Point(p_g3_lc_h5_t5) = {0, -rE, z3 - tD/2, mshSizeEtchHole};
p_g3_lc_h5_b1 = newp; Point(p_g3_lc_h5_b1) = {0, 0,   z3 - tD/2 - tC, mshSizeHoleCenter};
p_g3_lc_h5_b2 = newp; Point(p_g3_lc_h5_b2) = {rE, 0,  z3 - tD/2 - tC, mshSizeEtchHole};
p_g3_lc_h5_b3 = newp; Point(p_g3_lc_h5_b3) = {0, rE,  z3 - tD/2 - tC, mshSizeEtchHole};
p_g3_lc_h5_b4 = newp; Point(p_g3_lc_h5_b4) = {-rE, 0, z3 - tD/2 - tC, mshSizeEtchHole};
p_g3_lc_h5_b5 = newp; Point(p_g3_lc_h5_b5) = {0, -rE, z3 - tD/2 - tC, mshSizeEtchHole};

// points on upper dielectric layer(top, middel and bottom plane)
// hole1 (bottom left quarter hole)
p_g3_diel_h1_t1 = p_g3_uc_h1_b1;
p_g3_diel_h1_t2 = newp; Point(p_g3_diel_h1_t2) = {-bdrHalfX, -bdrHalfY + rO,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h1_t3 = newp; Point(p_g3_diel_h1_t3) = {-bdrHalfX + rO, -bdrHalfY,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h1_m1 = newp; Point(p_g3_diel_h1_m1) = {-bdrHalfX, -bdrHalfY,       z3, mshSizeHoleCenter};
p_g3_diel_h1_m2 = newp; Point(p_g3_diel_h1_m2) = {-bdrHalfX, -bdrHalfY + rI,  z3, mshSizeDielInnerHole};
p_g3_diel_h1_m3 = newp; Point(p_g3_diel_h1_m3) = {-bdrHalfX + rI, -bdrHalfY,  z3, mshSizeDielInnerHole};
p_g3_diel_h1_b1 = p_g3_lc_h1_t1;
p_g3_diel_h1_b2 = newp; Point(p_g3_diel_h1_b2) = {-bdrHalfX, -bdrHalfY + rO,  z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h1_b3 = newp; Point(p_g3_diel_h1_b3) = {-bdrHalfX + rO, -bdrHalfY,  z3 - tD/2, mshSizeDielOuterHole};
// hole2 (bottom right quarter hole)
p_g3_diel_h2_t1 = p_g3_uc_h2_b1;
p_g3_diel_h2_t2 = newp; Point(p_g3_diel_h2_t2) = {bdrHalfX - rO, -bdrHalfY,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h2_t3 = newp; Point(p_g3_diel_h2_t3) = {bdrHalfX, -bdrHalfY + rO,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h2_m1 = newp; Point(p_g3_diel_h2_m1) = {bdrHalfX, -bdrHalfY,       z3, mshSizeHoleCenter};
p_g3_diel_h2_m2 = newp; Point(p_g3_diel_h2_m2) = {bdrHalfX - rI, -bdrHalfY,  z3, mshSizeDielInnerHole};
p_g3_diel_h2_m3 = newp; Point(p_g3_diel_h2_m3) = {bdrHalfX, -bdrHalfY + rI,  z3, mshSizeDielInnerHole};
p_g3_diel_h2_b1 = p_g3_lc_h2_t1;
p_g3_diel_h2_b2 = newp; Point(p_g3_diel_h2_b2) = {bdrHalfX - rO, -bdrHalfY,  z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h2_b3 = newp; Point(p_g3_diel_h2_b3) = {bdrHalfX, -bdrHalfY + rO,  z3 - tD/2, mshSizeDielOuterHole};
// hole3 (top right quarter hole)
p_g3_diel_h3_t1 = p_g3_uc_h3_b1;
p_g3_diel_h3_t2 = newp; Point(p_g3_diel_h3_t2) = {bdrHalfX, bdrHalfY - rO,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h3_t3 = newp; Point(p_g3_diel_h3_t3) = {bdrHalfX - rO, bdrHalfY,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h3_m1 = newp; Point(p_g3_diel_h3_m1) = {bdrHalfX, bdrHalfY,       z3, mshSizeHoleCenter};
p_g3_diel_h3_m2 = newp; Point(p_g3_diel_h3_m2) = {bdrHalfX, bdrHalfY - rI,  z3, mshSizeDielOuterHole};
p_g3_diel_h3_m3 = newp; Point(p_g3_diel_h3_m3) = {bdrHalfX - rI, bdrHalfY,  z3, mshSizeDielOuterHole};
p_g3_diel_h3_b1 = p_g3_lc_h3_t1;
p_g3_diel_h3_b2 = newp; Point(p_g3_diel_h3_b2) = {bdrHalfX, bdrHalfY - rO,  z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h3_b3 = newp; Point(p_g3_diel_h3_b3) = {bdrHalfX - rO, bdrHalfY,  z3 - tD/2, mshSizeDielOuterHole};
// hole4 (top left quarter hole)
p_g3_diel_h4_t1 = p_g3_uc_h4_b1;
p_g3_diel_h4_t2 = newp; Point(p_g3_diel_h4_t2) = {-bdrHalfX + rO, bdrHalfY,    z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h4_t3 = newp; Point(p_g3_diel_h4_t3) = {-bdrHalfX, bdrHalfY - rO,    z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h4_m1 = newp; Point(p_g3_diel_h4_m1) = {-bdrHalfX, bdrHalfY,         z3, mshSizeHoleCenter};
p_g3_diel_h4_m2 = newp; Point(p_g3_diel_h4_m2) = {-bdrHalfX + rI, bdrHalfY,    z3, mshSizeDielInnerHole};
p_g3_diel_h4_m3 = newp; Point(p_g3_diel_h4_m3) = {-bdrHalfX, bdrHalfY - rI ,   z3, mshSizeDielInnerHole};
p_g3_diel_h4_b1 = p_g3_lc_h4_t1;
p_g3_diel_h4_b2 = newp; Point(p_g3_diel_h4_b2) = {-bdrHalfX + rO, bdrHalfY,    z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h4_b3 = newp; Point(p_g3_diel_h4_b3) = {-bdrHalfX, bdrHalfY - rO,    z3 - tD/2, mshSizeDielOuterHole};
// hole5 (center circular quarter hole)
p_g3_diel_h5_t1 = p_g3_uc_h5_b1;
p_g3_diel_h5_t2 = newp; Point(p_g3_diel_h5_t2) = {rO, 0,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_t3 = newp; Point(p_g3_diel_h5_t3) = {0, rO,  z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_t4 = newp; Point(p_g3_diel_h5_t4) = {-rO, 0, z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_t5 = newp; Point(p_g3_diel_h5_t5) = {0, -rO, z3 + tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_m1 = newp; Point(p_g3_diel_h5_m1) = {0, 0,   z3, mshSizeHoleCenter};
p_g3_diel_h5_m2 = newp; Point(p_g3_diel_h5_m2) = {rI, 0,  z3, mshSizeDielInnerHole};
p_g3_diel_h5_m3 = newp; Point(p_g3_diel_h5_m3) = {0, rI,  z3, mshSizeDielInnerHole};
p_g3_diel_h5_m4 = newp; Point(p_g3_diel_h5_m4) = {-rI, 0, z3, mshSizeDielInnerHole};
p_g3_diel_h5_m5 = newp; Point(p_g3_diel_h5_m5) = {0, -rI, z3, mshSizeDielInnerHole};
p_g3_diel_h5_b1 = p_g3_lc_h5_t1;
p_g3_diel_h5_b2 = newp; Point(p_g3_diel_h5_b2) = {rO, 0,  z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_b3 = newp; Point(p_g3_diel_h5_b3) = {0, rO,  z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_b4 = newp; Point(p_g3_diel_h5_b4) = {-rO, 0, z3 - tD/2, mshSizeDielOuterHole};
p_g3_diel_h5_b5 = newp; Point(p_g3_diel_h5_b5) = {0, -rO, z3 - tD/2, mshSizeDielOuterHole};
/*---------------------------------- Line Definitions -------------------------------------*/
// lines on upper plane
// on the top plane
l_g3_uc_h1_t1 = newl; Circle(l_g3_uc_h1_t1) = {p_g3_uc_h1_t2, p_g3_uc_h1_t1, p_g3_uc_h1_t3}; // hole 1
l_g3_uc_h12_t1 = newl; Line(l_g3_uc_h12_t1) = {p_g3_uc_h1_t3, p_g3_uc_h2_t2}; // between hole 1 and hole 2
l_g3_uc_h2_t1 = newl; Circle(l_g3_uc_h2_t1) = {p_g3_uc_h2_t2, p_g3_uc_h2_t1, p_g3_uc_h2_t3}; // hole 2
l_g3_uc_h23_t1 = newl; Line(l_g3_uc_h23_t1) = {p_g3_uc_h2_t3, p_g3_uc_h3_t2}; // between hole 2 and hole 3
l_g3_uc_h3_t1 = newl; Circle(l_g3_uc_h3_t1) = {p_g3_uc_h3_t2, p_g3_uc_h3_t1, p_g3_uc_h3_t3}; // hole 3
l_g3_uc_h34_t1 = newl; Line(l_g3_uc_h34_t1) = {p_g3_uc_h3_t3, p_g3_uc_h4_t2}; // between hole 3 and hole 4
l_g3_uc_h4_t1 = newl; Circle(l_g3_uc_h4_t1) = {p_g3_uc_h4_t2, p_g3_uc_h4_t1, p_g3_uc_h4_t3}; // hole 4
l_g3_uc_h41_t1 = newl; Line(l_g3_uc_h41_t1) = {p_g3_uc_h4_t3, p_g3_uc_h1_t2}; // between hole 4 and hole 1
l_g3_uc_h5_t1 = newl; Circle(l_g3_uc_h5_t1) = {p_g3_uc_h5_t2, p_g3_uc_h5_t1, p_g3_uc_h5_t3}; // hole 5
l_g3_uc_h5_t2 = newl; Circle(l_g3_uc_h5_t2) = {p_g3_uc_h5_t3, p_g3_uc_h5_t1, p_g3_uc_h5_t4};
l_g3_uc_h5_t3 = newl; Circle(l_g3_uc_h5_t3) = {p_g3_uc_h5_t4, p_g3_uc_h5_t1, p_g3_uc_h5_t5};
l_g3_uc_h5_t4 = newl; Circle(l_g3_uc_h5_t4) = {p_g3_uc_h5_t5, p_g3_uc_h5_t1, p_g3_uc_h5_t2};
// on the bottom plane
l_g3_uc_h1_b1 = newl; Circle(l_g3_uc_h1_b1) = {p_g3_uc_h1_b2, p_g3_uc_h1_b1, p_g3_uc_h1_b3}; // hole 1
l_g3_uc_h12_b1 = newl; Line(l_g3_uc_h12_b1) = {p_g3_uc_h1_b3, p_g3_uc_h2_b2}; // between hole 1 and hole 2
l_g3_uc_h2_b1 = newl; Circle(l_g3_uc_h2_b1) = {p_g3_uc_h2_b2, p_g3_uc_h2_b1, p_g3_uc_h2_b3}; // hole 2
l_g3_uc_h23_b1 = newl; Line(l_g3_uc_h23_b1) = {p_g3_uc_h2_b3, p_g3_uc_h3_b2}; // between hole 2 and hole 3
l_g3_uc_h3_b1 = newl; Circle(l_g3_uc_h3_b1) = {p_g3_uc_h3_b2, p_g3_uc_h3_b1, p_g3_uc_h3_b3}; // hole 3
l_g3_uc_h34_b1 = newl; Line(l_g3_uc_h34_b1) = {p_g3_uc_h3_b3, p_g3_uc_h4_b2}; // between hole 3 and hole 4
l_g3_uc_h4_b1 = newl; Circle(l_g3_uc_h4_b1) = {p_g3_uc_h4_b2, p_g3_uc_h4_b1, p_g3_uc_h4_b3}; // hole 4
l_g3_uc_h41_b1 = newl; Line(l_g3_uc_h41_b1) = {p_g3_uc_h4_b3, p_g3_uc_h1_b2}; // between hole 4 and hole 1
l_g3_uc_h5_b1 = newl; Circle(l_g3_uc_h5_b1) = {p_g3_uc_h5_b2, p_g3_uc_h5_b1, p_g3_uc_h5_b3}; // hole 5
l_g3_uc_h5_b2 = newl; Circle(l_g3_uc_h5_b2) = {p_g3_uc_h5_b3, p_g3_uc_h5_b1, p_g3_uc_h5_b4};
l_g3_uc_h5_b3 = newl; Circle(l_g3_uc_h5_b3) = {p_g3_uc_h5_b4, p_g3_uc_h5_b1, p_g3_uc_h5_b5};
l_g3_uc_h5_b4 = newl; Circle(l_g3_uc_h5_b4) = {p_g3_uc_h5_b5, p_g3_uc_h5_b1, p_g3_uc_h5_b2};
// connecting top and bottom
l_g3_uc_h1_tb1 = newl; Line(l_g3_uc_h1_tb1) = {p_g3_uc_h1_t1, p_g3_uc_h1_b1}; // hole 1
l_g3_uc_h1_tb2 = newl; Line(l_g3_uc_h1_tb2) = {p_g3_uc_h1_t2, p_g3_uc_h1_b2}; 
l_g3_uc_h1_tb3 = newl; Line(l_g3_uc_h1_tb3) = {p_g3_uc_h1_t3, p_g3_uc_h1_b3};
l_g3_uc_h2_tb1 = newl; Line(l_g3_uc_h2_tb1) = {p_g3_uc_h2_t1, p_g3_uc_h2_b1}; // hole 2
l_g3_uc_h2_tb2 = newl; Line(l_g3_uc_h2_tb2) = {p_g3_uc_h2_t2, p_g3_uc_h2_b2};
l_g3_uc_h2_tb3 = newl; Line(l_g3_uc_h2_tb3) = {p_g3_uc_h2_t3, p_g3_uc_h2_b3};
l_g3_uc_h3_tb1 = newl; Line(l_g3_uc_h3_tb1) = {p_g3_uc_h3_t1, p_g3_uc_h3_b1}; // hole 3
l_g3_uc_h3_tb2 = newl; Line(l_g3_uc_h3_tb2) = {p_g3_uc_h3_t2, p_g3_uc_h3_b2};
l_g3_uc_h3_tb3 = newl; Line(l_g3_uc_h3_tb3) = {p_g3_uc_h3_t3, p_g3_uc_h3_b3};
l_g3_uc_h4_tb1 = newl; Line(l_g3_uc_h4_tb1) = {p_g3_uc_h4_t1, p_g3_uc_h4_b1}; // hole 4
l_g3_uc_h4_tb2 = newl; Line(l_g3_uc_h4_tb2) = {p_g3_uc_h4_t2, p_g3_uc_h4_b2};
l_g3_uc_h4_tb3 = newl; Line(l_g3_uc_h4_tb3) = {p_g3_uc_h4_t3, p_g3_uc_h4_b3};
l_g3_uc_h5_tb1 = newl; Line(l_g3_uc_h5_tb1) = {p_g3_uc_h5_t1, p_g3_uc_h5_b1}; // hole 5
l_g3_uc_h5_tb2 = newl; Line(l_g3_uc_h5_tb2) = {p_g3_uc_h5_t2, p_g3_uc_h5_b2};
l_g3_uc_h5_tb3 = newl; Line(l_g3_uc_h5_tb3) = {p_g3_uc_h5_t3, p_g3_uc_h5_b3};
l_g3_uc_h5_tb4 = newl; Line(l_g3_uc_h5_tb4) = {p_g3_uc_h5_t4, p_g3_uc_h5_b4};
l_g3_uc_h5_tb5 = newl; Line(l_g3_uc_h5_tb5) = {p_g3_uc_h5_t5, p_g3_uc_h5_b5};

// lines on lower plane
// on the top plane
l_g3_lc_h1_t1 = newl; Circle(l_g3_lc_h1_t1) = {p_g3_lc_h1_t2, p_g3_lc_h1_t1, p_g3_lc_h1_t3}; // hole 1
l_g3_lc_h12_t1 = newl; Line(l_g3_lc_h12_t1) = {p_g3_lc_h1_t3, p_g3_lc_h2_t2}; // between hole 1 and hole 2
l_g3_lc_h2_t1 = newl; Circle(l_g3_lc_h2_t1) = {p_g3_lc_h2_t2, p_g3_lc_h2_t1, p_g3_lc_h2_t3}; // hole 2
l_g3_lc_h23_t1 = newl; Line(l_g3_lc_h23_t1) = {p_g3_lc_h2_t3, p_g3_lc_h3_t2}; // between hole 2 and hole 3
l_g3_lc_h3_t1 = newl; Circle(l_g3_lc_h3_t1) = {p_g3_lc_h3_t2, p_g3_lc_h3_t1, p_g3_lc_h3_t3}; // hole 3
l_g3_lc_h34_t1 = newl; Line(l_g3_lc_h34_t1) = {p_g3_lc_h3_t3, p_g3_lc_h4_t2}; // between hole 3 and hole 4
l_g3_lc_h4_t1 = newl; Circle(l_g3_lc_h4_t1) = {p_g3_lc_h4_t2, p_g3_lc_h4_t1, p_g3_lc_h4_t3}; // hole 4
l_g3_lc_h41_t1 = newl; Line(l_g3_lc_h41_t1) = {p_g3_lc_h4_t3, p_g3_lc_h1_t2}; // between hole 4 and hole 1
l_g3_lc_h5_t1 = newl; Circle(l_g3_lc_h5_t1) = {p_g3_lc_h5_t2, p_g3_lc_h5_t1, p_g3_lc_h5_t3}; // hole 5
l_g3_lc_h5_t2 = newl; Circle(l_g3_lc_h5_t2) = {p_g3_lc_h5_t3, p_g3_lc_h5_t1, p_g3_lc_h5_t4};
l_g3_lc_h5_t3 = newl; Circle(l_g3_lc_h5_t3) = {p_g3_lc_h5_t4, p_g3_lc_h5_t1, p_g3_lc_h5_t5};
l_g3_lc_h5_t4 = newl; Circle(l_g3_lc_h5_t4) = {p_g3_lc_h5_t5, p_g3_lc_h5_t1, p_g3_lc_h5_t2};
// on the bottom plane
l_g3_lc_h1_b1 = newl; Circle(l_g3_lc_h1_b1) = {p_g3_lc_h1_b2, p_g3_lc_h1_b1, p_g3_lc_h1_b3}; // hole 1
l_g3_lc_h12_b1 = newl; Line(l_g3_lc_h12_b1) = {p_g3_lc_h1_b3, p_g3_lc_h2_b2}; // between hole 1 and hole 2
l_g3_lc_h2_b1 = newl; Circle(l_g3_lc_h2_b1) = {p_g3_lc_h2_b2, p_g3_lc_h2_b1, p_g3_lc_h2_b3}; // hole 2
l_g3_lc_h23_b1 = newl; Line(l_g3_lc_h23_b1) = {p_g3_lc_h2_b3, p_g3_lc_h3_b2}; // between hole 2 and hole 3
l_g3_lc_h3_b1 = newl; Circle(l_g3_lc_h3_b1) = {p_g3_lc_h3_b2, p_g3_lc_h3_b1, p_g3_lc_h3_b3}; // hole 3
l_g3_lc_h34_b1 = newl; Line(l_g3_lc_h34_b1) = {p_g3_lc_h3_b3, p_g3_lc_h4_b2}; // between hole 3 and hole 4
l_g3_lc_h4_b1 = newl; Circle(l_g3_lc_h4_b1) = {p_g3_lc_h4_b2, p_g3_lc_h4_b1, p_g3_lc_h4_b3}; // hole 4
l_g3_lc_h41_b1 = newl; Line(l_g3_lc_h41_b1) = {p_g3_lc_h4_b3, p_g3_lc_h1_b2}; // between hole 4 and hole 1
l_g3_lc_h5_b1 = newl; Circle(l_g3_lc_h5_b1) = {p_g3_lc_h5_b2, p_g3_lc_h5_b1, p_g3_lc_h5_b3}; // hole 5
l_g3_lc_h5_b2 = newl; Circle(l_g3_lc_h5_b2) = {p_g3_lc_h5_b3, p_g3_lc_h5_b1, p_g3_lc_h5_b4};
l_g3_lc_h5_b3 = newl; Circle(l_g3_lc_h5_b3) = {p_g3_lc_h5_b4, p_g3_lc_h5_b1, p_g3_lc_h5_b5};
l_g3_lc_h5_b4 = newl; Circle(l_g3_lc_h5_b4) = {p_g3_lc_h5_b5, p_g3_lc_h5_b1, p_g3_lc_h5_b2};
// connecting top and bottom
l_g3_lc_h1_tb1 = newl; Line(l_g3_lc_h1_tb1) = {p_g3_lc_h1_t1, p_g3_lc_h1_b1}; // hole 1
l_g3_lc_h1_tb2 = newl; Line(l_g3_lc_h1_tb2) = {p_g3_lc_h1_t2, p_g3_lc_h1_b2}; 
l_g3_lc_h1_tb3 = newl; Line(l_g3_lc_h1_tb3) = {p_g3_lc_h1_t3, p_g3_lc_h1_b3};
l_g3_lc_h2_tb1 = newl; Line(l_g3_lc_h2_tb1) = {p_g3_lc_h2_t1, p_g3_lc_h2_b1}; // hole 2
l_g3_lc_h2_tb2 = newl; Line(l_g3_lc_h2_tb2) = {p_g3_lc_h2_t2, p_g3_lc_h2_b2};
l_g3_lc_h2_tb3 = newl; Line(l_g3_lc_h2_tb3) = {p_g3_lc_h2_t3, p_g3_lc_h2_b3};
l_g3_lc_h3_tb1 = newl; Line(l_g3_lc_h3_tb1) = {p_g3_lc_h3_t1, p_g3_lc_h3_b1}; // hole 3
l_g3_lc_h3_tb2 = newl; Line(l_g3_lc_h3_tb2) = {p_g3_lc_h3_t2, p_g3_lc_h3_b2};
l_g3_lc_h3_tb3 = newl; Line(l_g3_lc_h3_tb3) = {p_g3_lc_h3_t3, p_g3_lc_h3_b3};
l_g3_lc_h4_tb1 = newl; Line(l_g3_lc_h4_tb1) = {p_g3_lc_h4_t1, p_g3_lc_h4_b1}; // hole 4
l_g3_lc_h4_tb2 = newl; Line(l_g3_lc_h4_tb2) = {p_g3_lc_h4_t2, p_g3_lc_h4_b2};
l_g3_lc_h4_tb3 = newl; Line(l_g3_lc_h4_tb3) = {p_g3_lc_h4_t3, p_g3_lc_h4_b3};
l_g3_lc_h5_tb1 = newl; Line(l_g3_lc_h5_tb1) = {p_g3_lc_h5_t1, p_g3_lc_h5_b1}; // hole 5
l_g3_lc_h5_tb2 = newl; Line(l_g3_lc_h5_tb2) = {p_g3_lc_h5_t2, p_g3_lc_h5_b2};
l_g3_lc_h5_tb3 = newl; Line(l_g3_lc_h5_tb3) = {p_g3_lc_h5_t3, p_g3_lc_h5_b3};
l_g3_lc_h5_tb4 = newl; Line(l_g3_lc_h5_tb4) = {p_g3_lc_h5_t4, p_g3_lc_h5_b4};
l_g3_lc_h5_tb5 = newl; Line(l_g3_lc_h5_tb5) = {p_g3_lc_h5_t5, p_g3_lc_h5_b5};

// lines on the dielectric layer
// on the top plane
l_g3_diel_h1_t1 = newl; Circle(l_g3_diel_h1_t1) = {p_g3_diel_h1_t2, p_g3_diel_h1_t1, p_g3_diel_h1_t3}; // hole 1
l_g3_diel_h1_t2 = newl; Line(l_g3_diel_h1_t2) = {p_g3_diel_h1_t2, p_g3_uc_h1_b2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h1_t3 = newl; Line(l_g3_diel_h1_t3) = {p_g3_diel_h1_t3, p_g3_uc_h1_b3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h2_t1 = newl; Circle(l_g3_diel_h2_t1) = {p_g3_diel_h2_t2, p_g3_diel_h2_t1, p_g3_diel_h2_t3}; // hole 2
l_g3_diel_h2_t2 = newl; Line(l_g3_diel_h2_t2) = {p_g3_diel_h2_t2, p_g3_uc_h2_b2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h2_t3 = newl; Line(l_g3_diel_h2_t3) = {p_g3_diel_h2_t3, p_g3_uc_h2_b3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h3_t1 = newl; Circle(l_g3_diel_h3_t1) = {p_g3_diel_h3_t2, p_g3_diel_h3_t1, p_g3_diel_h3_t3}; // hole 3
l_g3_diel_h3_t2 = newl; Line(l_g3_diel_h3_t2) = {p_g3_diel_h3_t2, p_g3_uc_h3_b2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h3_t3 = newl; Line(l_g3_diel_h3_t3) = {p_g3_diel_h3_t3, p_g3_uc_h3_b3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h4_t1 = newl; Circle(l_g3_diel_h4_t1) = {p_g3_diel_h4_t2, p_g3_diel_h4_t1, p_g3_diel_h4_t3}; // hole 4
l_g3_diel_h4_t2 = newl; Line(l_g3_diel_h4_t2) = {p_g3_diel_h4_t2, p_g3_uc_h4_b2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h4_t3 = newl; Line(l_g3_diel_h4_t3) = {p_g3_diel_h4_t3, p_g3_uc_h4_b3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_t1 = newl; Circle(l_g3_diel_h5_t1) = {p_g3_diel_h5_t2, p_g3_diel_h5_t1, p_g3_diel_h5_t3}; // hole 5
l_g3_diel_h5_t2 = newl; Circle(l_g3_diel_h5_t2) = {p_g3_diel_h5_t3, p_g3_diel_h5_t1, p_g3_diel_h5_t4};
l_g3_diel_h5_t3 = newl; Circle(l_g3_diel_h5_t3) = {p_g3_diel_h5_t4, p_g3_diel_h5_t1, p_g3_diel_h5_t5};
l_g3_diel_h5_t4 = newl; Circle(l_g3_diel_h5_t4) = {p_g3_diel_h5_t5, p_g3_diel_h5_t1, p_g3_diel_h5_t2};
l_g3_diel_h5_t5 = newl; Line(l_g3_diel_h5_t5) = {p_g3_diel_h5_t2, p_g3_uc_h5_b2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_t6 = newl; Line(l_g3_diel_h5_t6) = {p_g3_diel_h5_t3, p_g3_uc_h5_b3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_t7 = newl; Line(l_g3_diel_h5_t7) = {p_g3_diel_h5_t4, p_g3_uc_h5_b4}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_t8 = newl; Line(l_g3_diel_h5_t8) = {p_g3_diel_h5_t5, p_g3_uc_h5_b5}; // this line is on the surface exposed to gas by etching
// connecting the top and middel plane
l_g3_diel_h1_tm1 = newl; Line(l_g3_diel_h1_tm1) = {p_g3_diel_h1_t1, p_g3_diel_h1_m1};
l_g3_diel_h1_tm2 = newl; Line(l_g3_diel_h1_tm2) = {p_g3_diel_h1_t2, p_g3_diel_h1_m2};
l_g3_diel_h1_tm3 = newl; Line(l_g3_diel_h1_tm3) = {p_g3_diel_h1_t3, p_g3_diel_h1_m3};
l_g3_diel_h2_tm1 = newl; Line(l_g3_diel_h2_tm1) = {p_g3_diel_h2_t1, p_g3_diel_h2_m1};
l_g3_diel_h2_tm2 = newl; Line(l_g3_diel_h2_tm2) = {p_g3_diel_h2_t2, p_g3_diel_h2_m2};
l_g3_diel_h2_tm3 = newl; Line(l_g3_diel_h2_tm3) = {p_g3_diel_h2_t3, p_g3_diel_h2_m3};
l_g3_diel_h3_tm1 = newl; Line(l_g3_diel_h3_tm1) = {p_g3_diel_h3_t1, p_g3_diel_h3_m1};
l_g3_diel_h3_tm2 = newl; Line(l_g3_diel_h3_tm2) = {p_g3_diel_h3_t2, p_g3_diel_h3_m2};
l_g3_diel_h3_tm3 = newl; Line(l_g3_diel_h3_tm3) = {p_g3_diel_h3_t3, p_g3_diel_h3_m3};
l_g3_diel_h4_tm1 = newl; Line(l_g3_diel_h4_tm1) = {p_g3_diel_h4_t1, p_g3_diel_h4_m1};
l_g3_diel_h4_tm2 = newl; Line(l_g3_diel_h4_tm2) = {p_g3_diel_h4_t2, p_g3_diel_h4_m2};
l_g3_diel_h4_tm3 = newl; Line(l_g3_diel_h4_tm3) = {p_g3_diel_h4_t3, p_g3_diel_h4_m3};
l_g3_diel_h5_tm1 = newl; Line(l_g3_diel_h5_tm1) = {p_g3_diel_h5_t1, p_g3_diel_h5_m1};
l_g3_diel_h5_tm2 = newl; Line(l_g3_diel_h5_tm2) = {p_g3_diel_h5_t2, p_g3_diel_h5_m2};
l_g3_diel_h5_tm3 = newl; Line(l_g3_diel_h5_tm3) = {p_g3_diel_h5_t3, p_g3_diel_h5_m3};
l_g3_diel_h5_tm4 = newl; Line(l_g3_diel_h5_tm4) = {p_g3_diel_h5_t4, p_g3_diel_h5_m4};
l_g3_diel_h5_tm5 = newl; Line(l_g3_diel_h5_tm5) = {p_g3_diel_h5_t5, p_g3_diel_h5_m5};
// on the middle plane
l_g3_diel_h1_m1 = newl; Circle(l_g3_diel_h1_m1) = {p_g3_diel_h1_m2, p_g3_diel_h1_m1, p_g3_diel_h1_m3}; // hole 1
l_g3_diel_h12_m1 = newl; Line(l_g3_diel_h12_m1) = {p_g3_diel_h1_m3, p_g3_diel_h2_m2}; // between hole 1 and hole 2
l_g3_diel_h2_m1 = newl; Circle(l_g3_diel_h2_m1) = {p_g3_diel_h2_m2, p_g3_diel_h2_m1, p_g3_diel_h2_m3}; // hole 2
l_g3_diel_h23_m1 = newl; Line(l_g3_diel_h23_m1) = {p_g3_diel_h2_m3, p_g3_diel_h3_m2}; // between hole 2 and hole 3
l_g3_diel_h3_m1 = newl; Circle(l_g3_diel_h3_m1) = {p_g3_diel_h3_m2, p_g3_diel_h3_m1, p_g3_diel_h3_m3}; // hole 3
l_g3_diel_h34_m1 = newl; Line(l_g3_diel_h34_m1) = {p_g3_diel_h3_m3, p_g3_diel_h4_m2}; // between hole 3 and hole 4
l_g3_diel_h4_m1 = newl; Circle(l_g3_diel_h4_m1) = {p_g3_diel_h4_m2, p_g3_diel_h4_m1, p_g3_diel_h4_m3}; // hole 4
l_g3_diel_h41_m1 = newl; Line(l_g3_diel_h41_m1) = {p_g3_diel_h4_m3, p_g3_diel_h1_m2}; // between hole 4 and hole 1
l_g3_diel_h5_m1 = newl; Circle(l_g3_diel_h5_m1) = {p_g3_diel_h5_m2, p_g3_diel_h5_m1, p_g3_diel_h5_m3}; // hole 5
l_g3_diel_h5_m2 = newl; Circle(l_g3_diel_h5_m2) = {p_g3_diel_h5_m3, p_g3_diel_h5_m1, p_g3_diel_h5_m4};
l_g3_diel_h5_m3 = newl; Circle(l_g3_diel_h5_m3) = {p_g3_diel_h5_m4, p_g3_diel_h5_m1, p_g3_diel_h5_m5};
l_g3_diel_h5_m4 = newl; Circle(l_g3_diel_h5_m4) = {p_g3_diel_h5_m5, p_g3_diel_h5_m1, p_g3_diel_h5_m2};
// connecting the middle and bottom plane
l_g3_diel_h1_mb1 = newl; Line(l_g3_diel_h1_mb1) = {p_g3_diel_h1_m1, p_g3_diel_h1_b1};
l_g3_diel_h1_mb2 = newl; Line(l_g3_diel_h1_mb2) = {p_g3_diel_h1_m2, p_g3_diel_h1_b2};
l_g3_diel_h1_mb3 = newl; Line(l_g3_diel_h1_mb3) = {p_g3_diel_h1_m3, p_g3_diel_h1_b3};
l_g3_diel_h2_mb1 = newl; Line(l_g3_diel_h2_mb1) = {p_g3_diel_h2_m1, p_g3_diel_h2_b1};
l_g3_diel_h2_mb2 = newl; Line(l_g3_diel_h2_mb2) = {p_g3_diel_h2_m2, p_g3_diel_h2_b2};
l_g3_diel_h2_mb3 = newl; Line(l_g3_diel_h2_mb3) = {p_g3_diel_h2_m3, p_g3_diel_h2_b3};
l_g3_diel_h3_mb1 = newl; Line(l_g3_diel_h3_mb1) = {p_g3_diel_h3_m1, p_g3_diel_h3_b1};
l_g3_diel_h3_mb2 = newl; Line(l_g3_diel_h3_mb2) = {p_g3_diel_h3_m2, p_g3_diel_h3_b2};
l_g3_diel_h3_mb3 = newl; Line(l_g3_diel_h3_mb3) = {p_g3_diel_h3_m3, p_g3_diel_h3_b3};
l_g3_diel_h4_mb1 = newl; Line(l_g3_diel_h4_mb1) = {p_g3_diel_h4_m1, p_g3_diel_h4_b1};
l_g3_diel_h4_mb2 = newl; Line(l_g3_diel_h4_mb2) = {p_g3_diel_h4_m2, p_g3_diel_h4_b2};
l_g3_diel_h4_mb3 = newl; Line(l_g3_diel_h4_mb3) = {p_g3_diel_h4_m3, p_g3_diel_h4_b3};
l_g3_diel_h5_mb1 = newl; Line(l_g3_diel_h5_mb1) = {p_g3_diel_h5_m1, p_g3_diel_h5_b1};
l_g3_diel_h5_mb2 = newl; Line(l_g3_diel_h5_mb2) = {p_g3_diel_h5_m2, p_g3_diel_h5_b2};
l_g3_diel_h5_mb3 = newl; Line(l_g3_diel_h5_mb3) = {p_g3_diel_h5_m3, p_g3_diel_h5_b3};
l_g3_diel_h5_mb4 = newl; Line(l_g3_diel_h5_mb4) = {p_g3_diel_h5_m4, p_g3_diel_h5_b4};
l_g3_diel_h5_mb5 = newl; Line(l_g3_diel_h5_mb5) = {p_g3_diel_h5_m5, p_g3_diel_h5_b5};
// on the bottom plane
l_g3_diel_h1_b1 = newl; Circle(l_g3_diel_h1_b1) = {p_g3_diel_h1_b2, p_g3_diel_h1_b1, p_g3_diel_h1_b3}; // hole 1
l_g3_diel_h1_b2 = newl; Line(l_g3_diel_h1_b2) = {p_g3_diel_h1_b2, p_g3_lc_h1_t2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h1_b3 = newl; Line(l_g3_diel_h1_b3) = {p_g3_diel_h1_b3, p_g3_lc_h1_t3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h2_b1 = newl; Circle(l_g3_diel_h2_b1) = {p_g3_diel_h2_b2, p_g3_diel_h2_b1, p_g3_diel_h2_b3}; // hole 2
l_g3_diel_h2_b2 = newl; Line(l_g3_diel_h2_b2) = {p_g3_diel_h2_b2, p_g3_lc_h2_t2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h2_b3 = newl; Line(l_g3_diel_h2_b3) = {p_g3_diel_h2_b3, p_g3_lc_h2_t3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h3_b1 = newl; Circle(l_g3_diel_h3_b1) = {p_g3_diel_h3_b2, p_g3_diel_h3_b1, p_g3_diel_h3_b3}; // hole 3
l_g3_diel_h3_b2 = newl; Line(l_g3_diel_h3_b2) = {p_g3_diel_h3_b2, p_g3_lc_h3_t2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h3_b3 = newl; Line(l_g3_diel_h3_b3) = {p_g3_diel_h3_b3, p_g3_lc_h3_t3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h4_b1 = newl; Circle(l_g3_diel_h4_b1) = {p_g3_diel_h4_b2, p_g3_diel_h4_b1, p_g3_diel_h4_b3}; // hole 4
l_g3_diel_h4_b2 = newl; Line(l_g3_diel_h4_b2) = {p_g3_diel_h4_b2, p_g3_lc_h4_t2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h4_b3 = newl; Line(l_g3_diel_h4_b3) = {p_g3_diel_h4_b3, p_g3_lc_h4_t3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_b1 = newl; Circle(l_g3_diel_h5_b1) = {p_g3_diel_h5_b2, p_g3_diel_h5_b1, p_g3_diel_h5_b3}; // hole 5
l_g3_diel_h5_b2 = newl; Circle(l_g3_diel_h5_b2) = {p_g3_diel_h5_b3, p_g3_diel_h5_b1, p_g3_diel_h5_b4};
l_g3_diel_h5_b3 = newl; Circle(l_g3_diel_h5_b3) = {p_g3_diel_h5_b4, p_g3_diel_h5_b1, p_g3_diel_h5_b5};
l_g3_diel_h5_b4 = newl; Circle(l_g3_diel_h5_b4) = {p_g3_diel_h5_b5, p_g3_diel_h5_b1, p_g3_diel_h5_b2};
l_g3_diel_h5_b5 = newl; Line(l_g3_diel_h5_b5) = {p_g3_diel_h5_b2, p_g3_lc_h5_t2}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_b6 = newl; Line(l_g3_diel_h5_b6) = {p_g3_diel_h5_b3, p_g3_lc_h5_t3}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_b7 = newl; Line(l_g3_diel_h5_b7) = {p_g3_diel_h5_b4, p_g3_lc_h5_t4}; // this line is on the surface exposed to gas by etching
l_g3_diel_h5_b8 = newl; Line(l_g3_diel_h5_b8) = {p_g3_diel_h5_b5, p_g3_lc_h5_t5}; // this line is on the surface exposed to gas by etching

/*---------------------------------- Curve Loop and Surface Definitions -------------------------------------*/
// top plane of upper copper plate
cl_g3_uc_t1 = newll; Curve Loop(cl_g3_uc_t1) = {l_g3_uc_h1_t1, l_g3_uc_h12_t1, l_g3_uc_h2_t1, l_g3_uc_h23_t1,
    l_g3_uc_h3_t1, l_g3_uc_h34_t1, l_g3_uc_h4_t1, l_g3_uc_h41_t1};
cl_g3_uc_t2 = newll; Curve Loop(cl_g3_uc_t2) = {l_g3_uc_h5_t1, l_g3_uc_h5_t2, l_g3_uc_h5_t3, l_g3_uc_h5_t4}; // curve loop for hole 5
// planes between top and bottom plane
cl_g3_uc_h1_tb1 = newll; Curve Loop(cl_g3_uc_h1_tb1) = {-l_g3_uc_h1_tb2, l_g3_uc_h1_t1, l_g3_uc_h1_tb3, -l_g3_uc_h1_b1}; // cylinderical surface on hole 1
cl_g3_uc_h12_tb1 = newll; Curve Loop(cl_g3_uc_h12_tb1) = {-l_g3_uc_h1_tb3, l_g3_uc_h12_t1, l_g3_uc_h2_tb2, -l_g3_uc_h12_b1}; // plane surface bwt hole 1 and hole  2
cl_g3_uc_h2_tb1 = newll; Curve Loop(cl_g3_uc_h2_tb1) = {-l_g3_uc_h2_tb2, l_g3_uc_h2_t1, l_g3_uc_h2_tb3, -l_g3_uc_h2_b1}; // cylinderical surface on hole 2
cl_g3_uc_h23_tb1 = newll; Curve Loop(cl_g3_uc_h23_tb1) = {-l_g3_uc_h2_tb3, l_g3_uc_h23_t1, l_g3_uc_h3_tb2, -l_g3_uc_h23_b1}; // plane surface bwt hole 2 and hole  3
cl_g3_uc_h3_tb1 = newll; Curve Loop(cl_g3_uc_h3_tb1) = {-l_g3_uc_h3_tb2, l_g3_uc_h3_t1, l_g3_uc_h3_tb3, -l_g3_uc_h3_b1}; // cylinderical surface on hole 3
cl_g3_uc_h34_tb1 = newll; Curve Loop(cl_g3_uc_h34_tb1) = {-l_g3_uc_h3_tb3, l_g3_uc_h34_t1, l_g3_uc_h4_tb2, -l_g3_uc_h34_b1}; // plane surface bwt hole 3 and hole  4
cl_g3_uc_h4_tb1 = newll; Curve Loop(cl_g3_uc_h4_tb1) = {-l_g3_uc_h4_tb2, l_g3_uc_h4_t1, l_g3_uc_h4_tb3, -l_g3_uc_h4_b1}; // cylinderical surface on hole 4
cl_g3_uc_h41_tb1 = newll; Curve Loop(cl_g3_uc_h41_tb1) = {-l_g3_uc_h4_tb3, l_g3_uc_h41_t1, l_g3_uc_h1_tb2, -l_g3_uc_h41_b1}; // plane surface bwt hole 4 and hole  1
cl_g3_uc_h5_tb1 = newll; Curve Loop(cl_g3_uc_h5_tb1) = {-l_g3_uc_h5_tb2, l_g3_uc_h5_t1, l_g3_uc_h5_tb3, -l_g3_uc_h5_b1}; // cylinderical surface on hole 5
cl_g3_uc_h5_tb2 = newll; Curve Loop(cl_g3_uc_h5_tb2) = {-l_g3_uc_h5_tb3, l_g3_uc_h5_t2, l_g3_uc_h5_tb4, -l_g3_uc_h5_b2};
cl_g3_uc_h5_tb3 = newll; Curve Loop(cl_g3_uc_h5_tb3) = {-l_g3_uc_h5_tb4, l_g3_uc_h5_t3, l_g3_uc_h5_tb5, -l_g3_uc_h5_b3};
cl_g3_uc_h5_tb4 = newll; Curve Loop(cl_g3_uc_h5_tb4) = {-l_g3_uc_h5_tb5, l_g3_uc_h5_t4, l_g3_uc_h5_tb2, -l_g3_uc_h5_b4};
// bottom plane of upper copper plate
cl_g3_uc_b1 = newll; Curve Loop(cl_g3_uc_b1) = {l_g3_uc_h1_b1, l_g3_uc_h12_b1, l_g3_uc_h2_b1, l_g3_uc_h23_b1,
    l_g3_uc_h3_b1, l_g3_uc_h34_b1, l_g3_uc_h4_b1, l_g3_uc_h41_b1};
cl_g3_uc_b2 = newll; Curve Loop(cl_g3_uc_b2) = {l_g3_uc_h5_b1, l_g3_uc_h5_b2, l_g3_uc_h5_b3, l_g3_uc_h5_b4}; // curve loop for hole 5
// surface definitions of upper plate
s_g3_uc_t1 = news; Plane Surface(s_g3_uc_t1) = {cl_g3_uc_t1, cl_g3_uc_t2};
s_g3_uc_h1_tb1 = news; Surface(s_g3_uc_h1_tb1) = {cl_g3_uc_h1_tb1}; // cylinderical surface on hole 1
s_g3_uc_h12_tb1 = news; Plane Surface(s_g3_uc_h12_tb1) = {cl_g3_uc_h12_tb1}; // plane surface bwt hole 1 and hole  2
s_g3_uc_h2_tb1 = news; Surface(s_g3_uc_h2_tb1) = {cl_g3_uc_h2_tb1}; // cylinderical surface on hole 2
s_g3_uc_h23_tb1 = news; Plane Surface(s_g3_uc_h23_tb1) = {cl_g3_uc_h23_tb1}; // plane surface bwt hole 2 and hole  3
s_g3_uc_h3_tb1 = news; Surface(s_g3_uc_h3_tb1) = {cl_g3_uc_h3_tb1}; // cylinderical surface on hole 3
s_g3_uc_h34_tb1 = news; Plane Surface(s_g3_uc_h34_tb1) = {cl_g3_uc_h34_tb1}; // plane surface bwt hole 3 and hole  4
s_g3_uc_h4_tb1 = news; Surface(s_g3_uc_h4_tb1) = {cl_g3_uc_h4_tb1}; // cylinderical surface on hole 4
s_g3_uc_h41_tb1 = news; Plane Surface(s_g3_uc_h41_tb1) = {cl_g3_uc_h41_tb1}; // plane surface bwt hole 4 and hole  1
s_g3_uc_h5_tb1 = news; Surface(s_g3_uc_h5_tb1) = {cl_g3_uc_h5_tb1}; // cylinderical surface on hole 5
s_g3_uc_h5_tb2 = news; Surface(s_g3_uc_h5_tb2) = {cl_g3_uc_h5_tb2};
s_g3_uc_h5_tb3 = news; Surface(s_g3_uc_h5_tb3) = {cl_g3_uc_h5_tb3};
s_g3_uc_h5_tb4 = news; Surface(s_g3_uc_h5_tb4) = {cl_g3_uc_h5_tb4};
s_g3_uc_b1 = news; Plane Surface(s_g3_uc_b1) = {cl_g3_uc_b1, cl_g3_uc_b2};

// curve loop definitions of upper copper plate
// top plane of upper copper plate
cl_g3_lc_t1 = newll; Curve Loop(cl_g3_lc_t1) = {l_g3_lc_h1_t1, l_g3_lc_h12_t1, l_g3_lc_h2_t1, l_g3_lc_h23_t1,
    l_g3_lc_h3_t1, l_g3_lc_h34_t1, l_g3_lc_h4_t1, l_g3_lc_h41_t1};
cl_g3_lc_t2 = newll; Curve Loop(cl_g3_lc_t2) = {l_g3_lc_h5_t1, l_g3_lc_h5_t2, l_g3_lc_h5_t3, l_g3_lc_h5_t4}; // curve loop for hole 5
// planes between top and bottom plane
cl_g3_lc_h1_tb1 = newll; Curve Loop(cl_g3_lc_h1_tb1) = {-l_g3_lc_h1_tb2, l_g3_lc_h1_t1, l_g3_lc_h1_tb3, -l_g3_lc_h1_b1}; // cylinderical surface on hole 1
cl_g3_lc_h12_tb1 = newll; Curve Loop(cl_g3_lc_h12_tb1) = {-l_g3_lc_h1_tb3, l_g3_lc_h12_t1, l_g3_lc_h2_tb2, -l_g3_lc_h12_b1}; // plane surface bwt hole 1 and hole  2
cl_g3_lc_h2_tb1 = newll; Curve Loop(cl_g3_lc_h2_tb1) = {-l_g3_lc_h2_tb2, l_g3_lc_h2_t1, l_g3_lc_h2_tb3, -l_g3_lc_h2_b1}; // cylinderical surface on hole 2
cl_g3_lc_h23_tb1 = newll; Curve Loop(cl_g3_lc_h23_tb1) = {-l_g3_lc_h2_tb3, l_g3_lc_h23_t1, l_g3_lc_h3_tb2, -l_g3_lc_h23_b1}; // plane surface bwt hole 2 and hole  3
cl_g3_lc_h3_tb1 = newll; Curve Loop(cl_g3_lc_h3_tb1) = {-l_g3_lc_h3_tb2, l_g3_lc_h3_t1, l_g3_lc_h3_tb3, -l_g3_lc_h3_b1}; // cylinderical surface on hole 3
cl_g3_lc_h34_tb1 = newll; Curve Loop(cl_g3_lc_h34_tb1) = {-l_g3_lc_h3_tb3, l_g3_lc_h34_t1, l_g3_lc_h4_tb2, -l_g3_lc_h34_b1}; // plane surface bwt hole 3 and hole  4
cl_g3_lc_h4_tb1 = newll; Curve Loop(cl_g3_lc_h4_tb1) = {-l_g3_lc_h4_tb2, l_g3_lc_h4_t1, l_g3_lc_h4_tb3, -l_g3_lc_h4_b1}; // cylinderical surface on hole 4
cl_g3_lc_h41_tb1 = newll; Curve Loop(cl_g3_lc_h41_tb1) = {-l_g3_lc_h4_tb3, l_g3_lc_h41_t1, l_g3_lc_h1_tb2, -l_g3_lc_h41_b1}; // plane surface bwt hole 4 and hole  1
cl_g3_lc_h5_tb1 = newll; Curve Loop(cl_g3_lc_h5_tb1) = {-l_g3_lc_h5_tb2, l_g3_lc_h5_t1, l_g3_lc_h5_tb3, -l_g3_lc_h5_b1}; // cylinderical surface on hole 5
cl_g3_lc_h5_tb2 = newll; Curve Loop(cl_g3_lc_h5_tb2) = {-l_g3_lc_h5_tb3, l_g3_lc_h5_t2, l_g3_lc_h5_tb4, -l_g3_lc_h5_b2};
cl_g3_lc_h5_tb3 = newll; Curve Loop(cl_g3_lc_h5_tb3) = {-l_g3_lc_h5_tb4, l_g3_lc_h5_t3, l_g3_lc_h5_tb5, -l_g3_lc_h5_b3};
cl_g3_lc_h5_tb4 = newll; Curve Loop(cl_g3_lc_h5_tb4) = {-l_g3_lc_h5_tb5, l_g3_lc_h5_t4, l_g3_lc_h5_tb2, -l_g3_lc_h5_b4};
// bottom plane of upper copper plate
cl_g3_lc_b1 = newll; Curve Loop(cl_g3_lc_b1) = {l_g3_lc_h1_b1, l_g3_lc_h12_b1, l_g3_lc_h2_b1, l_g3_lc_h23_b1,
    l_g3_lc_h3_b1, l_g3_lc_h34_b1, l_g3_lc_h4_b1, l_g3_lc_h41_b1};
cl_g3_lc_b2 = newll; Curve Loop(cl_g3_lc_b2) = {l_g3_lc_h5_b1, l_g3_lc_h5_b2, l_g3_lc_h5_b3, l_g3_lc_h5_b4}; // curve loop for hole 5
// surface definitions of upper plate
s_g3_lc_t1 = news; Plane Surface(s_g3_lc_t1) = {cl_g3_lc_t1, cl_g3_lc_t2};
s_g3_lc_h1_tb1 = news; Surface(s_g3_lc_h1_tb1) = {cl_g3_lc_h1_tb1}; // cylinderical surface on hole 1
s_g3_lc_h12_tb1 = news; Plane Surface(s_g3_lc_h12_tb1) = {cl_g3_lc_h12_tb1}; // plane surface bwt hole 1 and hole  2
s_g3_lc_h2_tb1 = news; Surface(s_g3_lc_h2_tb1) = {cl_g3_lc_h2_tb1}; // cylinderical surface on hole 2
s_g3_lc_h23_tb1 = news; Plane Surface(s_g3_lc_h23_tb1) = {cl_g3_lc_h23_tb1}; // plane surface bwt hole 2 and hole  3
s_g3_lc_h3_tb1 = news; Surface(s_g3_lc_h3_tb1) = {cl_g3_lc_h3_tb1}; // cylinderical surface on hole 3
s_g3_lc_h34_tb1 = news; Plane Surface(s_g3_lc_h34_tb1) = {cl_g3_lc_h34_tb1}; // plane surface bwt hole 3 and hole  4
s_g3_lc_h4_tb1 = news; Surface(s_g3_lc_h4_tb1) = {cl_g3_lc_h4_tb1}; // cylinderical surface on hole 4
s_g3_lc_h41_tb1 = news; Plane Surface(s_g3_lc_h41_tb1) = {cl_g3_lc_h41_tb1}; // plane surface bwt hole 4 and hole  1
s_g3_lc_h5_tb1 = news; Surface(s_g3_lc_h5_tb1) = {cl_g3_lc_h5_tb1}; // cylinderical surface on hole 5
s_g3_lc_h5_tb2 = news; Surface(s_g3_lc_h5_tb2) = {cl_g3_lc_h5_tb2};
s_g3_lc_h5_tb3 = news; Surface(s_g3_lc_h5_tb3) = {cl_g3_lc_h5_tb3};
s_g3_lc_h5_tb4 = news; Surface(s_g3_lc_h5_tb4) = {cl_g3_lc_h5_tb4};
s_g3_lc_b1 = news; Plane Surface(s_g3_lc_b1) = {cl_g3_lc_b1, cl_g3_lc_b2};

// curve loop definitions of dielectric layer
// top plane surfaces
// a surface exposed to gas by etching
cl_g3_diel_h1_t1 = newll; Curve Loop(cl_g3_diel_h1_t1) = {-l_g3_diel_h1_t2, l_g3_diel_h1_t1, l_g3_diel_h1_t3, -l_g3_uc_h1_b1}; 
cl_g3_diel_h2_t1 = newll; Curve Loop(cl_g3_diel_h2_t1) = {-l_g3_diel_h2_t2, l_g3_diel_h2_t1, l_g3_diel_h2_t3, -l_g3_uc_h2_b1};
cl_g3_diel_h3_t1 = newll; Curve Loop(cl_g3_diel_h3_t1) = {-l_g3_diel_h3_t2, l_g3_diel_h3_t1, l_g3_diel_h3_t3, -l_g3_uc_h3_b1};
cl_g3_diel_h4_t1 = newll; Curve Loop(cl_g3_diel_h4_t1) = {-l_g3_diel_h4_t2, l_g3_diel_h4_t1, l_g3_diel_h4_t3, -l_g3_uc_h4_b1};
cl_g3_diel_h5_t1 = newll; Curve Loop(cl_g3_diel_h5_t1) = {-l_g3_diel_h5_t5, l_g3_diel_h5_t1, l_g3_diel_h5_t6, -l_g3_uc_h5_b1};
cl_g3_diel_h5_t2 = newll; Curve Loop(cl_g3_diel_h5_t2) = {-l_g3_diel_h5_t6, l_g3_diel_h5_t2, l_g3_diel_h5_t7, -l_g3_uc_h5_b2};
cl_g3_diel_h5_t3 = newll; Curve Loop(cl_g3_diel_h5_t3) = {-l_g3_diel_h5_t7, l_g3_diel_h5_t3, l_g3_diel_h5_t8, -l_g3_uc_h5_b3};
cl_g3_diel_h5_t4 = newll; Curve Loop(cl_g3_diel_h5_t4) = {-l_g3_diel_h5_t8, l_g3_diel_h5_t4, l_g3_diel_h5_t5, -l_g3_uc_h5_b4};
// planes between top and middle plane
cl_g3_diel_h1_tm1 = newll; Curve Loop(cl_g3_diel_h1_tm1) = {l_g3_diel_h1_tm2, l_g3_diel_h1_m1, -l_g3_diel_h1_tm3, -l_g3_diel_h1_t1};
cl_g3_diel_h1_tm1 = newll; Curve Loop(cl_g3_diel_h1_tm1) = {l_g3_diel_h1_tm2, l_g3_diel_h1_m1, -l_g3_diel_h1_tm3, -l_g3_diel_h1_t1};
cl_g3_diel_h2_tm1 = newll; Curve Loop(cl_g3_diel_h2_tm1) = {l_g3_diel_h2_tm2, l_g3_diel_h2_m1, -l_g3_diel_h2_tm3, -l_g3_diel_h2_t1};
cl_g3_diel_h3_tm1 = newll; Curve Loop(cl_g3_diel_h3_tm1) = {l_g3_diel_h3_tm2, l_g3_diel_h3_m1, -l_g3_diel_h3_tm3, -l_g3_diel_h3_t1};
cl_g3_diel_h4_tm1 = newll; Curve Loop(cl_g3_diel_h4_tm1) = {l_g3_diel_h4_tm2, l_g3_diel_h4_m1, -l_g3_diel_h4_tm3, -l_g3_diel_h4_t1};
cl_g3_diel_h5_tm1 = newll; Curve Loop(cl_g3_diel_h5_tm1) = {l_g3_diel_h5_tm2, l_g3_diel_h5_m1, -l_g3_diel_h5_tm3, -l_g3_diel_h5_t1};
cl_g3_diel_h5_tm2 = newll; Curve Loop(cl_g3_diel_h5_tm2) = {l_g3_diel_h5_tm3, l_g3_diel_h5_m2, -l_g3_diel_h5_tm4, -l_g3_diel_h5_t2};
cl_g3_diel_h5_tm3 = newll; Curve Loop(cl_g3_diel_h5_tm3) = {l_g3_diel_h5_tm4, l_g3_diel_h5_m3, -l_g3_diel_h5_tm5, -l_g3_diel_h5_t3};
cl_g3_diel_h5_tm4 = newll; Curve Loop(cl_g3_diel_h5_tm4) = {l_g3_diel_h5_tm5, l_g3_diel_h5_m4, -l_g3_diel_h5_tm2, -l_g3_diel_h5_t4};
// planes between middle and bottom plane
cl_g3_diel_h1_mb1 = newll; Curve Loop(cl_g3_diel_h1_mb1) = {l_g3_diel_h1_mb2, l_g3_diel_h1_b1, -l_g3_diel_h1_mb3, -l_g3_diel_h1_m1};
cl_g3_diel_h2_mb1 = newll; Curve Loop(cl_g3_diel_h2_mb1) = {l_g3_diel_h2_mb2, l_g3_diel_h2_b1, -l_g3_diel_h2_mb3, -l_g3_diel_h2_m1};
cl_g3_diel_h3_mb1 = newll; Curve Loop(cl_g3_diel_h3_mb1) = {l_g3_diel_h3_mb2, l_g3_diel_h3_b1, -l_g3_diel_h3_mb3, -l_g3_diel_h3_m1};
cl_g3_diel_h4_mb1 = newll; Curve Loop(cl_g3_diel_h4_mb1) = {l_g3_diel_h4_mb2, l_g3_diel_h4_b1, -l_g3_diel_h4_mb3, -l_g3_diel_h4_m1};
cl_g3_diel_h5_mb1 = newll; Curve Loop(cl_g3_diel_h5_mb1) = {l_g3_diel_h5_mb2, l_g3_diel_h5_b1, -l_g3_diel_h5_mb3, -l_g3_diel_h5_m1};
cl_g3_diel_h5_mb2 = newll; Curve Loop(cl_g3_diel_h5_mb2) = {l_g3_diel_h5_mb3, l_g3_diel_h5_b2, -l_g3_diel_h5_mb4, -l_g3_diel_h5_m2};
cl_g3_diel_h5_mb3 = newll; Curve Loop(cl_g3_diel_h5_mb3) = {l_g3_diel_h5_mb4, l_g3_diel_h5_b3, -l_g3_diel_h5_mb5, -l_g3_diel_h5_m3};
cl_g3_diel_h5_mb4 = newll; Curve Loop(cl_g3_diel_h5_mb4) = {l_g3_diel_h5_mb5, l_g3_diel_h5_b4, -l_g3_diel_h5_mb2, -l_g3_diel_h5_m4};
// planes between top and bottom plane
cl_g3_diel_h12_tb1 = newll; Curve Loop(cl_g3_diel_h12_tb1) = {-l_g3_lc_h12_t1, -l_g3_diel_h1_b3, -l_g3_diel_h1_mb3, -l_g3_diel_h1_tm3, l_g3_diel_h1_t3,
    l_g3_uc_h12_b1, -l_g3_diel_h2_t2, l_g3_diel_h2_tm2, l_g3_diel_h2_mb2, l_g3_diel_h2_b2};
cl_g3_diel_h23_tb1 = newll; Curve Loop(cl_g3_diel_h23_tb1) = {-l_g3_lc_h23_t1, -l_g3_diel_h2_b3, -l_g3_diel_h2_mb3, -l_g3_diel_h2_tm3, l_g3_diel_h2_t3,
    l_g3_uc_h23_b1, -l_g3_diel_h3_t2, l_g3_diel_h3_tm2, l_g3_diel_h3_mb2, l_g3_diel_h3_b2};
cl_g3_diel_h34_tb1 = newll; Curve Loop(cl_g3_diel_h34_tb1) = {-l_g3_lc_h34_t1, -l_g3_diel_h3_b3, -l_g3_diel_h3_mb3, -l_g3_diel_h3_tm3, l_g3_diel_h3_t3,
    l_g3_uc_h34_b1, -l_g3_diel_h4_t2, l_g3_diel_h4_tm2, l_g3_diel_h4_mb2, l_g3_diel_h4_b2};
cl_g3_diel_h41_tb1 = newll; Curve Loop(cl_g3_diel_h41_tb1) = {-l_g3_lc_h41_t1, -l_g3_diel_h4_b3, -l_g3_diel_h4_mb3, -l_g3_diel_h4_tm3, l_g3_diel_h4_t3,
    l_g3_uc_h41_b1, -l_g3_diel_h1_t2, l_g3_diel_h1_tm2, l_g3_diel_h1_mb2, l_g3_diel_h1_b2};
// bottom plane surfaces
// a surface exposed to gas by etching
cl_g3_diel_h1_b1 = newll; Curve Loop(cl_g3_diel_h1_b1) = {-l_g3_diel_h1_b2, l_g3_diel_h1_b1, l_g3_diel_h1_b3, -l_g3_lc_h1_t1}; 
cl_g3_diel_h2_b1 = newll; Curve Loop(cl_g3_diel_h2_b1) = {-l_g3_diel_h2_b2, l_g3_diel_h2_b1, l_g3_diel_h2_b3, -l_g3_lc_h2_t1};
cl_g3_diel_h3_b1 = newll; Curve Loop(cl_g3_diel_h3_b1) = {-l_g3_diel_h3_b2, l_g3_diel_h3_b1, l_g3_diel_h3_b3, -l_g3_lc_h3_t1};
cl_g3_diel_h4_b1 = newll; Curve Loop(cl_g3_diel_h4_b1) = {-l_g3_diel_h4_b2, l_g3_diel_h4_b1, l_g3_diel_h4_b3, -l_g3_lc_h4_t1};
cl_g3_diel_h5_b1 = newll; Curve Loop(cl_g3_diel_h5_b1) = {-l_g3_diel_h5_b5, l_g3_diel_h5_b1, l_g3_diel_h5_b6, -l_g3_lc_h5_t1};
cl_g3_diel_h5_b2 = newll; Curve Loop(cl_g3_diel_h5_b2) = {-l_g3_diel_h5_b6, l_g3_diel_h5_b2, l_g3_diel_h5_b7, -l_g3_lc_h5_t2};
cl_g3_diel_h5_b3 = newll; Curve Loop(cl_g3_diel_h5_b3) = {-l_g3_diel_h5_b7, l_g3_diel_h5_b3, l_g3_diel_h5_b8, -l_g3_lc_h5_t3};
cl_g3_diel_h5_b4 = newll; Curve Loop(cl_g3_diel_h5_b4) = {-l_g3_diel_h5_b8, l_g3_diel_h5_b4, l_g3_diel_h5_b5, -l_g3_lc_h5_t4};

// surface definitions of dielectric layer
s_g3_diel_h1_t1 = news; Plane Surface(s_g3_diel_h1_t1) = {cl_g3_diel_h1_t1};
s_g3_diel_h2_t1 = news; Plane Surface(s_g3_diel_h2_t1) = {cl_g3_diel_h2_t1};
s_g3_diel_h3_t1 = news; Plane Surface(s_g3_diel_h3_t1) = {cl_g3_diel_h3_t1};
s_g3_diel_h4_t1 = news; Plane Surface(s_g3_diel_h4_t1) = {cl_g3_diel_h4_t1};
s_g3_diel_h5_t1 = news; Plane Surface(s_g3_diel_h5_t1) = {cl_g3_diel_h5_t1};
s_g3_diel_h5_t2 = news; Plane Surface(s_g3_diel_h5_t2) = {cl_g3_diel_h5_t2};
s_g3_diel_h5_t3 = news; Plane Surface(s_g3_diel_h5_t3) = {cl_g3_diel_h5_t3};
s_g3_diel_h5_t4 = news; Plane Surface(s_g3_diel_h5_t4) = {cl_g3_diel_h5_t4};
s_g3_diel_h1_tm1 = news; Surface(s_g3_diel_h1_tm1) = {cl_g3_diel_h1_tm1};
s_g3_diel_h1_tm1 = news; Surface(s_g3_diel_h1_tm1) = {cl_g3_diel_h1_tm1};
s_g3_diel_h2_tm1 = news; Surface(s_g3_diel_h2_tm1) = {cl_g3_diel_h2_tm1};
s_g3_diel_h3_tm1 = news; Surface(s_g3_diel_h3_tm1) = {cl_g3_diel_h3_tm1};
s_g3_diel_h4_tm1 = news; Surface(s_g3_diel_h4_tm1) = {cl_g3_diel_h4_tm1};
s_g3_diel_h5_tm1 = news; Surface(s_g3_diel_h5_tm1) = {cl_g3_diel_h5_tm1};
s_g3_diel_h5_tm2 = news; Surface(s_g3_diel_h5_tm2) = {cl_g3_diel_h5_tm2};
s_g3_diel_h5_tm3 = news; Surface(s_g3_diel_h5_tm3) = {cl_g3_diel_h5_tm3};
s_g3_diel_h5_tm4 = news; Surface(s_g3_diel_h5_tm4) = {cl_g3_diel_h5_tm4};
s_g3_diel_h1_mb1 = news; Surface(s_g3_diel_h1_mb1) = {cl_g3_diel_h1_mb1};
s_g3_diel_h2_mb1 = news; Surface(s_g3_diel_h2_mb1) = {cl_g3_diel_h2_mb1};
s_g3_diel_h3_mb1 = news; Surface(s_g3_diel_h3_mb1) = {cl_g3_diel_h3_mb1};
s_g3_diel_h4_mb1 = news; Surface(s_g3_diel_h4_mb1) = {cl_g3_diel_h4_mb1};
s_g3_diel_h5_mb1 = news; Surface(s_g3_diel_h5_mb1) = {cl_g3_diel_h5_mb1};
s_g3_diel_h5_mb2 = news; Surface(s_g3_diel_h5_mb2) = {cl_g3_diel_h5_mb2};
s_g3_diel_h5_mb3 = news; Surface(s_g3_diel_h5_mb3) = {cl_g3_diel_h5_mb3};
s_g3_diel_h5_mb4 = news; Surface(s_g3_diel_h5_mb4) = {cl_g3_diel_h5_mb4};
s_g3_diel_h12_tb1 = news; Plane Surface(s_g3_diel_h12_tb1) = {cl_g3_diel_h12_tb1};
s_g3_diel_h23_tb1 = news; Plane Surface(s_g3_diel_h23_tb1) = {cl_g3_diel_h23_tb1};
s_g3_diel_h34_tb1 = news; Plane Surface(s_g3_diel_h34_tb1) = {cl_g3_diel_h34_tb1};
s_g3_diel_h41_tb1 = news; Plane Surface(s_g3_diel_h41_tb1) = {cl_g3_diel_h41_tb1};
s_g3_diel_h1_b1 = news; Plane Surface(s_g3_diel_h1_b1) = {cl_g3_diel_h1_b1};
s_g3_diel_h2_b1 = news; Plane Surface(s_g3_diel_h2_b1) = {cl_g3_diel_h2_b1};
s_g3_diel_h3_b1 = news; Plane Surface(s_g3_diel_h3_b1) = {cl_g3_diel_h3_b1};
s_g3_diel_h4_b1 = news; Plane Surface(s_g3_diel_h4_b1) = {cl_g3_diel_h4_b1};
s_g3_diel_h5_b1 = news; Plane Surface(s_g3_diel_h5_b1) = {cl_g3_diel_h5_b1};
s_g3_diel_h5_b2 = news; Plane Surface(s_g3_diel_h5_b2) = {cl_g3_diel_h5_b2};
s_g3_diel_h5_b3 = news; Plane Surface(s_g3_diel_h5_b3) = {cl_g3_diel_h5_b3};
s_g3_diel_h5_b4 = news; Plane Surface(s_g3_diel_h5_b4) = {cl_g3_diel_h5_b4};
/*---------------------------------- Surface Loop and Volume Definitions -------------------------------------*/
// surface loop definitions of upper plate
sl_g3_uc = newsl; Surface Loop(sl_g3_uc) = {s_g3_uc_t1,
    s_g3_uc_h1_tb1, s_g3_uc_h2_tb1, s_g3_uc_h3_tb1, s_g3_uc_h4_tb1,
    s_g3_uc_h12_tb1, s_g3_uc_h23_tb1, s_g3_uc_h34_tb1, s_g3_uc_h41_tb1,
    s_g3_uc_h5_tb1, s_g3_uc_h5_tb2, s_g3_uc_h5_tb3, s_g3_uc_h5_tb4,
    s_g3_uc_b1};
// surface loop definitions of upper plate
sl_g3_lc = newsl; Surface Loop(sl_g3_lc) = {s_g3_lc_t1,
    s_g3_lc_h1_tb1, s_g3_lc_h2_tb1, s_g3_lc_h3_tb1, s_g3_lc_h4_tb1,
    s_g3_lc_h12_tb1, s_g3_lc_h23_tb1, s_g3_lc_h34_tb1, s_g3_lc_h41_tb1,
    s_g3_lc_h5_tb1, s_g3_lc_h5_tb2, s_g3_lc_h5_tb3, s_g3_lc_h5_tb4,
    s_g3_lc_b1};
// surface loop definitions of dielectric layer
sl_g3_diel = newsl; Surface Loop(sl_g3_diel) = {
    s_g3_uc_b1,
    s_g3_diel_h1_t1, s_g3_diel_h2_t1, s_g3_diel_h3_t1, s_g3_diel_h4_t1,
    s_g3_diel_h5_t1, s_g3_diel_h5_t2, s_g3_diel_h5_t3, s_g3_diel_h5_t4,
    s_g3_diel_h1_tm1, s_g3_diel_h2_tm1, s_g3_diel_h3_tm1, s_g3_diel_h4_tm1,
    s_g3_diel_h5_tm1, s_g3_diel_h5_tm2, s_g3_diel_h5_tm3, s_g3_diel_h5_tm4,
    s_g3_diel_h1_mb1, s_g3_diel_h2_mb1, s_g3_diel_h3_mb1, s_g3_diel_h4_mb1,
    s_g3_diel_h5_mb1, s_g3_diel_h5_mb2, s_g3_diel_h5_mb3, s_g3_diel_h5_mb4,
    s_g3_diel_h12_tb1, s_g3_diel_h23_tb1, s_g3_diel_h34_tb1, s_g3_diel_h41_tb1,
    s_g3_diel_h1_b1, s_g3_diel_h2_b1, s_g3_diel_h3_b1, s_g3_diel_h4_b1,
    s_g3_diel_h5_b1, s_g3_diel_h5_b2, s_g3_diel_h5_b3, s_g3_diel_h5_b4,
    s_g3_lc_t1
};
v_g3_uc = newv; Volume(v_g3_uc) = {sl_g3_uc};
v_g3_lc = newv; Volume(v_g3_lc) = {sl_g3_lc};
v_g3_diel = newv; Volume(v_g3_diel) = {sl_g3_diel};
// ---------------------------------------------------------------------------------------------------------------
// Gas
// ---------------------------------------------------------------------------------------------------------------
/*---------------------------------- Point Definitions -------------------------------------*/
// upper boundary plane
p_u1 = newp; Point(p_u1) = {-bdrHalfX, -bdrHalfY,   zU, mshSizeDrift};
p_u2 = newp; Point(p_u2) = {bdrHalfX, -bdrHalfY,    zU, mshSizeDrift};
p_u3 = newp; Point(p_u3) = {bdrHalfX, bdrHalfY,     zU, mshSizeDrift};
p_u4 = newp; Point(p_u4) = {-bdrHalfX, bdrHalfY,    zU, mshSizeDrift};
// pad plane
p_p1 = newp; Point(p_p1) = {-bdrHalfX, -bdrHalfY,   zP, mshSizeDrift};
p_p2 = newp; Point(p_p2) = {bdrHalfX, -bdrHalfY,    zP, mshSizeDrift};
p_p3 = newp; Point(p_p3) = {bdrHalfX, bdrHalfY,     zP, mshSizeDrift};
p_p4 = newp; Point(p_p4) = {-bdrHalfX, bdrHalfY,    zP, mshSizeDrift};
// point between plates and upper boudnary, pad
// this point is just for setting mesh size to prevent over refinement in drifing regions.
// between upper elecgrode and GEM2
p_ug1_1 = newp; Point(p_ug1_1) = {-bdrHalfX, -bdrHalfY, (zU + z1)/2, mshSizeDrift};
p_ug1_2 = newp; Point(p_ug1_2) = {bdrHalfX, -bdrHalfY,  (zU + z1)/2, mshSizeDrift};
p_ug1_3 = newp; Point(p_ug1_3) = {bdrHalfX, bdrHalfY,   (zU + z1)/2, mshSizeDrift};
p_ug1_4 = newp; Point(p_ug1_4) = {-bdrHalfX, bdrHalfY,  (zU + z1)/2, mshSizeDrift};
// between GEM2 and GEM1
p_g1g2_1 = newp; Point(p_g1g2_1) = {-bdrHalfX, -bdrHalfY, (z1 + z2)/2, mshSizeDrift};
p_g1g2_2 = newp; Point(p_g1g2_2) = {bdrHalfX, -bdrHalfY,  (z1 + z2)/2, mshSizeDrift};
p_g1g2_3 = newp; Point(p_g1g2_3) = {bdrHalfX, bdrHalfY,   (z1 + z2)/2, mshSizeDrift};
p_g1g2_4 = newp; Point(p_g1g2_4) = {-bdrHalfX, bdrHalfY,  (z1 + z2)/2, mshSizeDrift};
// between GEM1 and GEM3
p_g2g3_1 = newp; Point(p_g2g3_1) = {-bdrHalfX, -bdrHalfY, (z2 + z3)/2, mshSizeDrift};
p_g2g3_2 = newp; Point(p_g2g3_2) = {bdrHalfX, -bdrHalfY,  (z2 + z3)/2, mshSizeDrift};
p_g2g3_3 = newp; Point(p_g2g3_3) = {bdrHalfX, bdrHalfY,   (z2 + z3)/2, mshSizeDrift};
p_g2g3_4 = newp; Point(p_g2g3_4) = {-bdrHalfX, bdrHalfY,  (z2 + z3)/2, mshSizeDrift};
// between GEM3 and pad plane
p_g3p_1 = newp; Point(p_g3p_1) = {-bdrHalfX, -bdrHalfY, (z3 + zP)/2, mshSizeDrift};
p_g3p_2 = newp; Point(p_g3p_2) = {bdrHalfX, -bdrHalfY,  (z3 + zP)/2, mshSizeDrift};
p_g3p_3 = newp; Point(p_g3p_3) = {bdrHalfX, bdrHalfY,   (z3 + zP)/2, mshSizeDrift};
p_g3p_4 = newp; Point(p_g3p_4) = {-bdrHalfX, bdrHalfY,  (z3 + zP)/2, mshSizeDrift};
/*---------------------------------- Line Definitions -------------------------------------*/
// upper boundary
l_u1 = newl; Line(l_u1) = {p_u1, p_u2};
l_u2 = newl; Line(l_u2) = {p_u2, p_u3};
l_u3 = newl; Line(l_u3) = {p_u3, p_u4};
l_u4 = newl; Line(l_u4) = {p_u4, p_u1};
// between external eletrode and GEM1
l_ug1_u1 = newl; Line(l_ug1_u1) = {p_u1, p_ug1_1};
l_ug1_u2 = newl; Line(l_ug1_u2) = {p_u2, p_ug1_2};
l_ug1_u3 = newl; Line(l_ug1_u3) = {p_u3, p_ug1_3};
l_ug1_u4 = newl; Line(l_ug1_u4) = {p_u4, p_ug1_4};
l_ug1_l1 = newl; Line(l_ug1_l1) = {p_ug1_1, p_g1_uc_h1_t1};
l_ug1_l2 = newl; Line(l_ug1_l2) = {p_ug1_2, p_g1_uc_h2_t1};
l_ug1_l3 = newl; Line(l_ug1_l3) = {p_ug1_3, p_g1_uc_h3_t1};
l_ug1_l4 = newl; Line(l_ug1_l4) = {p_ug1_4, p_g1_uc_h4_t1};
// between GEM1 and GEM2
l_g1g2_u1 = newl; Line(l_g1g2_u1) = {p_g1_lc_h1_b1, p_g1g2_1};
l_g1g2_u2 = newl; Line(l_g1g2_u2) = {p_g1_lc_h2_b1, p_g1g2_2};
l_g1g2_u3 = newl; Line(l_g1g2_u3) = {p_g1_lc_h3_b1, p_g1g2_3};
l_g1g2_u4 = newl; Line(l_g1g2_u4) = {p_g1_lc_h4_b1, p_g1g2_4};
l_g1g2_l1 = newl; Line(l_g1g2_l1) = {p_g1g2_1, p_g2_uc_h1_t1};
l_g1g2_l2 = newl; Line(l_g1g2_l2) = {p_g1g2_2, p_g2_uc_h2_t1};
l_g1g2_l3 = newl; Line(l_g1g2_l3) = {p_g1g2_3, p_g2_uc_h3_t1};
l_g1g2_l4 = newl; Line(l_g1g2_l4) = {p_g1g2_4, p_g2_uc_h4_t1};
// between GEM2 and GEM3
l_g2g3_u1 = newl; Line(l_g2g3_u1) = {p_g2_lc_h1_b1, p_g2g3_1};
l_g2g3_u2 = newl; Line(l_g2g3_u2) = {p_g2_lc_h2_b1, p_g2g3_2};
l_g2g3_u3 = newl; Line(l_g2g3_u3) = {p_g2_lc_h3_b1, p_g2g3_3};
l_g2g3_u4 = newl; Line(l_g2g3_u4) = {p_g2_lc_h4_b1, p_g2g3_4};
l_g2g3_l1 = newl; Line(l_g2g3_l1) = {p_g2g3_1, p_g3_uc_h1_t1};
l_g2g3_l2 = newl; Line(l_g2g3_l2) = {p_g2g3_2, p_g3_uc_h2_t1};
l_g2g3_l3 = newl; Line(l_g2g3_l3) = {p_g2g3_3, p_g3_uc_h3_t1};
l_g2g3_l4 = newl; Line(l_g2g3_l4) = {p_g2g3_4, p_g3_uc_h4_t1};
// between GEM3 and pad plane
l_g3p_u1 = newl; Line(l_g3p_u1) = {p_g3_lc_h1_b1, p_g3p_1};
l_g3p_u2 = newl; Line(l_g3p_u2) = {p_g3_lc_h2_b1, p_g3p_2};
l_g3p_u3 = newl; Line(l_g3p_u3) = {p_g3_lc_h3_b1, p_g3p_3};
l_g3p_u4 = newl; Line(l_g3p_u4) = {p_g3_lc_h4_b1, p_g3p_4};
l_g3p_l1 = newl; Line(l_g3p_l1) = {p_g3p_1, p_p1};
l_g3p_l2 = newl; Line(l_g3p_l2) = {p_g3p_2, p_p2};
l_g3p_l3 = newl; Line(l_g3p_l3) = {p_g3p_3, p_p3};
l_g3p_l4 = newl; Line(l_g3p_l4) = {p_g3p_4, p_p4};
// pad plane
l_p1 = newl; Line(l_p1) = {p_p1, p_p2};
l_p2 = newl; Line(l_p2) = {p_p2, p_p3};
l_p3 = newl; Line(l_p3) = {p_p3, p_p4};
l_p4 = newl; Line(l_p4) = {p_p4, p_p1};
/*---------------------------------- Curve Loop and Surface Definitions -------------------------------------*/
// upper boundary plane
cl_u = newll; Curve Loop(cl_u) = {l_u1, l_u2, l_u3, l_u4};
// pad plane
cl_p = newll; Curve Loop(cl_p) = {l_p1, l_p2, l_p3, l_p4};
// side planes of gas volume
cl_up1 = newll; Curve Loop(cl_up1) = {
    l_u1,
    l_ug1_u2, l_ug1_l2, l_g1_uc_h2_tb1, l_g1_diel_h2_tm1, l_g1_diel_h2_mb1, l_g1_lc_h2_tb1,
    l_g1g2_u2, l_g1g2_l2, l_g2_uc_h2_tb1, l_g2_diel_h2_tm1, l_g2_diel_h2_mb1, l_g2_lc_h2_tb1,
    l_g2g3_u2, l_g2g3_l2, l_g3_uc_h2_tb1, l_g3_diel_h2_tm1, l_g3_diel_h2_mb1, l_g3_lc_h2_tb1,
    l_g3p_u2, l_g3p_l2,
    -l_g3p_l1, -l_g3p_u1, -l_g3_lc_h1_tb1, -l_g3_diel_h1_mb1, -l_g3_diel_h1_tm1, -l_g3_uc_h1_tb1,
    -l_g2g3_l1, -l_g2g3_u1, -l_g2_lc_h1_tb1, -l_g2_diel_h1_mb1, -l_g2_diel_h1_tm1, -l_g2_uc_h1_tb1,
    -l_g1g2_l1, -l_g1g2_u1, -l_g1_lc_h1_tb1, -l_g1_diel_h1_mb1, -l_g1_diel_h1_tm1, -l_g1_uc_h1_tb1,
    -l_ug1_l1, -l_ug1_u1,
    -l_p1};
cl_up2 = newll; Curve Loop(cl_up2) = {
    l_u2,
    l_ug1_u3, l_ug1_l3, l_g1_uc_h3_tb1, l_g1_diel_h3_tm1, l_g1_diel_h3_mb1, l_g1_lc_h3_tb1,
    l_g1g2_u3, l_g1g2_l3, l_g2_uc_h3_tb1, l_g2_diel_h3_tm1, l_g2_diel_h3_mb1, l_g2_lc_h3_tb1,
    l_g2g3_u3, l_g2g3_l3, l_g3_uc_h3_tb1, l_g3_diel_h3_tm1, l_g3_diel_h3_mb1, l_g3_lc_h3_tb1,
    l_g3p_u3, l_g3p_l3,
    -l_g3p_l2, -l_g3p_u2, -l_g3_lc_h2_tb1, -l_g3_diel_h2_mb1, -l_g3_diel_h2_tm1, -l_g3_uc_h2_tb1,
    -l_g2g3_l2, -l_g2g3_u2, -l_g2_lc_h2_tb1, -l_g2_diel_h2_mb1, -l_g2_diel_h2_tm1, -l_g2_uc_h2_tb1,
    -l_g1g2_l2, -l_g1g2_u2, -l_g1_lc_h2_tb1, -l_g1_diel_h2_mb1, -l_g1_diel_h2_tm1, -l_g1_uc_h2_tb1,
    -l_ug1_l2, -l_ug1_u2,
    -l_p2};
cl_up3 = newll; Curve Loop(cl_up3) = {
    l_u3,
    l_ug1_u4, l_ug1_l4, l_g1_uc_h4_tb1, l_g1_diel_h4_tm1, l_g1_diel_h4_mb1, l_g1_lc_h4_tb1,
    l_g1g2_u4, l_g1g2_l4, l_g2_uc_h4_tb1, l_g2_diel_h4_tm1, l_g2_diel_h4_mb1, l_g2_lc_h4_tb1,
    l_g2g3_u4, l_g2g3_l4, l_g3_uc_h4_tb1, l_g3_diel_h4_tm1, l_g3_diel_h4_mb1, l_g3_lc_h4_tb1,
    l_g3p_u4, l_g3p_l4,
    -l_g3p_l3, -l_g3p_u3, -l_g3_lc_h3_tb1, -l_g3_diel_h3_mb1, -l_g3_diel_h3_tm1, -l_g3_uc_h3_tb1,
    -l_g2g3_l3, -l_g2g3_u3, -l_g2_lc_h3_tb1, -l_g2_diel_h3_mb1, -l_g2_diel_h3_tm1, -l_g2_uc_h3_tb1,
    -l_g1g2_l3, -l_g1g2_u3, -l_g1_lc_h3_tb1, -l_g1_diel_h3_mb1, -l_g1_diel_h3_tm1, -l_g1_uc_h3_tb1,
    -l_ug1_l3, -l_ug1_u3,
    -l_p3};
cl_up4 = newll; Curve Loop(cl_up4) = {
    l_u4,
    l_ug1_u1, l_ug1_l1, l_g1_uc_h1_tb1, l_g1_diel_h1_tm1, l_g1_diel_h1_mb1, l_g1_lc_h1_tb1,
    l_g1g2_u1, l_g1g2_l1, l_g2_uc_h1_tb1, l_g2_diel_h1_tm1, l_g2_diel_h1_mb1, l_g2_lc_h1_tb1,
    l_g2g3_u1, l_g2g3_l1, l_g3_uc_h1_tb1, l_g3_diel_h1_tm1, l_g3_diel_h1_mb1, l_g3_lc_h1_tb1,
    l_g3p_u1, l_g3p_l1,
    -l_g3p_l4, -l_g3p_u4, -l_g3_lc_h4_tb1, -l_g3_diel_h4_mb1, -l_g3_diel_h4_tm1, -l_g3_uc_h4_tb1,
    -l_g2g3_l4, -l_g2g3_u4, -l_g2_lc_h4_tb1, -l_g2_diel_h4_mb1, -l_g2_diel_h4_tm1, -l_g2_uc_h4_tb1,
    -l_g1g2_l4, -l_g1g2_u4, -l_g1_lc_h4_tb1, -l_g1_diel_h4_mb1, -l_g1_diel_h4_tm1, -l_g1_uc_h4_tb1,
    -l_ug1_l4, -l_ug1_u4,
    -l_p4};
// hole in surface by GEM1
cl_up1_g1 = newll;Curve Loop(cl_up1_g1) = {-l_g1_lc_h1_tb3,
    -l_g1_diel_h1_b3, -l_g1_diel_h1_mb3, -l_g1_diel_h1_tm3,
    l_g1_diel_h1_t3, -l_g1_uc_h1_tb3,
    l_g1_uc_h12_t1,
    l_g1_uc_h2_tb2, -l_g1_diel_h2_t2,
    l_g1_diel_h2_tm2, l_g1_diel_h2_mb2, l_g1_diel_h2_b2,
    l_g1_lc_h2_tb2,
    -l_g1_lc_h12_b1};
cl_up2_g1 = newll;Curve Loop(cl_up2_g1) = {-l_g1_lc_h2_tb3,
    -l_g1_diel_h2_b3, -l_g1_diel_h2_mb3, -l_g1_diel_h2_tm3,
    l_g1_diel_h2_t3, -l_g1_uc_h2_tb3,
    l_g1_uc_h23_t1,
    l_g1_uc_h3_tb2, -l_g1_diel_h3_t2,
    l_g1_diel_h3_tm2, l_g1_diel_h3_mb2, l_g1_diel_h3_b2,
    l_g1_lc_h3_tb2,
    -l_g1_lc_h23_b1};
cl_up3_g1 = newll;Curve Loop(cl_up3_g1) = {-l_g1_lc_h3_tb3,
    -l_g1_diel_h3_b3, -l_g1_diel_h3_mb3, -l_g1_diel_h3_tm3,
    l_g1_diel_h3_t3, -l_g1_uc_h3_tb3,
    l_g1_uc_h34_t1,
    l_g1_uc_h4_tb2, -l_g1_diel_h4_t2,
    l_g1_diel_h4_tm2, l_g1_diel_h4_mb2, l_g1_diel_h4_b2,
    l_g1_lc_h4_tb2,
    -l_g1_lc_h34_b1};
cl_up4_g1 = newll;Curve Loop(cl_up4_g1) = {-l_g1_lc_h4_tb3,
    -l_g1_diel_h4_b3, -l_g1_diel_h4_mb3, -l_g1_diel_h4_tm3,
    l_g1_diel_h4_t3, -l_g1_uc_h4_tb3,
    l_g1_uc_h41_t1,
    l_g1_uc_h1_tb2, -l_g1_diel_h1_t2,
    l_g1_diel_h1_tm2, l_g1_diel_h1_mb2, l_g1_diel_h1_b2,
    l_g1_lc_h1_tb2,
    -l_g1_lc_h41_b1};
// hole in surface by GEM2
cl_up1_g2 = newll;Curve Loop(cl_up1_g2) = {-l_g2_lc_h1_tb3,
    -l_g2_diel_h1_b3, -l_g2_diel_h1_mb3, -l_g2_diel_h1_tm3,
    l_g2_diel_h1_t3, -l_g2_uc_h1_tb3,
    l_g2_uc_h12_t1,
    l_g2_uc_h2_tb2, -l_g2_diel_h2_t2,
    l_g2_diel_h2_tm2, l_g2_diel_h2_mb2, l_g2_diel_h2_b2,
    l_g2_lc_h2_tb2,
    -l_g2_lc_h12_b1};
cl_up2_g2 = newll;Curve Loop(cl_up2_g2) = {-l_g2_lc_h2_tb3,
    -l_g2_diel_h2_b3, -l_g2_diel_h2_mb3, -l_g2_diel_h2_tm3,
    l_g2_diel_h2_t3, -l_g2_uc_h2_tb3,
    l_g2_uc_h23_t1,
    l_g2_uc_h3_tb2, -l_g2_diel_h3_t2,
    l_g2_diel_h3_tm2, l_g2_diel_h3_mb2, l_g2_diel_h3_b2,
    l_g2_lc_h3_tb2,
    -l_g2_lc_h23_b1};
cl_up3_g2 = newll;Curve Loop(cl_up3_g2) = {-l_g2_lc_h3_tb3,
    -l_g2_diel_h3_b3, -l_g2_diel_h3_mb3, -l_g2_diel_h3_tm3,
    l_g2_diel_h3_t3, -l_g2_uc_h3_tb3,
    l_g2_uc_h34_t1,
    l_g2_uc_h4_tb2, -l_g2_diel_h4_t2,
    l_g2_diel_h4_tm2, l_g2_diel_h4_mb2, l_g2_diel_h4_b2,
    l_g2_lc_h4_tb2,
    -l_g2_lc_h34_b1};
cl_up4_g2 = newll;Curve Loop(cl_up4_g2) = {-l_g2_lc_h4_tb3,
    -l_g2_diel_h4_b3, -l_g2_diel_h4_mb3, -l_g2_diel_h4_tm3,
    l_g2_diel_h4_t3, -l_g2_uc_h4_tb3,
    l_g2_uc_h41_t1,
    l_g2_uc_h1_tb2, -l_g2_diel_h1_t2,
    l_g2_diel_h1_tm2, l_g2_diel_h1_mb2, l_g2_diel_h1_b2,
    l_g2_lc_h1_tb2,
    -l_g2_lc_h41_b1};
// hole in surface by GEM2
cl_up1_g3 = newll;Curve Loop(cl_up1_g3) = {-l_g3_lc_h1_tb3,
    -l_g3_diel_h1_b3, -l_g3_diel_h1_mb3, -l_g3_diel_h1_tm3,
    l_g3_diel_h1_t3, -l_g3_uc_h1_tb3,
    l_g3_uc_h12_t1,
    l_g3_uc_h2_tb2, -l_g3_diel_h2_t2,
    l_g3_diel_h2_tm2, l_g3_diel_h2_mb2, l_g3_diel_h2_b2,
    l_g3_lc_h2_tb2,
    -l_g3_lc_h12_b1};
cl_up2_g3 = newll;Curve Loop(cl_up2_g3) = {-l_g3_lc_h2_tb3,
    -l_g3_diel_h2_b3, -l_g3_diel_h2_mb3, -l_g3_diel_h2_tm3,
    l_g3_diel_h2_t3, -l_g3_uc_h2_tb3,
    l_g3_uc_h23_t1,
    l_g3_uc_h3_tb2, -l_g3_diel_h3_t2,
    l_g3_diel_h3_tm2, l_g3_diel_h3_mb2, l_g3_diel_h3_b2,
    l_g3_lc_h3_tb2,
    -l_g3_lc_h23_b1};
cl_up3_g3 = newll;Curve Loop(cl_up3_g3) = {-l_g3_lc_h3_tb3,
    -l_g3_diel_h3_b3, -l_g3_diel_h3_mb3, -l_g3_diel_h3_tm3,
    l_g3_diel_h3_t3, -l_g3_uc_h3_tb3,
    l_g3_uc_h34_t1,
    l_g3_uc_h4_tb2, -l_g3_diel_h4_t2,
    l_g3_diel_h4_tm2, l_g3_diel_h4_mb2, l_g3_diel_h4_b2,
    l_g3_lc_h4_tb2,
    -l_g3_lc_h34_b1};
cl_up4_g3 = newll;Curve Loop(cl_up4_g3) = {-l_g3_lc_h4_tb3,
    -l_g3_diel_h4_b3, -l_g3_diel_h4_mb3, -l_g3_diel_h4_tm3,
    l_g3_diel_h4_t3, -l_g3_uc_h4_tb3,
    l_g3_uc_h41_t1,
    l_g3_uc_h1_tb2, -l_g3_diel_h1_t2,
    l_g3_diel_h1_tm2, l_g3_diel_h1_mb2, l_g3_diel_h1_b2,
    l_g3_lc_h1_tb2,
    -l_g3_lc_h41_b1};
s_u = news; Plane Surface(s_u) = {cl_u};
s_p = news; Plane Surface(s_p) = {cl_p};
s_up1 = news; Plane Surface(s_up1) = {cl_up1, cl_up1_g1, cl_up1_g2, cl_up1_g3};
s_up2 = news; Plane Surface(s_up2) = {cl_up2, cl_up2_g1, cl_up2_g2, cl_up2_g3};
s_up3 = news; Plane Surface(s_up3) = {cl_up3, cl_up3_g1, cl_up3_g2, cl_up3_g3};
s_up4 = news; Plane Surface(s_up4) = {cl_up4, cl_up4_g1, cl_up4_g2, cl_up4_g3};
/*---------------------------------- Surface Loop and Volume Definitions -------------------------------------*/
sl_up = newsl; Surface Loop(sl_up) = {
    s_u, s_p,
    s_up1, s_up2, s_up3, s_up4,
    s_g1_uc_t1,
    s_g1_uc_h1_tb1, s_g1_uc_h2_tb1, s_g1_uc_h3_tb1, s_g1_uc_h4_tb1,
    s_g1_uc_h5_tb1, s_g1_uc_h5_tb2, s_g1_uc_h5_tb3, s_g1_uc_h5_tb4,
    s_g1_lc_h1_tb1, s_g1_lc_h2_tb1, s_g1_lc_h3_tb1, s_g1_lc_h4_tb1,
    s_g1_lc_h5_tb1, s_g1_lc_h5_tb2, s_g1_lc_h5_tb3, s_g1_lc_h5_tb4,
    s_g1_lc_b1,
    s_g1_diel_h1_t1, s_g1_diel_h2_t1, s_g1_diel_h3_t1, s_g1_diel_h4_t1,
    s_g1_diel_h5_t1, s_g1_diel_h5_t2, s_g1_diel_h5_t3, s_g1_diel_h5_t4,
    s_g1_diel_h1_tm1, s_g1_diel_h2_tm1, s_g1_diel_h3_tm1, s_g1_diel_h4_tm1,
    s_g1_diel_h5_tm1, s_g1_diel_h5_tm2, s_g1_diel_h5_tm3, s_g1_diel_h5_tm4,
    s_g1_diel_h1_mb1, s_g1_diel_h2_mb1, s_g1_diel_h3_mb1, s_g1_diel_h4_mb1,
    s_g1_diel_h5_mb1, s_g1_diel_h5_mb2, s_g1_diel_h5_mb3, s_g1_diel_h5_mb4,
    s_g1_diel_h1_b1, s_g1_diel_h2_b1, s_g1_diel_h3_b1, s_g1_diel_h4_b1,
    s_g1_diel_h5_b1, s_g1_diel_h5_b2, s_g1_diel_h5_b3, s_g1_diel_h5_b4,
    s_g2_uc_t1,
    s_g2_uc_h1_tb1, s_g2_uc_h2_tb1, s_g2_uc_h3_tb1, s_g2_uc_h4_tb1,
    s_g2_uc_h5_tb1, s_g2_uc_h5_tb2, s_g2_uc_h5_tb3, s_g2_uc_h5_tb4,
    s_g2_lc_h1_tb1, s_g2_lc_h2_tb1, s_g2_lc_h3_tb1, s_g2_lc_h4_tb1,
    s_g2_lc_h5_tb1, s_g2_lc_h5_tb2, s_g2_lc_h5_tb3, s_g2_lc_h5_tb4,
    s_g2_lc_b1,
    s_g2_diel_h1_t1, s_g2_diel_h2_t1, s_g2_diel_h3_t1, s_g2_diel_h4_t1,
    s_g2_diel_h5_t1, s_g2_diel_h5_t2, s_g2_diel_h5_t3, s_g2_diel_h5_t4,
    s_g2_diel_h1_tm1, s_g2_diel_h2_tm1, s_g2_diel_h3_tm1, s_g2_diel_h4_tm1,
    s_g2_diel_h5_tm1, s_g2_diel_h5_tm2, s_g2_diel_h5_tm3, s_g2_diel_h5_tm4,
    s_g2_diel_h1_mb1, s_g2_diel_h2_mb1, s_g2_diel_h3_mb1, s_g2_diel_h4_mb1,
    s_g2_diel_h5_mb1, s_g2_diel_h5_mb2, s_g2_diel_h5_mb3, s_g2_diel_h5_mb4,
    s_g2_diel_h1_b1, s_g2_diel_h2_b1, s_g2_diel_h3_b1, s_g2_diel_h4_b1,
    s_g2_diel_h5_b1, s_g2_diel_h5_b2, s_g2_diel_h5_b3, s_g2_diel_h5_b4,
    s_g3_uc_t1,
    s_g3_uc_h1_tb1, s_g3_uc_h2_tb1, s_g3_uc_h3_tb1, s_g3_uc_h4_tb1,
    s_g3_uc_h5_tb1, s_g3_uc_h5_tb2, s_g3_uc_h5_tb3, s_g3_uc_h5_tb4,
    s_g3_lc_h1_tb1, s_g3_lc_h2_tb1, s_g3_lc_h3_tb1, s_g3_lc_h4_tb1,
    s_g3_lc_h5_tb1, s_g3_lc_h5_tb2, s_g3_lc_h5_tb3, s_g3_lc_h5_tb4,
    s_g3_lc_b1,
    s_g3_diel_h1_t1, s_g3_diel_h2_t1, s_g3_diel_h3_t1, s_g3_diel_h4_t1,
    s_g3_diel_h5_t1, s_g3_diel_h5_t2, s_g3_diel_h5_t3, s_g3_diel_h5_t4,
    s_g3_diel_h1_tm1, s_g3_diel_h2_tm1, s_g3_diel_h3_tm1, s_g3_diel_h4_tm1,
    s_g3_diel_h5_tm1, s_g3_diel_h5_tm2, s_g3_diel_h5_tm3, s_g3_diel_h5_tm4,
    s_g3_diel_h1_mb1, s_g3_diel_h2_mb1, s_g3_diel_h3_mb1, s_g3_diel_h4_mb1,
    s_g3_diel_h5_mb1, s_g3_diel_h5_mb2, s_g3_diel_h5_mb3, s_g3_diel_h5_mb4,
    s_g3_diel_h1_b1, s_g3_diel_h2_b1, s_g3_diel_h3_b1, s_g3_diel_h4_b1,
    s_g3_diel_h5_b1, s_g3_diel_h5_b2, s_g3_diel_h5_b3, s_g3_diel_h5_b4};

v_up = newv; Volume(v_up) = {sl_up};
// ---------------------------------------------------------------------------------------------------------------
// Physical volumes and surfaces
// ---------------------------------------------------------------------------------------------------------------
/*---------------------------------- Physical volume definitions -------------------------------------*/
// gas
pv_up = newv; Physical Volume(pv_up) = {v_up}; // Body 1
// upper, lower copper and dielectric layers of GEM1, 2, 3
pv_g1_uc = newv; Physical Volume(pv_g1_uc) = {v_g1_uc}; // Body 2
pv_g1_lc = newv; Physical Volume(pv_g1_lc) = {v_g1_lc}; // Body 3
pv_g1_diel = newv; Physical Volume(pv_g1_diel) = {v_g1_diel}; // Body 4
pv_g2_uc = newv; Physical Volume(pv_g2_uc) = {v_g2_uc}; // Body 5
pv_g2_lc = newv; Physical Volume(pv_g2_lc) = {v_g2_lc}; // Body 6
pv_g2_diel = newv; Physical Volume(pv_g2_diel) = {v_g2_diel}; // Body 7
pv_g3_uc = newv; Physical Volume(pv_g3_uc) = {v_g3_uc}; // Body 8
pv_g3_lc = newv; Physical Volume(pv_g3_lc) = {v_g3_lc}; // Body 9
pv_g3_diel = newv; Physical Volume(pv_g3_diel) = {v_g3_diel}; // Body 10
/*---------------------------------- Physical surface definitions for boundary conditions -------------------------------------*/
// upper boundary plane
// Boundary 1
ps_u = news; Physical Surface(ps_u) = {s_u};
// Boundary 2
ps_g1_uc = news; Physical Surface(ps_g1_uc) = {s_g1_uc_t1,
    s_g1_uc_h1_tb1, s_g1_uc_h2_tb1, s_g1_uc_h3_tb1, s_g1_uc_h4_tb1,
    s_g1_uc_h12_tb1, s_g1_uc_h23_tb1, s_g1_uc_h34_tb1, s_g1_uc_h41_tb1,
    s_g1_uc_h5_tb1, s_g1_uc_h5_tb2, s_g1_uc_h5_tb3, s_g1_uc_h5_tb4,
    s_g1_uc_b1};
// Boundary 3
ps_g1_lc = news; Physical Surface(ps_g1_lc) = {s_g1_lc_t1,
    s_g1_lc_h1_tb1, s_g1_lc_h2_tb1, s_g1_lc_h3_tb1, s_g1_lc_h4_tb1,
    s_g1_lc_h12_tb1, s_g1_lc_h23_tb1, s_g1_lc_h34_tb1, s_g1_lc_h41_tb1,
    s_g1_lc_h5_tb1, s_g1_lc_h5_tb2, s_g1_lc_h5_tb3, s_g1_lc_h5_tb4,
    s_g1_lc_b1};
// Boundary 4
ps_g2_uc = news; Physical Surface(ps_g2_uc) = {s_g2_uc_t1,
    s_g2_uc_h1_tb1, s_g2_uc_h2_tb1, s_g2_uc_h3_tb1, s_g2_uc_h4_tb1,
    s_g2_uc_h12_tb1, s_g2_uc_h23_tb1, s_g2_uc_h34_tb1, s_g2_uc_h41_tb1,
    s_g2_uc_h5_tb1, s_g2_uc_h5_tb2, s_g2_uc_h5_tb3, s_g2_uc_h5_tb4,
    s_g2_uc_b1};
// Boundary 5
ps_g2_lc = news; Physical Surface(ps_g2_lc) = {s_g2_lc_t1,
    s_g2_lc_h1_tb1, s_g2_lc_h2_tb1, s_g2_lc_h3_tb1, s_g2_lc_h4_tb1,
    s_g2_lc_h12_tb1, s_g2_lc_h23_tb1, s_g2_lc_h34_tb1, s_g2_lc_h41_tb1,
    s_g2_lc_h5_tb1, s_g2_lc_h5_tb2, s_g2_lc_h5_tb3, s_g2_lc_h5_tb4,
    s_g2_lc_b1};
// Boundary 6
ps_g3_uc = news; Physical Surface(ps_g3_uc) = {s_g3_uc_t1,
    s_g3_uc_h1_tb1, s_g3_uc_h2_tb1, s_g3_uc_h3_tb1, s_g3_uc_h4_tb1,
    s_g3_uc_h12_tb1, s_g3_uc_h23_tb1, s_g3_uc_h34_tb1, s_g3_uc_h41_tb1,
    s_g3_uc_h5_tb1, s_g3_uc_h5_tb2, s_g3_uc_h5_tb3, s_g3_uc_h5_tb4,
    s_g3_uc_b1};
// Boundary 7
ps_g3_lc = news; Physical Surface(ps_g3_lc) = {s_g3_lc_t1,
    s_g3_lc_h1_tb1, s_g3_lc_h2_tb1, s_g3_lc_h3_tb1, s_g3_lc_h4_tb1,
    s_g3_lc_h12_tb1, s_g3_lc_h23_tb1, s_g3_lc_h34_tb1, s_g3_lc_h41_tb1,
    s_g3_lc_h5_tb1, s_g3_lc_h5_tb2, s_g3_lc_h5_tb3, s_g3_lc_h5_tb4,
    s_g3_lc_b1};
// Boundary 8
// pad plane
ps_p = news; Physical Surface(ps_p) = {s_p}; 

// surfaces for periodic boundary conditions
// Boundary 9
ps_per_xmin = news; Physical Surface(ps_per_xmin) = {
    s_up1,
    s_g1_uc_h12_tb1, s_g1_diel_h12_tb1, s_g1_lc_h12_tb1,
    s_g2_uc_h12_tb1, s_g2_diel_h12_tb1, s_g2_lc_h12_tb1,
    s_g3_uc_h12_tb1, s_g3_diel_h12_tb1, s_g3_lc_h12_tb1};  
// Boundary 10
ps_per_xmax = news; Physical Surface(ps_per_xmax) = {
    s_up3,
    s_g1_uc_h34_tb1, s_g1_diel_h34_tb1, s_g1_lc_h34_tb1,
    s_g2_uc_h34_tb1, s_g2_diel_h34_tb1, s_g2_lc_h34_tb1,
    s_g3_uc_h34_tb1, s_g3_diel_h34_tb1, s_g3_lc_h34_tb1};
// Boundary 11
ps_per_ymin = news; Physical Surface(ps_per_ymin) = {
    s_up4,
    s_g1_uc_h41_tb1, s_g1_diel_h41_tb1, s_g1_lc_h41_tb1,
    s_g2_uc_h41_tb1, s_g2_diel_h41_tb1, s_g2_lc_h41_tb1,
    s_g3_uc_h41_tb1, s_g3_diel_h41_tb1, s_g3_lc_h41_tb1};
// Boundary 12
ps_per_ymax = news; Physical Surface(ps_per_ymax) = {
    s_up2,
    s_g1_uc_h23_tb1, s_g1_diel_h23_tb1, s_g1_lc_h23_tb1,
    s_g2_uc_h23_tb1, s_g2_diel_h23_tb1, s_g2_lc_h23_tb1,
    s_g3_uc_h23_tb1, s_g3_diel_h23_tb1, s_g3_lc_h23_tb1};  