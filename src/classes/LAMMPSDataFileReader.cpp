#include <iostream>
#include <algorithm> // for isspace
#include <sstream>
#include "LAMMPSDataFileReader.h"
#include "error_code_defines.h"

LAMMPSDataFileReader::LAMMPSDataFileReader(std::string file) : _file(file) {
  _initKeywords();
  readFile();
}

void LAMMPSDataFileReader::_initKeywords() {
  _header_keywords.push_back("atoms"); // indicates # of atoms (2)
  _header_keywords.push_back("bonds"); // indicates # of bonds (2)
  _header_keywords.push_back("angles"); // indicates # of angles (2)
  _header_keywords.push_back("dihedrals"); // indicates # of dihedrals (2)
  _header_keywords.push_back("impropers"); // indicates # of impropers (2)
  _header_keywords.push_back("atom types"); // indicates # of atom types (3)
  _header_keywords.push_back("bond types"); // indicates # of bond types (3)
  _header_keywords.push_back("angle types"); // indicates # of angle types (3)
  _header_keywords.push_back("dihedral types"); // indicates # of dihedral types (3)
  _header_keywords.push_back("improper types"); // indicates # of improper types (3)
  _header_keywords.push_back("ellipsoids"); // indicates # of ellipsoids (2)
  _header_keywords.push_back("lines"); // indicates # of line segments (2)
  _header_keywords.push_back("triangles"); // indicates # of triangles (2)
  _header_keywords.push_back("bodies"); // indicates # of bodies (2)
  _header_keywords.push_back("xlo xhi"); // indicates box boundaries in x (4)
  _header_keywords.push_back("ylo yhi"); // indicates box boundaries in y (4)
  _header_keywords.push_back("zlo zhi"); // indicates box boundaries in z (4)
  _header_keywords.push_back("xy xz yz"); // indicates tilt factors for box (6)

  // Atom property sections
  _section_keywords.push_back("Atoms"); // looks for style comment
  _section_keywords.push_back("Velocities");
  _section_keywords.push_back("Masses");
  _section_keywords.push_back("Ellipsoids");
  _section_keywords.push_back("Lines");
  _section_keywords.push_back("Triangles");
  _section_keywords.push_back("Bodies");
  // Molecular topology sections
  _section_keywords.push_back("Bonds");
  _section_keywords.push_back("Angles");
  _section_keywords.push_back("Dihedrals");
  _section_keywords.push_back("Impropers");
  // Force field sections
  _section_keywords.push_back("Pair Coeffs"); // looks for style comment
  _section_keywords.push_back("PairIJ Coeffs"); // looks for style comment
  _section_keywords.push_back("Bond Coeffs"); // looks for style comment
  _section_keywords.push_back("Angle Coeffs"); // looks for style comment
  _section_keywords.push_back("Dihedral Coeffs"); // looks for style comment
  _section_keywords.push_back("Improper Coeffs"); // looks for style comment
  // class 2 force field sections
  _section_keywords.push_back("BondBond Coeffs");
  _section_keywords.push_back("BondAngle Coeffs");
  _section_keywords.push_back("MiddleBondTorsion Coeffs");
  _section_keywords.push_back("EndBondTorsion Coeffs");
  _section_keywords.push_back("AngleTorsion Coeffs");
  _section_keywords.push_back("AngleAngleTorsion Coeffs");
  _section_keywords.push_back("BondBond13 Coeffs");
  _section_keywords.push_back("AngleAngle Coeffs");

  _atom_style_keywords.push_back("angle");
  _atom_style_keywords.push_back("atomic");
  _atom_style_keywords.push_back("body");
  _atom_style_keywords.push_back("bond");
  _atom_style_keywords.push_back("charge");
  _atom_style_keywords.push_back("dipole");
  _atom_style_keywords.push_back("dpd");
  _atom_style_keywords.push_back("edpd");
  _atom_style_keywords.push_back("electron");
  _atom_style_keywords.push_back("ellipsoid");
  _atom_style_keywords.push_back("full");
  _atom_style_keywords.push_back("line");
  _atom_style_keywords.push_back("mdpd");
  _atom_style_keywords.push_back("molecular");
  _atom_style_keywords.push_back("peri");
  _atom_style_keywords.push_back("smd");
  _atom_style_keywords.push_back("sph");
  _atom_style_keywords.push_back("sphere");
  _atom_style_keywords.push_back("spin");
  _atom_style_keywords.push_back("tdpd");
  _atom_style_keywords.push_back("tri");
  _atom_style_keywords.push_back("template");
  _atom_style_keywords.push_back("hybrid");

  _pair_style_keywords.push_back("none");
  _pair_style_keywords.push_back("hybrid");
  _pair_style_keywords.push_back("hybrid/overlay");
  _pair_style_keywords.push_back("zero");
  _pair_style_keywords.push_back("adp");
  _pair_style_keywords.push_back("agni");
  _pair_style_keywords.push_back("airebo");
  _pair_style_keywords.push_back("airebo/morse");
  _pair_style_keywords.push_back("atm");
  _pair_style_keywords.push_back("awpmd/cut");
  _pair_style_keywords.push_back("beck");
  _pair_style_keywords.push_back("body/nparticle");
  _pair_style_keywords.push_back("body/rounded/polygon");
  _pair_style_keywords.push_back("body/rounded/polyhedron");
  _pair_style_keywords.push_back("bop");
  _pair_style_keywords.push_back("born");
  _pair_style_keywords.push_back("born/coul/dsf");
  _pair_style_keywords.push_back("born/coul/dsf/cs");
  _pair_style_keywords.push_back("born/coul/long");
  _pair_style_keywords.push_back("born/coul/long/cs");
  _pair_style_keywords.push_back("born/coul/msm");
  _pair_style_keywords.push_back("born/coul/wolf");
  _pair_style_keywords.push_back("born/coul/wolf/cs");
  _pair_style_keywords.push_back("brownian");
  _pair_style_keywords.push_back("brownian/poly");
  _pair_style_keywords.push_back("buck");
  _pair_style_keywords.push_back("buck/coul/cut");
  _pair_style_keywords.push_back("buck/coul/long");
  _pair_style_keywords.push_back("buck/coul/long/cs");
  _pair_style_keywords.push_back("buck/coul/msm");
  _pair_style_keywords.push_back("buck/long/coul/long");
  _pair_style_keywords.push_back("buck/mdf");
  _pair_style_keywords.push_back("buck6d/coul/gauss/dsf");
  _pair_style_keywords.push_back("buck6d/coul/gauss/long");
  _pair_style_keywords.push_back("colloid");
  _pair_style_keywords.push_back("comb");
  _pair_style_keywords.push_back("comb3");
  _pair_style_keywords.push_back("cosine/squared");
  _pair_style_keywords.push_back("coul/cut");
  _pair_style_keywords.push_back("coult/cut/soft");
  _pair_style_keywords.push_back("coul/debye");
  _pair_style_keywords.push_back("coul/diel");
  _pair_style_keywords.push_back("coul/dsf");
  _pair_style_keywords.push_back("coul/long");
  _pair_style_keywords.push_back("coul/long/cs");
  _pair_style_keywords.push_back("coul/long/soft");
  _pair_style_keywords.push_back("coul/msm");
  _pair_style_keywords.push_back("coul/slater/cut");
  _pair_style_keywords.push_back("coul/slater/long");
  _pair_style_keywords.push_back("coul/shield");
  _pair_style_keywords.push_back("coul/streitz");
  _pair_style_keywords.push_back("coul/wolf");
  _pair_style_keywords.push_back("coul/wolf/cs");
  _pair_style_keywords.push_back("dpd");
  _pair_style_keywords.push_back("dpd/fdt");
  _pair_style_keywords.push_back("dpd/fdt/energy");
  _pair_style_keywords.push_back("dpd/tstat");
  _pair_style_keywords.push_back("dsmc");
  _pair_style_keywords.push_back("e3b");
  _pair_style_keywords.push_back("drip");
  _pair_style_keywords.push_back("eam");
  _pair_style_keywords.push_back("eam/alloy");
  _pair_style_keywords.push_back("eam/cd");
  _pair_style_keywords.push_back("eam/cd/old");
  _pair_style_keywords.push_back("eam/fs");
  _pair_style_keywords.push_back("edip");
  _pair_style_keywords.push_back("eff/cut");
  _pair_style_keywords.push_back("eim");
  _pair_style_keywords.push_back("exp6/rx");
  _pair_style_keywords.push_back("extep");
  _pair_style_keywords.push_back("gauss");
  _pair_style_keywords.push_back("gauss/cut");
  _pair_style_keywords.push_back("gayberne");
  _pair_style_keywords.push_back("granular");
  _pair_style_keywords.push_back("gran/hertz/history");
  _pair_style_keywords.push_back("gran/hooke");
  _pair_style_keywords.push_back("gran/hooke/history");
  _pair_style_keywords.push_back("gw");
  _pair_style_keywords.push_back("gw/zbl");
  _pair_style_keywords.push_back("hbond/dreiding/lj");
  _pair_style_keywords.push_back("hbond/dreiding/morse");
  _pair_style_keywords.push_back("ilp/graphene/hbn");
  _pair_style_keywords.push_back("kim");
  _pair_style_keywords.push_back("kolmogorov/crespi/full");
  _pair_style_keywords.push_back("kolmogorov/crespi/z");
  _pair_style_keywords.push_back("lcbop");
  _pair_style_keywords.push_back("lebedeva/z");
  _pair_style_keywords.push_back("lennard/mdf");
  _pair_style_keywords.push_back("line/lj");
  _pair_style_keywords.push_back("list");
  _pair_style_keywords.push_back("lj/charmm/coul/charmm");
  _pair_style_keywords.push_back("lj/charmm/coul/charmm/implicit");
  _pair_style_keywords.push_back("lj/charmm/coul/long");
  _pair_style_keywords.push_back("lj/charmm/coul/long/soft");
  _pair_style_keywords.push_back("lj/charmm/coul/msm");
  _pair_style_keywords.push_back("lj/charmmfsw/coul/charmmfsh");
  _pair_style_keywords.push_back("lj/charmmfsw/coul/long");
  _pair_style_keywords.push_back("lj/class2");
  _pair_style_keywords.push_back("lj/class2/coul/cut");
  _pair_style_keywords.push_back("lj/class2/coul/cut/soft");
  _pair_style_keywords.push_back("lj/class2/coul/long");
  _pair_style_keywords.push_back("lj/class2/coul/long/cs");
  _pair_style_keywords.push_back("lj/class2/coul/long/soft");
  _pair_style_keywords.push_back("lj/class2/soft");
  _pair_style_keywords.push_back("lj/cubic");
  _pair_style_keywords.push_back("lj/cut");
  _pair_style_keywords.push_back("lj/cut/coul/cut");
  _pair_style_keywords.push_back("lj/cut/coul/cut/soft");
  _pair_style_keywords.push_back("lj/cut/coul/debye");
  _pair_style_keywords.push_back("lj/cut/coul/dsf");
  _pair_style_keywords.push_back("lj/cut/coul/long");
  _pair_style_keywords.push_back("lj/cut/coul/long/cs");
  _pair_style_keywords.push_back("lj/cut/coul/long/soft");
  _pair_style_keywords.push_back("lj/cut/coul/msm");
  _pair_style_keywords.push_back("lj/cut/coul/wolf");
  _pair_style_keywords.push_back("lj/cut/dipole/cut");
  _pair_style_keywords.push_back("lj/cut/dipole/long");
  _pair_style_keywords.push_back("lj/cut/soft");
  _pair_style_keywords.push_back("lj/cut/thole/long");
  _pair_style_keywords.push_back("lj/cut/tip4p/cut");
  _pair_style_keywords.push_back("lj/cut/tip4p/long");
  _pair_style_keywords.push_back("lj/cut/tip4p/long/soft");
  _pair_style_keywords.push_back("lj/expand");
  _pair_style_keywords.push_back("lj/expand/coul/long");
  _pair_style_keywords.push_back("lj/gromacs");
  _pair_style_keywords.push_back("lj/gromacs/coul/gromacs");
  _pair_style_keywords.push_back("lj/long/coul/long");
  _pair_style_keywords.push_back("lj/long/dipole/long");
  _pair_style_keywords.push_back("lj/long/tip4p/long");
  _pair_style_keywords.push_back("lj/mdf");
  _pair_style_keywords.push_back("lj/sdk");
  _pair_style_keywords.push_back("lj/sdk/coul/long");
  _pair_style_keywords.push_back("lj/sdk/coul/msm");
  _pair_style_keywords.push_back("lj/sf/dipole/sf");
  _pair_style_keywords.push_back("lj/smooth");
  _pair_style_keywords.push_back("lj/smooth/linear");
  _pair_style_keywords.push_back("lj/switch3/coulgauss/long");
  _pair_style_keywords.push_back("lj96/cut");
  _pair_style_keywords.push_back("local/density");
  _pair_style_keywords.push_back("lubricate");
  _pair_style_keywords.push_back("lubricate/poly");
  _pair_style_keywords.push_back("lubricateU");
  _pair_style_keywords.push_back("lubricateU/poly");
  _pair_style_keywords.push_back("mdpd");
  _pair_style_keywords.push_back("mdpd/rhosum");
  _pair_style_keywords.push_back("meam/c");
  _pair_style_keywords.push_back("meam/spline");
  _pair_style_keywords.push_back("meam/sw/spline");
  _pair_style_keywords.push_back("mesocnt");
  _pair_style_keywords.push_back("mgpt");
  _pair_style_keywords.push_back("mesont/tpm");
  _pair_style_keywords.push_back("mie/cut");
  _pair_style_keywords.push_back("mliap");
  _pair_style_keywords.push_back("mm3/switch3/coulgass/long");
  _pair_style_keywords.push_back("momb");
  _pair_style_keywords.push_back("morse");
  _pair_style_keywords.push_back("morse/smooth/linear");
  _pair_style_keywords.push_back("morse/soft");
  _pair_style_keywords.push_back("multi/lucy");
  _pair_style_keywords.push_back("multi/lucy/rx");
  _pair_style_keywords.push_back("nb3b/harmonic");
  _pair_style_keywords.push_back("nm/cut");
  _pair_style_keywords.push_back("nm/cut/coul/cut");
  _pair_style_keywords.push_back("nm/cut/coul/long");
  _pair_style_keywords.push_back("oxdna/coaxstk");
  _pair_style_keywords.push_back("oxdna/excv");
  _pair_style_keywords.push_back("oxdna/hbond");
  _pair_style_keywords.push_back("oxdna/stk");
  _pair_style_keywords.push_back("oxdna/xstk");
  _pair_style_keywords.push_back("oxdna2/coaxstk");
  _pair_style_keywords.push_back("oxdna2/dh");
  _pair_style_keywords.push_back("oxdna2/excv");
  _pair_style_keywords.push_back("oxdna2/hbond");
  _pair_style_keywords.push_back("oxdna2/stk");
  _pair_style_keywords.push_back("oxdna2/xstk");
  _pair_style_keywords.push_back("oxrna2/coaxstk");
  _pair_style_keywords.push_back("oxrna2/dh");
  _pair_style_keywords.push_back("oxrna2/excv");
  _pair_style_keywords.push_back("oxrna2/hbond");
  _pair_style_keywords.push_back("oxrna2/stk");
  _pair_style_keywords.push_back("oxrna2/xstk");
  _pair_style_keywords.push_back("peri/eps");
  _pair_style_keywords.push_back("peri/lps");
  _pair_style_keywords.push_back("peri/pmb");
  _pair_style_keywords.push_back("peri/ves");
  _pair_style_keywords.push_back("polymorphic");
  _pair_style_keywords.push_back("python");
  _pair_style_keywords.push_back("quip");
  _pair_style_keywords.push_back("reax/c");
  _pair_style_keywords.push_back("rebo");
  _pair_style_keywords.push_back("resquared");
  _pair_style_keywords.push_back("sdpd/taitwaiter/isothermal");
  _pair_style_keywords.push_back("smd/hertz");
  _pair_style_keywords.push_back("smd/tlsph");
  _pair_style_keywords.push_back("smd/tri_surface");
  _pair_style_keywords.push_back("smd/ulsph");
  _pair_style_keywords.push_back("smtbq");
  _pair_style_keywords.push_back("snap");
  _pair_style_keywords.push_back("soft");
  _pair_style_keywords.push_back("sph/heatconduction");
  _pair_style_keywords.push_back("sph/ideal/gas");
  _pair_style_keywords.push_back("sph/lj");
  _pair_style_keywords.push_back("sph/rhosum");
  _pair_style_keywords.push_back("sph/taitwaiter");
  _pair_style_keywords.push_back("sph/taitwaiter/morris");
  _pair_style_keywords.push_back("spin/dipole/cut");
  _pair_style_keywords.push_back("spin/dipole/long");
  _pair_style_keywords.push_back("spin/dmi");
  _pair_style_keywords.push_back("spin/exchange");
  _pair_style_keywords.push_back("spin/magelec");
  _pair_style_keywords.push_back("spin/neel");
  _pair_style_keywords.push_back("srp");
  _pair_style_keywords.push_back("sw");
  _pair_style_keywords.push_back("table");
  _pair_style_keywords.push_back("table/rx");
  _pair_style_keywords.push_back("tdpd");
  _pair_style_keywords.push_back("tersoff");
  _pair_style_keywords.push_back("tersoff/mod");
  _pair_style_keywords.push_back("tersoff/mod/c");
  _pair_style_keywords.push_back("tersoff/table");
  _pair_style_keywords.push_back("tersoff/zbl");
  _pair_style_keywords.push_back("thole");
  _pair_style_keywords.push_back("tip4p/cut");
  _pair_style_keywords.push_back("tip4p/long");
  _pair_style_keywords.push_back("tip4p/long/soft");
  _pair_style_keywords.push_back("tri/lj");
  _pair_style_keywords.push_back("ufm");
  _pair_style_keywords.push_back("vashishta");
  _pair_style_keywords.push_back("vashishta/table");
  _pair_style_keywords.push_back("yukawa");
  _pair_style_keywords.push_back("yukawa/colloid");
  _pair_style_keywords.push_back("zbl");

  _bond_style_keywords.push_back("class2");
  _bond_style_keywords.push_back("fene");
  _bond_style_keywords.push_back("fene/expand");
  _bond_style_keywords.push_back("harmoic");
  _bond_style_keywords.push_back("hybrid");
  _bond_style_keywords.push_back("morse");
  _bond_style_keywords.push_back("none");
  _bond_style_keywords.push_back("nonlinear");
  _bond_style_keywords.push_back("quartic");

  _angle_style_keywords.push_back("charmm");
  _angle_style_keywords.push_back("class2");
  _angle_style_keywords.push_back("cosine");
  _angle_style_keywords.push_back("cosine/squared");
  _angle_style_keywords.push_back("harmonic");
  _angle_style_keywords.push_back("hybrid");
  _angle_style_keywords.push_back("none");

  _dihedral_style_keywords.push_back("charmm");
  _dihedral_style_keywords.push_back("class2");
  _dihedral_style_keywords.push_back("harmonic");
  _dihedral_style_keywords.push_back("helix");
  _dihedral_style_keywords.push_back("hybrid");
  _dihedral_style_keywords.push_back("multi/harmonic");
  _dihedral_style_keywords.push_back("none");
  _dihedral_style_keywords.push_back("opls");

  _improper_style_keywords.push_back("class2");
  _improper_style_keywords.push_back("cvff");
  _improper_style_keywords.push_back("harmonic");
  _improper_style_keywords.push_back("hybrid");
  _improper_style_keywords.push_back("none");
}

void LAMMPSDataFileReader::readFile() {
  std::string str; // holds each line
  int i1; // for holding the integer values
  double d1, d2, d3; // for holding the double values
  std::vector <int> vals (14, -1); // for holding the integer values in the header region
  std::ifstream fin(_file.c_str());

  if (fin.fail()) {
    std::cerr << "Error opening file \"" << _file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }

  getline(fin, _comment); // First line is always a comment line

  while (getline(fin, str)) {
    std::size_t first_char = str.find_first_not_of(" \n\t\r\f\v");
    if (first_char == std::string::npos || first_char == '#') {
      continue; // ignore comments and empty lines
    }

    // check for header keywords
    for (std::size_t i = 0; i < _header_keywords.size(); ++i) {
      if (str.find(_header_keywords[i]) != std::string::npos) {
        std::stringstream ss(str);
        switch(i) {
          case 14: // xlo xhi
            ss >> d1 >> d2;
            _box.setXLow(d1); _box.setXHigh(d2);
          case 15: // ylo yhi
            ss >> d1 >> d2;
            _box.setYLow(d1); _box.setYHigh(d2);
          case 16: // zlo zhi
            ss >> d1 >> d2;
            _box.setZLow(d1); _box.setZHigh(d2);
          case 17: // xy xz yz
            ss >> d1 >> d2 >> d3;
            _box.setXYTilt(d1); _box.setXZTilt(d2); _box.setYZTilt(d3);
          default: // all the integer values
            ss >> i1;
            vals[i] = i1;
            break;
        }
      }
    }

    // check for section keywords
    for (std::size_t i = 0; i < _section_keywords.size(); ++i) {
      if (str.find(_section_keywords[i]) != std::string::npos) {

      }
    }
  }
}
