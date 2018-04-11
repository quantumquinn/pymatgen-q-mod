# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import re
import six

from os import listdir
from monty.io import zopen

from monty.re import regrep
from collections import defaultdict

from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.io_utils import clean_lines

"""
This module implements input and output processing from PWSCF.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "3/27/15"


class PWInput(object):
    """
    Base input file class. Right now, only supports no symmetry and is
    very basic.
    """

    def __init__(self, structure, pseudo_lib = 'SSSP',pseudo=None, control=None, system=None,
                 electrons=None, ions=None, cell=None, kpoints_mode="automatic",
                 kpoints_grid=(1, 1, 1),kpoints_shift=(0, 0, 0)):
        """
        Initializes a PWSCF input file.

        Args:
            structure (Structure): Input structure. For spin-polarized calculation,
                properties (e.g. {"starting_magnetization": -0.5, 
                "pseudo": "Mn.pbe-sp-van.UPF"}) on each site is needed instead of 
                pseudo (dict).
            pseudo (dict): A dict of the pseudopotentials to use. Default to None.
            control (dict): Control parameters. Refer to official PWSCF doc
                on supported parameters. Default to {"calculation": "scf"}
            system (dict): System parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            electrons (dict): Electron parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            ions (dict): Ions parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            cell (dict): Cell parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            kpoints_mode (str): Kpoints generation mode. Default to automatic.
            kpoints_grid (sequence): The kpoint grid. Default to (1, 1, 1).
            kpoints_shift (sequence): The shift for the kpoints. Defaults to
                (0, 0, 0).
        """
        self.structure = structure
        sections = {}
        sections["control"] = control or {"calculation": "scf"}
        sections["system"] = system or {}
        sections["electrons"] = electrons or {}
        sections["ions"] = ions or {}
        sections["cell"] = cell or {}

        sssp_pseudos = {"Ag" : "ag_pbe_v1.4.uspp.F.UPF",
          "Al" : "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
          "Ar" : "Ar.pbe-n-rrkjus_psl.1.0.0.UPF",
          "As" : "As.pbe-n-rrkjus_psl.0.2.UPF",
          "Au" : "Au_ONCV_PBE-1.0.upf",
          "Ba" : "Ba_ONCV_PBE-1.0.upf",
          "Be" : "Be_ONCV_PBE-1.0.upf",
          "Bi" : "Bi.pbe-dn-kjpaw_psl.0.2.2.UPF",
          "B"  : "B.pbe-n-kjpaw_psl.0.1.UPF",
          "Br" : "br_pbe_v1.4.uspp.F.UPF",
          "Ca" : "Ca_pbe_v1.uspp.F.UPF",
          "Cd" : "Cd.pbe-dn-rrkjus_psl.0.3.1.UPF",
          "Ce" : "Ce.GGA-PBE-paw-v1.0.UPF",
          "Cl" : "Cl.pbe-n-rrkjus_psl.1.0.0.UPF",
          "Co" : "Co_pbe_v1.2.uspp.F.UPF",
          "C"  : "C_pbe_v1.2.uspp.F.UPF",
          "Cr" : "cr_pbe_v1.5.uspp.F.UPF",
          "Cs" : "Cs_pbe_v1.uspp.F.UPF",
          "Cu" : "Cu_pbe_v1.2.uspp.F.UPF",
          "Dy" : "Dy.GGA-PBE-paw-v1.0.UPF",
          "Er" : "Er.GGA-PBE-paw-v1.0.UPF",
          "Eu" : "Eu.GGA-PBE-paw-v1.0.UPF",
          "Fe" : "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
          "F"  : "f_pbe_v1.4.uspp.F.UPF",
          "Ga" : "Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
          "Gd" : "Gd.GGA-PBE-paw-v1.0.UPF",
          "Ge" : "Ge.pbe-dn-kjpaw_psl.1.0.0.UPF",
          "He" : "He_ONCV_PBE-1.0.upf",
          "Hf" : "Hf.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
          "Hg" : "Hg_pbe_v1.uspp.F.UPF",
          "Ho" : "Ho.GGA-PBE-paw-v1.0.UPF",
          "H"  : "H.pbe-rrkjus_psl.0.1.UPF",
          "In" : "In.pbe-dn-rrkjus_psl.0.2.2.UPF",
          "I"  : "I_pbe_v1.uspp.F.UPF",
          "Ir" : "Ir_pbe_v1.2.uspp.F.UPF",
          "K"  : "K.pbe-spn-rrkjus_psl.1.0.0.UPF",
          "Kr" : "Kr.pbe-n-rrkjus_psl.0.2.3.UPF",
          "La" : "La.GGA-PBE-paw-v1.0.UPF",
          "Li" : "li_pbe_v1.4.uspp.F.UPF",
          "Lu" : "Lu.GGA-PBE-paw-v1.0.UPF",
          "Mg" : "mg_pbe_v1.4.uspp.F.UPF",
          "Mn" : "Mn.pbe-spn-kjpaw_psl.0.3.1.UPF",
	  "Mo" : "Mo_ONCV_PBE-1.0.upf",
          "Na" : "Na_ONCV_PBE-1.0.upf",
          "Nb" : "Nb.pbe-spn-kjpaw_psl.0.3.0.UPF",
          "Nd" : "Nd.GGA-PBE-paw-v1.0.UPF",
          "Ne" : "Ne.pbe-n-kjpaw_psl.1.0.0.UPF",
          "Ni" : "ni_pbe_v1.4.uspp.F.UPF",
	  "N"  : "N.pbe.theos.UPF",
          "O"  : "O.pbe-n-kjpaw_psl.0.1.UPF",
          "Os" : "Os.pbe-spfn-rrkjus_psl.1.0.0.UPF",
          "Pb" : "Pb_ONCV_PBE-1.0.upf",
          "Pd" : "Pd_ONCV_PBE-1.0.upf",
          "Pm" : "Pm.GGA-PBE-paw-v1.0.UPF",
          "Po" : "Po.pbe-dn-rrkjus_psl.1.0.0.UPF",
          "P"  : "P.pbe-n-rrkjus_psl.1.0.0.UPF",
          "Pr" : "Pr.GGA-PBE-paw-v1.0.UPF",
          "Pt" : "Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF",
          "Rb" : "Rb_ONCV_PBE-1.0.upf",
          "Re" : "Re_pbe_v1.2.uspp.F.UPF",
          "Rh" : "Rh.pbe-spn-kjpaw_psl.1.0.0.UPF",
          "Rn" : "Rn.pbe-dn-rrkjus_psl.1.0.0.UPF",
          "Ru" : "Ru_ONCV_PBE-1.0.upf",
          "Sb" : "sb_pbe_v1.4.uspp.F.UPF",
          "Sc" : "Sc_pbe_v1.uspp.F.UPF",
          "Se" : "Se_pbe_v1.uspp.F.UPF",
          "Si" : "Si.pbe-n-rrkjus_psl.1.0.0.UPF",
          "Sm" : "Sm.GGA-PBE-paw-v1.0.UPF",
          "Sn" : "Sn_pbe_v1.uspp.F.UPF",
          "S"  : "S_pbe_v1.2.uspp.F.UPF",
          "Sr" : "Sr.pbe-spn-rrkjus_psl.1.0.0.UPF",
          "Ta" : "Ta.pbe-spfn-rrkjus_psl.1.0.0.UPF",
          "Tb" : "Tb.GGA-PBE-paw-v1.0.UPF",
          "Tc" : "Tc_ONCV_PBE-1.0.upf",
          "Te" : "Te_pbe_v1.uspp.F.UPF",
          "Ti" : "ti_pbe_v1.4.uspp.F.UPF",
          "Tl" : "Tl.pbe-dn-rrkjus_psl.1.0.0.UPF",
          "Tm" : "Tm.GGA-PBE-paw-v1.0.UPF",
          "V"  : "V_pbe_v1.uspp.F.UPF",
          "W"  : "W_pbe_v1.2.uspp.F.UPF",
          "Xe" : "Xe.pbe-dn-rrkjus_psl.1.0.0.UPF",
          "Yb" : "Yb.GGA-PBE-paw-v1.0.UPF",
          "Y"  : "Y_pbe_v1.uspp.F.UPF",
          "Zn" : "Zn_pbe_v1.uspp.F.UPF",
          "Zr" : "Zr_pbe_v1.uspp.F.UPF"}

        nc_pseudos = {'Ac' : 'Ac.pbe-n-nc.UPF',
          'Ag' : 'Ag.pbe-n-nc.UPF',
          'Al' : 'Al.pbe-n-nc.UPF',
          'Ar' : 'Ar.pbe-n-nc.UPF',
          'As' : 'As.pbe-n-nc.UPF',
          'At' : 'At.pbe-n-nc.UPF',
          'Au' : 'Au.pbe-n-nc.UPF',
          'Ba' : 'Ba.pbe-n-nc.UPF',
          'Be' : 'Be.pbe-n-nc.UPF',
          'Bi' : 'Bi.pbe-n-nc.UPF',
          'B'  : 'B.pbe-nc.UPF',
          'Br' : 'Br.pbe-n-nc.UPF',
          'Ca' : 'Ca.pbe-n-nc.UPF',
          'Cd' : 'Cd.pbe-n-nc.UPF',
          'Ce' : 'Ce.pbe-n-nc.UPF',
          'Cl' : 'Cl.pbe-n-nc.UPF',
          'Co' : 'Co.pbe-n-nc.UPF',
          'C'  : 'C.pbe-nc.UPF',
          'Cr' : 'Cr.pbe-n-nc.UPF',
          'Cs' : 'Cs.pbe-n-nc.UPF',
          'Cu' : 'Cu.pbe-spn-rrkjus_psl.0.3.1.UPF',
          'Dy' : 'Dy.pbe-n-nc.UPF',
          'Er' : 'Er.pbe-n-nc.UPF',
          'Eu' : 'Eu.pbe-n-nc.UPF',
          'Fe' : 'Fe.pbe-n-nc.UPF',
          'F'  : 'F.pbe-n-nc.UPF',
          'Fr' : 'Fr.pbe-n-nc.UPF',
          'Ga' : 'Ga.pbe-n-nc.UPF',
          'Gd' : 'Gd.pbe-n-nc.UPF',
          'Ge' : 'Ge.pbe-n-nc.UPF',
          'He' : 'He.pbe-n-nc.UPF',
          'Hf' : 'Hf.pbe-n-nc.UPF',
          'Hg' : 'Hg.pbe-n-nc.UPF',
          'Ho' : 'Ho.pbe-n-nc.UPF',
          'H'  : 'H.pbe-n-nc.UPF',
          'In' : 'In.pbe-n-nc.UPF',
          'I'  : 'I.pbe-n-nc.UPF',
          'Ir' : 'Ir.pbe-n-nc.UPF',
          'K'  : 'K.pbe-n-nc.UPF',
          'Kr' : 'Kr.pbe-n-nc.UPF',
          'La' : 'La.pbe-n-nc.UPF',
          'Li' : 'Li.pbe-n-nc.UPF',
          'Lu' : 'Lu.pbe-n-nc.UPF',
          'Mg' : 'Mg.pbe-n-nc.UPF',
          'Mn' : 'Mn.pbe-n-nc.UPF',
          'Mo' : 'Mo.pbe-n-nc.UPF',
          'Na' : 'Na.pbe-n-nc.UPF',
          'Nb' : 'Nb.pbe-n-nc.UPF',
          'Nd' : 'Nd.pbe-n-nc.UPF',
          'Ne' : 'Ne.pbe-nc.UPF',
          'Ni' : 'Ni.pbe-n-nc.UPF',
          'N'  : 'N.pbe-nc.UPF',
          'Np' : 'Np.pbe-n-nc.UPF',
          'O'  : 'O.pbe-nc.UPF',
          'Os' : 'Os.pbe-n-nc.UPF',
          'Pa' : 'Pa.pbe-n-nc.UPF',
          'Pb' : 'Pb.pbe-n-nc.UPF',
          'Pd' : 'Pd.pbe-n-nc.UPF',
          'Pm' : 'Pm.pbe-n-nc.UPF',
          'Po' : 'Po.pbe-n-nc.UPF',
          'P'  : 'P.pbe-n-nc.UPF',
          'Pr' : 'Pr.pbe-n-nc.UPF',
          'Pt' : 'Pt.pbe-n-nc.UPF',
          'Pu' : 'Pu.pbe-n-nc.UPF',
          'Ra' : 'Ra.pbe-n-nc.UPF',
          'Rb' : 'Rb.pbe-n-nc.UPF',
          'Rh' : 'Rh.pbe-n-nc.UPF',
          'Rn' : 'Rn.pbe-n-nc.UPF',
          'Re' : 'Re.pbe-n-nc.UPF',
          'Ru' : 'Ru.pbe-n-nc.UPF',
          'Sb' : 'Sb.pbe-n-nc.UPF',
          'Sc' : 'Sc.pbe-n-nc.UPF',
          'Se' : 'Se.pbe-n-nc.UPF',
          'Si' : 'Si.pbe-n-nc.UPF',
          'Sm' : 'Sm.pbe-n-nc.UPF',
          'Sn' : 'Sn.pbe-n-nc.UPF',
          'S'  : 'S.pbe-n-nc.UPF',
          'Sr' : 'Sr.pbe-n-nc.UPF',
          'Ta' : 'Ta.pbe-n-nc.UPF',
          'Tb' : 'Tb.pbe-n-nc.UPF',
          'Tc' : 'Tc.pbe-n-nc.UPF',
          'Te' : 'Te.pbe-n-nc.UPF',
          'Th' : 'Th.pbe-n-nc.UPF',
          'Ti' : 'Ti.pbe-n-nc.UPF',
          'Tl' : 'Tl.pbe-n-nc.UPF',
          'Tm' : 'Tm.pbe-n-nc.UPF',
          'U'  : 'U.pbe-n-nc.UPF',
          'V'  : 'V.pbe-n-nc.UPF',
          'Xe' : 'Xe.pbe-n-nc.UPF',
          'Yb' : 'Yb.pbe-n-nc.UPF',
          'Y'  : 'Y.pbe-n-nc.UPF',
          'Zn' : 'Zn.pbe-nc.UPF',
          'Zr' : 'Zr.pbe-n-nc.UPF'}
          

        

        if pseudo == None:
             pseudo = {}
             for species in self.structure.composition.keys():
                 try:
                     if pseudo_lib == 'SSSP':
                         pseudo[species.symbol]=sssp_pseudos[species.symbol]
                     elif pseudo_lib == 'NC':
                         pseudo[species.symbol]=nc_pseudos[species.symbol]
                     elif pseudo_lib == 'GBRV_US':
                         #This obviously will only work on my computer, directory can easily be changed for others
                         possible_pseudos = listdir('/home/quinn/Documents/pymatgen_slab_tests/dft+U/GBRV_USPP_PBE_UPF_format')
                         for pseudo_name in possible_pseudos:
                             #print( pseudo_name )
                             if species.symbol.lower() in pseudo_name:
                                 #This gets around the issue of C matching to both C and Cr
                                 #There's almost certainly a better regular expression way to do this, but oh well ¯\_(ツ)_/¯
                                 if len(species.symbol) ==1:
                                     if pseudo_name[0] == species.symbol.lower() and pseudo_name[1] == '_':
                                         pseudo[species.symbol] = pseudo_name
                                         break
                                 else:
                                     pseudo[species.symbol] = pseudo_name
                                     break
                         if species.symbol not in pseudo.keys():
                             raise PWInputError("Missing %s pseudo"%species.symbol)          
                     else:
                         raise PWInputError("pseudo lib %s not supported yet. Have fun adding that yourself"%pseudo_lib)
                 except:
                     raise PWInputError("Missing %s pseudo"%species.symbol)
#            for site in structure:
 #               try:
  #                  site.properties['pseudo']
   #             except KeyError:
    #                raise PWInputError("Missing %s in pseudo specification!" 
     #                                  % site)
        else:
            for species in self.structure.composition.keys():
                if species.symbol not in pseudo:
                    raise PWInputError("Missing %s in pseudo specification!" 
                                       % species.symbol)
        self.pseudo = pseudo

        self.sections = sections
        self.kpoints_mode = kpoints_mode
        self.kpoints_grid = kpoints_grid
        self.kpoints_shift = kpoints_shift

    def __str__(self):
        out = []
        site_descriptions = {}

        if self.pseudo != None:
            site_descriptions = self.pseudo
        else:
            c = 1
            for site in self.structure:
                name = None
                for k, v in site_descriptions.items():
                    if site.properties == v:
                        name = k

                if name == None:
                    name = site.specie.symbol+str(c)
                    site_descriptions[name] = site.properties
                    c += 1

        def to_str(v):
            if isinstance(v, six.string_types):
                return "'%s'" % v
            elif isinstance(v, float):
                return "%s" % str(v).replace("e", "d")
            elif isinstance(v, bool):
                if v:
                    return ".TRUE."
                else:
                    return ".FALSE."
            return v


        elements_w_dftu = ['Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 
             'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
             'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'La', 'Sr','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', 
             'Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','C', 'N', 'O', 'Al', 'Si', 'As' , 'Ga', 'In']


        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                if isinstance(v1[k2], list):
                    n = 1
                    for l in v1[k2][:len(site_descriptions)]:
                        sub.append("  %s(%d) = %s" % (k2, n, to_str(v1[k2][n-1])))
                        n += 1
                else:
                    sub.append("  %s = %s" % (k2, to_str(v1[k2])))
            if k1 == "system":
                if 'ibrav' not in self.sections[k1]:
                    sub.append("  ibrav = 0")
                if 'nat' not in self.sections[k1]:
                    sub.append("  nat = %d" % len(self.structure))
                if 'ntyp' not in self.sections[k1]:
                    sub.append("  ntyp = %d" % len(site_descriptions))
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k).group(0)
            if self.pseudo is not None:
                p = v
            else:
                p = v['pseudo']
            out.append("  %s  %.4f %s" % (k, Element(e).atomic_mass, p)) 

        out.append("ATOMIC_POSITIONS crystal")
        #Making sure that elements with DFT+U are listed first. 
        #It's dumb, but that's how the auto DFT+U software works
        if self.pseudo is not None:
            for site in self.structure:
                if site.specie.symbol in elements_w_dftu:
                    out.append("  %s %.6f %.6f %.6f" % (site.specie.symbol, site.a,
                                                    site.b, site.c))
            for site in self.structure:
                if site.specie.symbol not in elements_w_dftu:
                    out.append("  %s %.6f %.6f %.6f" % (site.specie.symbol, site.a,
                                                    site.b, site.c))
        else:
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(),
                                   key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                if name in elements_w_dftu:
                    out.append("  %s %.6f %.6f %.6f" % (name, site.a, site.b, site.c))
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(),
                                   key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                if name not in elements_w_dftu:
                    out.append("  %s %.6f %.6f %.6f" % (name, site.a, site.b, site.c))

        out.append("K_POINTS %s" % self.kpoints_mode)
        kpt_str = ["%s" % i for i in self.kpoints_grid]
        kpt_str.extend(["%s" % i for i in self.kpoints_shift])
        out.append("  %s" % " ".join(kpt_str))
        out.append("CELL_PARAMETERS angstrom")
        for vec in self.structure.lattice.matrix:
            out.append("  %f %f %f" % (vec[0], vec[1], vec[2]))
        return "\n".join(out)

    def write_file(self, filename):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads an PWInput object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            PWInput object
        """
        with zopen(filename, "rt") as f:
            return PWInput.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an PWInput object from a string.

        Args:
            string (str): PWInput string

        Returns:
            PWInput object
        """
        lines = list(clean_lines(string.splitlines()))

        def input_mode(line):
            if line[0] == "&":
                return ("sections", line[1:].lower())
            elif "ATOMIC_SPECIES" in line:
                return ("pseudo", )
            elif "K_POINTS" in line:
                return ("kpoints", line.split("{")[1][:-1])
            elif "CELL_PARAMETERS" in line or "ATOMIC_POSITIONS" in line:
                return ("structure", line.split("{")[1][:-1])
            elif line == "/":
                return None
            else:
                return mode

        sections = {"control": {}, "system": {}, "electrons": {}, 
                    "ions": {}, "cell":{}}
        pseudo = {}
        pseudo_index = 0
        lattice = []
        species = []
        coords = []
        structure = None
        site_properties = {"pseudo":[]}
        mode = None
        for line in lines:
            mode = input_mode(line)
            if mode == None:
                pass
            elif mode[0] == "sections":
                section = mode[1]
                m = re.match(r'(\w+)\(?(\d*?)\)?\s*=\s*(.*)', line)
                if m:
                    key = m.group(1).strip()
                    key_ = m.group(2).strip()
                    val = m.group(3).strip()
                    if key_ != "":
                        if sections[section].get(key, None) == None:
                            val_ = [0.0]*20 # MAX NTYP DEFINITION
                            val_[int(key_)-1] = PWInput.proc_val(key, val)
                            sections[section][key] = val_

                            site_properties[key] = []
                        else:
                            sections[section][key][int(key_)-1] = PWInput.proc_val(key, val) 
                    else:
                        sections[section][key] = PWInput.proc_val(key, val)

            elif mode[0] == "pseudo":
                m = re.match(r'(\w+)\s+(\d*.\d*d?\d*)\s+(.*)', line)
                if m:
                    
                    pseudo[m.group(1).strip()] = {}
                    pseudo[m.group(1).strip()]["index"] = pseudo_index
                    pseudo[m.group(1).strip()]["pseudopot"] = m.group(3).strip()
                    pseudo_index += 1
                #print( pseudo)
            elif mode[0] == "kpoints":
                m = re.match(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
                if m:
                    kpoints_grid = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                    kpoints_shift = (int(m.group(4)), int(m.group(5)), int(m.group(6)))
                else:
                    kpoints_mode = mode[1]
            elif mode[0] == "structure":
                m_l = re.match(r'(-?\d+\.?\d*d?\d*)\s+(-?\d+\.?\d*d?\d*)\s+(-?\d+\.?\d*d?\d*)', line)
                m_p = re.match(r'(\w+)\s+(-?\d+\.\d*d?\d*)\s+(-?\d+\.?\d*d?\d*)\s+(-?\d+\.?\d*d?\d*)', line)
                alat = sections['system']['celldm'][0]*0.529177
                #print( pseudo)
                if m_l:
                    #print( m_l.groups())
                    if alat_multiply:
                        lattice += [ float(m_l.group(1).replace('d','e'))*alat, float(m_l.group(2).replace('d','e'))*alat, float(m_l.group(3).replace('d','e'))*alat ]
                    else:
                        lattice += [ float(m_l.group(1).replace('d','e')), float(m_l.group(2).replace('d','e')), float(m_l.group(3).replace('d','e')) ]
                elif m_p:
                    site_properties["pseudo"].append(pseudo[m_p.group(1)]["pseudopot"])
                    species += [m_p.group(1)]
                    if alat_multiply:
                        coords += [[float(m_p.group(2).replace('d','e'))*alat, float(m_p.group(3).replace('d','e'))*alat, float(m_p.group(4).replace('d','e'))*alat]]
                    else:
                        coords += [[float(m_p.group(2).replace('d','e')), float(m_p.group(3).replace('d','e')), float(m_p.group(4).replace('d','e'))]]
                    for k, v in site_properties.items():
                        if k != "pseudo":
                            site_properties[k].append(sections['system'][k][pseudo[m_p.group(1)]["index"]])
                if mode[1] == "angstrom":
                    coords_are_cartesian = True
                    alat_multiply = False
                elif mode[1] == "crystal":
                    coords_are_cartesian = False
                    alat_multiply = False
                elif mode[1] == "alat":
                    coords_are_cartesian = True
                    alat_multiply = True

        #print( site_properties)

        structure = Structure(Lattice(lattice), species, coords, 
                              coords_are_cartesian=coords_are_cartesian,
                              site_properties=site_properties)
        return PWInput(structure=structure, control=sections["control"],
                       system=sections["system"], electrons=sections["electrons"], 
                       ions=sections["ions"], cell=sections["cell"], kpoints_mode=kpoints_mode,
                       kpoints_grid=kpoints_grid, kpoints_shift=kpoints_shift)

    def proc_val(key, val):
        """
        Static helper method to convert PWINPUT parameters to proper type, e.g.,
        integers, floats, etc.

        Args:
            key: PWINPUT parameter key
            val: Actual value of PWINPUT parameter.
        """
        float_keys = ('etot_conv_thr','forc_conv_thr','conv_thr','Hubbard_U','Hubbard_J0','defauss',
                      'starting_magnetization','celldm')

        int_keys = ('nstep','iprint','nberrycyc','gdir','nppstr','ibrav','nat','ntyp','nbnd','nr1',
                    'nr2','nr3','nr1s','nr2s','nr3s','nspin','nqx1','nqx2','nqx3','lda_plus_u_kind',
                    'edir','report','esm_nfit','space_group','origin_choice','electron_maxstep',
                    'mixing_ndim','mixing_fixed_ns','ortho_para','diago_cg_maxiter','diago_david_ndim',
                    'nraise','bfgs_ndim','if_pos','nks','nk1','nk2','nk3','sk1','sk2','sk3','nconstr')

        bool_keys = ('wf_collect','tstress','tprnfor','lkpoint_dir','tefield','dipfield','lelfield',
                     'lorbm','lberry','lfcpopt','monopole','nosym','nosym_evc','noinv','no_t_rev',
                     'force_symmorphic','use_all_frac','one_atom_occupations','starting_spin_angle',
                     'noncolin','x_gamma_extrapolation','lda_plus_u','lspinorb','london',
                     'ts_vdw_isolated','xdm','uniqueb','rhombohedral','realxz','block',
                     'scf_must_converge','adaptive_thr','diago_full_acc','tqr','remove_rigid_rot',
                     'refold_pos')

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in bool_keys:
                if val.lower() == ".true.":
                    return True
                elif val.lower() == ".false.":
                    return False
                else:
                    raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*d?-?\d*", val.lower()).group(0).replace("d", "e"))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        try:
            val = val.replace("d","e")
            return smart_int_or_float(val)
        except ValueError:
            pass

        if "true" in val.lower():
            return True
        if "false" in val.lower():
            return False

        m = re.match(r"^[\"|'](.+)[\"|']$", val)
        if m:
            return m.group(1)



class PWInputError(BaseException):
    pass


class PWOutput(object):

    patterns = {
        "energies": r'total energy\s+=\s+([\d\.\-]+)\sRy',
        "ecut": r'kinetic\-energy cutoff\s+=\s+([\d\.\-]+)\s+Ry',
        "lattice_type": r'bravais\-lattice index\s+=\s+(\d+)',
        "celldm1": r"celldm\(1\)=\s+([\d\.]+)\s",
        "celldm2": r"celldm\(2\)=\s+([\d\.]+)\s",
        "celldm3": r"celldm\(3\)=\s+([\d\.]+)\s",
        "celldm4": r"celldm\(4\)=\s+([\d\.]+)\s",
        "celldm5": r"celldm\(5\)=\s+([\d\.]+)\s",
        "celldm6": r"celldm\(6\)=\s+([\d\.]+)\s",
        "nkpts": r"number of k points=\s+([\d]+)"
    }

    def __init__(self, filename):
        self.filename = filename
        self.data = defaultdict(list)
        self.read_pattern(PWOutput.patterns)
        for k, v in self.data.items():
            if k == "energies":
                self.data[k] = [float(i[0][0]) for i in v]
            elif k in ["lattice_type", "nkpts"]:
                self.data[k] = int(v[0][0][0])
            else:
                self.data[k] = float(v[0][0][0])

    def read_pattern(self, patterns, reverse=False,
                     terminate_on_match=False, postprocess=str):
        """
        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned
            values are lists of lists, because you can grep multiple
            items on one line.
        """
        matches = regrep(self.filename, patterns, reverse=reverse,
                         terminate_on_match=terminate_on_match,
                         postprocess=postprocess)
        self.data.update(matches)

    def get_celldm(self, i):
        return self.data["celldm%d" % i]

    @property
    def final_energy(self):
        return self.data["energies"][-1]

    @property
    def lattice_type(self):
        return self.data["lattice_type"]
