ITP_TEMPLATE = """\
;----------------------------TITLE -----------------------------------------------------------------------------------------
;   None
;
; This file was generated at {time} on {date} using the ATB API, using information generated by
;
;                  Automatic Topology Builder  
;
;                   REVISION {revision}
;---------------------------------------------------------------------------------------------------------------------------
; Authors     : Martin Stroet, Bertrand Caron, Alpeshkumar K. Malde, Thomas Lee, Alan E. Mark
;
; Institute   : Molecular Dynamics group, 
;               School of Chemistry and Molecular Biosciences (SCMB),
;               The University of Queensland, QLD 4072, Australia
; URL         : https://atb.uq.edu.au
; Citations   : 1. Malde AK, Zuo L, Breeze M, Stroet M, Poger D, Nair PC, Oostenbrink C, Mark AE.
;                  An Automated force field Topology Builder (ATB) and repository: version 1.0.
;                  Journal of Chemical Theory and Computation, 2011, 7, 4026-4037.
;               2. Stroet M, Caron B, Visscher K, Geerke D, Malde AK, Mark AE.
;                  Automated Topology Builder version 3.0: Prediction of solvation free enthalpies in water and hexane.
;                  DOI:10.1021/acs.jctc.8b00768
;
; Disclaimer  : 
;      While every effort has been made to ensure the accuracy and validity of parameters provided below
;      the assignment of parameters is being based on an automated procedure combining data provided by a
;      given user as well as calculations performed using third party software. They are provided as a guide.
;      The authors of the ATB cannot guarantee that the parameters are complete or that the parameters provided
;      are appropriate for use in any specific application. Users are advised to treat these parameters with discretion
;      and to perform additional validation tests for their specific application if required. Neither the authors
;      of the ATB or The University of Queensland except any responsibly for how the parameters may be used.
;
; Release notes and warnings: 
;  (1) The topology is based on a set of atomic coordinates and other data provided by the user after
;      after quantum mechanical optimization of the structure using different levels of theory depending on
;      the nature of the molecule.
;  (2) In some cases the automatic bond, bond angle and dihedral type assignment is ambiguous.
;      In these cases alternative type codes are provided at the end of the line.
;  (3) While bonded parameters are taken where possible from the nominated force field non-standard bond, angle and dihedral
;      type code may be incorporated in cases where an exact match could not be found. These are marked as "non-standard"
;      or "uncertain" in comments.
;  (4) In some cases it is not possible to assign an appropriate parameter automatically. "%%" is used as a place holder
;      for those fields that could not be determined automatically. The parameters in these fields must be assigned manually
;      before the file can be used.
;---------------------------------------------------------------------------------------------------------------------------
; Input Structure : {residue_name}
; Output          : {resolution_upper} topology
;	Use in conjunction with the corresponding {resolution} PDB file.
;---------------------------------------------------------------------------------------------------------------------------
; Citing this topology file
; ATB molid: {molecule_molid}
; ATB Topology Hash: {molecule_hash}
;---------------------------------------------------------------------------------------------------------------------------
; Intermediate Topology Generation was performed using:
; A B3LYP/6-31G* optimized geometry.
; Bonded and van der Waals parameters were taken from the GROMOS 54A7 parameter set.
; Initial charges were estimated using the ESP method of Merz-Kollman.
; Final charges and charge groups were generated by method described in the ATB paper.
;---------------------------------------------------------------------------------------------------------------------------
;
;
[ moleculetype ]
; Name   nrexcl
{residue_name}     3
[ atoms ]
;  nr  type  resnr  resid  atom  cgnr  charge    mass
{atoms}
; total charge of the molecule:  {total_charge:.3f}
[ bonds ]
;  ai   aj  funct   c0         c1
{bonds}
[ pairs ]
;  ai   aj  funct  ;  all 1-4 pairs but the ones excluded in GROMOS itp
{pairs}
[ angles ]
;  ai   aj   ak  funct   angle     fc
{angles}
[ dihedrals ]
; GROMOS improper dihedrals
;  ai   aj   ak   al  funct   angle     fc
{impropers}
[ dihedrals ]
;  ai   aj   ak   al  funct    ph0      cp     mult
{dihedrals}
[ exclusions ]
;  ai   aj  funct  ;  GROMOS 1-4 exclusions
{exclusions}"""