let
    vdw = PDB.vanderwaalsradius
    hyd = PDB._hydrophobic
    aro = PDB._aromatic
    cat = PDB._cationic
    ani = PDB._anionic
    SUITE["PDB"]["_generate_interaction_keys"]["defaults"] =
        @benchmarkable PDB._generate_interaction_keys($vdw, $hyd, $aro, $cat, $ani)
end
