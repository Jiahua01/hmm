import logging
from code_generation.quantities import NanoAODQuantity

log = logging.getLogger(__name__)

Tau_pt = NanoAODQuantity("Tau_pt")
Tau_eta = NanoAODQuantity("Tau_eta")
Tau_phi = NanoAODQuantity("Tau_phi")
Tau_mass = NanoAODQuantity("Tau_mass")
Tau_dz = NanoAODQuantity("Tau_dz")
Tau_dxy = NanoAODQuantity("Tau_dxy")
Tau_charge = NanoAODQuantity("Tau_charge")
Tau_decayMode = NanoAODQuantity("Tau_decayMode")
Tau_genMatch = NanoAODQuantity("Tau_genPartFlav")
Tau_IDraw = NanoAODQuantity("Tau_rawDeepTau2017v2p1VSjet")

Muon_pt = NanoAODQuantity("Muon_pt")
Muon_eta = NanoAODQuantity("Muon_eta")
Muon_phi = NanoAODQuantity("Muon_phi")
Muon_mass = NanoAODQuantity("Muon_mass")
Muon_iso = NanoAODQuantity("Muon_pfRelIso04_all")
Muon_dz = NanoAODQuantity("Muon_dz")
Muon_dxy = NanoAODQuantity("Muon_dxy")
Muon_charge = NanoAODQuantity("Muon_charge")

Electron_pt = NanoAODQuantity("Electron_pt")
Electron_eta = NanoAODQuantity("Electron_eta")
Electron_phi = NanoAODQuantity("Electron_phi")
Electron_mass = NanoAODQuantity("Electron_mass")
Electron_iso = NanoAODQuantity("Electron_pfRelIso03_all")

Jet_eta = NanoAODQuantity("Jet_eta")
Jet_phi = NanoAODQuantity("Jet_phi")
Jet_pt = NanoAODQuantity("Jet_pt")
Jet_mass = NanoAODQuantity("Jet_mass")