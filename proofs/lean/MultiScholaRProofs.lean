module

public import MultiScholaRProofs.Spec.FilterCount
public import MultiScholaRProofs.Spec.ReplayRule
public import MultiScholaRProofs.Spec.GoClosure
public import MultiScholaRProofs.Spec.QvalueBh
public import MultiScholaRProofs.Spec.Permutation
public import MultiScholaRProofs.Spec.RuvAlgebra
public import MultiScholaRProofs.Spec.RuvDownstream
public import MultiScholaRProofs.Spec.RuvIIIC
public import MultiScholaRProofs.Spec.SharperWitnesses
public import MultiScholaRProofs.Spec.RemovedFactors
public import MultiScholaRProofs.Spec.RuvIIICPattern
public import MultiScholaRProofs.Spec.RuvControlRules
public import MultiScholaRProofs.Oracle.GoClosure
public import MultiScholaRProofs.Oracle.FilterCount
public import MultiScholaRProofs.Oracle.ReplayRule
public import MultiScholaRProofs.Proof.GoClosure
public import MultiScholaRProofs.Proof.ReplayRule
public import MultiScholaRProofs.Proof.FilterCount
public import MultiScholaRProofs.Proof.QvalueBh
public import MultiScholaRProofs.Proof.Permutation
public import MultiScholaRProofs.Proof.RuvAlgebra
public import MultiScholaRProofs.Proof.RuvDownstream
public import MultiScholaRProofs.Proof.SharperWitnesses
public import MultiScholaRProofs.Proof.RuvAlgebra2
public import MultiScholaRProofs.Proof.RuvDownstream2
public import MultiScholaRProofs.Proof.RuvIIIC
public import MultiScholaRProofs.Proof.KeptFamily
public import MultiScholaRProofs.Proof.RemovedFactors
public import MultiScholaRProofs.Proof.RuvIIICPattern

/-!
# KIT-193: machine-checked statements about exact models of parts of MultiScholaR's methods

Root module: imports every Spec module (models and statements, `def …_stmt : Prop`, frozen: `SPEC.sha256`), every
Proof module and every Oracle module (the tables of section 5.5, written by `Kit193Oracles.lean`).
-/
