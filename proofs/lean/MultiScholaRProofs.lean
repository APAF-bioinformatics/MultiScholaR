module

public import MultiScholaRProofs.Spec.FilterCount
public import MultiScholaRProofs.Spec.ReplayRule
public import MultiScholaRProofs.Spec.GoClosure
public import MultiScholaRProofs.Spec.QvalueBh
public import MultiScholaRProofs.Spec.Permutation
public import MultiScholaRProofs.Spec.RuvAlgebra
public import MultiScholaRProofs.Spec.RuvDownstream
public import MultiScholaRProofs.Spec.RuvIIIC

/-!
# KIT-193: machine-checked statements about exact models of parts of MultiScholaR's methods

Root module: imports every Spec module (models and statements, `def …_stmt : Prop`) and, once written, every Proof
module. At this stage (the statements gate, KIT-193 section 5.3) only the Spec modules exist.
-/
