# KIT-056: RUV-III correction mode on the DIA_LIMPA route, alongside RUV as covariates

Status 2026-10-25: **in progress**. Plan approved 2026-10-22 (user decisions: position after KIT-053; the step-map unit confirmed whole matrix; default mode `covariates`; no counted member). Base build `kit_057` = `core/browser-smoke.1Nrzlk` (77cbd2c6: the route on the kit's mechanisms, the fastlimpadpc engine, the bounded DPC worker).

## Mechanism (one)

`ruv_mode` is a declared parameter of `protein_ruv_covariates` (`covariates`, the script's method and the default; `correction`, the DIA route's form). One RUV-III fit per run: the same replicate rule (KIT-035), sweep, k and control set; the mode changes only where the unwanted variation is accounted for. Correction commits `ruv_corrected`: the adjusted matrix (Y − Wα, RUV-III's newY) as the state's protein table and as the owned DPC EList's `E`; the DPC standard errors and observation counts kept as quantified; no W terms in the design; dpcDE fits the plain formula; the DA refuses a corrected object whose design still carries W terms (double correction). The record, replay, report overlay (`ruvcov_mode`, `da_ruv_mode`), the tab (three choices), the kit's step list and driver, the recorder (`KIT_LIMPA_TICKET`) and the bindings (ticket `056`, oracles `limpa_corr_screen`) name the mode.

## Evidence

(filled as the checks run)
