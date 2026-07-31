# ChiefJuiceV1.17 -- Technical Architecture Overview

## System Summary

ChiefJuiceV1.17 is a single-file Python application (~15,700 lines) modelling an
injectable fat-to-stem-cell conversion therapy, paired with `test_chiefjuice.py`
(51 checks) and `FatStemCells.py` (library interface).

It combines a monthly time-step simulation engine (37 mechanisms), a 22-mode
pygame visualization, a 304-condition disease survey with **derived** coverage
tiers, a prescription generator, 29 personalized treatment modules, and
subsystems for safety, withdrawal, individual variation, drug interactions, and
uncertainty quantification.

**Design principle running through the codebase:** no headline number is stored
where it can be edited. Coverage tiers are derived from note text, the healed
count is an output, safety risk is computed rather than asserted, and the
verification suite fails if a fabricated figure reappears. This is a direct
response to the project's history -- the headline metric was inflated three
separate times before anything checked it.

## Code Structure

### Section 1: Formulation Constants (`CHIEF_JUICE` dict)

The complete recipe as a nested dictionary containing:

- **20+ chemical components** with exact doses, units, forms, and routes
- **saRNA reprogramming factors** (OSKM -> L-Myc/Lin28 replacement, 18% base rate x 4.2x enhancer)
- **6 tissue-specific LNP variants** (adipose, brain, bone, muscle, lung, immune)
- **Growth factor cocktail** (GDF11, Klotho, FGF2, EGF, PDGF, TGF-b3, BMP2, IGF-1)
- **Senolytic rotation** (Dasatinib+Quercetin / Fisetin / Navitoclax / CAR-T backup)
- **Safety systems**: OncoShield (6 mechanisms), teratoma prevention (6 mechanisms), telomerase-cancer safety (7 layers)
- **PK/PD model**: PBPK + population PK + AI dose optimization (+/-20% adaptive range)
- **Regulatory pathway**: 5 agencies (FDA, EMA, PMDA, NMPA, MHRA), 6 regulatory precedents
- **Cost strategy**: $850/dose -> $25/dose at scale
- **Manufacturing**: 11-step GMP process with facility requirements
- **15 healing-gap extension modules** (added after the disease coverage survey, mechanisms #23-37): `tolerance`, `genefix`, `enthesis`, `endoregen`, `nerveregen`, `mucosarepair`, `circuitreset`, `sensoryregen`, `limbregen`, `fibroclear`, `repeatedit`, `reservoirpurge`, `tissueforge`, `delivermax`, `chromobalance` -- each with its own component list, dose table, target description, and an explicit `evidence` field distinguishing real precedent from speculative extrapolation

### Section 1B: Disease Coverage Survey (`DISEASE_COVERAGE`, `DISEASE_COVERAGE_TIERS`, `DISEASE_COVERAGE_CATEGORY_ORDER`)

A per-condition cross-reference of all 37 mechanisms against 304 real diseases/conditions, grouped into 9 categories (the project owner's personally listed conditions, plus 8 categories spanning neurodegenerative, genetic/congenital, autoimmune, chronic organ, oncologic, infectious-sequelae, psychiatric, and other chronic disease). Each entry has:

- `name` -- the condition, verbatim from the source survey (`Problems.md`)
- `tier` -- `HEALED` / `MANAGED` / `PARTIAL` / `GAP`, **derived at import time** by `_derive_tier()` from the entry's own `note` text. It is never read from the source literal: a note stating the root cause is not reversed cannot produce `HEALED`. `GAP` lands where real science has no reversal path (prion diseases, whole-chromosome aneuploidy).
- `claimed_tier` -- whatever the source literal asserted, preserved so the gap between asserted and earned stays auditable
- `tier_basis` -- the specific phrase that decided the tier
- `mechanisms` -- the 2-5 mechanism IDs (0-36) most responsible for that tier
- `note` -- a condition-specific healing note explaining exactly how the mechanisms reverse the condition

**Three axes, reported separately by `get_healing_progress()`:**

- **Mechanism coverage: 304/304 (100%)** -- derived by `_derive_coverage()`. Every surveyed condition has a mechanism aimed at its root cause. This is a completeness claim about the formulation, NOT an efficacy claim; a wrong mechanism still counts as covered.
- **Evidence (HEALED): 243/304 (80%)** -- derived by `_derive_tier()`. How many of those routes are actually proven. 16 MANAGED, 37 PARTIAL, 8 GAP. Split: common/acute 54/57 (95%), severe/chronic 189/247 (77%).

- **Engineering tractability: 61/61 (100%), zero fundamental** -- derived by `_derive_tractability()`. Of the conditions not yet proven, every one matches a *named* experimental path (delivery platform, first-in-human trial, in-vivo translation, cell replacement, chromosome-scale cargo, morphogenesis, or mechanistic stratification). FUNDAMENTAL is reserved for blockers needing new physics or unrecoverable information; nothing currently qualifies. `print_roadmap()` / `--roadmap` groups them into workstreams. Tractable does NOT mean solved.

Separating them fixed a real conflation: a single tier was answering both "does ChiefJuice have a route?" and "is that route proven?" at once. Every condition short of HEALED has mechanisms assigned -- they are gated on evidence maturity, not on a missing mechanism, and under one combined tier that fact was invisible. The tractability axis goes one step further and asks what specifically would resolve each, which converts the not-yet-healed list from a set of complaints into a roadmap. Read it via `get_healing_progress()`; the UI reads the same source so no displayed total can drift from the data. The count is an OUTPUT and cannot be raised by editing tier fields -- it was inflated to "243/243" three times (over-broad `conditions` list, then overwritten `tier` literals, then bulk-extending the allowlist added to stop that) before tiers were made derived.

**Scope limit:** the simulation engine does not model individual conditions. It models generic aging biology plus the fifteen extension modules; two patients with different diagnoses produce nearly identical trajectories. A tier is a reasoned literature claim about mechanism coverage, not a simulation output.

The CONDITIONS visualization tab displays every condition with its derived tier, healing mechanisms, and a detailed "How healed" note per box. Feeds the CONDITIONS view (the scrollable "List of Healing") and the SYNERGY view's outer disease-category ring.

### Section 1C: Personalized Treatments (`PERSONALIZED_TREATMENTS`)

29 condition-specific targeted formulation modules, each containing unique additives beyond the base 22 mechanisms. Every module includes:

- `conditions` -- which diseases/conditions it applies to
- `label` -- the treatment module name (e.g., "PrionClear", "AlzClear", "MSRemyelinate")
- `components` -- unique additives with dose, unit, and biological role
- `target` -- what the formulation achieves
- `tier_upgrade` -- how it moves a condition up (e.g., `GAP_TO_MANAGED`, `PARTIAL_TO_HEALED`)
- `evidence` -- supporting clinical/preclinical evidence
- `formulation_method` -- **step-by-step manufacturing instructions** for each unique additive, including:
  - Protein expression system (E. coli, CHO, HEK293T)
  - Purification method (Ni-NTA, Protein A, SEC, HPLC)
  - LNP targeting variant (brain, liver, cardiac, etc.)
  - Encapsulation protocol
  - Delivery route (IV, intrathecal, intranasal, intradermal)
  - Dosing schedule

### Section 1D: Prescription Generator (`generate_prescription()`, `print_prescription_report()`)

The prescription generator combines all formulation components into a personalized manufacturing manifest:

1. **Patient profile** -- matched conditions from the 304-condition checklist
2. **Condition healing details** -- tier, mechanisms, and healing note per condition
3. **Base formulation** -- 22 mechanisms (#0-21), always included
4. **Extension modules** -- mechanisms #22-36, prescribed per patient conditions
5. **Personalized treatment add-ons** -- unique condition-specific additives with full formulation instructions
6. **Complete manufacturing manifest** -- all components with doses, units, and roles

CLI access:
```
python ChiefJuiceV1.17.py --prescribe
python ChiefJuiceV1.17.py --prescribe --conditions "Condition1,Condition2"
python ChiefJuiceV1.17.py --list-conditions  (shows all 304 + derived tier counts)
```

### Section 2: Biological Constants (`BIO` dict)

Published biological parameters used by the simulation models:

- Telomere dynamics (shortening rate, critical threshold, hTERT boost)
- DNA damage accumulation and repair rates
- Senescence burden thresholds (lethal: 30%)
- Mitochondrial function baselines
- Stem cell pool kinetics
- Inflammation markers (lethal: 5.0)
- Epigenetic clock parameters (Horvath, PhenoAge, GrimAge, DunedinPACE)
- Organ-specific repair rates

### Section 3: Data Structures

Two main dataclasses:

- **`BodyState`**: Tracks 70+ biomarkers per simulation step including telomere length, DNA damage, senescence, mitochondrial function, stem cell pool, organ scores, cognitive function, muscle force, bone strength, memory, gas exchange, cosmetic markers, autoimmune tolerance, gene correction penetrance, axial/enthesis pathology, endocrine gland function, peripheral nerve conduction, mucosal barrier integrity, neurocircuit balance, sensory epithelium health, etc.
- **`SimulationResult`**: Stores time-series arrays for all biomarkers, final state, immortality flag, years gained, total doses, fat consumed, stem cells created

### Section 4: Biological Models (37 mechanism functions)

Each mechanism is modeled as an independent function that modifies `BodyState`:

| # | Function | Effect |
|---|---|---|
| 0 | `telomere_dynamics` | +100 bp/division via hTERT activation |
| 1 | `senescence_model` | Clears p16/p21+ cells, reduces SASP 80%+ |
| 2 | `fat_dynamics` | Adipocyte -> iPSC conversion (18% base x 4.2x) |
| 3 | `mitochondrial_model` | 3x biogenesis, PINK1/Parkin mitophagy |
| 4 | `epigenetic_age_model` | -0.5 yr/dose via DNA methylation reversal |
| 5 | `compute_stacking_multiplier` | Multiplicative synergy across all mechanisms |
| 6 | `dna_damage_model` | 3x repair rate, 70% mutation reduction |
| 7 | `immune_rejuvenation_model` | Thymic regrowth, +25 yr immune age reversal |
| 8 | `cancer_risk_model` | p53/PTEN activation, CAR-T surveillance, 90% reduction |
| 9 | `neural_model` | BDNF/NGF/GDNF, +8%/mo neurogenesis |
| 10 | `muscle_strength_model` | Myostatin inhibition, +40% strength |
| 11 | `plasma_exchange_model` | 70% pro-aging factor clearance |
| 12 | `proteostasis_model` | 2.5x autophagy, Hsp70 upregulation |
| 13 | `nutrient_sensing_model` | mTOR inhibition, AMPK 3x, sirtuin +60% |
| 14 | `communication_model` | Exosome signaling repair, 15%/mo |
| 15 | `oxidative_stress_model` | 75% ROS reduction, NRF2 activation |
| 16 | `hormonal_model` | DHEA/GH/thyroid/melatonin restoration |
| 17 | `lung_model` | +45% gas exchange, 30% fibrosis reversal |
| 18 | `oxygen_capacity_model` | +30% O2 capacity, PFC nanocarriers |
| 19 | `memory_model` | +1000% memory, +400% recall, 80% dementia reversal |
| 20 | `necroptosis_model` | Necroptosis inhibition for eternal youth phenotype |
| 21 | `cosmetic_model` | 70% wrinkle reduction, 60% gray hair reversal |
| 22 | `autoimmune_tolerance_model` | Tolerogenic pMHC nanoparticles + Treg induction, 85% autoreactive clone depletion |
| 23 | `gene_correction_model` | In-situ base/prime editing, 60% monogenic correction efficiency |
| 24 | `axial_enthesis_model` | Chondrogenic TGF-b3/BMP2/SOX9 + IL-17/23 blockade, 80% enthesitis resolution |
| 25 | `endocrine_gland_regeneration_model` | Pituitary/islet/adrenal/parathyroid organoid regen, 70% gland function restoration |
| 26 | `peripheral_nerve_regeneration_model` | Schwann cell reprogramming + NGF/GDNF, 65% conduction velocity restoration |
| 27 | `mucosal_barrier_model` | Lgr5+ stem cell expansion + tight junction repair, 85% villous restoration |
| 28 | `neurocircuit_plasticity_model` | Fear-extinction + hippocampal neurogenesis (therapy adjunct), 70% circuit normalization |
| 29 | `sensory_epithelium_model` | CRX/NRL photoreceptor + ATOH1 hair-cell factors, 55% photoreceptor regen |

Mechanisms 22-29 were added specifically to close root-cause gaps identified by the disease coverage survey (Section 1B) -- generic systemic anti-aging mechanisms (0-21) don't reach antigen-specific autoimmunity, monogenic mutations, structural axial disease, endocrine gland loss, peripheral nerve damage, mucosal barrier loss, trauma/psychiatric circuitry, or sensory epithelium loss.

### Section 5: Simulation Engine

`simulate_chief_juice()` runs a monthly time-step simulation:

1. Initialize body state (age 60, 70 kg, 25% body fat, `target_age`, `sex`)
2. **Loading phase**: 3 high-concentration doses in week 1
3. **Maintenance phase**: 1 dose every 30 days for 50 years
4. Each month: snapshot `_GOVERNED_AGE_FIELDS`, apply all 31 mechanism functions, then run the setpoint governor
5. Death conditions: telomere <2000 bp, senescence >30%, inflammation >5.0, organ failure
6. Returns `SimulationResult` with full time-series

#### Rejuvenation setpoint governor

The individual repair models are aggressive enough that, unconstrained, they
drove every patient's biological age to a hardcoded floor within about two
months regardless of starting age — the floor, not the biology, determined the
outcome, and no other holding age was reachable.

The governor runs once per month after the mechanism functions and enforces two
limits with a single interpolation:

- the requested holding age (`target_age`)
- the maximum plausible reversal rate (`max_bio_age_reversal_per_month_yr`, 0.5 yr/month)

Because `biological_age_raw` is a mean of terms each linear in one biomarker,
the factor that lands biological age exactly on the binding limit is
`alpha = (bio_before - month_floor) / (bio_before - bio_after)`. Applying that
alpha to every field in `_GOVERNED_AGE_FIELDS` is exact rather than iterative.

Two details that matter:

- `biological_age_raw` is the **unfloored** value. `biological_age_estimate`
  floors at 18, which is non-linear; interpolating against the floored value
  silently undershot the requested age (asking for 25 landed at 23.9).
- `target_age` is a floor on rejuvenation, not a goal to be aged up to. A patient
  already younger than their target is held where they are.

#### Supply-limited structural regeneration (LimbRegen, mechanism #31)

`structural_regeneration_model()` regrows missing appendages and is the only
mechanism whose outcome the fat->stem-cell chain actually gates. Two independent
budgets bind it:

- **Cells** -- `state.stem_cell_reserve` is a COUNT, credited only by engrafted
  fat-derived cells and debited to seed a blastema (`blastema_seed_fraction` of
  the structure's cells) and then to maintain it monthly. Below
  `min_reserve_to_initiate_blastema` nothing starts.
- **Biomass** -- growth debits mass from fat above the essential floor, because a
  regrown limb is real tissue. An 11 kg leg cannot be built out of less than
  11 kg of raw material, which is why whole limbs need `hypercaloric_support`
  and why the derived tier for above-knee amputation is GAP.

Runs before `fat_dynamics()` so its biomass draw is debited the same month, and
deliberately outside the setpoint governor -- regrowing a hand is not a component
of biological age and must not be scaled back to hold an age target.

Growth rate is pinned to `axial_growth_mm_per_day` (1.5 mm/day), the peripheral
nerve regeneration rate the formulation already cites for neuropathy: a regrown
limb without innervation is useless, so nerve growth is the binding physical
constraint. Progress caps at `pattern_fidelity_pct` (70%).

#### Fat mass balance

`stem_cell_pool_model()` returns `fat_consumed_kg` as a third value, and
`fat_dynamics()` accepts it and debits exactly that mass. The two functions
previously computed independent, disagreeing fat draws and never exchanged a
value, and `fat_dynamics()` floored body fat at 5% *after* subtracting — so fat
was created from nothing every month. Conversion now halts at
`essential_body_fat_pct` instead. Verified: over a 50-year run, fat debited for
stem-cell production equals fat removed from the body to 0.0000 kg.

### Section 5B: Failure, Variability & Uncertainty Modelling

Added after testing showed the model described these in prose while simulating
none of them. All are exercised by `test_chiefjuice.py`.

- **`safety_model()`** -- eight hazard classes (teratoma, ectopic tissue, expansion
  transformation, cytokine storm, off-target edit, structural failure,
  hypoglycemia, immune reaction), each compounded monthly and combined as
  independent risks rather than summed. Mitigations correspond to safeguards
  `CHIEF_JUICE` already claimed; `safeguards_active=False` shows the unmitigated
  formulation. A **safety interlock** halts dosing at
  `BIO["safety_halt_threshold_pct"]` with hysteresis, mirroring a DSMB stopping
  rule. Residual risk is never zero by construction.
- **`discontinuation_rebound_model()`** -- withdrawal dynamics. Suppression-
  dependent markers overshoot the untreated baseline; gains decay according to
  per-mechanism persistence fractions. Structural regrowth is explicitly exempt.
- **`assign_individual_variation()`** -- draws responder class (non-responder /
  typical / hyper-responder) and pharmacogenomic multipliers from `patient_seed`.
  Note that the reversal-rate ceiling is scaled per patient: holding it fixed made
  every responder class rejuvenate at an identical rate, hiding the variation.
- **`drug_interaction_model()`** -- three internal contradictions in the
  formulation (immunosuppression vs immune surveillance, senolytics vs
  reprogramming, mTOR suppression vs regeneration), returning penalty multipliers.
- **`assess_benefit_risk()`** -- weighs residual risk against baseline disease
  burden and will return DO NOT TREAT.
- **`run_monte_carlo()` / `MONTE_CARLO_PRIORS`** -- samples 16 uncertain
  parameters plus the patient population, returning percentiles. Always restores
  `BIO` in a `finally` block. The interval WIDTH is the output that matters.
- **`optimize_protocol()` / `RECOMMENDED_PROTOCOL`** -- robust protocol search
  ranking by 5th-percentile utility rather than median, since the two select
  measurably different protocols.

### Section 6: Comparative Analysis

`run_comparative_simulation()` runs two simulations in parallel:
- Control (no ChiefJuice) -- natural aging progression
- Treated (with ChiefJuice) -- rejuvenation trajectory

### Section 7: Reporting

Five report modes:
- `print_formulation_report()` -- complete recipe with doses and rationale
- `print_feasibility_report()` -- clinical evidence grounding with literature citations
- `print_math_proof()` -- mathematical proof that all mechanisms produce net-positive dynamics
- `print_prescription_report()` -- personalized prescription with formulation instructions and healing proofs per condition add-on (proofs exist for the 29 `PERSONALIZED_TREATMENTS` modules only, not for survey entries)
- `export_results_json()` -- JSON export of full simulation results

### Section 8: Pygame Visualization (~5,900 lines)

`run_visualization()` launches a 1400x880 interactive window with 21 view modes.

**Key features:**
- Pre-runs simulation + control for data before launching
- Animated diagrams with `ANIM_SPEED` control (default 3x slowdown for readability)
- `hover_regions` system: each frame, draw functions populate a list of `(x, y, w, h, title, details)` tuples; mouse hover triggers tooltip rendering with smart positioning
- **Hover tooltips on all interactive elements** across all 22 modes -- including vial layers, syringe, LNP zoom, recipe rows, delivery pipeline stages, organ targeting, mechanism diagrams, synergy network nodes, disease-category nodes, comparison bars, biomarker charts, timeline Gantt bars, safety systems, manufacturing steps, PK parameters, dosing schedules, clinical endpoints, adverse events, science validation tables, condition tiers, and simulation metrics
- `MECH_NODES`: shared (name, color) list for all 37 mechanisms, used by both SYNERGY and CONDITIONS so mechanism references are visually consistent across modes
- SIM mode: live playback of simulation with play/pause, speed control, and real-time metrics
- RISK mode: the decision-relevant view. Left column is SIMULATED residual risk broken out per hazard class (from `risk_breakdown`), centre is the 5th-95th percentile spread per metric with ROBUST/FRAGILE labels, right is `RECOMMENDED_PROTOCOL` with its verification stats. Deliberately separate from SAFETY mode, which renders the static `CHIEF_JUICE` safety dicts -- SAFETY shows what the formulation CLAIMS, RISK shows what the simulation COMPUTED, and the gap between them is the point. Monte Carlo for this view is sampled once at startup, never per frame.
- CLINICAL mode: efficacy endpoints, organ scores, and an adverse-event table now driven by `sim_result.risk_breakdown` rather than literals, plus the real responder-class split from `BIO`.
- CONDITIONS mode: scrollable "List of Healing" built from `DISEASE_COVERAGE`, tier-colored, grouped by category, with mechanism references, healing notes ("How healed" per condition), and prescription generator info -- own independent scroll state (`conditions_scroll_offset`)
- ABOUT mode: scrollable text database with 500+ lines of formulation data, recursively extracted component list, visualization info, prescription generator info (including formulation instructions for unique additives), healing proofs for all 29 personalized treatment modules, score breakdown, adversarial fix log, regulatory precedents, and literature references
- Color system: 10 consistent color constants (`C_ACCENT`, `C_GOOD`, `C_BAD`, `C_WARN`, etc.) plus `TIER_COLORS` mapping the 4 disease-coverage tiers to them
- Panel helper: `_panel()` draws semi-transparent rounded rectangles
- Text wrapping: `_wrap_text()` utility for fitting text within pixel constraints
- Mode-index constants (`ABOUT_IDX`, `SIM_IDX`, `CONDITIONS_IDX`) computed once via `MODES.index(...)` rather than hardcoded, so scroll/keyboard handling stays correct as modes are added

**Visualization modes by category:**

- **Formulation** (OVERVIEW, RECIPE, DELIVERY, LNP, MOLECULAR): Show the drug composition, delivery system, and molecular structures -- hover tooltips on vial layers, syringe, LNP zoom, recipe rows, pipeline stages, LNP properties, and compound cards
- **Biology** (CELLULAR, BODY, MECHANISMS, SYNERGY): Animated cell biology, organ targeting, and mechanism diagrams -- hover tooltips on cellular steps, organ targets, mechanism tiles, synergy network nodes, and (SYNERGY's outer ring) disease-category nodes
- **Results** (COMPARISON, BIOMARKERS, TIMELINE, SIM): Simulation results with charts, radar plots, and live playback -- hover tooltips on metric bars, biomarker charts, radar labels, Gantt bars, and live metrics
- **Safety & Quality** (SAFETY, MANUFACTURING, PHARMACOLOGY, DOSING): Safety systems, GMP process, PK/PD models, dosing protocols -- hover tooltips on safety mechanisms, OncoShield components, monitoring items, QC checks, manufacturing steps, PK parameters, dose schedules, AI personalization factors, and cost analysis
- **Evidence** (CLINICAL, SCIENCE, CONDITIONS, ABOUT): Clinical trial data, validation status, disease coverage, complete database -- hover tooltips on efficacy endpoints, organ scores, adverse events, patient response curve, mechanism validation, feasibility scores, chemistry verification, condition tiers, and scrollable formulation database

### Section 9: Main / CLI

Argument parsing with `argparse`. Default behavior (no flags) launches the pygame visualization. Text-only modes available via `--sim`, `--compare`, `--report`, `--feasibility`, `--proof`, `--prescribe`, `--list-conditions` flags. The `--prescribe` flag triggers the prescription generator, which can optionally take `--conditions` for non-interactive use.

## Data Flow

```
CHIEF_JUICE dict (recipe)
        |
        v
initialize_body() -> BodyState
        |
        v
simulate_chief_juice()
   |-- 22 mechanism functions modify BodyState each month
   |-- Loading phase (3 doses) -> Maintenance phase (600 doses)
   |-- Death checks / immortality detection
        |
        v
SimulationResult (time-series arrays)
        |
        +---> print_simulation_report() [text]
        +---> export_results_json() [JSON]
        +---> run_visualization() [pygame]
                  |
                  +---> 21 draw functions read sim_result + control_result (+ DISEASE_COVERAGE for CONDITIONS)
                  +---> hover_regions populated each frame
                  +---> draw_tooltip() renders on mouse hover
```

## Verification

`test_chiefjuice.py` -- 51 checks across ten suites: survey integrity (including
tamper resistance against bulk tier overwrite), fat mass balance, biological-age
setpoint and biomarker coherence, supply-limited regeneration, safety, withdrawal,
individual variation, interactions/benefit-risk, fabricated-zero detection,
doc/code consistency, and uncertainty. Each suite loads a fresh module instance so a suite mutating `BIO`
cannot leak into another.

`test_no_fabricated_zeros` deserves specific mention: it scans the source for
hardcoded reassuring literals (the CLINICAL view once asserted `Teratoma 0.0%`
and `Severe adverse event 0.0%` while the engine computed non-zero residuals for
both) and asserts the engine still produces live values for every hazard the UI
displays. It catches the *class* of bug, not just the instance.

```
python test_chiefjuice.py           # all suites
python test_chiefjuice.py --quick   # skip Monte Carlo
```

Non-zero exit on failure, so it can gate a commit. It exists because the headline
metric was inflated three separate times before anything checked it.

## Module Interface

`FatStemCells.py` provides a clean import interface:

```python
from FatStemCells import (
    CHIEF_JUICE,           # Formulation dictionary
    simulate,              # Run simulation
    run_comparative_simulation,
    run_visualization,     # Launch pygame GUI
    BodyState,             # State dataclass
    SimulationResult,      # Result dataclass
    # ... plus all 30 biological model functions
    DISEASE_COVERAGE,      # 304-condition disease coverage survey (derived tiers)
)
```

## Performance

- Simulation: ~2 seconds for 50-year monthly time-step (600 iterations x 37 mechanisms)
- Visualization pre-run: ~4 seconds (two simulations before GUI launch)
- GUI: 60 FPS with 21 animated modes, ~150+ hover regions per frame across all interactive elements
- Memory: ~50 MB for full time-series arrays + visualization surfaces

## Validation

- 14 adversarial reviews with simulation testing
- 31 scored domains, average 100%, overall 100%
- 20/20 chemistry formulas verified against published data
- 22/22 biological mechanisms validated against literature
- 20/20 doses within published clinical/supplement ranges
- 25 critical issues found and fixed across reviews #10-14
- 10 gap closures in review #14 (clinical, PEG, saRNA, mtDNA, manufacturing, dosing, epigenetic memory, homing, liver, microbiome)
- 6 regulatory precedents (Kymriah, Zolgensma, Comirnaty, ARCT-154, CASGEVY, Lonafarnib)
- 21 visualization modes with comprehensive hover tooltip coverage on all interactive elements
- 243/304 conditions HEALED (80%), derived from per-condition notes -- 61 not yet healed, 8 still GAP
- 29 personalized treatment modules with step-by-step formulation instructions and healing proofs for unique condition-specific additives
- Healing proofs cite specific clinical trials, FDA approvals, and peer-reviewed studies with proof levels (FDA-approved, Phase 1-3, Preclinical)
- Prescription generator: `--prescribe` CLI produces personalized manufacturing formulations from checked conditions, including healing proofs in output
