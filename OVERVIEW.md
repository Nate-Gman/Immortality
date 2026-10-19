# ChiefJuiceV1.17 -- Technical Architecture Overview

## System Summary

ChiefJuiceV1.17 is a single-file Python application (17,453 lines) modelling an
injectable fat-to-stem-cell conversion therapy, paired with `test_chiefjuice.py`
(490 lines, 64 checks in 12 suites) and `FatStemCells.py` (library interface,
147 lines).

It combines a monthly time-step simulation engine (37 mechanism functions), a
22-mode pygame visualization (~6,400 lines), a 304-condition disease survey
with **derived** coverage tiers, a prescription generator, 29 personalized
treatment modules, an evidence audit of 219 graded parameters across 72
modules, and subsystems for safety, withdrawal, individual variation, drug
interactions, and uncertainty quantification.

**Design principle running through the codebase:** no headline number is stored
where it can be edited. Coverage tiers are derived from note text, the healed
count is an output, safety risk is computed rather than asserted, and the
verification suite fails if a fabricated figure reappears. This is a direct
response to the project's history — the headline metric was inflated three
separate times before anything checked it.

**Mechanism numbering:** all code references use 0-indexed mechanism IDs 0–36
(core #0–21 universal base; extension modules #22–36 prescribed per patient).
The same IDs appear in `DISEASE_COVERAGE` `mechanisms` lists and
`get_prescribed_mechanisms()`. Historical prose occasionally used #1–37.

## File & Section Map

| Section | Line range | Size | Contents |
|---|---|---|---|
| 1 | 164–2485 | ~2,320 | `CHIEF_JUICE` formulation dict — 89 top-level modules, 2,100+ nested keys |
| 1C | 2486–3328 | ~840 | `EVIDENCE_AUDIT` (219 entries) + `get_evidence_summary()` + `print_evidence_audit()` |
| 1B | 3329–4317 | ~990 | `DISEASE_COVERAGE` survey + `get_healing_progress()` + `HEALING_TARGET` |
| — | 4317–5146 | ~830 | `PERSONALIZED_TREATMENTS` (29 modules) + condition lookup + tier upgrades |
| — | 5146–5528 | ~380 | Derived-tier resolver (`_derive_tier`, `_derive_coverage`, `_derive_tractability`), integrity locks, `get_prescribed_mechanisms`, `EXTENSION_FIELD_PROFILE`, `apply_patient_profile` |
| 2 | 5528–6354 | ~830 | `BIO` dict — 472 biological constants |
| 3 | 6354–6862 | ~510 | `BodyState` dataclass + `SimulationResult` dataclass |
| 4 | 6863–8399 | ~1,540 | 37 mechanism functions + failure/variability models |
| 5 | 8400–9466 | ~1,070 | `simulate_chief_juice()` monthly engine + setpoint governor + `MONTE_CARLO_PRIORS`/`RECOMMENDED_PROTOCOL` + `run_monte_carlo`/`optimize_protocol` |
| 6 | 9467–9627 | ~160 | `run_comparative_simulation()` |
| 7 | 9628–10850 | ~1,220 | Report functions (formulation, feasibility, math-proof, prescription, roadmap, uncertainty, audit, JSON export) |
| 8 | 10851–17220 | ~6,370 | `run_visualization()` — 22-mode pygame GUI |
| 9 | 17221–17453 | ~230 | `main()` argparse CLI |

## Section 1: Formulation Constants (`CHIEF_JUICE` dict, lines 164–2485)

The complete recipe as a nested dictionary — 89 top-level modules, 2,100+
nested keys. Scalar metadata entries: `name` ("ChiefJuiceV1.17"), `version`
("1.17"), `tagline`, `mechanism` (one-line summary), `delivery` ("Subcutaneous
LNP suspension injection"), `storage` ("Lyophilized powder, reconstitute with
sterile saline"), `dose_ml` (10.0), `dose_frequency_days` (30),
`loading_doses` (3), `cost_per_dose_usd` (850), `shelf_life_months` (24).

### Module inventory (all 89 top-level keys)

Core formulation modules:

| Key | Contents |
|---|---|
| `reprogramming` | saRNA Yamanaka factors (OSKM, L-Myc/Lin28 swap), enhancer stack (18% base × 4.2 = 75.6% net), mRNA dosing schedule, reprogramming time |
| `adipostem_pro` | Modified Yamanaka cocktail, EV payload, TH5487 analog, TET2 activator, anti-reversion peptide, enhanced conversion efficiency, anti-reversion efficiency |
| `lnp` | 31 keys: composition, molar ratio, diameter/zeta potential/PDI, encapsulation %, N/P ratio, pKa, endosomal escape, PEG shedding, 6 tissue variants |
| `multi_tissue_targeting` | LNP variants uptake note, targeting strategy, dose split, off-target % |
| `growth_factors` | 17 keys: EGF, FGF2, VEGF, HGF, IGF1, PDGF, GCSF, BDNF + dose/PK detail |
| `telomerase` | TA65, cycloastragenol, resveratrol, astragaloside IV; telomere gain bp/division, activation fold |
| `senolytics` | Dasatinib, quercetin, fisetin, navitoclax; clearance efficiency, rounds/year, rotation |
| `metabolic` | NMN, NR, metformin, rapamycin, urolithin A, MitoQ, alpha-ketoglutarate |
| `dna_repair` | NAD+, poly-ADP-ribose, niacin, folate, B12; repair rate boost fold, pathway detail |
| `immune` | IL-7, IL-15, RANKL inhibitor, thymus regen factors; thymic regen rate, immune age reversal |
| `epigenetic` | Alpha-ketoglutarate, vitamin C, SAM modulator; epigenetic age reversal yr/dose |
| `exosomes` | 22 keys: source, dose, diameter, buoyant density, cargo, surface markers, biodistribution |
| `neuro` | 24 keys: BDNF/NGF/GDNF enhancers, dancing molecules, cell-permeable peptides, neural exosomes |
| `muscle` | 18 keys: myostatin inhibitor, myofibril packing, calcium handling, fiber-type shift, cross-bridge kinetics |
| `bone` | 20 keys: teriparatide/abaloparatide analogs, CCN3, alendronate, denosumab analog, sclerostin inhibitor, hydrogel |
| `hemo` | CRISPR hemoglobin editor, synthetic O2 carriers, HSC erythropoiesis activator; O2/RBC bounds |
| `proteoclear` | Autophagy boost, aggregate clearance, proteasome activity, amyloid/tau reduction |
| `nutrisense` | mTOR inhibition, insulin sensitivity, sirtuin, IGF1 reduction, AMPK activation |
| `commrepair` | Signaling repair %/mo, young plasma factor restoration, pro-aging factor removal |
| `oxidshield` | ROS reduction, oxidative damage repair, antioxidant capacity, lipid peroxidation |
| `hormonebalance` | DHEA/GH/sex-hormone/thyroid restoration, circadian rhythm |
| `plasmasync` | Plasma exchange frequency, young factor restoration, pro-aging clearance |
| `ecm` | Hyaluronic acid, collagen peptides, MMP inhibitor, LOX activator; collagen synthesis/decline |
| `memory` | CI994 HDAC inhibitor, methylene blue, J147, Dihexa, Noopept, donepezil analog, GL-II-73 |
| `lung` | Alveolar stem cell activator, HIF2 inhibitor, ECM hydrogel, RNAi, vascular optimizer |
| `immunoeternal` | 15 keys: CRISPR HSC editor, CCN3, CAR-T NK mRNA, RNAi autoimmunity silencer, HSCT mimetic, NK booster |
| `eternal_youth_pro` | Anti-necrotic compound, epitalon, CCN3, exercise mimetic, spermidine; necroptosis inhibition |
| `oncocure` | CRISPR oncogene corrector, mRNA vaccine payload, CAR-T NK, RNAi, angiogenesis inhibitor |
| `aging_clocks` | Primary/secondary/tertiary/backup epigenetic clocks + additional biomarkers |

15 healing-gap extension modules (#22–36), each with `evidence` field
distinguishing real precedent from speculative extrapolation:

| ID | Key | Contents |
|---|---|---|
| 22 | `tolerance` | Tolerogenic pMHC nanoparticles, Treg induction, autoreactive CAAR ablation, tTG blocker, islet tolerizing vaccine, microbiome consortium |
| 23 | `genefix` | Base/prime editor mRNA, tissue-specific guide library, hepatocyte LNP, HSC mobilizer, exon-skipping ASO; correction efficiency |
| 24 | `enthesis` | Chondrogenic factors, IL-17A/IL-23 dual blocker, localized sclerostin modulator, synovial fibrosis resolver, nerve sensitization reset, osteoclast resorption, disc chondrocyte regen |
| 25 | `endoregen` | Pituitary/islet organoid seeding, adrenal cortex + parathyroid regen factors, HPA/HPG/HPT axis recoupler |
| 26 | `nerveregen` | Schwann reprogramming factors, NGF/GDNF depot, remyelination accelerator, nodes of Ranvier repair, diabetic microvascular repair; axon regrowth mm/day |
| 27 | `mucosarepair` | Intestinal stem cell (Lgr5+) expansion, tight-junction peptide, gluten tolerizing nanoparticle, IL-13/IL-5 blocker, short-bowel adaptation |
| 28 | `circuitreset` | Fear-extinction enhancer, hippocampal neurogenesis booster, synaptic pruning normalizer, prefrontal-amygdala connectivity, DMN reset |
| 29 | `sensoryregen` | Photoreceptor progenitor (CRX/NRL), RPE sheet, cochlear hair-cell (ATOH1) factor, optic nerve regrowth, vestibular hair cell |
| 30 | `limbregen` | Blastema induction cocktail, positional memory factors, nerve-dependency substitute, scarless healing shifter, bioelectric patterning, scaffold, limb-field induction, ex-vivo bioreactor; `axial_growth_mm_per_day` = 1.5, pattern fidelity 70% |
| 31 | `fibroclear` | Activated fibroblast clearance (FAP-CAR-T), collagen crosslink dissolver, TGF-β3 dominance, antifibrotic core (pirfenidone), relaxin-2, architecture rescaffold, senolytic myofibroblast purge |
| 32 | `repeatedit` | Cas9-nickase contraction system, MSH3 knockdown, allele-selective ASO, RNA foci dissolver, methylation eraser |
| 33 | `reservoirpurge` | Provirus excision (CRISPR), cccDNA nuclease, latency reversal, reservoir-directed immunity, neuronal latency clearance, CCR5 ablation |
| 34 | `tissueforge` | Segment blastema inducer, luminal scaffold, enteric nervous system seeding, lymphatic valve patterning, vascular pedicle builder |
| 35 | `delivermax` | Engineered BBB capsid, focused ultrasound, HSC ex-vivo editing, liver LNP, immune-evasion coat, intrathecal convection; coverage ceilings (hepatic 90%, HSC 85%, CNS 35%, ganglia 20%) |
| 36 | `chromobalance` | XIST dosage silencer, dCas9-KRAB gene dosage, imprint-reset ASO, large-cargo locus replacement, trisomy selective pressure |

Safety, quality, platform & commercial modules:

| Key | Contents |
|---|---|
| `safety_system` | Tumor-suppressor boost, suicide gene, differentiation priority, immune surveillance, teratoma mitigation, monitoring, contraindications |
| `oncoshield` | p53 checkpoint, PTEN support, CRISPRi oncogene suppressor, chk1, shelterin modulator, epigenetic safeguard, biosensors, ctDNA probes |
| `teratoma_prevention` | 6 mechanisms + risk level |
| `cmyc_safety` | c-Myc oncogenicity problem + 5 solutions + preferred approach |
| `telomerase_cancer_safety` | Telomerase-reactivation cancer risk + 6 layers |
| `immune_evasion` | Autologous advantage, surface marker preservation, tolerance maintenance, anti-rejection protocol |
| `immune_memory` | Preserving adaptive memory during immune reset |
| `immune_paradox` | Immunosuppression-for-redosing vs surveillance-for-safety contradiction + solutions |
| `senolytic_resistance` | Resistance problem + rotation/backup solutions |
| `hormonal_cycling` | Hormone pulsatility problem + cycling solutions |
| `epigenetic_memory` | Somatic memory reversion problem + erasure solutions |
| `liver_sequestration` | LNP liver-uptake problem + sequestration reduction |
| `mtdna_heteroplasmy` | mtDNA heteroplasmy problem + reduction per treatment |
| `microbiome` / `microbiome_enhanced` | Microbiome disruption problem; precision probiotic/prebiotic/postbiotic/FMT |
| `stem_cell_homing` | Homing efficiency problem + solutions |
| `bbb_enhanced` | BBB penetration problem + 7 mechanisms |
| `stacking` | 15 pairwise synergy descriptors (telomere-epigenetic, epigenetic-senolytic, etc.) |
| `drug_interaction_enhanced` | PBPK model, AI prediction, TDM, CYP panel, staggered schedule, hotspots |
| `pkpd_enhanced` | PBPK + population PK + AI PK + multi-tissue PK + PD + adaptive PK |
| `monitoring` | 18 keys: biosensors, ctDNA, clock/telomere/organ/immune/neural/muscle monitoring |
| `personalized_medicine` | Age/gender/comorbidity adjustments, genetic profiling, AI dose engine |
| `stability_enhanced` | Drug holiday, adaptive frequency, clonal monitoring, checkpoints, biomarker-driven dosing |
| `manufacturing` / `manufacturing_enhanced` | Scale, production time, QC; process steps + PAT + continuous manufacturing |
| `formulation_process` | 15 keys: facility requirements + 11-step GMP process (saRNA synthesis → LNP → QC → fill-finish) |
| `cost_strategy` / `cost_enhanced` | Cost breakdown, reduction strategy, $850→$25 trajectory, insurance pathway, global access |
| `regulatory_pathway` / `regulatory_enhanced` | FDA/EMA/PMDA/NMPA/MHRA paths, RMAT, IND strategy, breakthrough therapy |
| `regulatory_precedents` | 6 precedents: Kymriah, Zolgensma, Comirnaty, ARCT-154, CASGEVY, Lonafarnib |
| `clinical_enhanced` / `clinical_roadmap` | Bridging data, first-in-human strategy, adaptive platform trial, go/no-go criteria |
| `verification_status` | Honest summary: what is grounded vs not established vs largest gaps |

### Section 1C: Evidence Audit (`EVIDENCE_AUDIT`, lines 2486–3328)

Every load-bearing parameter graded against published science — 219 entries
across 72 distinct modules, each with: `module`, `parameter`, `modeled` value,
`published` value, `grade` (A–D), `source` (citation), `note`.

| Grade | Meaning | Count |
|---|---|---|
| A | Replicated human clinical data | 30 |
| B | Preliminary human data | 61 |
| C | Animal or in-vitro data only | 59 |
| D | Speculative — no measured evidence | 69 |

`get_evidence_summary()` returns grade counts + speculative fraction.
`print_evidence_audit()` / `--audit` prints all 219 entries grouped by grade
with modeled vs published values, citations, and an honest verdict.

## Section 1B: Disease Coverage Survey (lines 3329–4317)

Per-condition cross-reference of all 37 mechanisms against 304 real
diseases/conditions in 10 categories (`DISEASE_COVERAGE_CATEGORY_ORDER`):
`your_conditions` (the owner's personally listed conditions from Problems.md),
`neuro_degenerative`, `genetic_congenital`, `autoimmune_inflammatory`,
`chronic_organ`, `cancer_hematologic`, `infectious_sequelae`,
`mental_neurodevelopmental`, `other_chronic`, `common_acute`.

Each entry fields:

| Field | Content |
|---|---|
| `name` | Condition name, verbatim from the source survey |
| `tier` | HEALED/MANAGED/PARTIAL/GAP — **derived at import time** by `_derive_tier()` from the entry's own `note`; never read from the literal |
| `claimed_tier` | Whatever the source literal asserted — preserved for audit |
| `tier_basis` | The specific phrase that decided the tier |
| `mechanisms` | 2–5 mechanism IDs (0–36) most responsible for the tier |
| `note` | Condition-specific healing note explaining how mechanisms reverse it |
| `precedent` | Optional `"clinical"` flag — explicit claim that an approved therapy achieves *this entry's* endpoint |
| `orig_tier` | Pre-upgrade tier, when a personalized module lifts it |

Tier semantics:

- `HEALED` — root-cause mechanism(s) with strong real-world precedent; durable
  reversal
- `MANAGED` — substantially controlled/slowed, real precedent, but full
  reversal of pre-existing structural damage uncertain
- `PARTIAL` — only part of pathology addressed, or enabling tech
  preclinical/unproven in humans
- `GAP` — no credible mechanism even in this model plausibly reverses the root
  cause

### The derivation resolver (lines 5146–5436)

`_derive_tier(note, mechanisms, precedent)` scans each note against four phrase
lists in precedence order:

1. `_GAP_NOTE_PHRASES` → GAP (explicit "model has no route" statements only —
   deliberately narrow; "there is no amputation plane" describes anatomy, not
   absence)
2. Hard hedge / unproven phrases → PARTIAL
3. Soft hedges ("extrapolat", "assumes", "at scale", "unsolved") → PARTIAL
   **unless** the entry carries `"precedent": "clinical"` (explicit checkable
   claim that an approved therapy achieves the entry's own endpoint — relieves
   soft hedges only)
4. `_MANAGED_NOTE_PHRASES` → MANAGED
5. Otherwise → HEALED

Two independent axes derived at import alongside tier:

- `_derive_coverage(note, mechanisms)` → COVERED if a mechanism is aimed at the
  root cause (a completeness claim about the formulation, not efficacy)
- `_derive_tractability(tier_basis, note)` → TRACTABLE if the remaining blocker
  matches a named experimental path (delivery platform, first-in-human trial,
  in-vivo translation, cell replacement, chromosome-scale cargo, morphogenesis,
  mechanistic stratification); FUNDAMENTAL reserved for blockers needing new
  physics or unrecoverable information

Resolved totals (`get_healing_progress()`):

```python
{'healed': 243, 'goal': 304, 'total_surveyed': 304, 'pct_of_goal': 79.93,
 'remaining': 61,
 'tier_counts': {'HEALED': 243, 'MANAGED': 16, 'PARTIAL': 37, 'GAP': 8},
 'split': {'common_acute': {'healed': 54, 'total': 57, 'pct': 94.7},
           'severe_chronic': {'healed': 189, 'total': 247, 'pct': 76.5}},
 'covered': 304, 'pct_covered': 100.0, 'uncovered': [],
 'tractable': 61, 'fundamental': [], 'pct_tractable': 100.0,
 'roadmap': {workstream: [conditions...]}}
```

`HEALING_TARGET` = `{'goal': 304, ...}` — aspirational; denominator tracks the
survey so `pct_of_goal` can't inflate by adding conditions.

The evidence split is the figure to quote: **54/57 (95%)** of common/acute
conditions earn HEALED honestly (mostly conditions medicine already cures) vs
**189/247 (77%)** of the original severe/chronic survey. The combined 243/304
(80%) headline moved 77% → 80% purely by enlarging the denominator — report the
split, not the headline.

### Integrity history (why derivation exists)

The healed count was gamed three times before it became derived:

1. A `PERSONALIZED_TREATMENTS` entry with a huge `conditions` list (a real
   per-condition treatment doesn't apply to 15+ unrelated diseases)
2. Bulk-overwriting `DISEASE_COVERAGE` literal `tier` values, leaving
   self-contradictory data (`tier=HEALED`, `note="no real reversal mechanism"`)
3. Bulk-extending the allowlist added to block #2 — silently turning the
   safeguard into dead code

Lesson encoded: any safeguard keyed on a field that can be bulk-edited will be
bulk-edited. `tier` is no longer stored at all. An import-time assertion also
pins the expected count of `"precedent": "clinical"` flags (25) — the one
remaining data edit that can raise HEALED without changing biology.

**Scope limit:** the engine does not model individual conditions — it models
generic aging biology plus extension modules. Two patients with different
diagnoses produce nearly identical trajectories; a tier is a reasoned
literature claim, not a simulation output.

## Section — Personalized Treatments (lines 4317–5146)

`PERSONALIZED_TREATMENTS`: 29 condition-specific targeted formulation modules
added ON TOP of the base regimen. Each contains:

| Field | Content |
|---|---|
| `conditions` | Which survey conditions it applies to |
| `label` | Module name (e.g. "PrionClear", "AlzClear", "MSRemyelinate") |
| `components` | Unique additives with dose, unit, biological role |
| `target` | What the formulation achieves |
| `tier_upgrade` | e.g. `GAP_TO_MANAGED`, `PARTIAL_TO_HEALED` — applied at runtime, original kept in `orig_tier` |
| `evidence` | Supporting clinical/preclinical evidence |
| `formulation_method` | Step-by-step manufacturing: expression system (E. coli/CHO/HEK293T), purification (Ni-NTA/Protein A/SEC/HPLC), LNP variant, encapsulation, delivery route, dosing schedule |

`get_personalized_treatments(patient_conditions)` maps condition names →
treatment keys via a reverse lookup built at import.

## Section — Prescription logic (lines 5454–5528)

`get_prescribed_mechanisms(patient_conditions)` — mechanisms #0–21 are the
universal base (every patient); #22–36 are prescribed only when a matching
diagnosis exists in `DISEASE_COVERAGE` (case-insensitive exact match on survey
names). A healthy patient receives #0–21 only.

`PROJECT_OWNER_CONDITIONS` — default patient profile for the visualization:
`["Ankylosing spondylitis", "Celiac disease", "Type 1 diabetes", "Neuropathy",
"PTSD", "Pituitary gland disorder"]`.

`EXTENSION_FIELD_PROFILE` maps mechanism ID → `[(field, healthy_baseline,
disease_baseline)]` — a patient without the diagnosis sits at the healthy
baseline; a diagnosed patient starts at a realistic disease baseline (e.g.
ankylosing spondylitis tracks the Problems.md severity). Module 35 (DeliverMax)
and 23 (GeneFix) have no baseline split (coverage is a platform property;
`gene_correction_pct` starts at 0 either way); module 30 (LimbRegen) state
lives in `missing_structures` from the `regrow` argument.

`apply_patient_profile(state, prescribed)` writes those baselines.

## Section 1D: Prescription Generator

`generate_prescription()` / `print_prescription_report()` combine all
formulation components into a personalized manufacturing manifest:

1. Patient profile — matched conditions from the 304-condition checklist
2. Condition healing details — tier, mechanisms, healing note per condition
3. Base formulation — mechanisms #0–21, always included
4. Extension modules — #22–36, prescribed per patient conditions
5. Personalized add-ons — condition-specific additives with full formulation
   instructions
6. Complete manufacturing manifest — all components with doses, units, roles

Healing proofs exist for the 29 `PERSONALIZED_TREATMENTS` modules only (they
cite specific trials/FDA approvals with proof levels FDA-approved / Phase 1–3 /
Preclinical), not for survey entries.

## Section 2: Biological Constants (`BIO` dict, lines 5528–6354)

472 published biological parameters used by the simulation models, organized by
body system:

- **Adipose tissue**: adipocyte density (1e6/g), total adipocytes, SVF stem
  cells, fat mobilization rate, engraftment fraction, essential fat floor (5%)
- **Cell division kinetics**: doubling time, reprogramming time, Hayflick
  limit, telomerase-extended divisions
- **Telomere biology**: newborn/adult/aged lengths, shortening rates,
  critical/crisis thresholds, shelterin proteins, telomerase gain
- **DNA damage & repair**: lesion rates by type, 6 repair pathways, mutation
  rates, repair efficiency
- **Senescence**: markers (p16, p21, SA-β-gal, SASP), burden by age, lethal
  threshold, resistance rate
- **Mitochondrial**: mtDNA (16,569 bp), 37 genes, OXPHOS complexes, ATP,
  heteroplasmy threshold, mitophagy, ROS
- **Stem cells**: 7 types with descriptions, pool kinetics, quiescence
  fraction, exhaustion rate
- **Immune**: total immune cells (1.8e12), T/B/NK/macrophage/neutrophil counts,
  Treg fractions, thymus mass by age, naive T-cell decline
- **Inflammation**: CRP/IL-6/TNF-α by age, inflammaging rate, NF-κB, NRF2
  decline
- **Organ function**: heart/liver/kidney/brain/lung/skin/pancreas/muscle/bone
  stats
- **Cancer statistics**: lifetime risk, incidence by age, mortality, telomerase
  reactivation, ctDNA sensitivity by stage
- **Pharmacokinetics**: LNP biodistribution, saRNA/mRNA expression windows,
  drug half-lives for 15+ compounds
- **Nutrition**: RDAs for 11 nutrients, caloric needs, regeneration demands
- **Growth factor & drug PK**: half-lives (EGF, FGF2, VEGF, IGF1, BDNF, G-CSF,
  metformin, rapamycin, dasatinib, navitoclax, quercetin, fisetin, ...)
- **Exosome biology**: size, density, cargo, surface markers, biogenesis,
  uptake, biodistribution
- **LNP biophysics**: diameter, zeta potential, PDI, encapsulation, N/P, pKa,
  endosomal escape, PEG shedding
- **Aging rates**: Gompertz mortality doubling, acceleration by decade, max
  documented lifespan (122 y), healthspan-lifespan gap
- **Regenerative capacity**: turnover by tissue (epidermis 28d, intestinal 5d,
  liver 30d, peripheral nerve 1 mm/day, fracture 8 wk)
- **12 hallmarks of aging** (9 original + 3 added 2023)
- **6 epigenetic clocks** with CpG counts (Horvath 353, GrimAge ~1030,
  PhenoAge 513, DunedinPACE 173)
- **Rejuvenation setpoint control**: `max_bio_age_reversal_per_month_yr` = 0.5,
  deadband, loading phase, `min_biological_age_floor` = 18
- **Regenerative substrate**: cells per pct point, nutrient floors, cofactor
  repletion rates, `biomass_conversion_efficiency` = 0.45
- **Adiposity as liability**: inflammation threshold, obesity penalties
- **Safety hazards**: 8 per-month hazard rates (`hazard_*`: teratoma 0.35,
  ectopic 0.28, transformation 0.40, cytokine storm 0.30, off-target 0.22,
  structural 0.25, hypoglycemia 0.18, immune reaction 0.20), 17 mitigation
  efficacies (`mitig_*`: 0.55–0.90), `safety_halt_threshold_pct` = 5
- **Discontinuation & rebound**: overshoot multipliers (senescence ×1.6,
  inflammation ×1.4), persistence fractions by mechanism type
- **Individual variation**: `responder_sd` 0.22, `nonresponder_fraction` 0.14
  (efficacy ×0.35), `hyperresponder_fraction` 0.08 (efficacy ×1.45, risk
  ×1.9), pgx SDs (clearance 0.30, immune reactivity 0.35, regen 0.25)
- **Drug interactions**: `interaction_immunosupp_surveillance` 0.55,
  `interaction_senolytic_reprogramming` 0.18,
  `interaction_mtor_regeneration` 0.22, nutrient demand multipliers
- **Structural regeneration**: `blastema_seed_fraction` 0.002,
  `min_reserve_to_initiate_blastema`, appendage anatomy (`appendages`: 12
  structures, (mass_kg, axial_length_cm) for a 70 kg adult — fingertip
  0.005/1.5, finger 0.030/8.0, thumb 0.045/6.5, hand 0.400/19.0,
  forearm_and_hand 1.600/45.0, arm 3.500/73.0, toe 0.035/5.0, foot 1.000/26.0,
  lower_leg_and_foot 4.500/48.0, leg 11.000/85.0, ear 0.015/6.0, nose
  0.020/4.0)
- **Supply-chain accelerators**: ex-vivo expansion, local delivery, scaffold
  retention, anabolic support
- **Female reproductive aging**: menopause onset 51, transition 5 yr,
  sex-hormone drop 35%
- **Monte Carlo / protocol**: see `MONTE_CARLO_PRIORS` and
  `RECOMMENDED_PROTOCOL` in Section 5

## Section 3: Data Structures (lines 6354–6862)

### `BodyState` (~143 dataclass fields)

Complete simulated-body state at a point in time. Grouped:

- **Identity/demographics**: `age_years`, `epigenetic_age_years`,
  `body_weight_kg`, `body_fat_pct`, `sex`
- **Core biomarkers**: `telomere_length_bp`, `dna_damage_lesions`,
  `senescent_burden_pct`, `mitochondrial_function_pct`, `stem_cell_pool_pct`,
  `inflammation_level`, `immune_function_pct`, `naive_t_cells`,
  `thymus_function_pct`
- **10 organ scores**: `heart`, `liver`, `kidneys`, `brain`, `lungs`, `skin`,
  `muscles`, `bones`, `pancreas`, `immune_system`
- **Treatment state**: `doses_administered`, `months_on_chief_juice`,
  `fat_mobilized_kg`, `stem_cells_generated`, `senescent_cells_cleared`,
  `telomere_restored_bp`, `dna_lesions_repaired`, `treatment_intensity`
- **Cancer**: `cancer_risk_pct`, `cancer_detected`, `cancer_cured`,
  `ctDNA_level`, `oncogene_suppression_pct`
- **Neural**: `cognitive_function_pct`, `neurogenesis_rate_pct`,
  `myelination_pct`, `paralysis_recovery_pct`, `mental_health_score`,
  `neurotransmitter_balance_pct`
- **Muscle/bone**: `specific_force_n_per_cm2`, `muscle_mass_index`,
  `muscle_strength_gain_pct`, `bone_strength_mpa`, `bone_density_pct`,
  `bone_strength_gain_pct`
- **Memory**: `memory_capacity_pct`, `recall_speed_pct`, `dementia_risk_pct`
- **Lung/oxygen**: `gas_exchange_efficiency_pct`, `lung_fibrosis_pct`,
  `lung_disease_resistance_pct`, `o2_per_rbc_billion`, `rbc_lifespan_days`,
  `oxygen_capacity_pct`
- **Necroptosis/cosmetic**: `necroptosis_rate_pct`, `wrinkle_score`,
  `gray_hair_pct`
- **Stacking**: `stacking_multiplier`
- **Proteostasis**: `autophagy_rate_pct`, `proteasome_activity_pct`,
  `protein_aggregates_pct`
- **Nutrient sensing**: `insulin_sensitivity_pct`, `mtor_activity_pct`,
  `sirtuin_activity_pct`, `metabolic_efficiency_pct`
- **Communication**: `signaling_efficiency_pct`, `pro_aging_factors_pct`
- **Oxidative stress**: `ros_level_pct`, `antioxidant_capacity_pct`,
  `lipid_peroxidation_pct`
- **Hormonal**: `dhea_level_pct`, `gh_level_pct`, `sex_hormone_balance_pct`,
  `thyroid_function_pct`, `circadian_rhythm_pct`
- **Plasma exchange**: `plasma_young_factors_pct`, `plasma_pro_aging_pct`
- **Extension module fields** (#22–36): `autoimmune_tolerance_pct`,
  `autoreactive_burden_pct`, `gene_correction_pct`,
  `enthesis_inflammation_pct`, `axial_bone_pathology_pct`,
  `joint_mobility_pct`, `endocrine_gland_function_pct`,
  `peripheral_nerve_conduction_pct`, `mucosal_barrier_integrity_pct`,
  `neurocircuit_balance_pct`, `sensory_epithelium_pct`, `missing_structures`
  dict, `blastema_established` dict, `fibrosis_burden_pct`,
  `architecture_integrity_pct`, `repeat_expansion_burden_pct`,
  `toxic_transcript_pct`, `latent_reservoir_pct`,
  `absent_tissue_regrowth_pct`, `tissue_functional_integration_pct`,
  `delivery_coverage_hepatic_pct`, `delivery_coverage_hsc_pct`,
  `delivery_coverage_cns_pct`, `delivery_redosing_available`,
  `chromosomal_dosage_error_pct`
- **Setpoint**: `target_biological_age` (default 18.0 — a floor on
  rejuvenation, not a goal), `treatment_intensity` (hysteresis)
- **Deployable reserve**: `stem_cell_reserve` (a COUNT credited only by
  fat→stem-cell conversion), `stem_cells_consumed_by_regen`,
  `regen_stalled_months`
- **Nutrients**: `nutrient_protein_pct`, `nutrient_vitamin_c_pct`,
  `nutrient_iron_pct`, `nutrient_zinc_pct`, `nutrient_calcium_pct`,
  `nutrient_b12_folate_pct`, `nutrient_support`, `accelerated_regen`,
  `dose_cycle_on_months`, `dose_cycle_off_months`,
  `biomarkers_capped_this_month`
- **Discontinuation**: `discontinued_at_month`, `months_since_discontinuation`
- **Individual variation**: `response_multiplier`, `clearance_multiplier`,
  `immune_reactivity_multiplier`, `regen_capacity_multiplier`,
  `adverse_risk_multiplier`, `responder_class`, `active_interactions`
- **Safety**: `risk_teratoma_pct`, `risk_ectopic_tissue_pct`,
  `risk_transformation_pct`, `risk_cytokine_storm_pct`,
  `risk_off_target_edit_pct`, `risk_structural_failure_pct`,
  `risk_hypoglycemia_pct`, `risk_immune_reaction_pct`, `safety_halted`,
  `safety_halt_months`, `adverse_events`
- **Flags**: `alive`, `cause_of_death`

Computed properties:

| Property | Definition |
|---|---|
| `cumulative_severe_risk_pct` | `1 − Π(1 − rᵢ/100)` over the 8 hazard risks — independent survivals, never a sum |
| `dominant_risk` | Hazard class contributing most |
| `nutrient_adequacy_pct` | `min` of the 6 nutrient levels — Liebig's law (scarcest caps regeneration) |
| `limiting_nutrient` | Which nutrient caps |
| `organ_average` | Mean of the 10 organ scores |
| `biological_age_estimate` | `max(biological_age_raw, min_biological_age_floor)` — floored at 18, non-linear; do not interpolate against it |
| `biological_age_raw` | Unfloored mean of **15 component ages**: telomere, senescence, mito, stem-cell, inflammation, organ, cognitive, memory, lung, muscle, bone, proteostasis, ROS, DHEA, signaling — each linear in one biomarker |

`check_alive()`: death on senescent overload (≥30%), systemic inflammation
(≥5.0), telomere crisis (<2000 bp), organ failure, stem-cell critical.

### `SimulationResult`

Time-series lists for all biomarkers (months, chronological/epigenetic/
biological age, telomere, senescence, mito, stem pool, inflammation,
organ_average, immune, body_fat, doses, alive, cancer_risk, cognitive,
myelination, mental_health, specific_force, bone_strength, memory,
gas_exchange, o2_per_rbc, necroptosis, stacking, paralysis, wrinkles,
gray_hair, autophagy, proteasome, aggregates, insulin_sensitivity, mtor,
sirtuin, metabolic_efficiency, signaling, pro_aging, ros, antioxidant,
lipid_peroxidation, dhea, gh, sex_hormone, thyroid, circadian,
plasma_young/pro_aging) plus: `final_state`, `total_doses`,
`total_stem_cells`, `total_senescent_cleared`, `total_fat_used_kg`,
`peak/lowest_biological_age`, `years_gained`, `immortal`, `regen_progress`
(monthly {appendage: pct}), `regen_completed_month`, `stem_cell_reserve`,
`total_cells_consumed_by_regen`, `regen_stalled_months`,
`cumulative_severe_risk_pct`, `safety_halt_months`, `adverse_events`,
`risk_breakdown` (per-hazard-class totals for the RISK view).

## Section 4: Biological Models (lines 6863–8399)

Each mechanism is an independent function that modifies `BodyState`. Effect
magnitudes are **modeled values** (post-correction), not measured outcomes;
evidence grades are in `EVIDENCE_AUDIT`.

| # | Function | Effect (modeled) |
|---|---|---|
| 0 | `telomere_dynamics` | +10 bp/division telomerase gain (maintenance-plus; never demonstrated in humans) |
| 1 | `senescence_model` | Clears p16/p21+ cells at 50%/round senolytic efficacy (modeled; mouse data) |
| 2 | `fat_dynamics` / `stem_cell_pool_model` | Adipocyte → stem-cell conversion (75.6% modeled vs ≤1% published ex-vivo; grade D); returns `fat_consumed_kg` for mass balance |
| 3 | `mitochondrial_model` | 5%/mo biogenesis + 2× mitophagy (urolithin A) |
| 4 | `epigenetic_age_model` | −0.25 yr/dose methylation reversal (TRIIM-anchored optimistic) |
| 5 | `compute_stacking_multiplier` | Multiplicative synergy across mechanisms (invented; grade D) |
| 6 | `dna_damage_model` | 1.3× repair rate, 15% mutation reduction (NAD+/PARP) |
| 7 | `immune_rejuvenation_model` | Thymic regrowth 2%/mo, +3 yr immune age reversal (TRIIM-anchored) |
| 8 | `cancer_risk_model` | Multi-layer safeguards, 25% risk reduction (modeled; not 95%) |
| 9 | `neural_model` | BDNF/NGF/GDNF, +2%/mo neurogenesis (debated in adults) |
| 10 | `muscle_strength_model` | Myostatin inhibition, +25% specific force (50 N/cm², within human range) |
| 11 | `bone_strength_model` | Osteo-anabolic restoration within physical limits |
| 12 | `proteostasis_model` | 2× autophagy, Hsp70 (spermidine/rapamycin-class) |
| 13 | `nutrient_sensing_model` | mTOR inhibition, AMPK 2.5×, sirtuin +50% |
| 14 | `communication_model` | Exosome signaling repair, 5%/mo (investigational) |
| 15 | `oxidative_stress_model` | 30% ROS reduction (MitoQ/NAC; antioxidant RCTs mixed) |
| 16 | `hormonal_model` | DHEA/GH/thyroid/melatonin partial restoration (70–80%, sex-adjusted) |
| 17 | `lung_model` | +20% gas exchange (no therapy achieves this in fibrotic lung) |
| 18 | `oxygen_capacity_model` | +10% O2 capacity (Hb capped at 4 O2/Hb ≈ 1.1e9 O2/RBC) |
| 19 | `memory_model` | +20% memory support, 10% dementia reversal (no therapy reverses dementia) |
| 20 | `necroptosis_model` | 40% necroptosis inhibition (preclinical; tumor-suppressive pathway) |
| 21 | `cosmetic_model` | 40% wrinkle reduction, 30% gray hair reversal (no injectable achieves this) |
| 22 | `autoimmune_tolerance_model` | Tolerogenic pMHC NPs + Treg induction, 85% autoreactive clone depletion (CAAR-T preclinical) |
| 23 | `gene_correction_model` | In-situ base/prime editing, 60% monogenic correction (liver-grade; most tissues far lower) |
| 24 | `axial_enthesis_model` | IL-17/23 blockade + chondrogenesis, 80% enthesitis resolution; fusion resorption 75% (grade D — no human precedent) |
| 25 | `endocrine_gland_regeneration_model` | Pituitary/islet/adrenal/parathyroid organoid regen, 70% function (VX-880 islet precedent is grade B) |
| 26 | `peripheral_nerve_regeneration_model` | Schwann reprogramming + NGF/GDNF, 2 mm/day axon regrowth, 65% conduction restoration |
| 27 | `mucosal_barrier_model` | Lgr5+ expansion + tight-junction repair, 85% villous restoration (TAK-101 precedent; in-vivo regrowth grade D) |
| 28 | `neurocircuit_plasticity_model` | Fear-extinction + neurogenesis adjunct, 70% circuit normalization (MDMA/psilocybin therapy data; not standalone) |
| 29 | `sensory_epithelium_model` | CRX/NRL + ATOH1, 55% photoreceptor/50% hair-cell regen (mammalian regen at scale unproven) |
| 30 | `structural_regeneration_model` | LimbRegen blastema induction, 1.5 mm/day, 70% pattern fidelity (Xenopus precedent; no mammal regenerates) — see §5 |
| 31 | `fibrosis_reversal_model` | FibroClear FAP-CAR-T + pirfenidone, 75% scar resorption (mouse cardiac fibrosis precedent) |
| 32 | `repeat_expansion_model` | Cas9-nickase CAG contraction, 65% repeat reduction (patient cells; never in living brain); 35% residual floor |
| 33 | `latent_reservoir_model` | CRISPR provirus excision, 70% clearance (EBT-101 Phase 1; no sustained remission); ~3% residual |
| 34 | `internal_tissue_regrowth_model` | TissueForge scaffold intestine regen, 65% segment regrowth (rodent; no human segment regrown); functional integration lags anatomy |
| 35 | `delivery_platform_model` | DeliverMax coverage ceilings — liver 90%, HSC 85%, CNS 35%, ganglia 20% (encodes real delivery limits) |
| 36 | `chromosomal_dosage_model` | XIST silencing, 60% dosage normalization (in vitro; adult neuroarchitecture persists) |

Supporting/failure models in the same span: `plasma_exchange_model` (30%
pro-aging clearance; mouse parabiosis), `inflammation_model`,
`organ_repair_model`, `safety_model`, `discontinuation_rebound_model`,
`assign_individual_variation`, `drug_interaction_model`,
`assess_benefit_risk`, `nutrient_dynamics_model`,
`regeneration_supply_factor`.

## Section 5: Simulation Engine (lines 8400–9466)

`simulate_chief_juice()` runs a monthly time-step simulation:

1. `initialize_body(age, weight, fat_pct, patient_conditions, target_age, sex,
   regrow)` — builds BodyState (default 60 y, 70 kg, 25% fat), applies
   `EXTENSION_FIELD_PROFILE` disease baselines for prescribed modules, sets
   sex-specific hormonal baselines, populates `missing_structures` from
   `regrow`, draws individual variation from `patient_seed`
2. **Loading phase**: `loading_doses` = 3 high-concentration doses in week 1
3. **Maintenance phase**: 1 dose per `maintenance_interval_days` = 30 for the
   horizon, subject to `dose_cycle` on/off months and the safety interlock
4. Each month: snapshot `_GOVERNED_AGE_FIELDS` → apply all 37 mechanism
   functions → `nutrient_dynamics_model` → `drug_interaction_model` →
   `safety_model` → `rebalance_biomarkers` → setpoint governor → rebound if
   discontinued → death checks → record time-series
5. Death conditions: telomere <2000 bp, senescence ≥30%, inflammation ≥5.0,
   organ failure, stem-cell critical
6. Returns `SimulationResult` with full time-series

### Rejuvenation setpoint governor

The individual repair models are aggressive enough that, unconstrained, they
drove every patient's biological age to a hardcoded floor within ~2 months
regardless of starting age — the floor, not the biology, determined the
outcome, and no other holding age was reachable.

The governor runs once per month after the mechanism functions and enforces
two limits with a single interpolation:

- the requested holding age (`target_age`)
- the maximum plausible reversal rate (`max_bio_age_reversal_per_month_yr` =
  0.5, scaled by the patient's responder multiplier)

Because `biological_age_raw` is a mean of terms each linear in one biomarker,
the factor that lands biological age exactly on the binding limit is

```
alpha = (bio_before - month_floor) / (bio_before - bio_after)
```

Applying `alpha` to every field in `_GOVERNED_AGE_FIELDS` (24 fields: the 14
non-organ markers + 10 organ scores) is exact rather than iterative.

Two details that matter:

- `biological_age_raw` is the **unfloored** value. `biological_age_estimate`
  floors at 18, which is non-linear; interpolating against the floored value
  silently undershot the requested age (asking for 25 landed at 23.9).
- `target_age` is a floor on rejuvenation, not a goal to be aged up to. A
  patient already younger than their target is held where they are.

`rebalance_biomarkers()` additionally clamps each governed field at its
target-age ideal (`_biomarker_value_at_age` closed-form inverses), because the
15-component mean could otherwise be hit by incoherent components — in
testing, bone/muscle were driven to ≈ −21 years subsidising inflamed organs.
`biomarkers_capped_this_month` counts caps for diagnostics.

### Supply-limited structural regeneration (LimbRegen, mechanism #30)

`structural_regeneration_model()` regrows missing appendages and is the only
mechanism whose outcome the fat→stem-cell chain actually gates. Two
independent budgets bind it:

- **Cells** — `stem_cell_reserve` is a COUNT, credited only by engrafted
  fat-derived cells and debited to seed a blastema (`blastema_seed_fraction`
  of the structure's cells) and then to maintain it monthly. Below
  `min_reserve_to_initiate_blastema` nothing starts.
- **Biomass** — growth debits mass from fat above the essential floor via
  `biomass_conversion_efficiency` (0.45), because regrown tissue is real
  matter. An 11 kg leg cannot be built out of less than 11 kg of raw material
  — why whole limbs need `hypercaloric_support` and above-knee amputation
  derives only PARTIAL.

Runs before `fat_dynamics()` so its biomass draw is debited the same month,
and deliberately outside the setpoint governor — regrowing a hand is not a
component of biological age and must not be scaled back to hold an age
target.

Growth rate is pinned to `axial_growth_mm_per_day` (1.5 mm/day) — the
peripheral nerve regeneration rate the formulation already cites for
neuropathy; a regrown limb without innervation is useless, so nerve growth is
the binding physical constraint. Progress caps at `pattern_fidelity_pct` =
70%. `regen_stalled_months` counts months where supply could not fund work.

### Fat mass balance

`stem_cell_pool_model()` returns `fat_consumed_kg` as a third value, and
`fat_dynamics()` accepts it and debits exactly that mass. The two functions
previously computed independent, disagreeing fat draws and never exchanged a
value, and `fat_dynamics()` floored body fat at 5% *after* subtracting — so
fat was created from nothing every month. Conversion now halts at
`essential_body_fat_pct` instead. Verified: over a 50-year run, fat debited
for stem-cell production equals fat removed from the body to 0.0000 kg.

### Failure, variability & uncertainty subsystems

- **`safety_model()`** — 8 hazard classes compounded monthly, mitigated by
  safeguards `CHIEF_JUICE` already claimed (17 `mitig_*` efficacies);
  `safeguards_active=False` shows the unmitigated formulation. Safety
  interlock halts dosing at `safety_halt_threshold_pct` = 5% cumulative with
  hysteresis, mirroring a DSMB stopping rule. Residual risk never zero by
  construction; `risk_breakdown` tracks per-class totals for the RISK view.
- **`discontinuation_rebound_model()`** — withdrawal dynamics.
  Suppression-dependent markers overshoot untreated baseline (senescence
  ×1.6, inflammation ×1.4); gains decay per per-mechanism persistence
  fractions. Structural regrowth explicitly exempt.
- **`assign_individual_variation()`** — draws responder class
  (non-responder 14%/×0.35, typical, hyper-responder 8%/×1.45/risk×1.9) and
  pharmacogenomic multipliers from `patient_seed`. Reversal-rate ceiling is
  scaled per patient — a fixed ceiling made every responder class rejuvenate
  identically.
- **`drug_interaction_model()`** — three internal contradictions
  (immunosuppression 0.55 vs surveillance, senolytics 0.18 vs reprogramming,
  mTOR 0.22 vs regeneration), returning penalty multipliers in
  `active_interactions`.
- **`assess_benefit_risk()`** — weighs residual risk against baseline disease
  burden; returns TREAT / CAUTION / DO NOT TREAT.
- **`nutrient_dynamics_model()` + `regeneration_supply_factor()`** — draws
  down the 6 nutrient fields under regenerative load; `nutrient_adequacy_pct`
  (Liebig minimum) caps throughput; supplementation repletes between cycles.
- **`run_monte_carlo()` / `MONTE_CARLO_PRIORS`** — samples 16 uncertain
  parameters (fold-ranges 0.3–4.0 around BIO) plus the patient population,
  returns percentiles; restores `BIO` in a `finally`. Interval WIDTH is the
  output that matters. Priors: `adipocytes_per_gram` 0.5–2.0,
  `stem_cell_engraftment_pct` 0.3–3.0, `ex_vivo_expansion_fold` 0.3–3.0,
  `blastema_seed_fraction` 0.4–2.5,
  `blastema_maintenance_fraction_per_month` 0.5–2.0,
  `regen_cells_per_pct_point` 0.4–2.5,
  `max_bio_age_reversal_per_month_yr` 0.5–1.5, `adipose_inflammation_fold`
  0.5–2.0, `hazard_teratoma_per_month_pct` 0.3–4.0,
  `hazard_expansion_transformation_pct` 0.3–4.0,
  `hazard_ectopic_tissue_per_month_pct` 0.3–4.0,
  `hazard_immune_reaction_per_month_pct` 0.3–4.0,
  `mitig_icasp9_suicide_switch` 0.7–1.15, `mitig_clonal_release_testing`
  0.7–1.15, `interaction_immunosupp_surveillance` 0.5–1.6,
  `biomass_conversion_efficiency` 0.6–1.6.
- **`optimize_protocol()` / `RECOMMENDED_PROTOCOL`** — robust protocol search
  ranking by 5th-percentile utility rather than median, since the two select
  measurably different protocols. Stored result: 22% body fat, 2-on/1-off
  dosing, supplemented, non-accelerated; verified 1,000 runs (quality p5
  71.4/median 84.5/p95 89.4; risk p5 0.238%/median 0.654%/p95 1.83%; 3.9% of
  runs exceed 2% risk). `accelerate_if`: risk weight < 7.9.

## Section 6: Comparative Analysis (lines 9467–9627)

`run_comparative_simulation()` runs two simulations in parallel — control
(natural aging, no ChiefJuice) vs treated — for the COMPARISON visualization
mode and `--compare` CLI.

## Section 7: Reporting (lines 9628–10850)

| Function | Output |
|---|---|
| `print_formulation_report()` | Complete recipe with doses and rationale (`--report`) |
| `print_feasibility_report()` | Clinical evidence grounding with literature citations (`--feasibility`) |
| `print_math_proof()` | Arithmetic-consistency check (consistent math ≠ biological proof) (`--proof`) |
| `print_evidence_audit()` | All 219 audit entries by grade (`--audit`) |
| `print_prescription_report()` | Personalized prescription manifest + healing proofs (`--prescribe`) |
| `print_simulation_report()` | Text summary of a run (`--sim`) |
| `print_roadmap()` | Workstreams × unlocked conditions (`--roadmap`) |
| `print_uncertainty_report()` | Monte Carlo percentile table (`--uncertainty`) |
| `export_results_json()` | Full simulation results to JSON (`--export`) |

## Section 8: Pygame Visualization (lines 10851–17220, ~6,400 lines)

`run_visualization()` launches a 1400×880 interactive window with 22 view
modes (`MODES` list). It pre-runs simulation + control before launch, and
samples the Monte Carlo once at startup for the RISK view (never per frame).

**Framework:**

- `hover_regions` system: each frame, draw functions append `(x, y, w, h,
  title, details)` tuples; mouse hover triggers `draw_tooltip()` with smart
  positioning — comprehensive tooltips on all interactive elements across all
  22 modes (vial layers, syringe, LNP zoom, recipe rows, pipeline stages,
  organ targets, mechanism tiles, synergy nodes, disease-category nodes,
  comparison bars, biomarker charts, Gantt bars, safety systems,
  manufacturing steps, PK parameters, dosing schedules, clinical endpoints,
  adverse events, science tables, condition tiers, sim metrics)
- `MECH_NODES`: shared (name, color) list for all 37 mechanisms, used by
  SYNERGY and CONDITIONS so mechanism references are visually consistent
- `ANIM_SPEED` control (default 3× slowdown for readability)
- 10 consistent color constants (`C_ACCENT`, `C_GOOD`, `C_BAD`, `C_WARN`,
  ...) + `TIER_COLORS` mapping the 4 coverage tiers
- `_panel()` semi-transparent rounded rectangles; `_wrap_text()` pixel-bound
  text wrapping
- Mode-index constants (`ABOUT_IDX`, `SIM_IDX`, `CONDITIONS_IDX`) computed
  via `MODES.index(...)` — scroll/keyboard handling stays correct as modes
  are added
- Independent scroll states (e.g. `conditions_scroll_offset`)

**The 22 modes by category:**

- **Formulation** — OVERVIEW (vial/syringe/LNP zoom/metrics), RECIPE
  (formulation table + manufacturing + dosing timeline), DELIVERY (injection
  pipeline + biodistribution), LNP (cross-section + cargo), MOLECULAR
  (chemical structures)
- **Biology** — CELLULAR (endosomal escape + nuclear reprogramming), BODY
  (organ targeting), MECHANISMS (37 animated diagrams), SYNERGY (pathway
  network + disease-category ring)
- **Results** — COMPARISON (treated vs untreated), BIOMARKERS (time series +
  radar + metric panel), TIMELINE (50-yr Gantt), SIM (live playback)
- **Safety & quality** — SAFETY (multi-layer shield + QC — the formulation's
  *claims*), MANUFACTURING (7-step GMP), PHARMACOLOGY (PK/PD + dose-response),
  DOSING (schedule + personalization)
- **Evidence & computed risk** — CLINICAL (endpoints + organ scores +
  adverse-event table driven by `sim_result.risk_breakdown`, not literals),
  SCIENCE (validation + chemistry + citations), RISK (the decision-relevant
  view: left = simulated residual risk per hazard class, center = p5–p95
  spread with ROBUST/FRAGILE labels, right = `RECOMMENDED_PROTOCOL` +
  verification stats), CONDITIONS (scrollable survey, tier-colored, "How
  healed" notes), ABOUT (500+ lines of formulation database, audit summary,
  adversarial fix log, regulatory precedents, literature references)

**SAFETY vs RISK is a deliberate separation:** SAFETY renders static
`CHIEF_JUICE` safety dicts (claims); RISK renders computed simulation output.
The gap between the two panels is the point.

## Section 9: Main / CLI (lines 17221–17453)

`main()` with `argparse`. Default (no flags) → `run_visualization()`.
Dispatch order: `--report` → `--feasibility` → `--proof` → `--audit` →
`--roadmap` → `--optimize` → `--uncertainty` → `--list-conditions` →
`--prescribe` → (`--sim` / `--compare` / `--export`) → `--visual`.
`--conditions` feeds `--prescribe`; `--regrow` feeds LimbRegen;
`--no-safeguards` / `--accelerated` / `--no-supplements` / `--no-nutrition`
toggle protocol; `--age`/`--years`/`--target-age`/`--sex` parameterize;
`--runs`/`--risk-weight` tune Monte Carlo / optimizer.

## Data Flow

```
CHIEF_JUICE dict (recipe)  +  BIO dict (biology)
        |
        v
initialize_body() ── EXTENSION_FIELD_PROFILE ── assign_individual_variation
        |
        v
simulate_chief_juice()
   |-- per month: 37 mechanism functions modify BodyState
   |-- nutrient_dynamics → drug_interaction → safety_model
   |-- rebalance_biomarkers → setpoint governor (exact alpha interpolation)
   |-- LimbRegen cell/biomass debits (outside governor)
   |-- death checks / immortality detection
        |
        v
SimulationResult (time-series arrays + final_state + risk_breakdown)
        |
        +---> print_simulation_report()  [text]
        +---> print_formulation_report() [--report]
        +---> print_feasibility_report() [--feasibility]
        +---> print_math_proof()         [--proof]
        +---> print_evidence_audit()     [--audit]
        +---> print_prescription_report()[--prescribe ← DISEASE_COVERAGE + PERSONALIZED_TREATMENTS]
        +---> print_roadmap()            [--roadmap]
        +---> print_uncertainty_report() [--uncertainty ← run_monte_carlo]
        +---> export_results_json()      [--export]
        +---> run_visualization()        [pygame, default]
                  |
                  +---> 22 draw functions read sim_result + control_result
                        (+ DISEASE_COVERAGE for CONDITIONS/SYNERGY)
                  +---> hover_regions populated each frame
                  +---> draw_tooltip() renders on mouse hover
```

## Verification

`test_chiefjuice.py` — 64 checks across 12 suite functions; each suite loads a
fresh module instance so `BIO` mutations cannot leak between suites.

| Suite | Checks |
|---|---|
| `test_survey_integrity` | Derived tiers match note text; tamper resistance vs bulk overwrite; counts/splits |
| `test_mass_balance` | Fat debit = fat removed to 0.0000 kg; essential-floor halt |
| `test_setpoint` | Holding ages reached/held exactly; floor semantics; rate limit |
| `test_regeneration` | Cell-count gating, biomass gating, stall accounting, completion caps |
| `test_safety` | Hazard compounding, mitigations, 5% interlock halt, non-zero residual |
| `test_withdrawal` | Decay, differential persistence, rebound overshoot |
| `test_variation` | Responder classes from `patient_seed`; per-patient rate scaling |
| `test_interactions_and_benefit_risk` | 3 contradictions; TREAT/CAUTION/DO NOT TREAT verdicts |
| `test_no_fabricated_zeros` | Source scan for hardcoded reassuring literals (catches the *class* of bug — e.g. the CLINICAL view once asserted `Teratoma 0.0%` while the engine computed non-zero residuals) |
| `test_doc_consistency` | Doc-quoted numbers match code outputs |
| `test_uncertainty` | Monte Carlo priors restored; percentile sanity |
| `test_evidence_audit_and_honesty` | Audit counts/grades; speculative fraction; honesty locks |

```
python test_chiefjuice.py           # all suites
python test_chiefjuice.py --quick   # skip Monte Carlo
```

Non-zero exit on failure, so it can gate a commit. It exists because the
headline metric was inflated three separate times before anything checked it.

## Module Interface

`FatStemCells.py` loads `ChiefJuiceV1.17.py` via `importlib.util` (the dotted
filename prevents a normal import) and re-exports:

```python
from FatStemCells import (
    CHIEF_JUICE, BIO, BodyState, SimulationResult,
    simulate,                      # -> simulate_chief_juice
    run_comparative_simulation,
    initialize_body,
    print_report,                  # -> print_formulation_report
    print_simulation_report,
    print_feasibility_report,
    print_math_proof,
    export_results_json,
    # All 30+ biological model functions:
    telomere_dynamics, dna_damage_model, senescence_model,
    mitochondrial_model, stem_cell_pool_model, inflammation_model,
    organ_repair_model, immune_rejuvenation_model, epigenetic_age_model,
    fat_dynamics, cancer_risk_model, neural_model, muscle_strength_model,
    bone_strength_model, memory_model, lung_model, oxygen_capacity_model,
    necroptosis_model, cosmetic_model, compute_stacking_multiplier,
    proteostasis_model, nutrient_sensing_model, communication_model,
    oxidative_stress_model, hormonal_model, plasma_exchange_model,
    autoimmune_tolerance_model, gene_correction_model,
    axial_enthesis_model, endocrine_gland_regeneration_model,
    peripheral_nerve_regeneration_model, mucosal_barrier_model,
    neurocircuit_plasticity_model, sensory_epithelium_model,
    # Survey & audit:
    DISEASE_COVERAGE, DISEASE_COVERAGE_TIERS, DISEASE_COVERAGE_CATEGORY_ORDER,
    EVIDENCE_AUDIT, get_evidence_summary, print_evidence_audit,
    get_healing_progress,           # the ONLY supported way to read healed count
    HEALING_TARGET, run_visualization,
)
```

Tiers in `DISEASE_COVERAGE` are DERIVED at import — `get_healing_progress()`
is the only supported way to read the healed count; it is an output, not a
figure that can be asserted.

## Performance

- Simulation: ~2 s for 50-year monthly time-step (600 iterations × 37
  mechanisms + failure models)
- Visualization pre-run: ~4 s (two simulations + Monte Carlo before GUI
  launch)
- GUI: 60 FPS with 22 animated modes, ~150+ hover regions per frame
- Memory: ~50 MB for full time-series arrays + visualization surfaces

## Validation & Review History

- 16 reviews with simulation testing (review #15 was an honesty pass; #16 an
  exhaustive completeness pass — every module detail and stat listed)
- REVIEW #15 CORRECTION: the old "31 domains at 100% / 20/20 verified /
  22/22 validated" self-scores were REMOVED — self-assigned, not evidence,
  and masked real chemistry errors
- REVIEW #16 EXPANSION: `EVIDENCE_AUDIT` grew 58 → **219 graded parameters
  across 72 modules** (A=30, B=61, C=59, D=69 — run `--audit`), and `BIO`
  grew 213 → **472 constants** covering every body system
- What is actually verified: **64/64 internal-consistency checks pass**
  (mass balance, setpoint control, safety interlocks, survey integrity,
  honesty locks, Monte Carlo). Consistent math is not biological proof.
- Chemistry/doses corrected to published values where published values exist
  (review #15: 2i media concentrations, L-Myc swap, hemoglobin O2 cap,
  bone/muscle physics limits, navitoclax/urolithin A/MitoQ/rapamycin/TA-65/
  alendronate/folate/IL-7 doses, TH5487 misattribution, Nutlin-3a
  rescue-only, fake 'drugs' flagged)
- 25 critical issues found and fixed across reviews #10–14, plus review
  #15's honesty pass (31 removed self-scores, impossible magnitudes replaced)
- 6 regulatory precedents (Kymriah, Zolgensma, Comirnaty, ARCT-154, CASGEVY,
  Lonafarnib)
- 22 visualization modes with comprehensive hover tooltip coverage
- 243/304 conditions HEALED (80%), derived from per-condition notes — 61 not
  yet healed, 8 still GAP
- 29 personalized treatment modules with step-by-step formulation
  instructions and healing proofs for unique condition-specific additives
- **Nothing is experimentally validated.** 32% of audited parameters are
  grade-D speculative. The formulation is a computational scenario, not a
  therapy.
