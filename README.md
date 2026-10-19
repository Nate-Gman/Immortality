# ChiefJuiceV1.17 -- Injectable Fat-to-Stem-Cell Conversion & Longevity Engine

> **"Your fat is your cure. Your body is your pharmacy."**

A computational biology model (17,453 lines, single file) of an injectable
formulation that reprograms adipose tissue into stem cells in vivo, then deploys
them to reverse the hallmarks of aging, hold the patient at a chosen biological
age, regrow missing structures, and treat a 304-condition disease survey.
Includes a 22-mode interactive visualization, a prescription generator, safety
and uncertainty modelling, and a 64-check verification suite.

> ### Read this before quoting any number
>
> **Nothing here has been experimentally validated.** Not in a cell, an animal, or
> a person. Every parameter is an estimate. A simulation cannot show that a
> treatment heals anything — only that a set of assumptions is internally
> consistent. This model's own uncertainty analysis says its safety figure could
> be off by an order of magnitude.
>
> **Three different numbers, and they answer different questions:**
>
> | | Result | Means |
> |---|---|---|
> | Mechanism coverage | **304/304 — 100%** | Every condition has a mechanism aimed at its root cause. **Not a claim any of them work.** |
> | Engineering tractability | **61/61 — 100%**, 0 fundamental | Every unproven blocker has a *named experimental route*. Not a claim the experiments will succeed. |
> | **Evidence (HEALED)** | **243/304 — 80%** | How many routes are actually proven. **This is the number that matters for a real person.** |
>
> And the evidence figure mixes two populations that should not be averaged:
> **common/acute illness 54/57 (95%)** vs **severe/chronic 189/247 (77%)**. The
> first is mostly diseases medicine already cures.
>
> **The model does not endorse this for healthy people.** Its own benefit/risk
> assessment returns CAUTION for a healthy 30-year-old, because residual risk
> exceeds what a person with nothing wrong should accept.
>
> Earlier versions of this README claimed **243/243 HEALED (100%)**. That figure
> came from bulk-overwriting a `tier` field and was contradicted by the project's
> own per-condition notes. It has been retracted, and `test_chiefjuice.py` now
> fails if anything like it returns.

## Repository contents

```
Immortality.2/
  ChiefJuiceV1.17.py    # Main engine: 17,453 lines, sections 1-9 (see Architecture)
  test_chiefjuice.py    # Verification suite: 490 lines, 64 checks in 12 suites
  FatStemCells.py       # Thin library wrapper re-exporting the engine API
  Goal.md               # Original specification / target document
  Problems.md           # The owner's condition list — the survey's seed
  projectworkflow.md    # Working notes on the review/fix workflow
  OVERVIEW.md           # Technical architecture (companion to this file)
  ReferenceCode/        # Earlier iteration kept for reference
  README.md             # This file
```

## Concept

The human body stores vast reserves of adipose tissue. Each adipocyte retains
its full genome and can be reprogrammed. ChiefJuiceV1.17 models an injectable
cocktail of:

- **Reprogramming factors** — Yamanaka OSKM delivered as saRNA in lipid
  nanoparticles, with an L-Myc/Lin28 substitution and a 4.2× enhancer stack on an
  18% base conversion rate (net 75.6% in the modeled scenario — the single
  largest grade-D assumption in the file; ≤1% is the best published *ex-vivo*
  figure and in-vivo is unquantified)
- **Growth factor cocktail** — GDF11, Klotho, FGF2, EGF, PDGF, TGF-β3, BMP2,
  IGF-1, VEGF, HGF, G-CSF, BDNF, NGF, GDNF
- **Telomerase activators** — Epitalon, Astragaloside IV, Cycloastragenol,
  TA-65, resveratrol
- **Senolytic rotation** — Dasatinib + Quercetin / Fisetin / Navitoclax /
  CAR-T backup, on rotation to blunt senolytic resistance
- **Mitochondrial rejuvenators** — NMN, NR, MitoQ, Urolithin A
- **DNA repair enhancers** — NAD+ precursors, PARP1/2 support, niacin, folate,
  B12
- **Exosomes** — stem-cell-derived vesicles for BBB crossing and intercellular
  signaling repair
- **Nanocarrier delivery** — 6 tissue-specific LNP variants (adipose, brain,
  bone, muscle, lung, immune) plus the DeliverMax coverage model
- **15 healing-gap extension modules** (mechanisms #22–36) added after
  cross-referencing the formulation against the disease survey: `tolerance`
  (autoimmune tolerance reset), `genefix` (in-situ monogenic correction),
  `enthesis`/`SpondyRegen` (axial/enthesis regeneration), `endoregen` (endocrine
  organoid regrowth), `nerveregen` (peripheral nerve), `mucosarepair` (mucosal
  barrier), `circuitreset` (psychiatric neurocircuit plasticity),
  `sensoryregen` (sensory epithelium), `limbregen` (structural/appendage
  regrowth), `fibroclear` (established-fibrosis reversal), `repeatedit`
  (repeat-expansion contraction), `reservoirpurge` (latent-reservoir
  clearance), `tissueforge` (de novo internal tissue regrowth), `delivermax`
  (delivery-platform coverage ceilings), `chromobalance` (chromosomal dosage
  correction)
- **29 personalized treatment modules** (`PERSONALIZED_TREATMENTS`) —
  condition-specific targeted formulations with unique additives, each including
  step-by-step manufacturing instructions (protein expression, purification,
  LNP encapsulation, delivery route, dosing schedule)

*Numbering note: mechanism IDs in this file are the 0-indexed codes used by the
engine (core #0–21, extension modules #22–36, 37 total) — the same IDs stored in
`DISEASE_COVERAGE` entries and returned by `get_prescribed_mechanisms()`. Older
project prose sometimes cited them as #1–37.*

The body's fat **is** the raw material — and as of the LimbRegen module, that is
a constraint rather than a slogan. Fat is *finite*: it is mass-balanced (debit
equals credit to 0.0000 kg over a 50-year run), it runs out at the 5%
essential-fat floor, and structural regeneration draws on it hard enough to
stall when a patient does not have it. "Infinitely renewable" is not what the
model shows.

## Key Results

Default comparative run (`--compare`, start age 60, 50 years, holding age 18):

| Metric | Control | Treated |
|---|---|---|
| Final age | 64.3 (deceased) | 110 (alive) |
| Biological age | 52.5 | 18.0 (held exactly) |
| Telomere length | 1,919 bp (critical) | 7,253 bp |
| Organ function | 66.5% | 73.4% |
| Cancer risk | 32.4% | 0.1% |
| Fat converted | — | 20.0 kg (mass-balanced to 0.0000 kg) |

**Recommended operating point** (from `--optimize`), selected on 5th-percentile
utility rather than median, verified over 1,000 Monte Carlo runs under full
parameter and patient uncertainty:

**22% body fat · 2-on/1-off dosing · supplemented · non-accelerated**

| | p5 | median | p95 |
|---|---|---|---|
| Quality | 71.4 | 84.5 | 89.4 |
| Severe-event risk | 0.24% | 0.65% | 1.83% |

Headline facts:

- 37 mechanisms (22 core #0–21 + 15 extension modules #22–36)
- 304 conditions surveyed in 10 categories; tiers **derived** from each
  condition's note text, never asserted
- Biological-age reversal rate-limited to 0.5 yr/month, and varies per patient
- Fat is finite, mass-balanced, and conversion halts at the 5% essential floor
- 14% of simulated patients are non-responders (efficacy ×0.35), 8% are
  hyper-responders (efficacy ×1.45, risk ×1.9)
- Residual severe-event risk is never reported as zero

## Quick Start

### Prerequisites

- Python 3.8+
- `numpy`
- `pygame` (for visualization)

```bash
pip install numpy pygame
```

### Complete CLI reference

No flags → launches the interactive pygame visualization (22 modes).

| Flag | Type | Default | Effect |
|---|---|---|---|
| `--sim` | switch | off | Text-only single simulation (no GUI) |
| `--compare` | switch | off | Control vs treated comparative run |
| `--no-cj` | switch | off | Control only (no ChiefJuice) |
| `--report` | switch | off | Full formulation report (all components, doses, rationale) |
| `--feasibility` | switch | off | Clinical-evidence grounding report with citations |
| `--proof` | switch | off | Math-consistency check (arithmetic only — not biology) |
| `--audit` | switch | off | Evidence audit: 219 parameters graded A–D vs published science |
| `--prescribe` | switch | off | Prescription generator (interactive without `--conditions`) |
| `--conditions "A,B"` | str | — | Condition list for `--prescribe` (verbatim survey names) |
| `--list-conditions` | switch | off | Print all 304 survey conditions + derived tier counts |
| `--roadmap` | switch | off | Route to 100% proven: 7 workstreams × unlocked conditions |
| `--uncertainty` | switch | off | Monte Carlo confidence intervals + benefit/risk verdict |
| `--runs N` | int | 200 | Monte Carlo run count for `--uncertainty` |
| `--optimize` | switch | off | Robust protocol search (re-derives RECOMMENDED_PROTOCOL) |
| `--risk-weight W` | float | 12.0 | Utility weight on severe-event risk in `--optimize` |
| `--age A` | float | 60.0 | Starting chronological age |
| `--years Y` | float | 50.0 | Simulation horizon in years |
| `--target-age T` | float | 18.0 | Holding biological age (floor on rejuvenation) |
| `--sex S` | str | `unspecified` | `female`/`male`/`unspecified` — hormonal baselines |
| `--regrow "x,y"` | str | — | Missing structures for LimbRegen (see anatomy table) |
| `--no-safeguards` | switch | off | Disable all safety mitigations (shows raw hazard) |
| `--accelerated` | switch | off | Accelerated structural-regen protocol (ex-vivo expansion, local delivery, scaffold retention, anabolic support) |
| `--no-supplements` | switch | off | Turn off cofactor supplementation |
| `--no-nutrition` | switch | off | Turn off all nutritional support (adipose only) |
| `--export FILE` | str | — | Export simulation results to JSON |
| `--visual` | switch | default | Explicitly request the GUI |

### Command recipes

```bash
# Default: launches interactive pygame visualization (22 modes)
python ChiefJuiceV1.17.py

# Text-only simulation (no GUI)
python ChiefJuiceV1.17.py --sim

# Comparative simulation (control vs treated)
python ChiefJuiceV1.17.py --compare

# Formulation report (full recipe with doses)
python ChiefJuiceV1.17.py --report

# Feasibility analysis (clinical evidence grounding)
python ChiefJuiceV1.17.py --feasibility

# Math consistency check (arithmetic, not biology)
python ChiefJuiceV1.17.py --proof

# Evidence audit: every parameter graded A-D vs published science
python ChiefJuiceV1.17.py --audit

# Export results to JSON
python ChiefJuiceV1.17.py --sim --export results.json

# Prescription generator (interactive)
python ChiefJuiceV1.17.py --prescribe

# Prescription generator (with conditions)
python ChiefJuiceV1.17.py --prescribe --conditions "Alzheimer's disease and other dementias,Multiple sclerosis (MS)"

# List all surveyed conditions (with derived tier counts)
python ChiefJuiceV1.17.py --list-conditions

# Custom simulation parameters
python ChiefJuiceV1.17.py --sim --age 40 --years 70

# Rejuvenate to a chosen biological age and HOLD there
python ChiefJuiceV1.17.py --sim --age 70 --years 40 --target-age 30 --sex female

# Regrow a missing appendage (supply-limited by fat-derived stem cells)
python ChiefJuiceV1.17.py --sim --age 35 --years 12 --regrow "hand"

# Confidence intervals instead of point estimates, + benefit/risk verdict
python ChiefJuiceV1.17.py --uncertainty --runs 200

# Recommended operating point, re-derived by robust protocol search
python ChiefJuiceV1.17.py --optimize --risk-weight 12

# Route to 100% proven: which experiment unlocks which conditions
python ChiefJuiceV1.17.py --roadmap
```

### Verification suite

```bash
python test_chiefjuice.py           # 64 checks across 12 suites
python test_chiefjuice.py --quick   # skip the slow Monte Carlo suite
```

Run it after any change to `ChiefJuiceV1.17.py`. It exists because this project's
headline metric was inflated three separate times before anything checked it.

The 12 suites, in file order:

| # | Suite | What it checks |
|---|---|---|
| 1 | `test_survey_integrity` | Derived tiers consistent with note text; tamper resistance against bulk tier overwrite; category/split counts |
| 2 | `test_mass_balance` | Fat debited for conversion equals fat removed from the body (50-yr run, 0.0000 kg tolerance); essential-floor halt |
| 3 | `test_setpoint` | Holding ages 18/25/30/40 reached and held exactly; floor semantics (younger patients held where they are); reversal rate limit |
| 4 | `test_regeneration` | LimbRegen supply limits: cell-count gating, biomass gating, stall accounting, completion caps |
| 5 | `test_safety` | Hazard compounding, mitigation effects, interlock halt at 5%, non-zero residual |
| 6 | `test_withdrawal` | Discontinuation decay, differential persistence, rebound overshoot |
| 7 | `test_variation` | Responder classes drawn from `patient_seed`; per-patient reversal-rate scaling |
| 8 | `test_interactions_and_benefit_risk` | The 3 internal contradictions applied; TREAT/CAUTION/DO NOT TREAT verdicts |
| 9 | `test_no_fabricated_zeros` | Source scan for hardcoded reassuring literals that contradict live engine values |
| 10 | `test_doc_consistency` | Numbers quoted in docs match what the code computes |
| 11 | `test_uncertainty` | Monte Carlo priors restored after run; percentile sanity |
| 12 | `test_evidence_audit_and_honesty` | Audit entry counts/grades; speculative fraction; honesty locks |

Each suite loads a fresh module instance, so a suite that mutates `BIO` cannot
leak into another.

### Evidence audit (`--audit`)

Every load-bearing parameter in the formulation is graded against published
science in `EVIDENCE_AUDIT` — 219 entries across 72 distinct modules.

| Grade | Meaning | Count | Share |
|---|---|---|---|
| **A** | Replicated human clinical data | 30 | 14% |
| **B** | Preliminary human data | 61 | 28% |
| **C** | Animal or in-vitro data only | 59 | 27% |
| **D** | Speculative — no measured evidence | 69 | 32% |

The model's headline outcome — indefinite lifespan maintenance — rests openly
on its 69 grade-D assumptions. The three largest:

- **In-vivo conversion efficiency**: 75.6% modeled vs ≤1% published (ex-vivo)
  and unquantified in vivo
- **Endosomal escape**: 78% modeled vs ≤15% published
- **Multiplicative enhancer stacking**: invented

The audit covers all 89 `CHIEF_JUICE` modules (reprogramming, LNP, growth
factors, senolytics, metabolic, DNA repair, immune, epigenetic, neuro, muscle,
bone, hemo, proteoclear, nutrisense, commrepair, oxidshield, hormonebalance,
plasmasync, ecm, memory, lung, immunoeternal, eternal_youth_pro, all 15
extension modules, safety systems, monitoring, personalized medicine,
drug interactions, PK/PD, stability, and 70+ `BIO` constants) plus a `BIO`
section with 472 biological constants covering hematology, body composition,
cell biology, telomere dynamics, DNA damage/repair, senescence, mitochondrial
function, stem cell pools, immune system, inflammation, organ function, cancer
statistics, pharmacokinetics, nutritional requirements, growth factor PK, drug
half-lives, exosome biology, LNP biophysics, aging rates, regenerative capacity
by tissue, the 12 hallmarks of aging, and 6 epigenetic clocks.

`python ChiefJuiceV1.17.py --audit` prints the full table with modeled value,
published value, grade, source citation, and an honest verdict.

### Choosing a holding age

`--target-age` (default 18) is the biological age the patient is rejuvenated to
and then held at. Treatment is titrated against this setpoint every month by an
explicit governor, so the age is genuinely held rather than clamped:

| `--target-age` | Result for a 70-year-old woman over 40 years |
|---|---|
| 18 | reaches 18.0 at month 54, holds exactly |
| 25 | reaches 25.0 at month 40, holds exactly |
| 30 | reaches 30.0 at month 30, holds exactly |
| 40 | reaches 40.0 at month 10, holds exactly |

Two properties worth knowing:

- **It is a floor on rejuvenation, not a target to be aged up to.** A patient
  already biologically younger than the requested age is held where they are.
  Requesting 55 for a patient at biological age 44.7 holds them at 44.7.
- **Approach is rate-limited** to 0.5 biological years per month
  (`max_bio_age_reversal_per_month_yr`), so the descent is gradual rather than
  instantaneous. The ceiling is also scaled per patient by responder class.
- **18 is the physiological floor** (`min_biological_age_floor`), not a default
  goal — the governor's deadband sits above it.

`--sex` (`female` / `male` / `unspecified`) sets sex-specific hormonal
baselines, including the post-menopausal estradiol deficit for female patients
(menopause onset 51, transition 5 yr, additional sex-hormone drop 35%).
Before this existed, "old woman" and "old man" ran as the same simulation.

## Limb and Appendage Regrowth (LimbRegen, mechanism #30)

```bash
# Regrow a hand
python ChiefJuiceV1.17.py --sim --age 35 --years 12 --regrow "hand"

# Multiple structures, competing for one cell reserve
python ChiefJuiceV1.17.py --sim --years 25 --regrow "finger,hand,lower_leg_and_foot"

# Adipose reserve only — no nutritional support
python ChiefJuiceV1.17.py --sim --years 25 --regrow "leg" --no-nutrition
```

Known structures (`BIO["appendages"]`, mass/axial length for a 70 kg adult):

| Structure | Mass (kg) | Axial length (cm) |
|---|---|---|
| `fingertip` | 0.005 | 1.5 |
| `finger` | 0.030 | 8.0 |
| `thumb` | 0.045 | 6.5 |
| `hand` | 0.400 | 19.0 |
| `forearm_and_hand` | 1.600 | 45.0 |
| `arm` | 3.500 | 73.0 |
| `toe` | 0.035 | 5.0 |
| `foot` | 1.000 | 26.0 |
| `lower_leg_and_foot` | 4.500 | 48.0 |
| `leg` | 11.000 | 85.0 |
| `ear` | 0.015 | 6.0 |
| `nose` | 0.020 | 4.0 |

Mass drives the cell budget (`cells_per_gram_soft_tissue`); axial length drives
regrowth *time* via the nerve-limited growth rate (1.5 mm/day).

### This is the one mechanism the fat supply actually gates

Every other mechanism in the formulation runs off a percentage. LimbRegen spends
a **cell count** drawn from `stem_cell_reserve`, which is credited *only* by
fat→stem-cell conversion, and it spends **biomass**, because a regrown hand is
made of matter that has to come from somewhere. Both constraints bind
independently:

| Patient | Result (regrow a hand, adipose only) |
|---|---|
| 6% body fat | **0.0%** — cannot even seed a blastema |
| 12% body fat | complete at month 25 |
| 25% body fat | complete at month 8 |
| 45% body fat | complete at month 6 |

Cripple engraftment to 0.02% and regrowth stalls at 0.0% even with fat available.
Turn nutrition off and a leg never starts, because 11 kg of leg cannot be built
from 14 kg of usable fat while that same fat is also feeding conversion.

### What the model says is reachable

With sustained nutritional support, 30% starting body fat:

| Structure | Mass | Complete at | Binding constraint |
|---|---|---|---|
| fingertip | 0.005 kg | month 1 | none |
| finger | 0.03 kg | month 3 | none |
| hand | 0.40 kg | month 7 | cell supply (3 stalled months) |
| forearm + hand | 1.6 kg | month 28 | cell supply (20 stalled) |
| arm | 3.5 kg | month 82 | cell supply (69 stalled) |
| lower leg + foot | 4.5 kg | month 84 | cell supply (75 stalled) |
| leg | 11.0 kg | month 282 (23.5 y) | cell supply (267 stalled) |

The ordering falls out of the cell and mass budgets, not from tuning — and it
happens to match the real-world evidence gradient, which is why the derived
coverage tiers land where they do:

| Condition | Derived tier | Why |
|---|---|---|
| Fingertip amputation (distal to nail matrix) | **HEALED** | genuinely happens in humans, especially children |
| Digit amputation | PARTIAL | extends real fingertip biology past the nail matrix — never observed |
| Traumatic limb amputation (hand/forearm/below-knee) | PARTIAL | no mammal regenerates a limb |
| Whole-limb amputation (above-knee/shoulder) | PARTIAL | the model's own supply wall: 23 years, and impossible without feeding |
| Congenital limb difference / amelia | PARTIAL | no amputation plane — but a limb-field induction component addresses the template problem (unproven) |

Adding the module did **not** let it award itself HEALED tiers — `_derive_tier()`
read the honest notes and assigned PARTIAL at best.

### Honest limits

- **No mammal regenerates a limb.** The human wound response forms a fibrotic
  scar that actively forecloses blastema formation. Reversing that is a larger
  unproven step than anything else in this formulation, including syndesmophyte
  resorption.
- **Adult positional-memory re-specification has no experimental support at all.**
- Pattern fidelity caps at **70%** — the model does not claim a regrown limb
  matches the original in strength, sensation, or anatomy.
- Real precedent it does extend: human fingertip regrowth, MRL/MpJ mouse ear-hole
  closure, annual deer antler regeneration, and the Levin lab's 2022 induction of
  substantial leg-like regrowth in adult *Xenopus* with a 5-drug cocktail.

Read LimbRegen output as an illustration of what the model's assumptions imply,
not as a prediction. See `CHIEF_JUICE["limbregen"]["evidence"]`.

## Modelling the Things That Go Wrong

Five subsystems were added after testing showed the model was describing safety
and variability in prose while simulating neither.

### Safety (`safety_model`)

Eight hazard classes, each the known failure mode of its technology — residual
pluripotency after OSKM causes teratoma, ex vivo expanded cells accumulate
transforming mutations, "shock and kill" latency reversal causes cytokine
release syndrome, systemic collagenase weakens vessel walls. Mitigations are the
safeguards the formulation already claimed (iCasp9 suicide switch, p53
activation, forced differentiation, clonal release testing, fibrosis-activated
protease). A **safety interlock** halts dosing at 5% cumulative risk
(`safety_halt_threshold_pct`), as a DSMB stopping rule would, with hysteresis
before resumption.

Per-month hazard rates and mitigation efficacies (`BIO`):

| Hazard | Base rate/mo | Mapped mitigations (efficacy) |
|---|---|---|
| Teratoma | 0.35% | iCasp9 suicide switch 0.82 · p53 tumor suppressor 0.70 · forced differentiation 0.75 · immune surveillance 0.55 · episomal non-integrating 0.80 |
| Ectopic tissue | 0.28% | forced differentiation 0.75 · tissue-targeted delivery 0.68 · spatial confinement 0.75 |
| Expansion transformation | 0.40% | p53 0.70 · clonal release testing 0.85 |
| Cytokine storm | 0.30% | staged latency reversal 0.72 · transient immunosuppression 0.60 |
| Off-target edit | 0.22% | high-fidelity nuclease 0.78 · tissue-targeted delivery 0.68 |
| Structural failure | 0.25% | fibrosis-activated protease 0.88 · vascular surveillance 0.70 |
| Hypoglycemia | 0.18% | closed-loop glucose 0.85 · glucose-responsive graft 0.90 · islet mass titration 0.75 |
| Immune reaction | 0.20% | capsid rotation 0.72 · transient immunosuppression 0.60 |

Risks compound monthly and combine as independent survivals
(1 − Π(1 − rᵢ)), never by summing. Residual risk is never zero by construction.

| | 8-year course |
|---|---|
| Unmitigated (`--no-safeguards`) | halts after 6 doses |
| Full safeguards | **0.78%, zero halts, all 98 doses** |

### Discontinuation & rebound

Stopping at month 60 sends biological age 18 → 41 in 4.5 years, ending roughly
where the untreated control does. Gains decay **differentially**: a regrown hand
persists at 70% (structural regrowth is explicitly exempt from rebound), immune
tolerance falls to 50%, the epigenetic clock drifts back. Suppression-dependent
markers overshoot the untreated baseline — senescence ×1.6, inflammation ×1.4
(`rebound_senescence_overshoot_mult`, `rebound_inflammation_overshoot_mult`) —
which is a documented failure mode of anti-cytokine and senolytic withdrawal.

### Individual variation

Drawn once per patient from `patient_seed` (same seed → same individual):

| Parameter | Value |
|---|---|
| Responder SD | 0.22 |
| Non-responder fraction / efficacy | 14% / ×0.35 |
| Hyper-responder fraction / efficacy / risk | 8% / ×1.45 / ×1.9 |
| Pharmacogenomic clearance SD | 0.30 |
| Immune-reactivity SD | 0.35 |
| Regen-capacity SD | 0.25 |

The reversal-rate ceiling is scaled per patient — holding it fixed made every
responder class rejuvenate at an identical rate, hiding the variation.
Cohort time-to-target spans **50–151 months**; before this every patient was the
mean patient, which systematically overstates reliability.

### Internal drug interactions

Three are contradictions in the formulation itself, not external drugs. The
sharpest: it uses **immunosuppression to permit redosing while relying on immune
surveillance to catch transformed clones** — opposite instructions to the same
system.

| Conflict | Penalty | Effect |
|---|---|---|
| Immunosuppression vs surveillance | `interaction_immunosupp_surveillance` = 0.55 | surveillance runs at 67% of nominal |
| Senolytics vs reprogramming | `interaction_senolytic_reprogramming` = 0.18 | reprogramming efficiency reduced during senolytic rounds |
| mTOR suppression vs regeneration | `interaction_mtor_regeneration` = 0.22 | regenerative throughput penalized |

### Benefit/risk in context

Risk only means something next to what the patient faces untreated:

| Patient | Risk | Tolerance | Verdict |
|---|---|---|---|
| Healthy 30-year-old | 0.52% | 0.36% | CAUTION — marginal |
| 3 chronic conditions, 55 | 0.78% | 1.72% | TREAT |
| Healthy 30, safeguards off | 3.62% | 0.36% | **DO NOT TREAT** |

**The model does not endorse this for healthy people** — the population the
"immortality engine" framing implies.

## Uncertainty: How Much of This to Believe (`--uncertainty`)

Every figure here is computed from parameters that were *estimated, not
measured*. `run_monte_carlo()` samples 16 uncertain parameters
(`MONTE_CARLO_PRIORS`, uniform fold-ranges around the `BIO` base values) plus
the patient population; `BIO` is always restored in a `finally` block. Sampling
across ~11,500 runs:

| Metric | 5th | median | 95th | Spread | Verdict |
|---|---|---|---|---|---|
| biological age | 18.02 | 18.08 | 18.97 | 1.1× | **ROBUST** |
| organ average | 80.95 | 86.60 | 87.15 | 1.1× | **ROBUST** |
| fat used | 14.83 | 14.83 | 14.83 | 1.0× | **ROBUST** |
| regeneration | 59.1 | 81.2 | 95.5 | 1.6× | moderate |
| months to target | 28 | 44 | 90 | 3.2× | FRAGILE |
| **severe risk %** | **0.21** | **0.82** | **2.75** | **12.9×** | **FRAGILE** |

**The safety number is the least trustworthy figure in the model.** Quote it as
a range, never as a point estimate.

The 16 sampled priors (multiplicative fold ranges on `BIO` values):

| Parameter | Fold range | Why it matters |
|---|---|---|
| `adipocytes_per_gram` | 0.5–2.0 | Fat→cell feedstock density |
| `stem_cell_engraftment_pct` | 0.3–3.0 | How many converted cells engraft |
| `ex_vivo_expansion_fold` | 0.3–3.0 | Accelerated-protocol cell supply |
| `blastema_seed_fraction` | 0.4–2.5 | Cells to seed a blastema |
| `blastema_maintenance_fraction_per_month` | 0.5–2.0 | Monthly blastema upkeep |
| `regen_cells_per_pct_point` | 0.4–2.5 | Cell cost per progress point |
| `max_bio_age_reversal_per_month_yr` | 0.5–1.5 | Rejuvenation rate ceiling |
| `adipose_inflammation_fold` | 0.5–2.0 | Adiposity→inflammation coupling |
| `hazard_teratoma_per_month_pct` | 0.3–4.0 | |
| `hazard_expansion_transformation_pct` | 0.3–4.0 | |
| `hazard_ectopic_tissue_per_month_pct` | 0.3–4.0 | |
| `hazard_immune_reaction_per_month_pct` | 0.3–4.0 | **Dominant risk driver** |
| `mitig_icasp9_suicide_switch` | 0.7–1.15 | |
| `mitig_clonal_release_testing` | 0.7–1.15 | |
| `interaction_immunosupp_surveillance` | 0.5–1.6 | |
| `biomass_conversion_efficiency` | 0.6–1.6 | |

**Value of information:** one parameter — `hazard_immune_reaction_per_month_pct`
— dominates risk uncertainty (Spearman +0.39). Measuring it collapses most of
the spread. Nine of the sixteen sampled parameters move nothing and can stay
estimates.

### Recommended operating point (`--optimize`)

`optimize_protocol()` searches body-fat level × dosing schedule ×
supplementation × acceleration, ranking by **5th-percentile utility** rather
than median, because the two select measurably different protocols.
`RECOMMENDED_PROTOCOL` in code:

```python
{'start_fat_pct': 22.0, 'dose_cycle': (2, 1), 'nutrient_support': True,
 'accelerated_regen': False,
 'rationale': 'best 5th-percentile utility across parameter and patient uncertainty',
 'verified': {'runs': 1000, 'quality_p5': 71.4, 'quality_median': 84.5,
              'quality_p95': 89.4, 'risk_p5_pct': 0.238, 'risk_median_pct': 0.654,
              'risk_p95_pct': 1.83, 'frac_runs_risk_over_2pct': 3.9},
 'accelerate_if': 'risk weight W < 7.9 (high-burden indications such as limb loss)'}
```

Two findings retire the earlier "accelerated is strictly better" claim:

1. **Acceleration is mean-improving but tail-worsening** — median quality 86.3 →
   90.0, but risk p95 2.05% → 3.78%.
2. **Tuning the cheap levers beats adding the risky one** — with body fat and
   schedule free to vary, non-accelerated wins at *every* risk weight tested.
   Acceleration buys with risk what scheduling supplies for free.

`--risk-weight` is exposed because that trade is a value judgement, not a
measurement.

## Route to 100% Proven (`--roadmap`)

Treating every blocker as an engineering hurdle rather than an impossibility
turns out to be checkable, not just optimistic. `_derive_tractability()`
requires each remaining blocker to match a **named experimental path** — a
specific platform, trial, or scale-up — or it is classified FUNDAMENTAL.

**Result: 61/61 tractable, 0 fundamental.** Nothing in the survey is blocked by
missing physics or unrecoverable information. The 61 group into 7 workstreams:

| Conditions | Workstream | The step |
|---|---|---|
| **18** | Clinical trial | Mechanism has preclinical support and a delivery route; needs a Phase 1. Regulatory and funding, not scientific unknown |
| **12** | Stratification | Syndromic labels covering several distinct pathologies — resolve into biologically defined subtypes with measurable targets |
| **10** | Cell replacement | Tissue replacement rather than rescue: iPSC grafts + scaffold-guided circuit reconstruction (bemdaneprocel precedent) |
| **10** | Morphogenesis | Blastema induction, scaffold-guided growth, de novo field induction — embryo/rodent literature into adult human |
| **5** | Delivery platform | Directed-evolution capsids, focused-ultrasound BBB opening, convection infusion, PEG-shielded redosing |
| **3** | Cargo/scale | Targeted chromosome elimination, in vivo XIST, integrase multi-megabase cassettes |
| **3** | Translation | Demonstrated in vitro; package into a vehicle, large-animal efficacy, IND |

`python ChiefJuiceV1.17.py --roadmap` prints this with the conditions each
workstream unlocks. The delivery-platform workstream is the highest leverage per
unit of work; the clinical-trial bucket is the largest but is 18 separate
programs.

## Common & Acute Illness (57 conditions) — and an honest warning about the percentage

The survey was 247 conditions of almost entirely rare, severe, chronic disease.
Checking it against a list of 50 everyday illnesses found **40 were simply
absent** — no influenza, no UTI, no strep throat, no appendicitis, no acne, no
peptic ulcer. As a picture of "what people get sick with," it was badly skewed.

57 common conditions were added. Most genuinely resolve or are cured by existing
standard of care, so **54 of 57 earn HEALED honestly** — H. pylori eradication
really does cure peptic ulcer disease, penicillin really does cure strep throat,
appendectomy really is definitive. The three that don't: irritable bowel
syndrome (functional, no consistent lesion), migraine (pathophysiology
incompletely understood), and osteoarthritis (cartilage regeneration unproven
in humans).

### ⚠️ This raised the headline number without curing anything hard

| Population | Healed | Rate |
|---|---|---|
| **Common & acute illness** | 54/57 | **95%** — everyday sickness, mostly already curable |
| **Severe & chronic** (the original survey) | 189/247 | **77%** — where the hard problems are |
| Combined headline | 243/304 | 80% |

The overall figure moved 77% → 80% purely because easier conditions were added
to the denominator. **Not one severe or chronic condition improved.** Prion
disease, Down syndrome and whole-limb amputation are exactly where they were.

This is a real way to make a coverage number look better while accomplishing
nothing, so the split is computed in `get_healing_progress()["split"]` and
printed by `--list-conditions` rather than left to whoever writes the README.
Quote the split, not the headline.

## Closing the Remaining Gaps (mechanisms #31–34)

Four modules were added to attack the largest clusters of not-yet-healed
conditions. Each has real state in `BodyState` and real monthly dynamics —
under treatment the pathology reverses, and untreated it progresses:

| Module | Attacks | Real precedent it extends |
|---|---|---|
| **FibroClear** (#31) | Established scar, the "irreversible fibrosis" cluster | Pirfenidone + nintedanib (FDA-approved); FAP-CAR-T reversed established cardiac fibrosis in mice (Nature 2019) |
| **RepeatEdit** (#32) | CAG/CGG/CTG expansions GeneFix couldn't touch | Cas9-nickase CAG contraction in patient cells (Cinesi 2016); tofersen (FDA 2023); AMT-130 Phase 1/2 |
| **ReservoirPurge** (#33) | Latent provirus, cccDNA, ganglionic latency | Berlin/London/Düsseldorf patients are confirmed HIV cures; EBT-101 dosed in Phase 1 |
| **TissueForge** (#34) | Absent internal structures (bowel, lymphatics) | Tissue-engineered intestine in rodents; enteric neural crest restores motility; teduglutide approved |

Sample dynamics — a patient at 65% fibrosis / 40% architecture:

| | month 1 | month 12 | month 60 | untreated at 60 |
|---|---|---|---|---|
| fibrosis burden | 53.6% | 6.5% | 0.0% | 74.0% |
| architecture | 44.2% | 75.1% | 85.0% (capped) | 34.0% |

Every module keeps an honest residual: repeat tracts contract to a 35% floor
(partial, not normalised), reservoirs to ~3% (not sterilising), and
TissueForge's functional integration deliberately lags anatomical regrowth —
a regrown bowel segment without peristalsis is present but not working.

### What this actually moved

**GAP fell from 15 to 11.** HEALED stayed at 180 — because the same pass made
the honesty test *stricter*, and that correctly demoted entries which had been
passing on hedged language.

Two loopholes were closed in the resolver:

1. **Soft hedges are now caught.** Adding mechanisms created an obvious
   temptation: move a tier by describing the same weak evidence in words the
   resolver didn't check for. Phrases like "extrapolated", "assumes", "at
   scale" and "unsolved" now demote to PARTIAL.
2. **Precedent must be declared, not inferred.** Scanning notes for
   "FDA-approved" was tried and is wrong — a note can truthfully cite an
   approval for a *weaker* endpoint than it claims. Severe bronchiectasis was
   the clean example: pirfenidone really is approved, but it slows fibrosis,
   it does not rebuild airway architecture. Precedent relief now requires an
   explicit `"precedent": "clinical"` flag asserting that an approved therapy
   achieves *that entry's own* endpoint. 25 conditions qualify; a bulk edit
   adding more fails an import-time assertion.

### Why 100% is not reachable

The push toward 100% produced two more modules (**DeliverMax** #35,
**ChromoBalance** #36) and moved the severe/chronic count 180 → **189/247
(77%)**. Then it hit a wall that is not made of code. Those remaining 59 break
down as:

| # | Blocker | Can more mechanisms fix it? |
|---|---|---|
| **25** | Mechanism exists; **evidence is immature** — never done in a human | **No.** Needs trials, not code |
| **20** | **Intrinsic limit** — aneuploidy, destroyed tissue, no positional template | No |
| **12** | **Polygenic / circuit-level** — the root cause is not a lesion | No |
| 2 | Residual uncertainty only | Possibly |

The 25 in the first row are the honest core of it. Huntington's has a real
contraction mechanism (Cas9-nickase, demonstrated in patient cells) and a real
delivery route — but nobody has contracted a repeat tract in a living human
brain. No module written here changes that. Writing a note that omits it would
move the number and would be a lie.

**DeliverMax makes this concrete.** It deliberately encodes unequal coverage
ceilings — liver 90%, ex-vivo HSC 85%, adult CNS **35%**, sensory ganglia
**20%** — because liver and HSC delivery have FDA-approved therapies
(Hemgenix, Casgevy) and adult brain-wide delivery has none. Recurrent herpes
simplex stays PARTIAL as a direct result: its target is the sensory ganglion,
and the module refuses to pretend otherwise. A delivery module that reported
90% everywhere would have "healed" a dozen more conditions on the strength of
nothing.

### The irreducible 8 (GAP)

Prion disease, Down syndrome, Prader-Willi, Williams, corticobasal
degeneration, locked-in syndrome, severe autism, irreversible blindness/
deafness — the final GAP count reported by the derived resolver is 8
(`get_healing_progress()["tier_counts"]["GAP"]`).

ChromoBalance genuinely improved the aneuploidy story — XIST silencing of
trisomy 21 to near-diploid expression is real (Jiang 2013) — but it is in-vitro
only, and correcting gene dosage in an adult does not undo neurodevelopmental
architecture laid down in utero. Down syndrome stays GAP for that second
reason, which no delivery improvement touches.

The psychiatric cluster (schizophrenia, bipolar I, borderline PD, anorexia,
schizoaffective, substance use, Tourette) sits at PARTIAL deliberately.
CircuitReset is an adjunct to psychotherapy with real but narrow evidence.
Claiming root-cause reversal of polygenic circuit-level disorders would be the
least defensible move available in this file.

### What "100%" would actually require

Not more code. It would require: a repeat tract contracted in a human brain, a
blastema induced in a mammal, XIST delivered in vivo, an absent bowel segment
regrown with working peristalsis, and a fibrotic lung rebuilt to normal
architecture. Those are experiments. **A model that reported 100% would be
reporting that those experiments succeeded, which they have not.**

The 61 non-healed entries are the honest to-do list, and most of them are
waiting on laboratories rather than on this repository.

## Interactive Visualization

22 view modes, cycled with TAB. Hover for detail tooltips on all interactive
elements.

| # | Mode | Shows |
|---|---|---|
| 1 | OVERVIEW | Vial, syringe, LNP zoom, live metrics |
| 2 | RECIPE | Full formulation table, manufacturing, dosing timeline |
| 3 | DELIVERY | Injection pipeline, biodistribution |
| 4 | LNP | Lipid nanoparticle cross-section + cargo |
| 5 | MOLECULAR | Chemical structure diagrams |
| 6 | CELLULAR | Endosomal escape, nuclear reprogramming |
| 7 | BODY | Organ targeting map |
| 8 | MECHANISMS | All 37 animated mechanism diagrams |
| 9 | SYNERGY | Interconnected pathway network + disease-category ring |
| 10 | COMPARISON | Treated vs untreated aging |
| 11 | BIOMARKERS | Time-series charts, radar, metric panel |
| 12 | TIMELINE | 50-year Gantt with milestones |
| 13 | SAFETY | Multi-layer shield + QC — what the formulation **claims** |
| 14 | MANUFACTURING | 7-step GMP flowchart |
| 15 | PHARMACOLOGY | PK/PD, dose-response, biodistribution |
| 16 | DOSING | Schedule, component breakdown, personalization |
| 17 | CLINICAL | Efficacy endpoints, adverse events, organ scores |
| 18 | SCIENCE | Mechanism validation, chemistry, citations |
| **19** | **RISK** | **Simulated residual risk per hazard class, uncertainty spread, recommended protocol — what the model *computed*** |
| 20 | CONDITIONS | Scrollable survey with derived tiers |
| 21 | SIM | Live playback with bio-age trajectory |
| 22 | ABOUT | Full formulation database (scrollable) |

**SAFETY and RISK are deliberately separate.** SAFETY renders the static
`CHIEF_JUICE` safety dicts — the claims. RISK renders `risk_breakdown`,
percentile spreads, and `RECOMMENDED_PROTOCOL` — the simulated result. The gap
between the two panels is the most informative thing in the app.

CLINICAL's adverse-event table is also simulation-driven. It previously held
hardcoded literals asserting `Teratoma 0.0%` and `Severe adverse event 0.0%`
while `safety_model()` computed non-zero residuals for exactly those hazards —
the UI contradicting the engine. `test_no_fabricated_zeros` now fails if such
a literal returns.

### Controls

| Key | Action |
|---|---|
| TAB | Cycle view modes |
| SPACE | Play/pause (SIM mode) |
| `,` / `.` | Sim speed |
| `[` / `]` | Animation speed |
| UP/DOWN, PGUP/PGDN, HOME/END | Scroll (ABOUT / CONDITIONS / RISK) |
| Mouse wheel | Scroll |
| H | Help overlay |
| ESC / Q | Quit |

## Architecture

| Section | Lines | Contents |
|---|---|---|
| 1 | ~164–2485 | Formulation constants (`CHIEF_JUICE`) — 89 top-level modules, 2,100+ nested keys |
| 1C | ~2486–3328 | Evidence audit (`EVIDENCE_AUDIT`) — 219 graded parameters across 72 modules |
| 1B | ~3329–4317 | Disease coverage survey (`DISEASE_COVERAGE`) — tiers **derived**, never asserted |
| — | ~4317–5146 | Personalized treatments (29 modules with formulation instructions) |
| — | ~5146–5528 | Tier derivation resolver + integrity locks + prescription logic |
| 2 | ~5528–6354 | Biological constants (`BIO`) — 472 constants covering all body systems |
| 3 | ~6354–6862 | Data structures (`BodyState` ~143 fields, `SimulationResult`) |
| 4 | ~6863–8399 | Biological models (37 mechanism functions + failure models) |
| 5 | ~8400–9466 | Simulation engine (monthly step, setpoint governor, safety interlock) |
| 5B | inside 4–5 | Failure/variability modelling — safety, rebound, variation, interactions, Monte Carlo |
| 6 | ~9467–9627 | Comparative analysis |
| 7 | ~9628–10850 | Reporting (formulation, feasibility, roadmap, uncertainty, protocol, audit) |
| 8 | ~10851–17220 | Pygame visualization (22 modes, ~6,400 lines) |
| 9 | ~17221–17453 | CLI (argparse) |

See `OVERVIEW.md` for the full technical architecture document.

## Library Usage

`FatStemCells.py` loads `ChiefJuiceV1.17.py` via `importlib` (the dotted
filename prevents a normal import) and re-exports the public API:

```python
from FatStemCells import simulate, get_healing_progress, CHIEF_JUICE, BIO

r = simulate(start_age=70, duration_years=40, target_age=30, sex="female",
             dose_cycle=(2, 1), patient_seed=42)
print(r.biological_age[-1], r.cumulative_severe_risk_pct)

print(get_healing_progress()["split"])   # never quote the headline alone
```

Re-exported: `CHIEF_JUICE`, `BIO`, `BodyState`, `SimulationResult`, `simulate`,
`run_comparative_simulation`, `initialize_body`, all report/export functions,
all 37+ biological model functions, `DISEASE_COVERAGE` (+ tier/category
constants), `EVIDENCE_AUDIT`, `get_evidence_summary`, `get_healing_progress`,
`HEALING_TARGET`, `run_visualization`. See `OVERVIEW.md` §Module Interface for
the complete list.

## License

Computational biology model for research and educational purposes. The
formulation is hypothetical. **Nothing here is medical advice, nothing has been
experimentally validated, and the model does not endorse this for healthy
people.**
