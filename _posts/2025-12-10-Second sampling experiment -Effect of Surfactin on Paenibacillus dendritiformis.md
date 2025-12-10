# Effect of Surfactin on *Paenibacillus dendritiformis* - second sampling experiment

This document summarizes the experimental design, preparation workflow, and sampling logic for an experiment comparing *P. dendritiformis* growth and physiology with vs. without 6.25 µM surfactin. Each condition contains five biological replicates. Three negative-control tubes (LB only) are also included.

---

## Overview of the Experimental Structure

### Conditions and Replicates
- with Surfactin (+): 15 tubes (5 per time point)
- without Surfactin (-): 15 tubes (5 per time point)
- Negative control (LB only): 3 tubes (1 per time point)

### Time Points
- 2 h  
- 12 h  
- 24 h  

At each time point:
- 11 tubes are used (5 +surfactin, 5 –surfactactin, 1 control)
- A total of 33 tubes across the experiment

### Measurements per Tube
- OD600 measurement (technical triplicates)
- Metabolomics sample (for LC-MS; only +/– surfactin tubes)
- FACS sample (all tubes, including control)
- RNA sample (only +/– surfactin tubes)

---

# Part 1 — Preparation Phase

## Starter Cultures (Night Before)
- Three independent 5 mL LB cultures of *P. dendritiformis* were grown overnight, together with a negative-control LB-only tube.
- 200 RPM, 30 degrees.


## Morning — Pre-Growth to Ensure Log Phase
- 30 tubes were filled with 4500 µL fresh LB each.
- The 3 starter tubes were combined into one 50 mL falcon and vortexed.
- from the combined culture, 0.5 mL were transferred each of the 30 new tubes. basically we diluted the culture by 10.  
- The diluted (1:10) 30 cultures were grown for an aditional 3 hours to avoid stationary phase.
- Tt the end of the 3 hours OD readings averaged ~0.3.
- To reach an initial OD around 0.1, an 8-fold dilution was ultimately performed, by preparing another 30 tubes with 7 mL fresh LB each and transferring 1 mL of the "~0.3 OD-ed" tubes to the new ones.

Final culture stock used for inoculation:
- 30 tubes × 8 mL = 240 mL (sufficient volume for the experiment)

---

## Preparation of Surfactin and Control Tubes

Surfactin arrived concentrated and was diluted with ethanol(100%) to 1 mg/mL (we diluted 1:10).  
Because ethanol can affect cells, experimental tubes were pre-loaded, and ethanol was allowed to evaporate before adding bacteria. This step can be done during the 3 hours additional growth phase in order to save time (recomended), since the tubes need to sit open in the hood for a few hours in order for the ethanol to evaporate.


### 1. Surfactin Tubes (15 tubes)
Each tube contained:
- 32.5 µL of 1 mg/mL surfactin solution (in ethanol)
- 67.5 µL LB  
Total pre-load volume: 100 µL

### 2. No-surfactin Tubes (18 tubes)
Each tube contained:
- 30 µL 100% ethanol  
- 70 µL LB  
Total pre-load volume: 100 µL

These include the (–) surfactin tubes and the three LB-only negative controls.

All 33 tubes were left open in a hood for ethanol evaporation.



### Surfactin Quantity Calculation 

- Final concentration = 6.25 µM surfactin
- Final volume = 5 mL
- Molecular weight (MW) = 1036.3 g/mol

```

Volume required from a 1 mg/mL (1 µg/µL) surfactin solution:

Mass (g) = Concentration (mol/L) × Volume (L) × MW (g/mol)
Mass = (6.25 × 10⁻⁶ mol/L) × (5 × 10⁻³ L) × (1036.3 g/mol)
Mass ≈ 3.24 × 10⁻⁵ g
Volume (µL) = Mass (g) / Concentration (g/µL)
Volume ≈ 32.5 µL
```


Ethanol content:
- Stock solution is 90% ethanol
- 0.9 × 32.5 µL ≈ 29.2 µL ≈ 30 µL ethanol  
This matches the ethanol-only control volume.

---

## Inoculation of Experimental Tubes

After ethanol evaporation:
- Each tube contained approximately 100 µL.
- 4.9 mL of the "~0.1 OD-ed culture" was added to all experimental tubes.
- Negative controls received 4.9 mL LB instead of culture.

Total working volume per tube: **5 mL**

---

# Part 2 — Experimental Phase

All tubes were incubated under optimal growth conditions for *P. dendritiformis*. (200 RPM, 30 degrees.)

At each time point, tubes were vortexed and handled as follows:

---

## 1. OD Measurement
- 600 µL was aliquoted into 3 wells (technical triplicates, 200 µL in each well ) for each tube.
- 15 wells per treatment (3 wells × 5 replicates).
- Control tubes also sampled.
- Plate layout and readings recorded.

![]()
- the results of OD measurments can be found [here]()

---

## 2. Metabolomics Sampling
- 1 mL was collected from each +surfactin and –surfactin tube.
- Negative controls were not used for metabolomics.
- Samples stored at −80°C in a cryo tube.

---

## 3. FACS Sampling
- From each time point (11 tubes, including the negative control), 500 µL was collected and mixed with 2.5 µL glutaraldehyde.
- Stored at −80°C in a cryo tube for later flow cytometry.

---

## 4. RNA Sampling
- 10 test tubes (+/- surfactin, not the N.C) were cenrifuged in 7000g for 5 minuts.
- Supernatant discarded, pellets resuspended in 1 mL RNAsave (pipate thorothoroughly).
- Stored at −80°C.

---

# Labeling Strategy
- Tubes with *blue* writing contain samples of PD without surfactin. Tubes with *red* writing contain samples of PD with surfactin. Control samples are contained within tubes with *black* writing.
- Samples were labeled A–E within each condition and time point.
- OD irregularities were noted on tube labels.
- Cross-referencing allows linking metabolomics, FACS, and RNA results.
- Note: excspt of the 24 hour sampling - there is **NO** connectio between the location of samples on the 96 plate and the labeling.

