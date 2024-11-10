# Protocol: RNA Extraction for RNA-seq of *Paenibacillus dendritiformis* with and without Surfactin

## Objective
To extract RNA from *Paenibacillus dendritiformis* grown in exponential and stationary phases, with and without surfactin treatment, for RNA-seq sequencing.

## Background
In previous experiments ([see here](../_posts/2024-09-16-Take%20two%20Effect%20of%20Surfactin%20on%20Paenibacillus%20dendritiformis%20Growth.md)), growth curves of *P. dendritiformis* with varying concentrations of surfactin were analyzed, determining the exponential phase at 20 hr and stationary phase at 44 hr. For this experiment, we are using a surfactin concentration of 6.25 µM.

## Materials
- *Paenibacillus dendritiformis* culture
- LB media
- Surfactin solution (1 mg/mL, contains ethanol)
- Incubator (30°C, 200 RPM)
- Fume hood
- Plate reader (600 nm OD measurements)
- Breathing test tubes
- RNA extraction kit ([protocol link](../_posts/2024-07-24-RNA%20Extraction%20Protocol.md))

## Procedure

### Day 1: Starter Culture Preparation
1. Prepare 3 starter cultures of *P. dendritiformis* in LB medium.
2. Incubate overnight at 30°C, 200 RPM.

### Day 2: Experimental Setup

#### Step 1: Dilution of Starter Cultures
1. Combine all starter cultures into one test tube for uniformity.
2. Dilute 1:10 with fresh LB in breathing test tubes.
3. Incubate for 3 hours at 30°C, 200 RPM, then dilute again to reach an OD of ~0.1.

#### Step 2: Tube Preparation
1. Calculate the volume of surfactin solution needed to reach a final concentration of 6.25 µM in 5 mL.
   - **Calculation:**
     - **Mass (g) = Concentration (mol/L) × Volume (L) × MW (g/mol)**
     - Mass = (6.25 × 10⁻⁶ mol/L) × (5 × 10⁻³ L) × (1036.3 g/mol)
     - Mass ≈ 3.24 × 10⁻⁵ g
   - **Volume (µL) = Mass (g) / Concentration (g/µL)**
     - Volume ≈ 32.5 µL

2. Pipette **32.5 µL** of surfactin solution into **13 breathing tubes** (10 for extraction, 3 for growth monitoring) and add **67.5 µL of LB**.
3. Prepare **13 Pd tubes without surfactin** by adding 100 µL of LB.
4. Prepare **3 negative control tubes** with 100 µL of LB only.

5. Place all tubes in a fume hood for **3 hours** to allow ethanol evaporation, then add **4900 µL of *P. dendritiformis* culture at ~0.1 OD** to each tube except the negative controls.

*Note:* We will prepare a total of **29 test tubes** for RNA extraction and growth stage estimation as follows.

### Tube Allocation
- **Extraction Tubes:** 20 tubes
  - **Exponential Phase Extraction:** 5 with surfactin, 5 without surfactin
  - **Stationary Phase Extraction:** 5 with surfactin, 5 without surfactin
- **Growth Stage Estimation Tubes:** 9 tubes
  - 3 with surfactin, 3 without surfactin, 3 LB tubes (negative controls)

![Tube Allocation Diagram](../images/growth%20extraction/tube%20destribution.png)

### Step 3: Monitoring Growth
1. Every few hours (t=0, 8, 16, 20, 24, 28, 34, 38, 42, 44):
   - **Sample Collection:** Take **200 µL** from each of the nine growth monitoring tubes.
   - **OD Measurement:** Measure OD using the plate reader at 600 nm.
   
   *Purpose:* This allows tracking of bacterial growth to confirm extraction timing for the desired stages.
   
2. **FACS Sampling:**
   - After OD measurement, take **200 µL** from each tube and add **1 µL of glutaraldehyde** for fixation.
   - Store in a cryotube at **-80°C**.

### Step 4: RNA Extraction Points

At the desired growth phases (exponential and stationary), perform the following:

1. **OD Measurement:** Measure OD for each of the nine monitoring tubes using the plate reader.
   
2. **RNA Preservation:**
   - Centrifuge extraction tubes for **5 minutes at 7000g**.
   - Resuspend the pellet in **1 mL RNAsave** solution.
   - Transfer the resuspended sample to a cryotube and store at **-80°C** until RNA extraction.

3. **FACS Sampling:** Take **200 µL** from the extraction sample and add **1 µL of glutaraldehyde** for FACS. Store at **-80°C**.

![sampling time table](../images/growth%20extraction/exp%20planing%201.png)

## Results
moitoring tubes graph based on OD:
![sampling time table](../images/growth%20extraction/the%20nine%20graph.png)
