# Experimental Protocol: Evaluating the Effect of Surfactin on *Paenibacillus dendritiformis* Growth

## Goal
To evaluate how different concentrations of surfactin (6.25 µM, 12.5 µM, 25 µM) affect the growth of *Paenibacillus dendritiformis* using a plate reader over a 67-hour period with hourly measurements. Ethanol effects are to be eliminated by ensuring equal bacterial volumes in each well.

## Materials
- Surfactin (1 mg/mL stock solution in ethanol)
- 100% ethanol
- *Paenibacillus dendritiformis* bacterial culture
- 96-well plate
- Plate reader
- Micropipettes and sterile tips

## Procedure

### 1. Bacterial Growth
1. Thaw *P. dendritiformis* from stock and grow in a breathing test tube containing LB medium at 30°C with shaking at 200 RPM for 48 hours.
2. After incubation, vortex the test tube and mix thoroughly to disrupt the bacterial structure.
3. Measure the OD of the bacterial culture at 600 nm using a plate reader.
4. Adjust the culture to OD ~0.1 by diluting with fresh LB medium.

### 2. Preparation of Surfactin Dilutions

1. **Stock Preparation**:
   - Start with a surfactin stock solution at **1 mg/mL** in ethanol.
   - Prepare a **0.1 mg/mL** working solution by transferring **100 µL** of the 1 mg/mL stock into **900 µL** of 100% ethanol.

2. **Mass Calculation for 25 µM Surfactin**:
   - To calculate the mass of surfactin required for a concentration of **25 µM** in **200 µL**:
     - Molecular weight (MW) of surfactin = **1036.3 g/mol**.
     - Use the formula:
       ```
       Mass (g) = Concentration (mol/L) × Volume (L) × MW (g/mol)
       ```
     - For 25 µM:
       ```
       Mass = (25 × 10⁻⁶ mol/L) × (200 × 10⁻⁶ L) × (1036.3 g/mol)
       Mass ≈ 5.18 × 10⁻⁶ g
       ```

3. **Volume Calculation for 0.1 mg/mL Surfactin**:
   - To calculate how much volume of the **0.1 mg/mL** surfactin solution is needed to provide **5.18 × 10⁻⁶ g**:
     - Use the formula:
       ```
       Volume (µL) = Mass (g) / Concentration (g/µL)
       ```
     - For 0.1 mg/mL (which is **0.1 g/L = 0.1 × 10⁻³ g/mL**):
       ```
       Volume ≈ 51.8 µL
       ```

4. **Serial Dilutions**:
   - To prepare **12.5 µM** and **6.25 µM** solutions, perform 2-fold serial dilutions:
     - For **12.5 µM**, Dilute the ethanol-surfactin solution 0.1mg/mL 2 times (with 100% ethanol) to reach a solution with a concentration of 0.05 mg/mL - and evaporate it in the well to reach a concentration of 12.5 µM.
     - For **6.25 µM**, do the same with the ethanol-surfactin solution concenrtation of 0.05 mg/mL to get a concentration of 0.025 mg/mL.

![results](../images/growth%20curves/plate%20planing1.png)
     

### 3. Plate Loading

1. **Loading Surfactin**:
   - Load a treplicate of  **51.8 µL** of the appropriate surfactin+ethanol solution into each well for 25 µM, 12.5 µM, and 6.25 µM concentrations.

2. **Ethanol Evaporation**:
   - Wait for the ethanol to fully evaporate from the wells, leaving only the surfactin.

3. **Loading Bacteria**:
   - Add **200 µL** of *Paenibacillus dendritiformis* bacterial culture into each well after ethanol has evaporated.

### 4. Plate Reader Settings
- Measure OD at 600 nm every hour for 67 hours (68 measurements in total).
- Shake plate at high intensity for 5 minutes before each measurement.

## Notes

- Ensure ethanol completely evaporates before adding bacteria to avoid interference.
- Perform three technical replicates per condition and include control wells with no surfactin (bacteria only).
