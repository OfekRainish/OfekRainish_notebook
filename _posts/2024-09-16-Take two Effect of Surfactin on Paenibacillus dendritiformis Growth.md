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

## Results

![results](../images/growth%20curves/19.9%20growth%20curve%20evaporation.png)

[raw data](../exel%20files/growth%20curve/datacsv%20evapo%20data%20filtered%2019.9.csv)

## take two
We repeated the experiment when we changed the concentrations, instead of a concentration of 6.25 µM we put a concentration of 50 µM surfuctin.
 To do this we doubled the amount of the serfectin+ethanol solution to 103.36 µL.

 ## Results

![results](../images/growth%20curves/2.10%20growth%20curve%20evaporation.png)

[raw data](../exel%20files/growth%20curve/data%202.10%20sortedcvs.csv)

## take three

We repeated the growth curves experiment with modifications. This time, surfactin was diluted in LB, and 100 µL of the solution (containing LB, surfactin, and ethanol) was transferred into each well. Ethanol was allowed to evaporate for approximately 3 hours before the addition of bacterial cultures.

Additionally, the overnight bacterial culture was diluted 10 times and incubated for 3 hours to ensure that the bacteria were in the exponential growth phase before adding them to the wells.

### Surfactin Dilution and Concentration Calculation

To prepare the required concentrations of surfactin, we followed these steps:

1. **Calculate the volume needed from the initial surfactin stock (1 mg/mL) to create a solution with a concentration of 50 µM in 1200 µL**:

   **Formula**:
   ```markdown
   Mass (g) = Concentration (mol/L) × Volume (L) × MW (g/mol)

   Mass = (50 × 10⁻⁶ mol/L) × (1200 × 10⁻⁶ L) × 1036.3 g/mol
   Mass ≈ 6.2 × 10⁻⁵ g

   Determine the volume of surfactin stock (1 mg/mL) needed:
  
   Volume (µL) = Mass (g) / Concentration (g/µL)
   Volume ≈ 62 µL

We took **62 µL** of surfactin from the stock solution (1 mg/mL) and made up the volume to **1200 µL** with LB. This test tube contained a surfactin concentration of **50 µM**.

### 2-Fold Dilutions
A series of **2-fold dilutions** was prepared with LB to achieve the following concentrations:
- **25 µM**
- **12.5 µM**

### Preparing the Wells
Each test tube's diluted solution was transferred to the wells in **5 riplicate** axpect the loest concentration which we prepered 10 replicates. From each test tube, **100 µL** was transferred into three wells, as shown in the following image:

![](../images/growth%20curves/plate%20planning%203.png)

### Incubation and Bacterial Addition
After the ethanol evaporation period (3 hours), **100 µL** of bacterial culture (OD ~0.1) was added to each well, excluding the negative control wells (so another 1:2 dilution).

### Plate Reader Setup
- The plate was placed in the plate reader for **84 hours**, with a total of **85 measurements**.
- **Intense shaking** for **5 minutes** was performed before each measurement.

### Results:
![](../images/growth%20curves/28.10.24%205%20replicates%20res.png)
[raw data](../exel%20files/growth%20curve/28.10.24%205%20replicates%20res%20sorted.csv)



