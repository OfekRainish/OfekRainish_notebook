# PCR Calibration Protocol

This protocol describes the steps to calibrate a PCR reaction and determine the efficiency of amplification for a specific gene using a standard curve method.

## Objective

To establish a standard curve for quantifying DNA and calculate the PCR efficiency, ensuring accurate and reliable quantification of the gene of interest.

## Materials

- DNA samples of known concentrations (e.g., 4, 2, 0.5, 0.25, 0.125 ng/µL)
- qPCR machine
- qPCR reaction mix (including primers specific to the gene of interest, master mix, etc.)
- PCR tubes or 96-well plate
- Pipettes and tips
- Nuclease-free water

## Procedure

### 1. Prepare Serial Dilutions of DNA Template

1. Prepare a series of known DNA concentrations by serially diluting a stock solution.
   - Example concentrations: 4 ng/µL, 2 ng/µL, 0.5 ng/µL, 0.25 ng/µL, 0.125 ng/µL.
2. Label each dilution tube or well clearly.

### 2. Set Up qPCR Reactions

1. Prepare a qPCR reaction mix according to the manufacturer's instructions.
2. Pipette the appropriate volume of each DNA dilution into separate qPCR wells (e.g., 20 µL per well).
3. Include no-template controls (NTC) to check for contamination.

### 3. Run the qPCR Program (that suits your expiriment)



### 4. Record Ct Values

1. Monitor the qPCR run and record the threshold cycle (Ct) values for each dilution.
2. Higher concentrations will show lower Ct values (cross the threshold earlier).
3. Lower concentrations will show higher Ct values (cross the threshold later).

### 5. Plot the Standard Curve

1. Plot the logarithm (base 10) of the initial DNA concentration (x-axis) against the Ct values (y-axis).
2. Fit a linear regression line through the points.

### 6. Calculate PCR Efficiency

1. Determine the slope of the standard curve from the linear regression.
2. Calculate the PCR efficiency using the formula:

Efficiency = 10^(-1/slope)
- Example: If the slope is -3.32, calculate the efficiency:

   Efficiency = 10^(-1/-3.32) ≈ 2 (or 100% efficiency)





### 7. Interpretation of Results

- An ideal slope is approximately -3.32, indicating 100% efficiency (the DNA quantity doubles with each cycle).
- Acceptable PCR efficiency ranges from 90% to 110%, corresponding to a slope between -3.1 and -3.6.
- Consistency in efficiency across different concentrations confirms the reliability of the qPCR assay.

## Notes

- Ensure all reagents and equipment are contamination-free.
- Accurate pipetting is crucial for reliable standard curve data.
- Repeat the experiment if the slope deviates significantly from the expected range.

## Conclusion

By following this protocol, you can calibrate your qPCR reaction, determine the efficiency of your PCR amplification, and ensure accurate quantification of your target gene.

