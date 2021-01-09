# PhazeOpt

`PhazeOpt` is a set of tools for Reservoir Fluid Charachterization.

This is a work in progress. So, any contribution or advice is welcomed and highly appreciated.

I'm using this project to:

* Deepen my understaning of Reservoir Fluid Charachterization
* Increase and enhance my Programing Skills
* Refurbish my understaning of Numerical Methods and Math in general

## References

* SPE Monograph [PhaseBehavior] (https://store.spe.org/Phase-Behavior-P46.aspx) is the backbone of this work

## Ready tools for Testing

1. `C7+` Charachterization using `Gamma` function using Gas Chromotograph data

## `Gamma` Function

Using `Gamma.py`, it contains the `Gamma` Class.

The `Gamma` Class inuts are:

1. File name in `Excel` format as a string 'data.xlsx', make sure to include the file extension

2. Sheet name which cotains that data table as a string 'Data'

3. Column name containing the Single Carbon Number `SCN` as a string 'SCN'

4. Column name containing the laboratory measured `Molecular Weights` for each `SCN` as a string 'Lab Mi'

5. Column name containing the laboratory measured `Molar Composition` for each `SCN` as a string 'mol-%'

6. Optionally, you can specify the `h` value to be used for calculating the 'Charachterized Molecular Weights' as an integer format, default value is -4
