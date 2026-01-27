This repository accompanies the paper "Impurity-driven spontaneous and infrared-activated Fischer–Tropsch chemistry on commercial aluminum" by Seon Young Hwang, Seungwon Kim, Jihye Lee, Gi Beom Sim, Go Eun Park, Soohaeng Yoo Willow, Chang Woo Myung, and Youngku Sohn
    • Data/ contains initial input data of 1) geometry relaxation and 2) vibrational frequency calculation of reaction intermediates.
        ◦ Calculated data file FT_on_AlFe.xlsx for Gibbs free energy exists in Data/.
    • Figures/ contains figures of the paper that obtained from Data/.
    • Scripts/ contains Jupyter notebook and Python scripts that can reproduce figures.
        ◦ Zero-point energy, entropy, and enthalpy corrections can be calculated using the Python scripts Scripts/calculate_*.py.
If you want to check your calculation, please visit our output data of project. (Zenodo link: --- )


## Pseudopotentials The following PAW potentials (PBE) and BEEF-vdW functionals were used in this study. 
Note 1: Due to licensing restrictions, the POTCAR files are not included in this repository. 
| Element | Label | Date (Version) | 
| Al | PAW_PBE Al | 04Jan2001 |
| Fe | PAW_PBE Fe | 06Sep2000 | 
| C | PAW_PBE C | 08Apr2002 | 
| O | PAW_PBE O | 08Apr2002 |
| H | PAW_PBE H | 15Jun2001 |

Note 2: Due to licensing restrictions, the BEEF-vdW functionals (libbeef) are not included in this repository. Please visit this repository, and install it when you want to do your project.
https://github.com/vossjo/libbeef


Paper abstract
Fischer–Tropsch (F–T) synthesis conventionally requires high temperatures, elevated pressures, and externally supplied H2 to convert CO into hydrocarbons. Here, we report a fundamentally distinct ambient-temperature analogue: spontaneous F–T chemistry occurring on commercial aluminum (Al) in aqueous solution. This reaction is driven by Al corrosion, which supplies electrons and surface hydrogen (*H), while trace Fe impurities embedded in the Al matrix act as catalytic microdomains that stabilize *CO and promote C–C coupling. Under CO-saturated alkaline conditions, this self-sustained Al–Fe redox system produces H2, CH4, and long-chain C2–C7 hydrocarbons with product distributions following Anderson–Schulz–Flory kinetics, whereas ultrapure Al yields only H2. Infrared irradiation further enhances hydrocarbons via localized surface restructuring and transient photothermal heating. Density functional theory identifies the Al–Fe interface as the critical site for CO adsorption and *CH2-mediated chain propagation, establishing spontaneous redox interfaces as a platform for low-temperature CO valorization.
