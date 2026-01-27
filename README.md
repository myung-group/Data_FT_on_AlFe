# Impurity-driven spontaneous and infrared-activated Fischer–Tropsch chemistry on commercial aluminum

This repository accompanies the paper:

> **"Impurity-driven spontaneous and infrared-activated Fischer–Tropsch chemistry on commercial aluminum"**

**Authors:** Seon Young Hwang, Seungwon Kim, Jihye Lee, Gi Beom Sim, Go Eun Park, Soohaeng Yoo Willow, Chang Woo Myung, and Youngku Sohn

## 📂 Repository Structure

* **`Data/`**: Contains initial input data of 1) geometry relaxation and 2) vibrational frequency calculation of reaction intermediates.
    * Calculated data file `FT_on_AlFe.xlsx` for Gibbs free energy exists in `Data/`.
* **`Figures/`**: Contains figures of the paper that were obtained from `Data/`.
* **`Scripts/`**: Contains Jupyter notebook and Python scripts that can reproduce figures.
    * Zero-point energy, entropy, and enthalpy corrections can be calculated using the Python scripts `Scripts/calculate_*.py`.

---

### 💾 Output Data
If you want to check your calculation, please visit our output data of the project.
* **Zenodo Link:** [Insert Link Here](---)

---

## ⚙️ Computational Details

### Pseudopotentials
The following PAW potentials (PBE) and BEEF-vdW functionals were used in this study.

> **Note 1:** Due to licensing restrictions, the POTCAR files are not included in this repository.

| Element | Label | Date (Version) |
| :--- | :--- | :--- |
| **Al** | PAW_PBE Al | 04Jan2001 |
| **Fe** | PAW_PBE Fe | 06Sep2000 |
| **C** | PAW_PBE C | 08Apr2002 |
| **O** | PAW_PBE O | 08Apr2002 |
| **H** | PAW_PBE H | 15Jun2001 |

### Functionals
> **Note 2:** Due to licensing restrictions, the BEEF-vdW functionals (`libbeef`) are not included in this repository.

Please visit the repository below and install it when you want to proceed with your project:
* [https://github.com/vossjo/libbeef](https://github.com/vossjo/libbeef)
