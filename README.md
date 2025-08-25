# ABM-SIR Model for Modeling COVID-19 Spread

This repository contains an **Agent-Based Model (ABM)** and **Susceptible–Infected–Recovered (SIR)** framework to study the spread of COVID-19.  
The project was developed as part of my master’s thesis in mathematics, with a focus on **mathematical modeling of infectious diseases**.
The purpose of this model is to understand and extend the modeling process of a disease progression inspired by "H. De-Leon and D. Aran" work.
---

## 📌 ABM Model Functions
- Average :- To average the age specific sub groups in the population
- ComputeInfectious :- To generate unique infection days per variant type.
- HealthStatus :- To initalize the health status in the population.
- VaccineEffectiveness :- To create vaccine and natural induced immunity models.
- MAM :- To control the core process/algorithm of infection probability.
- Simulations :- The driver function to initialize, average, and generate plots.

## 📌 SIR Model Functions
- SIRbasic :- Traditional susceptiblle, infected, recoevered model.
- SIRMultiAge :- Multi age-structured SIR model
- SIRMultiAgeTest :- The driver function to initialize, analyze and generate plots.

---
## Additional Files Used
- VaccinationRate.xlsx contains data on the vaccination rate in Israel as a function of time, divided into age groups.
- VaccRate_3gps :- MATLAB file that plots the vaccination rate data for 4 vaccination doses into age groups.
- Rt_Values.xlsx provides the data to calculate the variant specific theoretical reproduction number values.
- CC_Data.xlsx contains data on the daily confirmed cases and total number of population in Israel. 

---

## 🛠️ Technologies Used
- **Language**:  MATLAB 

---

## Clone the repository
```bash
git clone https://github.com/yashleenss/ABM-SIR-Model-for-Modelling-COVID-19-Spread.git

---

This template is professional, clear, and easy for others to navigate.
You can modify the **Technologies**, **Data**, and **Results** sections depending on the requirement of the disease.  
