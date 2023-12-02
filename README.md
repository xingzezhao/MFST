# MFST READMD

- version: v 1.0
  
- update time: 12/02/2023
  

Welcome to MFST, where MFST stands for MolFilmStabTool. We provide an end-to-end automated solution for researching the stability of foam liquid films.

To ensure that MFST serves you well, please check if you have the following software or libraries installed:

- Python 3.7 and above
- Numpy
- MDAnalysis
- Subprocess
- Matplotlib
- Docx
- Scipy
- Datetime
- Glob
- Shutil
- Logging

To use MFST correctly, you need to configure your **GROMACS path**, **Packmol path**, and the name of your job management program in the corresponding locations in MFST.rc. These programs should be added to the PATH. Currently supported job management programs include Slurm and PBS. You can easily configure other job management programs by modifying the relevant locations in the program – it's straightforward, so no need to worry.

MFST preserves the option for users to customize structure and force field files for simulations. However, this feature is still under development. In the current version, you can customize structure data by replacing existing structures or force field names in the database or by adding customized structure data using the paths provided in the dictionary in the pyModeling module.

- Core Values

We believe in openness, sharing, and collaboration as outstanding human qualities. We remain open and willing to engage in discussions and improvements with all friends working on similar projects. We hope that these efforts can help more people enter the field of molecular simulation and foam stability research, contributing to the progress of chemical engineering and energy.

- About Me

### Hello, I'm [Xingze Zhao] 👋

- 🏢 I work at China University of Petroleum (East China).
  - College of Materials Science and Engineering
  - Advanced Chemical Engineering and Energy Materials Research Center
  - Lab of Materials Theory Simulation and Design
- 🧪 My research focuses on molecular simulations of nanomaterials and chemicals for enhanced oil recovery (*EOR*).
- ✉️ You can reach me at [xingzezhao@gmail.com](mailto:xingzezhao@gmail.com).

Feel free to explore my research work and contributions!
