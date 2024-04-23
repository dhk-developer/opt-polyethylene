# opt-polyethylene
  Project for In2ResearchUK - in collaboration with University of Cambridge.

</br>

## Polyethylene Structure Generator

This Python package generates random polyethylene (PE) structures, which are useful for initial configurations in molecular dynamics simulations or for exploring structural diversity in polymer research.
Features

  1. Generates random PE chains with conventional bonding angles.
  2. Randomly places chains in a randomised scaling box without overlaps.
  3. Applies random rotations to increase structural diversity.
  4. Relaxes structures using the AIREBO potential via LAMMPS.
  5. Outputs structures in various formats for visualization and further analysis.
</br>

## Installation

  Clone this repository:
    
    ``` git clone https://github.com/dhk-developer/opt-polyethylene.git ```

  Navigate to the cloned directory:
    
    ``` cd opt-polyethylene ```

  Install the required dependencies:
    
    ``` pip install numpy ase ```

  Note: You will need to have LAMMPS installed and properly configured for the relaxation step. A Python-native version of LAMMPS exists but this has not been thoroughly tested.
</br>

## Usage

  Run the main script, passing the pressure as the first command-line argument and the number of runs as the second argument. 
    
    ``` python run_generator.py 0.1 10 # Pressure in pascals``` 

  Generates both pre-converged .CIF files and converged .RES / .CASTEP for visualisation (VESTA recommended). Converged files contain information on structural enthalpy and most stable configuration based on the force-field 
  used. By default - the AIREBO potential is used, but this can be changed to different force-field inputs.
</br>

## Requirements

  NumPy  
  ASE  
  LAMMPS  
  AIRSS (Requires gcc and gfortran).
</br>

## License

  This project is licensed under the MIT License - see the LICENSE file for details.
</br>

## Acknowledgments
  
  This code was inspired by work done by Professor Chris Pickard from the Department of Materials Science and Metallurgy in the University of Cambridge. This work is based on the AIRSS formulation that generates a robust      first-principle derivation of crystal structures from randomised structure starting points.
