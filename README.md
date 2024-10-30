This work mainly discusses the aggregation characteristics of charged peptide chains in ionic environments. The relevant files are explained as follows:

The AIMD folder contains two initial inp files (2-Arg system and 2-Asp system), and a program for analyzing the distance between sodium ions and oxygen atoms.

The GaMD folder contains initial structure file, force field file, simulation control file, as well as a program for performing reweighting analysis.

The coordinates_and_ff folder contains the initial structural files and force field files for all two-chain systems.

The cMD folder contains minimization, equilibrium, and productivity control files for cMD simulations of all two-chain systems.

The analysis folder contains all the analysis programs, which are responsible for the following tasks:
1. bridge_ion.cpp is used to calculate the average number and standard deviation of ions shared between two side-chain charge centers;
2. lifetime.cpp is used to calculate the lifetime of complementary residues after forming salt bridge contacts;
3. main_Side_comr.cpp is used to calculate the motion correlation between backbone CA atoms and side-chain charge centers;
4. s_b_series.cpp is used to calculate the evolution of the number of salt bridges in a single trajectory over time;
5. stacking_env.cpp is used to analyze the stacking interactions formed between guanidine groups on Arg side-chains;
6. v_ion_env.cpp is used to analyze the coupling between the velocity of ions and their surrounding side-chain environment.
