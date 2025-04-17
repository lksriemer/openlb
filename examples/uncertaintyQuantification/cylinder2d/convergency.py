import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = pd.read_csv('drag_coefficients.txt', sep='\t')
mc = pd.read_csv('mc_drag_coefficients.txt', sep='\t')

N = 10  # Number used for naming output files

# Display the data
print(data)

# Extract the columns
orders = data['Order'].values
nqs = data['Nq'].values
mean_drag = data['Mean_Drag'].values
std_drag = data['Std_Drag'].values


orders_mc = mc['Order'].values
nqs_mc = mc['Nq'].values
mean_drag_mc = mc['Mean_Drag'].values
std_drag_mc = mc['Std_Drag'].values

# Reference values from the highest order (last element)
reference_mean = mean_drag[-1]
reference_std = std_drag[-1]

reference_mean_mc = mean_drag_mc[-1]
reference_std_mc = std_drag_mc[-1]

# Compute relative errors (exclude the last term)
relative_error_mean = np.abs((mean_drag[:-1] - reference_mean) / np.abs(reference_mean))
relative_error_std = np.abs((std_drag[:-1] - reference_std) / np.abs(reference_std))

relative_error_mean_mc = np.abs((mean_drag_mc[:-1] - reference_mean_mc) / np.abs(reference_mean_mc))
relative_error_std_mc = np.abs((std_drag_mc[:-1] - reference_std_mc) / np.abs(reference_std_mc))
# Create arrays excluding the last data point for relative errors
nqs_rel_error = nqs[:-1]
nqs_rel_error_mc = nqs_mc[:-1]

# Function to create tick labels in the form 2^n
def format_ticks_as_powers_of_two(values):
    return [f'$2^{{{int(np.log2(val))}}}$' for val in values]


# Plot Relative Error of Mean Drag vs. Nq (excluding last term)
plt.figure(figsize=(8, 6))
plt.plot(nqs_rel_error, relative_error_mean, marker='o', linestyle='-', color='b', label='SC LBM')
plt.plot(nqs_rel_error_mc, relative_error_mean_mc, marker='o', linestyle='-', color='r', label='MC LBM')
plt.xlabel('Nq', fontsize='large')
plt.ylabel(r'$\delta(\bar{C_D})$', fontsize='large')
plt.yscale("log")
plt.grid(True)
# plt.xticks(nqs_rel_error, format_ticks_as_powers_of_two(nqs_rel_error))
plt.legend(fontsize='large')
plt.savefig('cylinder_relative_error_mean_drag_vs_nq.png')
plt.close()

# Plot Relative Error of Std Drag vs. Nq (excluding last term)
plt.figure(figsize=(8, 6))
plt.semilogy(nqs_rel_error, relative_error_std, marker='o', linestyle='-', color='b', label='SC LBM')
plt.semilogy(nqs_rel_error_mc, relative_error_std_mc, marker='o', linestyle='-', color='r', label='MC LBM')
plt.title('Relative Error of Std Drag vs. Nq')
plt.xlabel('Nq', fontsize='large')
plt.ylabel(r'$\delta(\sigma(C_D))$', fontsize='large')
plt.yscale("log")
plt.grid(True)
# plt.xticks(nqs_rel_error, format_ticks_as_powers_of_two(nqs_rel_error))
plt.legend(fontsize='large', loc='lower right')
plt.savefig('cylinder_relative_error_std_drag_vs_nq.png')
plt.close()