import matplotlib.pyplot as plt
import numpy as np

with open("data.txt", 'r') as file:
    data = file.read()

# data = """0.069165 0.095521 0.843737 1.744045
# 0.083086 0.107675 0.87256 1.761941
# 0.065187 0.11778 0.949845 1.893914"""

# Splitting the string into lines, then each line into numbers, and converting to float
data_list = [[float(number) for number in line.split()] for line in data.strip().split('\n')]
# print(data_list)

# Taking transpose
data_list = [[row[i] for row in data_list] for i in range(len(data_list[0]))]

# Initialize a list to hold all performance data and another for the labels
labels = []

# Loop through each combination of n and stencil size
for n in [4096*4096, 8192*8192]:
    for P in [8, 12]:
        for leader in ["w l", "w/o l"]:
        # Create a label for this set of performance data indicating n and stencil size
            labels.append(f"N={n}, P={P}, {leader}")

# Plotting configurations
plt.figure(figsize=(30, 15))  
plt.boxplot(data_list, labels=labels)
plt.xlabel("Configuration")
plt.ylabel("Execution Time (s)")
# plt.yticks(np.arange(0, 40, 0.1))
# plt.xticks(rotation=0, ha="right")  # Rotate x-axis labels for better readability
plt.title("Halo Exchange Comparison Across Different Configurations")
plt.tight_layout()  
# plt.legend()  

# Save the plot as a PNG file
plt.savefig("halo_exchange_comparison.png")