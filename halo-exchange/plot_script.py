import matplotlib.pyplot as plt

times = {
    (512**2, 5): [0.087172, 0.096136, 0.083290],
    (512**2, 9): [0.128304, 0.136145, 0.139221],
    (2048**2, 5): [1.353963, 0.844899, 0.857837],
    (2048**2, 9): [1.950937, 1.444299, 1.472888]
}

N_values = sorted(set(N for N, _ in times.keys()))
stencil_values = sorted(set(stencil for _, stencil in times.keys()))

fig, axes = plt.subplots(len(stencil_values), figsize=(10, 8), sharex=True)

for i, stencil in enumerate(stencil_values):
    data = [times[(N, stencil)] for N in N_values]
    
    axes[i].boxplot(data)
    axes[i].set_title(f'Stencil {stencil}')
    axes[i].set_xlabel('Data Size (N)')
    axes[i].set_ylabel('Time (seconds)')
    axes[i].set_xticks(range(1, len(N_values) + 1))
    axes[i].set_xticklabels([str(N) for N in N_values])

plt.show()
