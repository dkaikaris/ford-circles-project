from sympy import *
import matplotlib.pyplot as plt
from math import sqrt

# Set upper bound for q
q_max = 30

''' GENERATING FRACTIONS '''
L = []

# Generate unique fractions with 1 <= p <= q <= q_max
for q in range(2, q_max + 1):  
    for p in range(1, q): 
        fraction = Rational(p, q)
        L.append((fraction, (fraction.p, fraction.q)))

L.append((Rational(1/1), (1, 1)))

# Remove duplicate fractions and sort fractions by ascending order
uniqueFractions = {item[0]: item for item in L} 
L = list(uniqueFractions.values())
L.sort(key=lambda x: x[0])


''' GENERATING FAREY SUMS '''
fareyCount = 0

for i in range(len(L) - 2):
    fraction1, (p, q) = L[i]
    fraction2, (farey_p, farey_q) = L[i+1]
    fraction3, (p1, q2) = L[i+2]
    if Rational(p + p1, q + q2) == Rational(farey_p, farey_q):
        fareyCount += 1


''' PLOTTING FORD CIRCLES '''
fig, ax = plt.subplots()

# Title + axis labels (added)
ax.set_title("Ford Circles for Farey Fractions up to q = 30")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Plot base circle at 0/1 since 0/1 not included in list of fractions at first
circle = plt.Circle((0, 1/2), (1/2), fill=False)
ax.add_patch(circle)

# Plot all Ford circles
for frac, (p, q) in L:
    circle = plt.Circle((p/q, 1/(2*q**2)), (1/(2*q**2)), fill=False)
    ax.add_patch(circle)

ax.set_aspect('equal', adjustable='box')

plt.show()
fig.savefig("farey_ford_circles.png", dpi=300, bbox_inches='tight')


''' GOLDEN RATIO (phi) WITH FORD CIRCLES '''

phi_minus_1 = (sqrt(5) - 1) / 2

# Generating a few Fibonacci numbers and Fibonacci fractions
fib = [1, 1]
for i in range(10):
    fib.append(fib[-1] + fib[-2])

fib_fracs = []
for i in range(1, 7):
    p = fib[i]
    q = fib[i+1]
    fib_fracs.append((Rational(p, q), (p, q)))

# Plot Ford circles that correspond to Fibonacci fractions
fig2, ax2 = plt.subplots(figsize=(11, 4))

ax2.axvline(phi_minus_1, color='gold', linestyle='--', linewidth=1)
ax2.text(phi_minus_1 + 0.0067, 0.11, "φ−1", color='gold', fontsize=12)

# Using colors and labels for ford circles for better visualisation
cmap = plt.cm.viridis
colors = [cmap(i / len(fib_fracs)) for i in range(len(fib_fracs))]

legend_handles = []
legend_labels = []

for i, ((frac, (p, q))) in enumerate(fib_fracs):
    x = p/q
    r = 1/(2*q**2)

    circle = plt.Circle((x, r), r, fill=False, color=colors[i], linewidth=2)
    ax2.add_patch(circle)

    handle = plt.Line2D([0], [0], color=colors[i], linewidth=3)
    label = f"{p}/{q} (radius = {r:.5f})"
    legend_handles.append(handle)
    legend_labels.append(label)

ax2.set_xlim(phi_minus_1 - 0.17, phi_minus_1 + 0.17)
ax2.set_ylim(0, 0.12)

ax2.set_title("Ford Circles Approaching the Conjugate Golden Ratio (φ−1)")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

ax2.set_aspect('equal', adjustable='box')
ax2.legend(legend_handles, legend_labels, fontsize=10, loc='upper right')

plt.show()
fig2.savefig("golden_ratio_ford_circles.png", dpi=300, bbox_inches='tight')
