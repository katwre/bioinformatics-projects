"""A solution to the Rosalind problem "Calculating Expected Offspring."
http://rosalind.info/problems/iev/
a. AA-AA
b. AA-Aa
c. AA-aa
d. Aa-Aa
e. Aa-aa
f. aa-aa
"""
import sys

with open(sys.argv[1]) as f:
    a,b,c,d,e,f = map(int, f.read().split())
    expected =	a * 1 + b * 1 + c * 1 + d * 0.75 + e * 0.5 + f * 0
    print(2 * expected)
