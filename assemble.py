import numpy as np
import math
import matplotlib.pyplot as plt
import networkx as nx
import sys

import dimod
import minorminer
# import dwave_networkx as dnx
# from dwave.cloud import Client
from dwave.embedding import embed_ising, unembed_sampleset
from dwave.embedding.utils import edgelist_to_adjacency
# from dwave.system.samplers import DWaveSampler
from dwave.embedding.chain_breaks import majority_vote


"""
Overlap between pair-wise reads
"""
def align(read1, read2, mm):
    l1, l2 = len(read1), len(read2)
    for shift in range(l1 - l2, l1):
        mmr = 0
        r2i = 0
        for r1i in range(shift, l1):
            if read1[r1i] != read2[r2i]:
                mmr += 1
                r2i += 1
            if mmr > mm:
                break
        if mmr <= mm:
            return l2 - shift
    return 0


"""
Convert set of reads to adjacency matrix of pair-wise overlap for TSP
"""


def reads_to_tspAdjM(reads, max_mismatch=0):
    n_reads = len(reads)
    O_matrix = np.zeros((n_reads, n_reads))  # edge directivity = (row id, col id)
    for r1 in range(0, n_reads):
        for r2 in range(0, n_reads):
            if r1 != r2:
                O_matrix[r1][r2] = align(reads[r1], reads[r2], max_mismatch)
                O_matrix = O_matrix / np.linalg.norm(O_matrix)
    return O_matrix


"""
Convert adjacency matrix of pair-wise overlap for TSP to QUBO matrix of TSP
"""


def tspAdjM_to_quboAdjM(tspAdjM, p0, p1, p2):
    n_reads = len(tspAdjM)
    # Initialize
    Q_matrix = np.memmap('large_array.dat', dtype=np.float32, mode='w+', shape=(n_reads**2, n_reads**2))
    # Q_matrix = np.zeros((n_reads**2, n_reads**2))
    # Qubit index semantics: {c(0)t(0) |..| c(i)-t(j) | c(i)t(j+1) |..| c(i)t(n-1) | c(i+1)t(0) |..| c(n-1)t(n-1)}
    # Assignment reward (self-bias)
    p0 = -1.6
    for ct in range(0, n_reads**2):
        Q_matrix[ct][ct] += p0
        # Multi-location penalty
        p1 = -p0  # fixed emperically by trail-and-error
    for c in range(0, n_reads):
        for t1 in range(0, n_reads):
            for t2 in range(0, n_reads):
                if t1 != t2:
                    Q_matrix[c * n_reads + t1][c * n_reads + t2] += p1
                    # Visit repetation penalty
                    p2 = p1
    for t in range(0, n_reads):
        for c1 in range(0, n_reads):
            for c2 in range(0, n_reads):
                if c1 != c2:
                    Q_matrix[c1 * n_reads + t][c2 * n_reads + t] += p2
                    # Path cost
                    # kron of tspAdjM and a shifted diagonal matrix
    for ci in range(0, n_reads):
        for cj in range(0, n_reads):
            for ti in range(0, n_reads):
                tj = (ti + 1) % n_reads
                Q_matrix[ci * n_reads + ti][cj * n_reads + tj] += -tspAdjM[ci][cj]
    print(Q_matrix)
    return Q_matrix


"""
Convert QUBO matrix of TSP to QUBO dictionary of weighted adjacency list
"""


def quboAdjM_to_quboDict(Q_matrix):
    n_reads = int(math.sqrt(len(Q_matrix)))
    Q = {}
    for i in range(0, n_reads**2):
        ni = "n" + str(int(i / n_reads)) + "t" + str(int(i % n_reads))
        for j in range(0, n_reads**2):
            nj = "n" + str(int(j / n_reads)) + "t" + str(int(j % n_reads))
            if Q_matrix[i][j] != 0:
                Q[(ni, nj)] = Q_matrix[i][j]
    print(Q)
    return Q


"""
Solve a QUBO model using dimod exact solver
"""


def solve_qubo_exact(Q, all=False):
    solver = dimod.ExactSolver()
    response = solver.sample_qubo(Q)
    minE = min(response.data(["sample", "energy"]), key=lambda x: x[1])
    for sample, energy in response.data(["sample", "energy"]):
        if all or energy == minE[1]:
            print(sample)


"""
Solve an Ising model using dimod exact solver
"""


def solve_ising_exact(hii, Jij, plotIt=False):
    solver = dimod.ExactSolver()
    response = solver.sample_ising(hii, Jij)
    print("Minimum Energy Configurations\t===>")
    minE = min(response.data(["sample", "energy"]), key=lambda x: x[1])
    for sample, energy in response.data(["sample", "energy"]):
        if energy == minE[1]:
            print(sample, energy)
    if plotIt:
        y = []
        for sample, energy in response.data(["sample", "energy"]):
            y.append(energy)
        plt.plot(y)
        plt.xlabel("Solution landscape")
        plt.ylabel("Energy")
        plt.savefig("ising.png")
        plt.show()
    # print(hii)
    # print(Jij)

def solve_ising_dwave(hii,Jij):
	config_file='QA_DeNovoAsb/dwcloud.conf'
	client = Client.from_config(config_file, profile='aritra')
	solver = client.get_solver() # Available QPUs: DW_2000Q_2_1 (2038 qubits), DW_2000Q_5 (2030 qubits)
	dwsampler = DWaveSampler(config_file=config_file)

	edgelist = solver.edges
	adjdict = edgelist_to_adjacency(edgelist)
	embed = minorminer.find_embedding(Jij.keys(),edgelist)
	[h_qpu, j_qpu] = embed_ising(hii, Jij, embed, adjdict)

	response_qpt = dwsampler.sample_ising(h_qpu, j_qpu, num_reads=solver.max_num_reads())
	client.close()

	bqm = dimod.BinaryQuadraticModel.from_ising(hii, Jij)
	unembedded = unembed_sampleset(response_qpt, embed, bqm, chain_break_method=majority_vote)
	print("Maximum Sampled Configurations from D-Wave\t===>")
	solnsMaxSample = sorted(unembedded.record,key=lambda x: -x[2])
	for i in range(0,10):
		print(solnsMaxSample[i])
	print("Minimum Energy Configurations from D-Wave\t===>")
	solnsMinEnergy = sorted(unembedded.record,key=lambda x: +x[1])
	for i in range(0,10):
		print(solnsMinEnergy[i])

def ising_solver(hii, Jij):
    solver = dimod.SimulatedAnnealingSampler()
    # edgelist = solver.edges
    # adjdict = edgelist_to_adjacency(edgelist)
    # embed = minorminer.find_embedding(Jij.keys(),edgelist)
    # [h_qpu, j_qpu] = embed_ising(hii, Jij, embed, adjdict)

    bqm = dimod.BinaryQuadraticModel.from_ising(hii, Jij)
    print(bqm)
    sampler = solver.sample(bqm)
    print(sampler)
    print(type(sampler))

    for sample, energy in sampler.data(fields=['sample', 'energy']):
        print(sample, energy)


# r1 = "NAACCTCTCTGTTTACTGATAAGTTCCAGATCCTCCTGGCAACTTGCACAAGTCCGACAACCCTGAACGACCAGGCGTCTTCGTTCATCTATCGGATCTCCACACTCACAACAATGAGTGGCAGATATAGCCTGGTGGTTCAGGCGGCGCA"
# r2 = "NGCACGGATGCTACACGAACCTGATGAACAAACTGGATACGATTGGATTCGACAACAAAAAAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAACACTTTTCAGCTAC"
# r3 = "NGGATTGTCGGGAGTATCGGCAGCGCCATTGGCGGGGCTGGTTGTGGGGGCGCCTCCCCGCCCCGCGGGACAACCCCTCAGGCCCCCGCGGCGGAAATTCCTTTTTTTAACCGAGGGGTTTACTGGACCCGGATGTGGGCTTTTCCACCAC"
# r4 = "NTGGACATGGATACCCCGTGAGTTACCCGGCGGGCGCGCCTCGTTCATTCACGTTTTTGAACCCGTGGAGGACGGGCAGACTCGCGGTGCAAATGTGTTTTACAGCGTGATGGAGCAGATGAAGATGCTCGACACGCTGCAGAACACGCAG"

r1 = "ATGCGTG"
r2 = "CGTAGCA"
r3 = "ACTTCAG"
r4 = "CAGCTAG"
r5 = "CGATCAG"
r6 = "TGTGCAA"
r7 = "CAGCTAG"

def deNovo_locally():
    reads = [ r1, r2, r3, r4 ]
    tspAdjM = reads_to_tspAdjM(reads)
    quboAdjM = tspAdjM_to_quboAdjM(tspAdjM, -1.6, 1.6, 1.6)
    quboDict = quboAdjM_to_quboDict(quboAdjM)
    hii, Jij, offset = dimod.qubo_to_ising(quboDict)
    solve_ising_exact(hii, Jij, plotIt=True)
    # solve_ising_dwave(hii,Jij)
    ising_solve(hii, Jij)
    print(f"Reads used:\n{reads}")

deNovo_locally()
