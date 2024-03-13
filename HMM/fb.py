from HMM import *

DRAW = False
if DRAW:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(nrows=3, figsize=(30, 15), dpi=150)

test_sequence = 'gataggattatcattcataagtttcagagcaatgtccttattctggaacttggatttatggctcttttggtttaatttcgcctgattcttgatctcctttagcttctcgacgtgggcctttttcttgccatatggatccgctgcacggtcctgttccctagcatgtacgtgagcgtatttccttttaaaccacgacgctttgtcttcattcaacgtttcccattgtttttttctactattgctttgctgtgggaaaaacttatcgaaagatgacgactttttcttaattctcgttttaagagcttggtgagcgctaggagtcactgccag'
test_sequence = test_sequence.upper()

def run_fb(model: HMM, sequence):
    res_fb = model.forward_backward(sequence)
    # to draw sites better
    # {'S0': [0.7, 0, 0], 'S1': [0, 0.7, 0], 'S2': [0, 0, 0.7]} -> [0.7 0.7, 0.7]
    site_probas = {i: j for i, j in res_fb.items() if i.startswith("S")}
    res_fb["S"] = list(map(max, zip(*site_probas.values())))
    return res_fb


if __name__ == "__main__":
    x = range(len(test_sequence))
    for i, p_site in enumerate([0.01, 0.05, 0.0001]):
        model = create_restriction_model("ACGT", p_site=p_site)
        res_fb = run_fb(model, test_sequence)
        if DRAW:
            ax[i].plot(x, res_fb["S"], label="Site")
            ax[i].plot(x, res_fb["B"], label="Background")
            ax[i].set_xticks(x, test_sequence, fontsize=3.6)
            ax[i].set_xlim([min(x), max(x)])
            ax[i].legend(prop={'size':10})
            ax[i].text(5, 0.5, f"p_site = {p_site}")
        else:
            print(f"p_site: {p_site}, FB probabilities: {res_fb}")
    if DRAW:
        plt.show()