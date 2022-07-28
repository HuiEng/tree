>>>data.fasta

pick random 42 seqs from silva dataset, 
run all-against-all needleman-wunsch with "needle", all score < 80

from the 42 seqs, generate toy dataset where each seq is introduced to substitution rate of:
0.01 to 0.05 (step = 0.01), each rate has 5 copies (excluding original)


tree pbuild -i data.fasta -o data.bin -k 9 -w 50 -s 5
tree tree -i data.bin

