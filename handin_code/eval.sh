echo "Question 1:" 
python3 exc2test.py exc2.cost eval.fasta seq1 seq2 --backtrack
echo
echo "Question 2:" 
python3 exc2test.py exc2.cost eval.fasta seq1 seq2 --backtrack --affine
echo
echo "Question 3:" 
python3 exc2test.py exc2.cost eval.fasta seq1 seq2 --backtrack --matrix
echo
echo "Question 4:" 
python3 exc2test.py exc2.cost eval.fasta seq1 seq2 --backtrack --affine --matrix

