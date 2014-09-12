for i in ''{500..750..25}
	do
	>&2 echo "Generating data for $i"
	rm -f tests.fasta
	python3 gen.py acgt $i > tests.fasta
	>&2 echo "Running linear with size $i"	
	perf stat -r 10 -x, python3 exc2test.py exc2.cost tests.fasta seq1 seq2 > /dev/null
	>&2 echo "Running affine with size $i"	
	perf stat -r 10 -x, python3 exc2test.py exc2.cost tests.fasta seq1 seq2 --affine > /dev/null
	done
