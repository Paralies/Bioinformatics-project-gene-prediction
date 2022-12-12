#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "HMM.cpp"
#include "read_file.cpp"
using namespace std;

void tracking_iterate(string& genome, int state, int** backtracking, int squence_len)
{
	int state_prev = state;
	for (int i = squence_len - 1; i > 0; i--)
	{
		genome += (char)(48 + backtracking[state_prev][i]);	//backtraking to generate genome
		state_prev = backtracking[state_prev][i];			//state backtracking
	}
	reverse(genome.begin(), genome.end());
	cout << "Genome: " << genome << endl << endl;
	for (int i = 0; i < genome.length() - 1; i++)
	{
		if (genome[i] == '0' && genome[i + 1] == '1')
			cout << "exon start : " << i + 1 << endl;		//print exon start site
		if (genome[i] == '1' && genome[i + 1] == '0')
			cout << "exon end : " << i + 1 << endl << endl;	//print exon end site
	}
}

class viterbi
{
public:
	long double s0, s1;
	double NtoE, EtoN, EtoE, NtoN;
	int DNA_len;

public:
	viterbi(int len, double* StateChange)
	{
		NtoE = StateChange[1];		//intron to exon probability
		NtoN = 1 - StateChange[1];	//intron to intron probability
		EtoN = StateChange[0];		//exon to intron probability
		EtoE = 1 - StateChange[0];	//exon to exon probability
		DNA_len = len;
	}

	int convert(char nucleotide)
	{
		if (nucleotide == 'A') return 0;
		if (nucleotide == 'T') return 1;
		if (nucleotide == 'C') return 2;
		if (nucleotide == 'G') return 3;
		return -1;
	}

	void calculate(string DNAseq, HMM& DNA_1)
	{

		string sequence = DNAseq;
		s0 = log(0.99 * DNA_1.HMMprob[0][convert(sequence[0])]); // Intron
		s1 = log(0.01 * DNA_1.HMMprob[1][convert(sequence[0])]); // Exon

		int* backtracking[2];
		backtracking[0] = new int[sequence.length()]; // Intron
		backtracking[1] = new int[sequence.length()]; // Exon

		for (int i = 0; i < sequence.length(); i++)
		{
			backtracking[0][i] = 0;
			backtracking[1][i] = 0;
		}
		
		long double s0_prev;
		long double s1_prev;

		// Making Viterbi graph
		for (int i = 1; i < sequence.length(); i++)
		{
			s0_prev = s0;
			s1_prev = s1;

			if (s0_prev + log(NtoN * DNA_1.HMMprob[0][convert(sequence[i])]) >= s1_prev + log(EtoN * DNA_1.HMMprob[0][convert(sequence[i])])) // pre-intron to post-intron
			{
				s0 = s0_prev + log(NtoN * DNA_1.HMMprob[0][convert(sequence[i])]);
				backtracking[0][i] = 0;
			}
			else if (s0_prev + log(NtoN * DNA_1.HMMprob[0][convert(sequence[i])]) < s1_prev + log(EtoN * DNA_1.HMMprob[0][convert(sequence[i])])) // pre-exon to post-intron
			{
				s0 = s1_prev + log(EtoN * DNA_1.HMMprob[0][convert(sequence[i])]);
				backtracking[0][i] = 1;
			}

			if (s0_prev + log(NtoE * DNA_1.HMMprob[1][convert(sequence[i])]) >= s1_prev + log(EtoE * DNA_1.HMMprob[1][convert(sequence[i])])) // pre intron to post-exon
			{
				s1 = s0_prev + log(NtoE * DNA_1.HMMprob[1][convert(sequence[i])]);
				backtracking[1][i] = 0;
			}
			else if (s0_prev + log(NtoE * DNA_1.HMMprob[1][convert(sequence[i])]) < s1_prev + log(EtoE * DNA_1.HMMprob[1][convert(sequence[i])])) // pre-exon to next-exon
			{
				s1 = s1_prev + log(EtoE * DNA_1.HMMprob[1][convert(sequence[i])]);
				backtracking[1][i] = 1;
			}			
		}
		


		// Show exon-intron part
		string genome = "";
		if (s0 > s1) // intron > exon
		{
			//cout << "Start s0" << endl;
			genome += '0';
			tracking_iterate(genome, 0, backtracking, sequence.length());
		}
		else if (s0 < s1) // exon > intron
		{
			//cout << "Start s1" << endl;
			genome += '1';
			tracking_iterate(genome, 1, backtracking, sequence.length());
		}
	}
};

int main()
{
	ReadGenomeFile CAGenome; 
	ReadGenomeFile yeastGenome; 
	
	ifstream CA_genome_file;
	CA_genome_file.open("CA_genome.fa"); // open the file of genome in Candida albican
	CAGenome.readLines(CA_genome_file);

	ifstream yeast_genome_file;
	yeast_genome_file.open("training.fa"); // open the file of genome in Baker's yeast
	yeastGenome.readLines(yeast_genome_file);

	ifstream yeast_CDS_file;
	yeast_CDS_file.open("training_CDS.txt");

	ReadGeneFile yeastCDS(yeastGenome.getDNALength());
	yeastCDS.getExonLocInFIle(yeast_CDS_file);

	HMM yeast_HMM(yeastGenome.getDNALength());	

	yeastGenome.syncInfoInHMM_Genome(yeast_HMM.DNAseq);
	yeastCDS.syncInfoInHMM_CDS(yeast_HMM.exonNum, yeast_HMM.exonCnt, yeast_HMM.codingRegion);
	
	yeast_HMM.makeHMMprob();	//calulate HMM probability
	yeast_HMM.printHMMProb();	//print HMM probability

	viterbi viter(yeastGenome.getDNALength(), yeast_HMM.StateChange);	
	viter.calculate(CAGenome.getDNASequence(), yeast_HMM);	//calculate viterbi
}