#include <iostream>
#include <vector>
#include <string>
using namespace std;

class HMM {

public:
	string DNAseq;
	int DNALen;
	vector<int> codingRegion;
	int exonNum;
	int exonCnt;
	double HMMprob[2][4]; // Exon/Intron, A/T/C/G
	double StateChange[2];

	HMM(int len) {
		DNALen = len;
		codingRegion.resize(len);
		exonCnt = 0;
	}

	void makeHMMprob() {

		int countNeucleo[2][4] = { 0, }; // Exon/Intron, A/T/C/G

		StateChange[0] = (exonNum) / (double)DNALen;
		StateChange[1] = (exonNum) / (double)DNALen;

		for (int i = 0; i < DNALen; i++) {
			if (codingRegion[i] == 0) {
				switch (DNAseq[i]) {
				case 'A':
					countNeucleo[0][0] += 1;	//probability for detecting A at intron
					break;
				case 'T':
					countNeucleo[0][1] += 1;	//probability for detecting T at intron
					break;
				case 'C':
					countNeucleo[0][2] += 1;	//probability for detecting C at intron
					break;
				case 'G':
					countNeucleo[0][3] += 1;	//probability for detecting G at intron
					break;
				default:
					cout << "Not a neucleotide expected!!" << endl;
					exit(1);
				}
			}
			else {
				switch (DNAseq[i]) {

				case 'A':
					countNeucleo[1][0] += 1;	//probability for detecting A at exon	
					break;
				case 'T':
					countNeucleo[1][1] += 1;	//probability for detecting T at exon
					break;
				case 'C':
					countNeucleo[1][2] += 1;	//probability for detecting C at exon
					break;
				case 'G':
					countNeucleo[1][3] += 1;	//probability for detecting AG at exon
					break;
				default:
					cout << "Not a neucleotide expected!!" << endl;
					exit(1);
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			HMMprob[0][i] = countNeucleo[0][i] / (double)(DNALen - exonCnt);
			HMMprob[1][i] = countNeucleo[1][i] / (double)(exonCnt);
		}
	}

	void printHMMProb() {
		cout << "==================< HMM Probability >==================" << endl;
		cout << "Genome length: " << DNALen << endl;
		cout << "The number of exon sites: " << exonCnt << endl;
		cout << "The number of exon: " << exonNum << endl;
		for (int i = 0; i < 2; i++) {
			if (i == 1) cout << "Exon A|T|C|G prob: \t";
			else cout << "Intron A|T|C|G prob: \t";

			for (int j = 0; j < 4; j++) {
				cout << HMMprob[i][j] << '/t' << "|";
			}
			cout << endl;
		}
		cout << "Probabilty of exon to intron: " << StateChange[0] << endl;
		cout << "Probabilty of intron to exon: " << StateChange[1] << endl;
		cout << "=======================================================" << endl;
	}
};