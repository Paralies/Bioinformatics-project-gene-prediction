#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class ReadGenomeFile {
	string DNASeq;
	int DNALen;

public:
	void readLines(ifstream& genome) {
		string tmp_line;
		while (getline(genome, tmp_line)) {
			DNASeq += tmp_line;
		}

		DNALen = DNASeq.size();
	}

	int getDNALength() {
		return DNALen;
	}

	string getDNASequence() {
		return DNASeq;
	}

	void syncInfoInHMM_Genome(string& DNAseqInHMM) {
		DNAseqInHMM = DNASeq;
	}
};

class ReadGeneFile {
	int exonNum = 0;
	int exonCnt = 0;
	vector<int> codingRegion;

public:
	ReadGeneFile(int DNALen) {
		codingRegion.resize(DNALen);
	}

	void getExonLocInFIle(ifstream& CDS) {

		int exonStart;
		int exonEnd;
		string tmp_line;
		string l_site;
		string r_site;
		int curLoc = 0;
		
		getline(CDS, tmp_line);
		exonNum = stoi(tmp_line);

		while (getline(CDS, tmp_line)) {

			curLoc = 0;
			l_site.clear();
			r_site.clear();

			while (tmp_line[curLoc] != '.') {
				l_site += tmp_line[curLoc];
				curLoc++;
			}
			exonStart = stoi(l_site) - 1;
			curLoc += 2;
			
			while (curLoc < tmp_line.length()) {
				r_site += tmp_line[curLoc];
				curLoc++;
			}
			exonEnd = stoi(r_site) - 1;

			exonCnt += (exonEnd - exonStart + 1); // Counting exon

			for (int j = exonStart; j <= exonEnd; j++) {
				codingRegion[j] = 1;
			}
		}
	}

	void syncInfoInHMM_CDS(int& exonNumInHMM, int& exonCntInnHMM, vector<int>& codingRegionInHMM) {
		exonNumInHMM = exonNum;
		exonCntInnHMM = exonCnt;
		codingRegionInHMM = codingRegion;
	}
};