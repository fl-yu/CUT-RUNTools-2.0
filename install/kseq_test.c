#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1 || argc==2 || argc==3) {
		fprintf(stderr, "Usage: %s <in.seq> <read.length> <outfile>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	int rd_len = atoi(argv[2]);
	seq = kseq_init(fp);
	int ind = 0;
	char seq2[200];
	char a[] = "AGATCG";

	//outfile
	gzFile fp2;
	fp2 = gzopen(argv[3], "wb");

	//master branch
	int default_cut_index = 6;
	if(rd_len>=149 && rd_len<=151){
		default_cut_index = 100;
	}else if(rd_len>=74 && rd_len<=76){
		default_cut_index = 35;
	}else{
		default_cut_index = 6;
	}


	while ((l = kseq_read(seq)) >= 0) {
		//printf("name: %s\n", seq->name.s);
		//if (seq->comment.l)
			//printf("comment: %s\n", seq->comment.s);
		//printf("seq: %s ", seq->seq.s);
		int cut_index = -1;
		if(strlen(seq->seq.s)==rd_len){
			cut_index = default_cut_index;			
/*			
			for(ind = 6; ind>=1; ind--){ //previously ind>=3
				a[ind] = '\0';
				if(!strcmp(seq->seq.s + strlen(seq->seq.s) - ind, a)){
					cut_index = ind;
					//printf("CUT");
					break;
				}
			}
*/
		}

		//printf("\n");
		//if (seq->qual.l)
		//	printf("qual: %s\n", seq->qual.s);
		gzprintf(fp2, "@%s %s\n", seq->name.s, seq->comment.s);
		if(cut_index==-1){
			gzprintf(fp2, "%s\n", seq->seq.s);
			gzprintf(fp2, "+\n");
			gzprintf(fp2, "%s\n", seq->qual.s);
		}else{
			//printf("seq: %s CUT %d\n", seq->seq.s, cut_index);
			strncpy(seq2, seq->seq.s, strlen(seq->seq.s) - cut_index);
			seq2[strlen(seq->seq.s) - cut_index] = '\0';
			gzprintf(fp2, "%s\n", seq2);
			gzprintf(fp2, "+\n");
			strncpy(seq2, seq->qual.s, strlen(seq->seq.s) - cut_index);
			seq2[strlen(seq->seq.s) - cut_index] = '\0';
			gzprintf(fp2, "%s\n", seq2);
		}		
	}

	//printf("return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	gzclose(fp2);
	return 0;
}
