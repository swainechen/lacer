#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "sam.h"
#include "sam_header.h"
#include "htslib/sam.h"
#include "htslib/kseq.h"
#include "lacepr.h"

const char *version = "0.2";
// from bam_import.c
// char *bam_nt16_rev_tabl            = "=ACMGRSVTWYHKDBN";
// the appropriate complement would be
// char *bam_nt16_rev_tabl_complement = "=TGKCYSBAWRDMHVN";
// we're only going to accept "real" nts
const char *nt16_table      = ".AC.G...T.......";
const char *nt16_table_comp = ".TG.C...A.......";
const char* const KNOWN_RGFIELD[] = { "PU", "LB", "SM" };

KSEQ_INIT(gzFile, gzread);
cache_t *cache;

void init_cache (int n) {
  int i, j, k, l;
  cache = realloc(cache, n * sizeof(cache_t));
  for (i=0; i < n; i++) {
    for (j=0; j < MAX_REASONABLE_Q+1; j++) {
      for (k=0; k < MAX_CONTEXT+1; k++) {
        for (l=0; l < MAX_CYCLE+1; l++) {
          cache[i][j][k][l] = 0;
        }
      }
    }
  }
}

long choose (int n, int k) {
  long result = 1;
  int j = 1;
  if (k == n || k == 0) {
    return 1;
  }
  if (k > n || k < 0 || n < 1) {
    return 0;
  }
  if (n - k < k) {
    k = n - k;
  }
  while (j < k) {
    result *= n/j;
    n--;
    j++;
  }
}

double log10binomial (double kd, int n, double p) {
  int k;
  double mean;
  double stdev;
  k = lround(kd);
//printf("kd %f, n %d, p %f\n", kd, n, p);
  if (k >= 0 && n > 0 && p > 0 && p < 1) {
    if (n > 100) {
      // use normal approximation
      mean = n * p;
      stdev = sqrt(mean * (1 - p));
      return (-0.5 * ( (k - mean) * (k - mean) / n / p / (1-p) ) / log(10) - 0.5 * log10(2 * M_PI * n * p * (1-p)));
    } else {
      return ( log10(choose(n, k)) + k * log10(p) + (n-k) * log10(1-p) );
    }
  }
  return 0;
}

int delta (int origq, int obs, double err) {
  int i;
  int maxi = -1;
  double temp;
  double max = -DBL_MAX;
  double p[MAX_REASONABLE_Q + 1];
  for (i = MIN_Q; i <= MAX_REASONABLE_Q; i++) {
    if (abs(origq - i) > MAX_USABLE_Q) {
      p[i] = prior[MAX_USABLE_Q];
    } else {
      p[i] = prior[abs(origq - i)];
    }
    temp = log10binomial(err, obs, pow(10, -i/10.0));
    if (temp != 0) {
      p[i] += temp;
    }
    if (p[i] > max) {
      maxi = i;
      max = p[i];
    }
//printf ("origq %d, i %d, diff %f temp %f pi %f\n", origq, i, prior[abs(origq - i)], temp, p[i]);
  }
  if (maxi >= MIN_Q && maxi <= MAX_REASONABLE_Q) {
    return (maxi - origq);
  } else {
    return 0;
  }
}

int newq (int8_t origq, int cycle, char preceding, char current, recal_t *recaldata, int recaldata_index) {
  int nq = 40;
  int running_q;
  int adjust_q = 0;
  int adjust_cyc = 0;
  int adjust_con = 0;
  char *rg_search = NULL;
  int i;
  char context[3];
  context[0] = preceding;
  context[1] = current;
  context[2] = '\0';
  int con_i = get_context_index(context);
  int cyc_i = get_cycle_index(cycle);

  if (cache[recaldata_index][origq][con_i][cyc_i] > 0) {
    return cache[recaldata_index][origq][con_i][cyc_i];
  }

  running_q = recaldata[recaldata_index].Quality;
  if (recaldata[recaldata_index].OrigQual[origq].Observations > 0) {
    adjust_q = delta(running_q, recaldata[recaldata_index].OrigQual[origq].Observations, recaldata[recaldata_index].OrigQual[origq].Errors);
//printf("base %d, obs %d, err %f, adjust_q %d\n", running_q, recaldata[recaldata_index].OrigQual[origq].Observations, recaldata[recaldata_index].OrigQual[origq].Errors, adjust_q);
  }
  running_q += adjust_q;
 
  if (recaldata[recaldata_index].Context[con_i].OrigQual[origq].Observations > 0) {
    adjust_con = delta(running_q, recaldata[recaldata_index].Context[con_i].OrigQual[origq].Observations, recaldata[recaldata_index].Context[con_i].OrigQual[origq].Errors);
//printf("base %d, context %s, index %d, obs %d, err %f, adjust_q %d\n", running_q, context, con_i, recaldata[recaldata_index].Context[con_i].OrigQual[origq].Observations, recaldata[recaldata_index].Context[con_i].OrigQual[origq].Errors, adjust_con);
  }
  if (recaldata[recaldata_index].Cycle[get_cycle_index(cycle)].OrigQual[origq].Observations > 0) {
    adjust_cyc = delta(running_q, recaldata[recaldata_index].Cycle[cyc_i].OrigQual[origq].Observations, recaldata[recaldata_index].Cycle[cyc_i].OrigQual[origq].Errors);
//printf("base %d, cycle %d, index %d, obs %d, err %f, adjust_q %d\n", running_q, cycle, cyc_i, recaldata[recaldata_index].Cycle[cyc_i].OrigQual[origq].Observations, recaldata[recaldata_index].Cycle[cyc_i].OrigQual[origq].Errors, adjust_cyc);
  }
  nq = running_q + adjust_con + adjust_cyc;
//  printf ("orig %d prior %d, cycle %d, context %s, adjust %d %d %d, new %d\n\n\n", origq, recaldata[recaldata_index].Quality, cycle, context, adjust_q, adjust_con, adjust_cyc, nq);
  if (nq == 0) {
    nq = origq;
  }
  if (nq < MIN_Q) {
    nq = MIN_Q;
  }
  if (nq > MAX_REASONABLE_Q) {
    nq = MAX_REASONABLE_Q;
  }
  cache[recaldata_index][origq][con_i][cyc_i] = nq;
  return nq;
}

int get_rg_index (char** rglist, char* rg, bool insert) {
  // we build the list of read groups (if insert) and help access them here
  // we always track the end of the array with a null in the next element
  int i;
  if (rglist[0] == NULL) {
    if (insert) {
      rglist[0] = malloc(MAX_FIELD * sizeof(char));
      strcpy(rglist[0], rg);
      rglist[1] = NULL;
      return(0);
    } else {
      return(-1);
    }
  } else {
    for (i=0; i < MAX_RG; i++) {
      if (rglist[i] == NULL) {
        if (insert) {
          rglist[i] = malloc(MAX_FIELD * sizeof(char));
          strcpy(rglist[i], rg);
	  if (i+1 < MAX_RG) { rglist[i+1] = NULL; }
	  return(i);
        } else {
          return(-1);
        }
      }
      if (strcmp(rg, rglist[i]) == 0) {
        return(i);
      }
    }
  }
  // we should never actually get here
  return(-1);
}

// convert bases to integer encoding as per htslib / samtools
int nt2int (char s) {
  switch (s) {
    case 'A': return(1);
    case 'C': return(2);
    case 'G': return(4);
    case 'T': return(8);
    default : return(0);
  }
}

// we're really only going to assume there is 2 base context
int get_context_index (char *s) {
  int i = 0;
  int sum = 0;
  while (s[i] != '\0') {
    sum *= 4;
    if (s[i] == 'A') {
      // this is 0
    } else if (s[i] == 'C') {
      sum += 1;
    } else if (s[i] == 'G') {
      sum += 2;
    } else if (s[i] == 'T') {
      sum += 3;
    } else {
      return 0;	// bad data
    }
    i++;
  }
  sum++;	// we want to give 1-16, not 0-15
  return sum;
}

// cycle here should be an integer, positive or negative
// encode negatives as even numbers
// 0 is bad data
int get_cycle_index (int i) {
  if (i > 0 && i < MAX_CYCLE) {
    return (i*2-1);
  } else if (i < 0 && -i < MAX_CYCLE) {
    return (-i*2);
  }
  return 0;	// bad data
}

// assume we just get an unallocated pointer, make memory and initialize
recal_t* init_recal (int num_rg) {
  int i, j, k;
  recal_t *data = malloc(sizeof(recal_t) * num_rg);
  for (i=0; i < num_rg; i++) {
    data[i].Quality = 0;
    data[i].Observations = 0;
    data[i].Errors = 0;
    for (j=0; j < MAX_Q+1; j++) {
      data[i].OrigQual[j].Quality = 0;
      data[i].OrigQual[j].Observations = 0;
      data[i].OrigQual[j].Errors = 0;
    }
    for (j=0; j < MAX_CONTEXT; j++) {
      for (k=0; k < MAX_Q+1; k++) {
        data[i].Context[j].OrigQual[k].Quality = 0;
        data[i].Context[j].OrigQual[k].Observations = 0;
        data[i].Context[j].OrigQual[k].Errors = 0;
      }
    }
    for (j=0; j < MAX_CYCLE; j++) {
      for (k=0; k < MAX_Q+1; k++) {
        data[i].Cycle[j].OrigQual[k].Quality = 0;
        data[i].Cycle[j].OrigQual[k].Observations = 0;
        data[i].Cycle[j].OrigQual[k].Errors = 0;
      }
    }
  }
  return data;
}

int read_recal (char* file, char** rglist, recal_t *data) {
  char *in = 0;
  size_t n = 0;
  int num_rg = 0;
  ssize_t bytes;
  char *tkn;
  int i;
  int parsed;
  char rg[MAX_FIELD];
  char covariate[MAX_FIELD];
  char covname[MAX_FIELD];
  double empqual;
  double estqual;
  int qual;
  char eventtype[4];
  int obs;
  double err;
  int rgindex;
  int covindex;
  point_recal_t *temp_table0;

  FILE *r = fopen(file, "r");
  temp_table0 = malloc(MAX_RG * sizeof(point_recal_t));

  while ( (bytes = getline(&in, &n, r)) != -1 ) {
    if (bytes > 0) {
      tkn = strtok(in, ":");
      if (tkn[0] == '#') {
        while (tkn != NULL) {
          tkn = strtok(NULL, ":");
          if (tkn != NULL) {
            if (strcmp(tkn, "RecalTable0") == 0) {
              // RecalTable0 should be:
              // 0 - ReadGroup
              // 1 - EventType
              // 2 - EmpiricalQuality
              // 3 - EstimatedQReported
              // 4 - Observations
              // 5 - Errors
              bytes = getline(&in, &n, r);
              while ((bytes = getline(&in, &n, r)) != -1 && !whitespaceline(in)) {
                parsed = sscanf(in, "%s %s %lf %lf %d %lf", rg, eventtype, &empqual, &estqual, &obs, &err);
                if (parsed == 6) {
                  rgindex = get_rg_index(rglist, rg, true);
                  temp_table0[rgindex].Quality = lround(empqual);
                  temp_table0[rgindex].Observations = obs;
                  temp_table0[rgindex].Errors = err;
//printf("rg %s index %d obs %d err %f\n", rg, rgindex, data[rgindex].Observations, data[rgindex].Errors);
//printf("rg %s index %d obs %d err %f\n", rg, rgindex, obs, err);
                } else {
                  fprintf(stderr, "Parse error for RecalTable0 on line:\n%s", in);
                }
              }
	      // after this first RecalTable0 we should know the read groups
	      // and we can allocate memory for the real data structure
	      for (i = 0; i < MAX_RG; i++) {
                if (rglist[i] == NULL) {
                  num_rg = i;
                  break;
                }
              }
	      if (num_rg == 0) {
                fprintf(stderr, "Couldn't find any read groups in the recal file\n");
		return(0);
              }
              data = realloc(data, sizeof(recal_t) * num_rg);
	      for (i = 0; i < num_rg; i++) {
                data[rgindex].Quality = temp_table0[rgindex].Quality;
                data[rgindex].Observations = temp_table0[rgindex].Observations;
                data[rgindex].Errors = temp_table0[rgindex].Errors;
              }
            } else if (strcmp(tkn, "RecalTable1") == 0) {
              // RecalTable1 should be:
              // 0 - ReadGroup
              // 1 - QualityScore
              // 2 - EventType
              // 3 - EmpiricalQuality
              // 4 - Observations
              // 5 - Errors
              bytes = getline(&in, &n, r);
              while ((bytes = getline(&in, &n, r)) != -1 && !whitespaceline(in)) {
                parsed = sscanf(in, "%s %d %s %lf %d %lf", rg, &qual, eventtype, &empqual, &obs, &err);
                if (parsed == 6 && strcmp(eventtype,"M") == 0) {
                  rgindex = get_rg_index(rglist, rg, false);
                  data[rgindex].OrigQual[qual].Quality = lround(empqual);
                  data[rgindex].OrigQual[qual].Observations = obs;
                  data[rgindex].OrigQual[qual].Errors = err;
                } else {
                  fprintf(stderr, "Parse error for RecalTable1 on line:\n%s", in);
                }
//                printf("Table 1: %d, %s, %d, %s, %f, %d, %f\n", parsed, rg, qual, eventtype, empqual, data[rgindex].OrigQual[qual].Observations, data[rgindex].OrigQual[qual].Errors);
              }
            } else if (strcmp(tkn, "RecalTable2") == 0) {
              // RecalTable2 should be:
              // 0 - ReadGroup
              // 1 - QualityScore
              // 2 - CovariateValue
              // 3 - CovariateName
              // 4 - EventType
              // 5 - EmpiricalQuality
              // 6 - Observations
              // 7 - Errors
              bytes = getline(&in, &n, r);
              while ((bytes = getline(&in, &n, r)) != -1 && !whitespaceline(in)) {
                parsed = sscanf(in, "%s %d %s %s %s %lf %d %lf", rg, &qual, covariate, covname, eventtype, &empqual, &obs, &err);
                if (parsed == 8 && strcmp(eventtype,"M") == 0) {
                  rgindex = get_rg_index(rglist, rg, false);
                  if (strcmp(covname, "Cycle") == 0) {
                    // covariate is a string but if it's cycle it's a number
                    if ((covindex = get_cycle_index(atoi(covariate))) != -1) {
                      data[rgindex].Cycle[covindex].OrigQual[qual].Quality = lround(empqual);
                      data[rgindex].Cycle[covindex].OrigQual[qual].Observations = obs;
                      data[rgindex].Cycle[covindex].OrigQual[qual].Errors = err;
//                  printf("Table 2: %d, %s, %d, %s, %s, %s, %f, %d, %f\n", parsed, rg, qual, covariate, covname, eventtype, empqual, data[rgindex].Cycle[covindex].OrigQual[qual].Observations, data[rgindex].Cycle[covindex].OrigQual[qual].Errors);
                    }
                  } else if (strcmp(covname, "Context") == 0) {
                    if ((covindex = get_context_index(covariate)) != -1) {
                      data[rgindex].Context[covindex].OrigQual[qual].Quality = lround(empqual);
                      data[rgindex].Context[covindex].OrigQual[qual].Observations = obs;
                      data[rgindex].Context[covindex].OrigQual[qual].Errors = err;
//                    printf("Table 2: %d, %s, %d, %s, %s, %s, %f, %d, %f\n", parsed, rg, qual, covariate, covname, eventtype, empqual, data[rgindex].Context[covindex].OrigQual[qual].Observations, data[rgindex].Context[covindex].OrigQual[qual].Errors);
                    }
                  }
                
                } else {
                  fprintf (stderr, "Parse error for RecalTable2 on line:\n%s", in);
                }
              }
            }
          }
        }
      }
    }
  }
  return(num_rg);
}

int read_group_check (samfile_t *fp, char** rglist, int num_rg, rg_item_t* rg_data[3][MAX_RG], char* rg_field) {
  // read in read groups and collect info to map back later
  // rg_data has dimensions 3 (for field, use KNOWN_RGFIELD) and MAX_RG
  // rg_item_t has ID, Value (strings)
  // we'll track size by making the last item rg_data[i][j] = NULL
  int rg_ok = -1; // this will return the KNOWN_RGFIELD index that we can use
  bool rg_check[3];
  void *iter;
  const char *key, *val;
  int i;
  int j;

  rg_ok = -1;
  for(i = 0; i < 3; i++) {
    iter = sam_header_parse2(fp->header->text);
    j = 0;
    while (iter = sam_header2key_val(iter, "RG", "ID", KNOWN_RGFIELD[i], &key, &val)) {
      rg_data[i][j] = malloc(MAX_FIELD * sizeof(rg_item_t));
      strcpy(rg_data[i][j]->ID, key);
      strcpy(rg_data[i][j]->Value, val);
      j++;
      if (j >= MAX_RG) { break; }
    }
    if (j < MAX_RG) {
      rg_data[i][j] = malloc(sizeof(rg_item_t));
      rg_data[i][j]->ID[0] = '\0';
      rg_data[i][j]->Value[0] = '\0';
    }
  }
  for(i = 0; i < 3; i++) {
    rg_check[i] = true;
    for(j = 0; j < MAX_RG; j++) {
      if (rg_data[i][j]->Value[0] == '\0') {
        if (j == 0) { rg_check[i] = false; }
        break;
      } else if (get_rg_index(rglist, rg_data[i][j]->Value, false) < 0) {
        rg_check[i] = false;
      }
    }
    if (rg_check[i] && strcmp(KNOWN_RGFIELD[i], rg_field) == 0) {
      rg_ok = i;
    }
  }
  if (rg_ok < 0) {
    for(i = 0; i < 3; i++) {
      if (rg_check[i]) {
        fprintf(stderr, "RG field specified (%s) doesn't seem to match between recal and bam files\nIt looks like using --field %s may work instead\n", rg_field, KNOWN_RGFIELD[i]);
        break;
      }
    }
    if (i == 3) {
      fprintf(stderr, "Can't seem to find any match for read groups.\n");
      fprintf(stderr, "In the recal file, I see:\n");
      for(j = 0; j < num_rg; j++) {
        if (rglist[j] == NULL) { break; }
        fprintf(stderr, "  %s\n", rglist[j]);
      }
      fprintf(stderr, "In the bam file, I see:\n");
      for(i = 0; i < 3; i++) {
        fprintf(stderr, "- RG Field %s\n", KNOWN_RGFIELD[i]);
        for(j = 0; j < MAX_RG; j++) {
          if (rg_data[i][j]->Value[0] == '\0') { break; }
          fprintf(stderr, "  %s\n", rg_data[i][j]->Value);
        }
      }
    }
    return(-1);
  } else {
    return(rg_ok);
  }
}

int main (int argc, char *argv[])
{

  // handle command line and options, print help if needed
  int getopt_c;
  char *inbam;
  char *infastq;
  char *outfile;
  char *recal_file;
  char *use_rg = NULL;
  char *rg_field = "PU";
  int read_pairnum = 1;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"field", required_argument, 0, 'f'},
      {"bam",   required_argument, 0, 'b'},
      {"fastq", required_argument, 0, 'q'},
      {"pair",  required_argument, 0, 'p'},
      {"rg",    required_argument, 0, 'g'},
      {"recal", required_argument, 0, 'r'},
      {"out",   required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };
    getopt_c = getopt_long_only(argc, argv, "f:i:r:o:", long_options, &option_index);
    if (getopt_c == -1)
      break;
    switch (getopt_c) {
      case 'f': rg_field = optarg; break;
      case 'b': inbam = optarg; break;
      case 'q': infastq = optarg; break;
      case 'p': read_pairnum = atoi(optarg); break;
      case 'g': use_rg = optarg; break;
      case 'r': recal_file = optarg; break;
      case 'o': outfile = optarg; break;
    }
  }
  if(strlen(rg_field) > 2) { rg_field[2] = '\0'; }

  if ( ( access(inbam, F_OK) != 0 && access(infastq, F_OK) != 0 ) ||
       access(recal_file, F_OK) != 0 ) {
    fprintf(stderr, "Lacepr version %s; headers %s\n", version, header_version);
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  lacepr [ options ] --bam <in.bam> --recal <recal file> --out <out.bam>\n");
    fprintf(stderr, "  lacepr [ options ] --fastq <in fastq.gz> --recal <recal file> --out <out fastq.gz>\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  --field PU|LB|SM  for bam input, header field with read group information\n");
    fprintf(stderr, "                    Default PU.\n");
    fprintf(stderr, "  --rg <string>     for bam/fastq input, force usage of data for this read group\n");
    fprintf(stderr, "                    from the recalibration file.\n");
    fprintf(stderr, "  --pair 1|2        for fastq input, specify R1/R2 read. Default 1.\n");
    fprintf(stderr, "\nTakes recalibration file from lacer or GATK, applies it to a bam or fastq file\n");
    fprintf(stderr, "Will output bam if the input is bam, fastq if the input is fastq\n");
    fprintf(stderr, "Can force use of a read group in the recalibration file with --rg.\n");
    fprintf(stderr, "When processing fastq files, will default to the first read group in the recalibation file if --rg is not specified.\n");
    fprintf(stderr, "By default, will output gzipped fastq if processing a fastq file\n");
    fprintf(stderr, "\nNOTE: For clarity, only long options (--field, --fastq, etc.) are allowed\n");
    return 1;
  }

  // main declarations
  bam1_t *b;
  uint8_t *rawdata;
  int8_t *sequence;
  uint8_t *quality;
  uint8_t *aux;
  int bytes;
  uint8_t *rg;
  char *rg_const = "foo";
  char **rglist;
  int num_rg = 0;
  char *temp;
  char preceding;
  char current;
  int cycle;
  int i;
  int j;
  int debug = 0;
  int qlen;
  int max_length = 0;
  int max_aux = 0;
  int lines = 0;
  int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
  recal_t *recaldata;
  char *q_pointer;
  uint8_t *q_temp;
  int good_rgfield_index;

  rg_item_t *rg_data[3][MAX_RG];
  char *forceable_rg;
  char *rgID;
  char *rg_search;
  int rg_index;
  int force_rg_index = -1;

  sequence = 0;
  quality = NULL;
  rg = 0;

  // read in the recal table
  rglist = malloc(MAX_RG * sizeof(char*));
  rglist[0] = NULL;
  recaldata = init_recal(1);
  num_rg = read_recal(recal_file, rglist, recaldata);
  if (use_rg != NULL) {
    force_rg_index = get_rg_index(rglist, use_rg, false);
    if (force_rg_index == -1) {
      fprintf(stderr, "Can't find data in recalibration file %s for read group %s, aborting.\n", recal_file, use_rg);
    }
  }
  init_cache(num_rg);

  if (access(inbam, F_OK) == 0) {	// we are processing a bam file
    samfile_t *fp = samopen(inbam, "rb", 0);
    samfile_t *outfp = samopen(outfile, "wb", fp->header);
    bam_header_t *bh = fp->header;
    // set up rg_data so we can map from IDs back to read groups, which then
    // will get matched in rglist to choose the right recalibration data
    if (force_rg_index == -1) {
      good_rgfield_index = read_group_check(fp, rglist, num_rg, rg_data, rg_field);
      if (good_rgfield_index == -1) {
        return(-1);
      }
    }

    b = bam_init1();
    while ((bytes = samread(fp, b)) > 0) {
      qlen = b->core.l_qseq;
      if (qlen + 1 > max_length) {
        max_length = qlen + 1;
        kroundup32(max_length);
        sequence = realloc(sequence, max_length);
//      quality = realloc(quality, max_length);
      }
      sequence[qlen] = '\0';
//    quality[qlen] = '\0';
      rawdata = bam_get_seq(b);
      uint8_t backup[b->l_data];
      for(i=0; i < b->l_data; i++) {
        backup[i] = b->data[i];
      }
      int offset = ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1);
      quality = backup + offset;

      // deal with read groups again
      if (force_rg_index == -1) {
        rg = bam_aux_get(b, "RG");
        rgID = malloc(sizeof(char) * strlen(rg));
        if (rg[0] == 90) {
          strcpy(rgID, rg + 1);
        }
        // rg comes from the bam alignment structure as an ID
        // the ID is listed in the @RG lines in the sam header
        // This should then go to one of the RG fields (PU, SM, LB)
        // Then we can pick up where this is in rglist
        for (i = 0; i < MAX_RG; i++) {
          if (rg_data[good_rgfield_index][i]->Value[0] == '\0') {
            break;
          }
          if (strcmp(rg_data[good_rgfield_index][i]->ID, rgID) == 0) {
            rg_search = malloc(sizeof(char) * (strlen(rg_data[good_rgfield_index][i]->Value)+1));
            strcpy(rg_search, rg_data[good_rgfield_index][i]->Value);
            break;
          }
        }

        if (rg_search == NULL) {
          fprintf(stderr, "Can't find read group for ID %s; not recalibrating this\n", rgID);
        } else {
          rg_index = get_rg_index(rglist, rg_search, false);
        }
      }

      for (i = 0; i < qlen; i++) {
        sequence[i] = bam_seqi(rawdata, i);
      }
      // this is how they access qualities in bam.h
      q_pointer = b->data + b->core.n_cigar * 4 + b->core.l_qname + ((b->core.l_qseq+1)>>1);
      for (i = 0; i < qlen; i++) {
        // get context and cycle
        if (b->core.flag & 16) {	// mapped to reverse strand
          current = nt16_table_comp[sequence[i]];
          if (i < qlen-1) {
            preceding = nt16_table_comp[sequence[i+1]];
          } else {
            preceding = '.';
          }
          cycle = qlen - i;
        } else {	// mapped to forward strand
          current = nt16_table[sequence[i]];
          if (i > 0) {
            preceding = nt16_table[sequence[i-1]];
          } else {
            preceding = '.';
          }
          cycle = i + 1;
        }
        // take care of cycle for first / second read
        if (b->core.flag & 128) {	// second of a pair
          cycle = -cycle;
        }
        debug && fprintf(stderr, "rgID %s, cycle %d, sequence %d, pre %c, cur %c, orig %d, new %d\n", rgID, cycle, sequence[i], preceding, current, quality[i], newq(quality[i], cycle, preceding, current, recaldata, rg_index));
	if (force_rg_index == -1) {
          q_pointer[i] = newq(quality[i], cycle, preceding, current, recaldata, rg_index);
        } else {
          q_pointer[i] = newq(quality[i], cycle, preceding, current, recaldata, force_rg_index);
        }
      }
      samwrite(outfp, b);
    }
    bam_destroy1(b);
    samclose(outfp);
    samclose(fp);
    return(0);

  } else if (access(infastq, F_OK) == 0) {	// processing a fastq file

    gzFile fp;
    gzFile outp;
    kseq_t *seq;
    int l;
    char *newquality;
    newquality = malloc(MAX_FIELD * sizeof(char));

    if (read_pairnum != 1 && read_pairnum != 2) {
      read_pairnum = 1;	// default value if invalid
    }

    // we will not have read groups by default - take the first in the recal table
    if (force_rg_index != -1) {
      rg_index = force_rg_index;
    } else {
      rg_index = 0;
    }

    // read fastq and recalibrate
    // we are going to only output gzipped files
    fp = gzopen(infastq, "r");
    if (outfile[strlen(outfile)-2] != 'g' || outfile[strlen(outfile)-1] != 'z') {
      outfile = strcat(outfile, ".gz");
    }
    outp = gzopen(outfile, "w");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      sequence = seq->seq.s;
      qlen = strlen(sequence);
      quality = seq->qual.s;
      if (qlen > strlen(newquality)) {
        newquality = realloc(newquality, qlen);
      }
      strcpy(newquality, quality);
      for (i = 0; i < qlen; i++) {
        current = sequence[i];
        cycle = i+1;
        if (read_pairnum == 2) { cycle = -cycle; }
        if (i > 0) {
          preceding = sequence[i-1];
        } else {
          preceding = '.';
        }
        debug && fprintf(stderr, "cycle %d, sequence %d, pre %c, cur %c, orig %d, new %d\n", cycle, sequence[i], preceding, current, quality[i] - 33, newq(quality[i] - 33, cycle, preceding, current, recaldata, rg_index));
        newquality[i] = newq(quality[i] - 33, cycle, preceding, current, recaldata, rg_index) + 33;
      }
      if (seq->comment.s == NULL || strlen(seq->comment.s) < 1) {
        gzprintf(outp, "@%s\n%s\n+\n%s\n", seq->name.s, sequence, newquality);
      } else {
        gzprintf(outp, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, sequence, newquality);
      }
    }
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(outp);
  }
  return(0);
}
