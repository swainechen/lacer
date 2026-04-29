#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
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
int num_recal_rg = 0;

void init_cache (int n) {
  int i, j, k, l;
  cache_t *tmp_cache = realloc(cache, n * sizeof(cache_t));
  if (!tmp_cache && n > 0) {
    fprintf(stderr, "Memory allocation failed for cache\n");
    return;
  }
  cache = tmp_cache;
  for (i=0; i < n; i++) {
    for (j=0; j < MAX_Q+1; j++) {
      for (k=0; k < MAX_CONTEXT+1; k++) {
        for (l=0; l < MAX_CYCLE_BINS; l++) {
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
  return result;
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

int newq (uint8_t origq, int cycle, char preceding, char current, recal_t *recaldata, int recaldata_index) {
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

  if (origq > MAX_Q || recaldata_index < 0 || recaldata_index >= num_recal_rg || con_i < 0 || con_i > MAX_CONTEXT || cyc_i < 0 || cyc_i >= MAX_CYCLE_BINS) {
    return origq;
  }

  if (cache == NULL) {
    return origq;
  }

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
  if (recaldata[recaldata_index].Cycle[cyc_i].OrigQual[origq].Observations > 0) {
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
      if (rglist[0]) {
        strncpy(rglist[0], rg, MAX_FIELD - 1);
        rglist[0][MAX_FIELD - 1] = '\0';
      } else {
        return(-1);
      }
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
          if (rglist[i]) {
            strncpy(rglist[i], rg, MAX_FIELD - 1);
            rglist[i][MAX_FIELD - 1] = '\0';
          } else {
            return(-1);
          }
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
  if (strlen(s) > CONTEXT_LENGTH) { return 0; }
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
  if (i > 0 && i <= MAX_CYCLE) {
    return (i*2-1);
  } else if (i < 0 && i != INT_MIN && -i <= MAX_CYCLE) {
    return (-i*2);
  }
  return 0;	// bad data
}

// assume we just get an unallocated pointer, make memory and initialize
recal_t* init_recal (int num_rg) {
  recal_t *data = calloc(num_rg, sizeof(recal_t));
  if (!data && num_rg > 0) {
    fprintf(stderr, "Memory allocation failed for recal data\n");
    return NULL;
  }
  return data;
}

int read_recal (char* file, char** rglist, recal_t **data_ptr) {
  recal_t *data = *data_ptr;
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
  if (!r) {
    fprintf(stderr, "Cannot open recal file %s\n", file);
    return 0;
  }
  temp_table0 = calloc(MAX_RG, sizeof(point_recal_t));
  if (!temp_table0) {
    fprintf(stderr, "Memory allocation failed for temp_table0\n");
    fclose(r);
    return 0;
  }

  while ( (bytes = getline(&in, &n, r)) != -1 ) {
    if (bytes > 0) {
      tkn = strtok(in, ":");
      if (tkn != NULL && tkn[0] == '#') {
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
                // Limit input lengths to prevent buffer overflows
                parsed = sscanf(in, "%255s %3s %lf %lf %d %lf", rg, eventtype, &empqual, &estqual, &obs, &err);
                if (parsed == 6) {
                  rgindex = get_rg_index(rglist, rg, true);
                  if (rgindex < 0 || rgindex >= MAX_RG) {
                    fprintf(stderr, "Invalid read group index %d\n", rgindex);
                    continue;
                  }
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
                free(temp_table0);
                free(in);
                fclose(r);
		return(0);
              }
              recal_t *tmp_data = realloc(data, sizeof(recal_t) * num_rg);
              if (!tmp_data) {
                fprintf(stderr, "Memory allocation failed for recal data\n");
                free(temp_table0);
                free(in);
                fclose(r);
                return 0;
              }
              data = tmp_data;
              // The original data had 1 element (initialized in main)
              // Ensure newly added elements are zero-initialized
              if (num_rg > 1) {
                memset(data + 1, 0, sizeof(recal_t) * (num_rg - 1));
              }
              *data_ptr = data;
	      for (i = 0; i < num_rg; i++) {
                data[i].Quality = temp_table0[i].Quality;
                data[i].Observations = temp_table0[i].Observations;
                data[i].Errors = temp_table0[i].Errors;
              }
              break;
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
                // Limit input lengths to prevent buffer overflows
                parsed = sscanf(in, "%255s %d %3s %lf %d %lf", rg, &qual, eventtype, &empqual, &obs, &err);
                if (parsed == 6 && strcmp(eventtype,"M") == 0) {
                  rgindex = get_rg_index(rglist, rg, false);
                  if (rgindex < 0 || rgindex >= num_rg || qual < 0 || qual > MAX_Q) {
                    fprintf(stderr, "Invalid RG index %d or Quality %d\n", rgindex, qual);
                    continue;
                  }
                  data[rgindex].OrigQual[qual].Quality = lround(empqual);
                  data[rgindex].OrigQual[qual].Observations = obs;
                  data[rgindex].OrigQual[qual].Errors = err;
                } else {
                  fprintf(stderr, "Parse error for RecalTable1 on line:\n%s", in);
                }
//                printf("Table 1: %d, %s, %d, %s, %f, %d, %f\n", parsed, rg, qual, eventtype, empqual, data[rgindex].OrigQual[qual].Observations, data[rgindex].OrigQual[qual].Errors);
              }
              break;
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
                // Limit input lengths to prevent buffer overflows
                parsed = sscanf(in, "%255s %d %255s %255s %3s %lf %d %lf", rg, &qual, covariate, covname, eventtype, &empqual, &obs, &err);
                if (parsed == 8 && strcmp(eventtype,"M") == 0) {
                  rgindex = get_rg_index(rglist, rg, false);
                  if (rgindex < 0 || rgindex >= num_rg || qual < 0 || qual > MAX_Q) {
                    fprintf(stderr, "Invalid RG index %d or Quality %d\n", rgindex, qual);
                    continue;
                  }
                  if (strcmp(covname, "Cycle") == 0) {
                    // covariate is a string but if it's cycle it's a number
                    if ((covindex = get_cycle_index(atoi(covariate))) != -1 && covindex < MAX_CYCLE_BINS) {
                      data[rgindex].Cycle[covindex].OrigQual[qual].Quality = lround(empqual);
                      data[rgindex].Cycle[covindex].OrigQual[qual].Observations = obs;
                      data[rgindex].Cycle[covindex].OrigQual[qual].Errors = err;
//                  printf("Table 2: %d, %s, %d, %s, %s, %s, %f, %d, %f\n", parsed, rg, qual, covariate, covname, eventtype, empqual, data[rgindex].Cycle[covindex].OrigQual[qual].Observations, data[rgindex].Cycle[covindex].OrigQual[qual].Errors);
                    }
                  } else if (strcmp(covname, "Context") == 0) {
                    if ((covindex = get_context_index(covariate)) != -1 && covindex <= MAX_CONTEXT) {
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
              break;
            }
          }
        }
      }
    }
  }
  free(temp_table0);
  free(in);
  fclose(r);
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
    void *header_ptr = sam_header_parse2(fp->header->text);
    if (!header_ptr) continue;
    iter = header_ptr;
    j = 0;
    while (iter = sam_header2key_val(iter, "RG", "ID", KNOWN_RGFIELD[i], &key, &val)) {
      rg_data[i][j] = malloc(sizeof(rg_item_t));
      if (!rg_data[i][j]) {
        fprintf(stderr, "Memory allocation failed for rg_data\n");
        continue;
      }
      strncpy(rg_data[i][j]->ID, key, MAX_FIELD - 1);
      rg_data[i][j]->ID[MAX_FIELD - 1] = '\0';
      strncpy(rg_data[i][j]->Value, val, MAX_FIELD - 1);
      rg_data[i][j]->Value[MAX_FIELD - 1] = '\0';
      j++;
      if (j >= MAX_RG) { break; }
    }
    if (j < MAX_RG) {
      rg_data[i][j] = malloc(sizeof(rg_item_t));
      if (rg_data[i][j]) {
        rg_data[i][j]->ID[0] = '\0';
        rg_data[i][j]->Value[0] = '\0';
      }
    }
    sam_header_free(header_ptr);
  }
  for(i = 0; i < 3; i++) {
    rg_check[i] = true;
    for(j = 0; j < MAX_RG; j++) {
      if (rg_data[i][j] == NULL || rg_data[i][j]->Value[0] == '\0') {
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
          if (rg_data[i][j] == NULL || rg_data[i][j]->Value[0] == '\0') { break; }
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
  char *inbam = NULL;
  char *infastq = NULL;
  char *outfile = NULL;
  char *recal_file = NULL;
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

  if ( !recal_file || !outfile || ( !inbam && !infastq ) ||
       ( (inbam ? access(inbam, F_OK) : -1) != 0 && (infastq ? access(infastq, F_OK) : -1) != 0 ) ||
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
  for(i=0; i<3; i++) {
    for(j=0; j<MAX_RG; j++) {
      rg_data[i][j] = NULL;
    }
  }
  char *forceable_rg;
  char *rgID;
  char *rg_search;
  int rg_index = -1;
  int force_rg_index = -1;

  sequence = 0;
  quality = NULL;
  rg = 0;

  // read in the recal table
  rglist = calloc(MAX_RG, sizeof(char*));
  if (!rglist) {
    fprintf(stderr, "Memory allocation failed for rglist\n");
    return 1;
  }
  recaldata = init_recal(1);
  if (!recaldata) {
    free(rglist);
    return 1;
  }
  num_rg = read_recal(recal_file, rglist, &recaldata);
  if (num_rg <= 0) {
    fprintf(stderr, "Failed to read recalibration file or no read groups found\n");
    for (i = 0; i < MAX_RG; i++) {
      if (rglist[i]) free(rglist[i]);
    }
    free(rglist);
    if (recaldata) free(recaldata);
    return 1;
  }
  num_recal_rg = num_rg;
  if (use_rg != NULL) {
    force_rg_index = get_rg_index(rglist, use_rg, false);
    if (force_rg_index == -1) {
      fprintf(stderr, "Can't find data in recalibration file %s for read group %s, aborting.\n", recal_file, use_rg);
    }
  }
  init_cache(num_rg);

  if (inbam && access(inbam, F_OK) == 0) {	// we are processing a bam file
    samfile_t *fp = samopen(inbam, "rb", 0);
    if (!fp) {
      fprintf(stderr, "Cannot open BAM file %s\n", inbam);
      return 1;
    }
    samfile_t *outfp = samopen(outfile, "wb", fp->header);
    if (!outfp) {
      fprintf(stderr, "Cannot open output BAM file %s\n", outfile);
      samclose(fp);
      return 1;
    }
    bam_header_t *bh = fp->header;
    // set up rg_data so we can map from IDs back to read groups, which then
    // will get matched in rglist to choose the right recalibration data
    if (force_rg_index == -1) {
      good_rgfield_index = read_group_check(fp, rglist, num_rg, rg_data, rg_field);
      if (good_rgfield_index == -1) {
        samclose(outfp);
        samclose(fp);
        return(-1);
      }
    }

    b = bam_init1();
    while ((bytes = samread(fp, b)) > 0) {
      rg_index = -1;
      qlen = b->core.l_qseq;
      if (qlen + 1 > max_length) {
        max_length = qlen + 1;
        kroundup32(max_length);
        int8_t *tmp_sequence = realloc(sequence, max_length);
        if (!tmp_sequence) {
          fprintf(stderr, "Memory allocation failed for sequence\n");
          break; // Exit the loop if we can't allocate memory
        } else {
          sequence = tmp_sequence;
        }
//      quality = realloc(quality, max_length);
      }
      sequence[qlen] = '\0';
//    quality[qlen] = '\0';
      rawdata = bam_get_seq(b);
      // Removed stack-allocated backup array to prevent stack overflow
      int offset = ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1);
      quality = b->data + offset;

      // deal with read groups again
      rg_search = NULL;
      rgID = NULL;
      if (force_rg_index == -1) {
        rg = bam_aux_get(b, "RG");
        if (rg && rg[0] == 'Z') {
          rgID = (char*)(rg + 1);
        }

        if (rgID) {
          for (i = 0; i < MAX_RG; i++) {
            if (rg_data[good_rgfield_index][i] == NULL || rg_data[good_rgfield_index][i]->Value[0] == '\0') {
              break;
            }
            if (strcmp(rg_data[good_rgfield_index][i]->ID, rgID) == 0) {
              rg_search = rg_data[good_rgfield_index][i]->Value;
              break;
            }
          }
        }

        if (rg_search == NULL) {
          fprintf(stderr, "Can't find read group for ID %s; not recalibrating this\n", rgID ? rgID : "NULL");
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

    // Cleanup dynamically allocated resources
    for (i = 0; i < MAX_RG; i++) {
      if (rglist[i]) free(rglist[i]);
    }
    free(rglist);
    free(recaldata);
    free(cache);
    for (i = 0; i < 3; i++) {
      for (j = 0; j < MAX_RG; j++) {
        if (rg_data[i][j]) free(rg_data[i][j]);
      }
    }
    if (sequence) free(sequence);

    return(0);

  } else if (infastq && access(infastq, F_OK) == 0) {	// processing a fastq file

    gzFile fp;
    gzFile outp;
    kseq_t *seq;
    int l;
    char *newquality;
    size_t newquality_size = MAX_FIELD;
    newquality = malloc(newquality_size);

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
    if (!fp) {
      fprintf(stderr, "Cannot open input FastQ file %s\n", infastq);
      return 1;
    }
    size_t outlen = strlen(outfile);
    if (outlen < 3 || strcmp(outfile + outlen - 3, ".gz") != 0) {
      char *new_outfile = malloc(outlen + 4);
      if (new_outfile == NULL) {
        fprintf(stderr, "Memory allocation failed for outfile name\n");
        return(1);
      }
      strcpy(new_outfile, outfile);
      strcat(new_outfile, ".gz");
      outfile = new_outfile;
    }
    outp = gzopen(outfile, "w");
    if (!outp) {
      fprintf(stderr, "Cannot open output FastQ file %s\n", outfile);
      gzclose(fp);
      return 1;
    }
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      char *seq_ptr = seq->seq.s;
      char *qual_ptr = seq->qual.s;
      if (seq->seq.l != seq->qual.l) {
        fprintf(stderr, "Sequence and quality lengths differ for %s; skipping\n", seq->name.s);
        continue;
      }
      qlen = (int)seq->seq.l;
      if (qlen + 1 > newquality_size) {
        newquality_size = qlen + 1;
        char *tmp_q = realloc(newquality, newquality_size);
        if (!tmp_q) {
          fprintf(stderr, "Memory allocation failed for newquality\n");
          free(newquality);
          kseq_destroy(seq);
          gzclose(fp);
          gzclose(outp);
          return 1;
        }
        newquality = tmp_q;
      }
      if (newquality) {
        memcpy(newquality, qual_ptr, qlen);
        newquality[qlen] = '\0';
      }
      for (i = 0; i < qlen; i++) {
        current = seq_ptr[i];
        cycle = i+1;
        if (read_pairnum == 2) { cycle = -cycle; }
        if (i > 0) {
          preceding = seq_ptr[i-1];
        } else {
          preceding = '.';
        }
        debug && fprintf(stderr, "cycle %d, sequence %d, pre %c, cur %c, orig %d, new %d\n", cycle, seq_ptr[i], preceding, current, qual_ptr[i] - 33, newq(qual_ptr[i] - 33, cycle, preceding, current, recaldata, rg_index));
        newquality[i] = newq(qual_ptr[i] - 33, cycle, preceding, current, recaldata, rg_index) + 33;
      }
      if (seq->comment.s == NULL || seq->comment.l < 1) {
        gzprintf(outp, "@%s\n%s\n+\n%s\n", seq->name.s, seq_ptr, newquality);
      } else {
        gzprintf(outp, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq_ptr, newquality);
      }
    }
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(outp);
    free(newquality);
  }

  // Cleanup dynamically allocated resources
  for (i = 0; i < MAX_RG; i++) {
    if (rglist[i]) free(rglist[i]);
  }
  free(rglist);
  free(recaldata);
  free(cache);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < MAX_RG; j++) {
      if (rg_data[i][j]) free(rg_data[i][j]);
    }
  }

  return(0);
}
